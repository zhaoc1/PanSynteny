# ------------------------------------------------------------------------------
# pipeline.R
#
# Strain-aware operon analysis pipeline.
#
# Author:   Chunyu Zhao <chunyu.zhao@gladstone.ucsf.edu>
# Created:  2025-07-14
# Updated:  2026-05-01
# ------------------------------------------------------------------------------

# --- Parse command-line arguments ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: pipeline.R <config.yaml>")
}

suppressPackageStartupMessages({
  library(pheatmap)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(gggenes)
  library(tidyverse)
  library(ggraph)
  library(tibble)
  library(gridExtra)
  library(viridis)
  library(igraph)
  library(data.table)   # also attached inside midas.R / neighbor.R / plot.R; preload here silences their startup conflict messages
})

source("R/config.R")
source("R/model.R")
source("R/graph.R")
source("R/path.R")
source("R/neighbor.R")
source("R/midas.R")
source("R/blocks.R")
source("R/plot.R")
source("R/parse.R")

load_job_config(args[1])

species_id <- cfg_get(job_config, "species_id")

# Pipeline input: focal_c80_df is produced by Step 0 (`prepare.R`): prepare.R
# reads the user-provided `focal_meta` TSV, optionally applies the |score_col|
# thresholds (prepare: section), and caches the derived focal table here. This
# driver consumes `is_focal` as already encoded. Required columns (minimum):
# focal_c80, focal_label, is_focal, gene_label. Optional carry-through columns:
# cor_to_b, beta, sample_prevalence, trait, genome_counts (consumed when present).
focal_fp <- get_target("focal_meta")
if (!file.exists(focal_fp)) {
  stop(sprintf("focal_c80_df not found at %s. Run `Rscript prepare.R %s` first.", focal_fp, args[1]))
}
focal_c80_df <- read_delim(focal_fp, delim = "\t", show_col_types = F)

# Re-check that every is_focal centroid still has a neighbor TSV. prepare.R
# is the canonical place that enumerates missing files; this guard exists so
# a user who runs pipeline.R directly without re-running prepare.R after a
# partial neighbor-TSV materialisation gets a loud, actionable error
# instead of a silently-incomplete gene_neighbors table.
genes_focal <- unique(focal_c80_df %>% filter(is_focal) %>% .$focal_c80)
genes_missing <- genes_focal[!file.exists(file.path(get_target("neighbor_list"), paste0(genes_focal, ".tsv")))]
if (length(genes_missing) > 0) {
  stop(sprintf(
    "%d/%d focal centroid(s) missing per-focal neighbor TSVs under %s. Run `Rscript prepare.R %s` to refresh the missing list and materialise them externally before re-running.",
    length(genes_missing), length(genes_focal),
    get_target("neighbor_list"), args[1]
  ))
}

c80_tables <- load_c80_tables(get_target("clusters_80_updated"), get_target("catalog_genes_info"))
cluster_80 <- c80_tables$cluster_80
gene_to_c80 <- c80_tables$gene_to_c80
rm(c80_tables)

cat("\n[1/6] Step 1 : per-focal neighborhood extraction...\n")
# run_step1_neighbor_extraction (neighbor.R) runs the three sub-stages:
# per-focal extraction, cross-genome assembly, and label attachment
# (synthetic small-ORF labels + length variants). Writes three RDS caches
# and returns gene_neighbors / short_gene_prevalence / c80_variants_mapping.
res <- run_step1_neighbor_extraction(focal_c80_df, gene_to_c80)
gene_neighbors <- res$gene_neighbors
short_gene_prevalence <- res$short_gene_prevalence
c80_variants_mapping <- res$c80_variants_mapping


cat("\n[2/6] Step 2 : per-genome operon graphs to maximal paths...\n")
# run_step2_path_stitching (graph.R) calls stitch_paths_across_focal_genes,
# derives the path_genome_comp join key. Writes path_df.rds + esupport_df.rds
# and returns path_df.
path_df <- run_step2_path_stitching(gene_neighbors)


cat("\n[3/6] Step 3 : cross-genome consolidation (three granularity levels)...\n")
# run_step3_consolidation (path.R) builds all six L1/L2/L3 frames, decorates
# with small-ORF + truncation/fragmentation flags, writes the five
# canonical-paths TSVs, and returns the three data frames. Step 6 (blocks)
# needs c_paths, collapsed_paths, c80s_coarse in memory.
res <- run_step3_consolidation(path_df, c80_variants_mapping, focal_c80_df, cluster_80, short_gene_prevalence)
c_paths <- res$c_paths
collapsed_paths <- res$collapsed_paths
c80s_coarse <- res$c80s_coarse


cat("\n[4/6] Step 4 : parse + summaries (selection, sampling, BLAST gene lists)...\n")
# run_step4_parse (parse.R) loads the three Step 3 TSVs at the call site to
# keep inputs explicit; itself it builds summaries, applies the
# fine-coverage isoform filter, attaches the isoform map to coarse_summary,
# samples exemplar genomes per surviving fine isoform, enriches with
# per-isoform context, and writes the BLAST gene-id lists. Returns
# fine_long for inspection.
c80s_coarse <- read_delim(get_target("canonical_paths_c80s"), delim = "\t", show_col_types = FALSE)
c80s_fine <- read_delim(get_target("canonical_paths_fine_c80s"), delim = "\t", show_col_types = FALSE)
per_genome <- read_delim(get_target("canonical_paths_per_genome"), delim = "\t", show_col_types = FALSE)
fine_long <- run_step4_parse(c80s_coarse, c80s_fine, per_genome, gene_neighbors)


cat("\n[5/6] Step 5 : gggenes figures (writes under step5_figures/)...\n")
# run_step5_figures (plot.R) : renders global + per-component PDFs for each
# fill mode. The selection sets are read here at the call site to keep the
# function inputs explicit; c80s_coarse and c80s_fine are reused from
# Step 4 above.
selected_coarse <- read_delim(get_target("parse_selected_coarse"), delim = "\t", show_col_types = FALSE)
selected_fine <- read_delim(get_target("parse_selected_fine"), delim = "\t", show_col_types = FALSE)
run_step5_figures(selected_coarse, selected_fine, c80s_coarse, c80s_fine)


cat("\n[6/6] Step 6 : focal block extraction (gated; writes under step6_blocks/)...\n")
# run_step6_blocks (blocks.R) : mines c80s_coarse for runs of focal genes,
# aggregates equivalent blocks across canonical paths, ranks per
# (component, path_type), drops subset-redundant blocks, walks reps back to
# per-genome attribution, prints the rep-overlap diagnostic, optionally
# renders rep_heatmap.pdf. Independent of Step 5; safe to skip.
if (isTRUE(cfg_get(job_config, "skip_block"))) {
  cat("    blocks.skip_block = true -> Step 6 skipped\n")
} else {
  rep_slim <- run_step6_blocks(c80s_coarse, c_paths, collapsed_paths, path_df)
}
