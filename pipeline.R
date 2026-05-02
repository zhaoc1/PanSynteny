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
#args <- c("example.yaml")
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
})

source("config.R")
source("model.R")
source("graph.R")
source("path.R")
source('neighbor.R')
source('midas.R')
source('blocks.R')
source('plot.R')
source("parse.R")

load_job_config(args[1])

species_id <- cfg_get(job_config, "species_id")
my_trait <- cfg_get(job_config, "trait")

c80_tables <- load_c80_tables(get_target("clusters_80_updated"), get_target("genes_info"))
cluster_80 <- c80_tables$cluster_80
gene_to_c80 <- c80_tables$gene_to_c80
rm(c80_tables)

# Pipeline input: focal_c80_df is produced by Step 0 (`prepare.R`): it filters
# `corrected_genes` to species/trait, applies the |cor_to_b| thresholds, and
# writes the result here. This driver consumes `is_focal` already encodes the 
# focal selection decision. Required columns: focal_c80, focal_label, is_focal, 
# cor_to_b, beta, sample_prevalence, [trait, genome_counts].
focal_fp <- get_target("gene_meta")
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


###### STEP 1: per-focal neighborhoods -> small-ORF + length-variant labels
# Step 1 is encapsulated in run_step1_neighbor_extraction (neighbor.R). It
# runs the three sub-stages: per-focal extraction, cross-genome assembly,
# and label attachment (synthetic small-ORF labels + length variants).
# Output: three RDS caches, and returns gene_neighbors, short_gene_prevalence, and c80_variants_mapping. 
res <- run_step1_neighbor_extraction(focal_c80_df, gene_to_c80)
gene_neighbors <- res$gene_neighbors
short_gene_prevalence <- res$short_gene_prevalence
c80_variants_mapping <- res$c80_variants_mapping


###### STEP 2: per-genome operon graphs -> maximal paths
# Step 2 is encapsulated in run_step2_path_stitching (graph.R). It calls
# stitch_paths_across_focal_genes, derives the path_genome_comp join key,
# Output: path_df.rds + esupport_df.rds, and returns path_df.
path_df <- run_step2_path_stitching(gene_neighbors)


###### STEP 3: cross-genome consolidation -> three granularity levels
# Step 3 is encapsulated in run_step3_consolidation (path.R). It builds all
# six L1/L2/L3 frames, decorates with small-ORF + truncation/fragmentation
# flags, writes the five canonical-paths TSVs, and returns the three data frames.
# Step 4 needs in memory (c_paths, collapsed_paths, c80s_coarse).
res <- run_step3_consolidation(path_df, c80_variants_mapping, focal_c80_df, cluster_80, short_gene_prevalence)
c_paths <- res$c_paths
collapsed_paths <- res$collapsed_paths
c80s_coarse <- res$c80s_coarse


###### STEP 4: focal block extraction + representative ranking
# Step 4 is encapsulated in run_step4_block_extraction (blocks.R). It mines
# c80s_coarse for runs of focal genes, aggregates equivalent blocks across
# canonical paths, picks a dominant block per (component, path_type) as
# reference, ranks the others, drops subset-redundant blocks, walks reps
# back to per-genome attribution, prints the rep-overlap diagnostic, and
# (optionally) renders rep_heatmap.pdf. Outputs the slim per-genome
# attribution table (rep_slim) for downstream flanking backbone analysis.
rep_slim <- run_step4_block_extraction(c80s_coarse, c_paths, collapsed_paths, path_df)


###### STEP 5 — summaries, selection, exemplar sampling, BLAST gene lists
# Step 5 is encapsulated in run_step5_parse (parse.R). The three Step 3
# TSVs are loaded here at the call site to keep the function inputs
# explicit; run_step5_parse itself builds summaries, applies the
# fine-coverage isoform filter, attaches the isoform map to
# coarse_summary, samples exemplar genomes per surviving fine isoform,
# enriches with per-isoform context, and writes the BLAST gene-id
# lists. Returns fine_long for inspection.
c80s_coarse <- read_delim(get_target("canonical_paths_c80s"), delim = "\t", show_col_types = FALSE)
c80s_fine <- read_delim(get_target("canonical_paths_fine_c80s"), delim = "\t", show_col_types = FALSE)
per_genome <- read_delim(get_target("canonical_paths_per_genome"), delim = "\t", show_col_types = FALSE)
fine_long <- run_step5_parse(c80s_coarse, c80s_fine, per_genome, gene_neighbors)


###### STEP 6 — gggenes figures
# Step 6 is encapsulated in run_step6_figures (plot.R). The Step 5
# selection sets are loaded here at the call site to keep the function
# inputs explicit; c80s_coarse and c80s_fine are reused from the
# loads above. Renders global + per-component PDFs for each fill mode.
selected_coarse <- read_delim(get_target("parse_selected_coarse"), delim = "\t", show_col_types = FALSE)
selected_fine <- read_delim(get_target("parse_selected_fine"), delim = "\t", show_col_types = FALSE)
run_step6_figures(selected_coarse, selected_fine, c80s_coarse, c80s_fine)
