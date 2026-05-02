# ------------------------------------------------------------------------------
# prepare.R
#
# Step 0 — input setup for the strain-aware operon pipeline.
#
# Two responsibilities:
#   (1) Build focal_c80_df from corrected_genes and write it as a TSV at
#       get_target("gene_meta").
#   (2) Discover which focal centroids still need their per-focal neighbor
#       TSVs (an external preprocessing job materialises these), and write
#       the missing list to get_target("gene_list").
#
# The main pipeline (pipeline.R) consumes both products as-is. If any
# neighbor TSV is still missing, pipeline.R aborts with a pointer back to
# this script.
#
# Usage:
#   Rscript prepare.R <config.yaml>
#
# Always overwrites the destination TSVs. Cheap to re-run.
# ------------------------------------------------------------------------------

# --- Parse command-line arguments ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: prepare.R <config.yaml>")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

source("config.R")
source("model.R")

load_job_config(args[1])

species_id       <- cfg_get(job_config, "species_id")
my_trait         <- cfg_get(job_config, "trait")
score_col        <- cfg_get(job_config, "score_col")
inclusion_cutoff <- cfg_get(job_config, "inclusion_cutoff")
focal_cutoff     <- cfg_get(job_config, "focal_cutoff")

stopifnot(focal_cutoff >= inclusion_cutoff)

# Build focal_c80_df from corrected_genes. The trait-stat column named by
# `score_col` (typically "cor_to_b" or "beta") drives three decisions:
#   * inclusion_cutoff — minimum |score| to retain a gene at all.
#   * focal_cutoff     — minimum |score| to flag is_focal = TRUE
#                        (drives downstream Step 1 neighbor extraction).
#   * focal_label sign — sign(score) determines "pos" vs "neg".
# Both cor_to_b and beta survive in the output as separate columns; only
# the gating logic above is parameterised.
gene_meta <- readRDS(get_target("corrected_genes")) %>%
  filter(species_id == !!species_id, trait_model == my_trait)

if (!score_col %in% names(gene_meta)) {
  stop(sprintf("score_col='%s' not found in corrected_genes columns: %s",
               score_col, paste(names(gene_meta), collapse = ", ")))
}

gene_meta <- gene_meta %>%
  filter(abs(.data[[score_col]]) >= inclusion_cutoff) %>%
  mutate(is_focal = abs(.data[[score_col]]) >= focal_cutoff)

focal_c80_df <- gene_meta %>%
  mutate(beta = round(beta, 3), cor_to_b = round(cor_to_b, 3)) %>%
  mutate(focal_label = ifelse(.data[[score_col]] > 0, "pos", "neg")) %>%
  select(gene_id, sample_prevalence, cor_to_b, focal_label, beta,
         trait_model, genome_counts, is_focal) %>%
  dplyr::rename(focal_c80 = gene_id, trait = trait_model)

out_fp <- get_target("gene_meta")
write.table(focal_c80_df, out_fp, sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf(
  "Wrote %s: %d rows total, %d flagged is_focal (|%s| >= %s)\n",
  out_fp, nrow(focal_c80_df), sum(focal_c80_df$is_focal), score_col, focal_cutoff
))

# (2) Discover which focal centroids still need their per-focal neighbor
# TSVs. These are produced by an external (non-R) preprocessing job; this
# script only enumerates what's still missing. pipeline.R aborts at startup
# if any are still absent.
genes_focal <- unique(focal_c80_df %>% filter(is_focal) %>% .$focal_c80)
genes_missing <- character()
for (gid in genes_focal) {
  in_fp <- file.path(get_target("neighbor_list"), paste0(gid, ".tsv"))
  if (!file.exists(in_fp)) {
    genes_missing <- c(genes_missing, gid)
  }
}

gene_list_fp <- get_target("gene_list")
if (length(genes_missing) > 0) {
  write.table(genes_missing, gene_list_fp, sep = "\n",
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  cat(sprintf(
    "Missing %d/%d neighbor TSVs under %s. List written to %s.\n",
    length(genes_missing), length(genes_focal),
    get_target("neighbor_list"), gene_list_fp
  ))
  cat("Materialise the missing TSVs externally, then re-run prepare.R or pipeline.R.\n")
} else {
  if (file.exists(gene_list_fp)) file.remove(gene_list_fp)
  cat(sprintf("All %d focal centroids have neighbor TSVs. Ready to run pipeline.R.\n",
              length(genes_focal)))
}
