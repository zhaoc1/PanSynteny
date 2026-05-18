# ------------------------------------------------------------------------------
# prepare.R
#
# Step 0: input setup for the strain-aware operon pipeline.
#
# Four responsibilities:
#   (1) Snapshot the input YAML to the project directory for provenance
#       (`get_target("run_config")` -> {proj_dir}/step1_setup/run_config.yaml).
#   (2) Seed the local copy of clusters_80_updated from `data.clusters_80_updated`
#       into `get_target("clusters_80_updated")` on first run (overwrite=FALSE);
#       subsequent runs preserve any hand edits to the local copy.
#   (3) Build focal_c80_df from the user-provided focal_meta TSV (path
#       comes from `data.focal_meta` in the YAML); optionally apply the
#       |score_col| thresholds; cache the result inside step1_setup/
#       at `get_target("focal_meta")`. From there pipeline.R reads it.
#   (4) Discover which focal centroids still need their per-focal neighbor
#       TSVs (those are materialised by build_genome_catalog.py +
#       run_species.sh), and write the missing list to `get_target("gene_list")`.
#
# Threshold handling (from the `prepare:` YAML section):
#   - If `score_col` is empty/missing, focal_meta is passed through as-is.
#     In that mode the input MUST already contain an `is_focal` column.
#   - If `score_col` is set (non-empty), filter to rows where
#     |focal_meta[[score_col]]| >= inclusion_cutoff and set
#     is_focal = (|focal_meta[[score_col]]| >= focal_cutoff). This is the
#     "derive from raw scores" mode.
#
# Usage:
#   Rscript prepare.R <config.yaml>
#
# Always overwrites the destination files. Cheap to re-run.
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

species_id <- cfg_get(job_config, "species_id")

# (1) Snapshot <config.yaml> to {proj_dir}/step1_setup/run_config.yaml so the run
# directory carries the exact config used. 
cat("[1/4] Snapshotting config...\n")
invisible(file.copy(args[1], get_target("run_config"), overwrite = TRUE))
cat(sprintf("  %s -> %s\n", args[1], get_target("run_config")))

# (1b) Seed clusters_80_updated from the YAML-declared source into the local
# data_dir copy (overwrite = FALSE). 
# To refresh from the source, `rm` the local copy and re-run.
cat("\n[2/4] Seeding clusters_80_updated...\n")
clusters_src <- cfg_get(job_config, "clusters_80_updated")
if (is.null(clusters_src) || !is.character(clusters_src) || !nzchar(clusters_src)) {
  stop("data.clusters_80_updated is not set in <config.yaml>. Add the absolute path under data:.")
}
clusters_src <- path.expand(clusters_src)
clusters_dst <- get_target("clusters_80_updated")
if (!file.exists(clusters_dst)) {
  if (!file.exists(clusters_src)) {
    stop(sprintf("clusters_80_updated source not found: %s (from data.clusters_80_updated)", clusters_src))
  }
  file.copy(clusters_src, clusters_dst, overwrite = FALSE)
  cat(sprintf("  seeded from %s -> %s\n", clusters_src, clusters_dst))
} else {
  cat(sprintf("  already present -> %s\n", clusters_dst))
}

# (2) Build focal_c80_df from focal_meta (input) and cache the (possibly
# threshold-derived) result back as focal_meta (output, in step1_setup/).
# The YAML key data.focal_meta names the *input* file - wherever the user
# keeps it - and get_target("focal_meta") names the *cached output* under
# step1_setup/. pipeline.R reads the cached output.
cat("\n[3/4] Processing focal_meta...\n")
focal_fp <- cfg_get(job_config, "focal_meta")
if (is.null(focal_fp) || !is.character(focal_fp) || !nzchar(focal_fp)) {
  stop("data.focal_meta is not set in <config.yaml>. Add `focal_meta: \"/path/to/focal_meta.tsv\"` under data:.")
}
focal_fp <- path.expand(focal_fp)
if (!file.exists(focal_fp)) {
  stop(sprintf(
    "focal_meta not found at %s (resolved from data.focal_meta in the YAML).",
    focal_fp
  ))
}
focal_c80_df <- read_delim(focal_fp, delim = "\t", show_col_types = FALSE)

score_col <- cfg_get(job_config, "score_col")
if (!is.null(score_col) && is.character(score_col) && nzchar(score_col)) {
  inclusion_cutoff <- cfg_get(job_config, "inclusion_cutoff")
  focal_cutoff     <- cfg_get(job_config, "focal_cutoff")
  if (!(score_col %in% names(focal_c80_df))) {
    stop(sprintf(
      "score_col '%s' not present in focal_meta columns: %s",
      score_col, paste(names(focal_c80_df), collapse = ", ")
    ))
  }
  # NB: if the input already carries an is_focal column, the threshold
  # derivation BELOW OVERWRITES it. Loud warning so a hand-curated mix of
  # focal + context rows doesn't get silently reflagged by `|score_col|`.
  if ("is_focal" %in% names(focal_c80_df)) {
    warning(sprintf(
      "is_focal column in focal_meta is being overwritten by the prepare.score_col threshold ('%s' >= %s). To preserve the input is_focal, set prepare.score_col: \"\".",
      score_col, focal_cutoff
    ))
  }
  focal_c80_df <- focal_c80_df %>%
    filter(abs(.data[[score_col]]) >= inclusion_cutoff) %>%
    mutate(is_focal = abs(.data[[score_col]]) >= focal_cutoff)
  cat(sprintf(
    "  applied thresholds on |%s|: inclusion>=%s, focal>=%s -> %d rows (%d is_focal)\n",
    score_col, inclusion_cutoff, focal_cutoff,
    nrow(focal_c80_df), sum(focal_c80_df$is_focal, na.rm = TRUE)
  ))
} else {
  if (!("is_focal" %in% names(focal_c80_df))) {
    stop(sprintf(
      "focal_meta has no `is_focal` column and prepare.score_col is unset - provide one or the other (focal_meta columns: %s)",
      paste(names(focal_c80_df), collapse = ", ")
    ))
  }
  cat(sprintf("  pass-through (no score_col threshold): %d rows (%d is_focal)\n",
              nrow(focal_c80_df), sum(focal_c80_df$is_focal, na.rm = TRUE)))
}

out_fp <- get_target("focal_meta")
write_delim(focal_c80_df, out_fp, delim = "\t")
cat(sprintf("  cached -> %s\n", out_fp))

# (3) Discover which focal centroids still need their per-focal neighbor TSVs.
# Those TSVs are produced by build_genome_catalog.py + run_species.sh;
# this script only enumerates what's still missing. pipeline.R aborts at
# startup if any are still absent.
cat("\n[4/4] Checking neighbor TSVs...\n")
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
    "  missing %d/%d under %s\n",
    length(genes_missing), length(genes_focal),
    get_target("neighbor_list")
  ))
  cat(sprintf("  list written -> %s\n", gene_list_fp))
  cat(sprintf("  next: bash run_species.sh %s\n", args[1]))
} else {
  if (file.exists(gene_list_fp)) file.remove(gene_list_fp)
  cat(sprintf("  all %d focal centroids have neighbor TSVs - ready to run pipeline.R\n",
              length(genes_focal)))
}
