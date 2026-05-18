# ------------------------------------------------------------------------------
# midas.R
#
# MIDAS Database utilities: centroid_80 metadata loading and c80 label assignment
# for short / length-variant genes.
#
# Author:   Chunyu Zhao <chunyu.zhao@gladstone.ucsf.edu>
# Created:  2025-07-14
# ------------------------------------------------------------------------------

library(tidyverse)
library(dplyr)
library(data.table)
library(readr)
library(purrr)
library(stringr)


#' Load cluster_80 metadata and gene-to-c80 mapping
#'
#' Read the per-cluster centroid_80 metadata TSV and the gene-level info TSV,
#' normalize column names, and produce both the per-c80 metadata table and a
#' per-gene mapping with the parent c80's reference length and genome
#' prevalence attached.
#'
#' @param c80_fp Path to the clusters_80 TSV (species-scoped midasdb file;
#'   `get_target("clusters_80_updated")`).
#' @param genes_fp Path to the unified genome-catalog `genes_info.tsv` built
#'   by `build_genome_catalog.py` - header `gene_id, centroid_80,
#'   gene_length`. Use `get_target("catalog_genes_info")`. (The legacy
#'   midasdb `genes_annotated.tsv` works too, since the function only reads
#'   those three columns.)
#' @return A list with two elements: `cluster_80` (per-c80 metadata) and
#'   `gene_to_c80` (per-gene mapping with parent-c80 attributes). For genes
#'   whose `centroid_80` is not present in `cluster_80` (e.g. ECOR genes
#'   mapping to a c80 from a different species), `neighbor_c80_length_coarse`
#'   and `genome_prevalence` come through as NA.
#' @export
load_c80_tables <- function(c80_fp, genes_fp) {
  cluster_80 <- read_delim(c80_fp, delim = "\t", show_col_types = FALSE) %>%
    dplyr::rename(genome_prevalence = centroid_80_genome_prevalence,
                  neighbor_c80_length_coarse = centroid_80_gene_length) %>%
    select(c80:centroid_80_genome_counts)

  gene_to_c80 <- read_delim(genes_fp, delim = "\t", show_col_types = FALSE) %>%
    select(gene_id, centroid_80, gene_length) %>%
    unique() %>%
    left_join(
      cluster_80 %>% select(c80, neighbor_c80_length_coarse, genome_prevalence),
      by = c("centroid_80" = "c80")
    )

  list(cluster_80 = cluster_80, gene_to_c80 = gene_to_c80)
}


#' Compute small ORFs prevalence and assign synthetic c80_label
#'
#' Replace missing `neighbor_c80_coarse` values with synthetic, focal-scoped cluster
#' labels so that small ORFs can be treated as first-class
#' clusters in downstream pattern extraction rather than collapsing into a
#' single undifferentiated `NA` category.
#'
#' For each `focal_c80`, the function takes the distinct
#' `(neighbor_gene_type, neighbor_gene_length)` combinations among NA-c80 rows
#' and assigns each a synthetic label of the form
#' `"_<focal_c80>-<neighbor_gene_type>_<rank>"`, where `<rank>` is the ordinal
#' rank of `neighbor_gene_length` within the focal (ties broken by input
#' order). Annotated rows (those with a non-NA `neighbor_c80_coarse`) are left
#' unchanged. The leading underscore distinguishes synthetic labels from real
#' centroid_80 identifiers so the two never collide.
#'
#' In addition to writing synthetic labels, the function emits a
#' `neighbor_c80_genome_prevalence` column for short-gene rows holding a
#' **within-focal intra-short-gene proportion** - the fraction of distinct
#' `neighbor_gene_id`s that share the synthetic label, relative to all short
#' genes for that focal, negated (`-1 *`) as a sign-based marker so
#' downstream code can distinguish short-gene prevalence from genome-wide
#' cluster prevalence. Note: this column replaces the earlier behavior of
#' overwriting `sample_prevalence`; `sample_prevalence` is no longer touched
#' for short-gene rows.
#'
#' @details
#' **Scope of synthetic labels.** Because labels contain `focal_c80`, the same
#' physical short gene appearing next to two different focals receives two
#' different synthetic labels - intentional, since analyses are focal-scoped.
#' Do not compare synthetic labels across focals.
#'
#' **Tie-breaking in rank.** When multiple `(gene_type, length)` pairs share
#' the same length within a focal, `rank(..., ties.method = "first")` breaks
#' ties by input order. Labels remain unique, but the numeric suffix is
#' assigned globally within the focal rather than per gene type, so numbering
#' can look interleaved (e.g., `_A-phage_1`, `_A-IS_2`, `_A-phage_3`).
#'
#' **`neighbor_c80_genome_prevalence` semantics shift.** For annotated rows, this column
#' carries its original genome-wide meaning; for synthetic short-gene rows it
#' becomes an intra-focal proportion. Downstream consumers that mix both should
#' be aware of this.
#'
#' @export
compute_short_gene_prevalence <- function(gene_neighbors) {
  # step 1: isolate the short-gene rows
  short_genes <- gene_neighbors %>% filter(is.na(neighbor_c80_coarse))
  
  # step 2: build a distinct-label table
  # for short genes without centroid_80, we used the focal_gene + gene_type + gene_length as the unique id
  # synthetic label: _<focal_c80>-<gene_type>_<rank>
  c80_for_short_genes <- gene_neighbors %>% 
    filter(is.na(neighbor_c80_coarse)) %>% 
    select(focal_c80, neighbor_gene_type, neighbor_gene_length, neighbor_c80_coarse) %>% unique() %>%
    group_by(focal_c80) %>%
    mutate(
      gene_length_count = n_distinct(neighbor_gene_length),
      rank = if_else(gene_length_count > 1, rank(neighbor_gene_length, ties.method = "first"), NA_real_),
      neighbor_c80_new = if_else(
        !is.na(rank),
        paste0("_", focal_c80, "-", neighbor_gene_type, "_", rank),
        paste0("_", focal_c80, "-", neighbor_gene_type, "_", 1)
      )
    ) %>%
    ungroup() %>%
    select(focal_c80, neighbor_gene_type, neighbor_gene_length, neighbor_c80_new)
  
  # step3: attach label back to the full short-gene row
  short_genes <- short_genes %>% 
    left_join(c80_for_short_genes, by = c("focal_c80", "neighbor_gene_type", "neighbor_gene_length"))
  
  # step4: compute synthetic prevalence
  prevalence_for_short_genes <- short_genes %>%
    group_by(focal_c80, neighbor_c80_new) %>%
    summarise(ng = n_distinct(neighbor_gene_id), .groups = "drop") %>%
    group_by(focal_c80) %>%
    mutate(prevalence = ng / sum(ng)) %>%
    ungroup() %>%
    select(focal_c80, neighbor_c80_new, prevalence)
  
  # step 5: write the synthetic label into neighbor_c80_coarse (was NA) and
  # store the within-focal short-gene proportion in neighbor_c80_genome_prevalence
  # as a negative number (sign-flagged so downstream code can distinguish
  # short-gene prevalence from genome-wide cluster prevalence).
  short_genes <- short_genes %>%
    left_join(prevalence_for_short_genes, by = c("focal_c80", "neighbor_c80_new")) %>%
    mutate(neighbor_c80_genome_prevalence = -1 * prevalence, neighbor_c80_coarse = neighbor_c80_new) %>%
    select(all_of(colnames(gene_neighbors)))
  short_genes
}


#' Add small ORFs and compute their genome prevalence
#'
#' Identifies short genes (small ORFs) from `gene_neighbors` using
#' [compute_short_gene_prevalence()], appends them to the input table,
#' and constructs a genome-prevalence summary at the `neighbor_c80_coarse` level.
#'
#' @details
#' Small ORFs are computed via [compute_short_gene_prevalence()], which
#' returns a subset of `gene_neighbors` with associated
#' `neighbor_c80_genome_prevalence`. This function reshapes that output
#' into a distinct per-c80 prevalence table and appends the small ORFs
#' back to the original data.
#'
#' @export
assign_c80_to_short_genes <- function(gene_neighbors) {
  short_genes <- compute_short_gene_prevalence(gene_neighbors)
  short_gene_prevalence <- short_genes %>% 
    select(neighbor_c80_coarse, neighbor_c80_genome_prevalence) %>% distinct() %>%
    dplyr::rename(neighbor_c80_fine = neighbor_c80_coarse, genome_prevalence = neighbor_c80_genome_prevalence)
  
  gene_neighbors <- bind_rows(gene_neighbors %>% filter(!is.na(neighbor_c80_coarse)), short_genes)
  
  list(
    gene_neighbors = gene_neighbors,
    short_gene_prevalence = short_gene_prevalence
  )
}


#' Label length variants of annotated c80 clusters
#'
#' Build a compact mapping that disambiguates length variants within each
#' annotated centroid_80 cluster by appending a numeric rank suffix to the
#' cluster label. Clusters observed at a single length keep their original
#' label unchanged; clusters observed at multiple lengths get per-variant
#' labels of the form `"<neighbor_c80_coarse>_<rank>"`, where `<rank>` is the
#' ordinal position of the length (smallest = 1) within the cluster.
#'
#' Unlike [compute_short_gene_prevalence()], this function is the sibling operation
#' for **annotated** (non-NA) `neighbor_c80_coarse` values. It is globally scoped -
#' not per-focal - because real c80 clusters are defined once at database
#' build time, so their length variants are intrinsic properties of the
#' cluster rather than focal-dependent artifacts.
#'
#' @details
#' **Ordering dependency.** This function should run **after**
#' [compute_short_gene_prevalence()]. If it runs first, any `NA` `neighbor_c80_coarse`
#' rows are grouped together and can receive surprising labels such as
#' `"NA_1"`, `"NA_2"`, because `dplyr::group_by` treats `NA` as a valid group.
#' In the standard pipeline flow, NAs have already been replaced with
#' focal-scoped synthetic labels by the time this function runs, so this edge
#' case does not arise.
#'
#' **Tie-breaking is defensive.** The upstream `unique()` step collapses
#' duplicate `(c80, length)` pairs before ranking, so `rank(...)` never sees
#' tied lengths within a cluster. `ties.method = "first"` is specified for
#' safety but is not reached in normal input.
#'
#' **Global vs focal scope.** The resulting `neighbor_c80_fine` is a global
#' identifier - the same length variant of the same cluster receives the same
#' label regardless of which focal neighborhood it appears in. Downstream
#' code can therefore use `neighbor_c80_coarse` for cluster-level coarse grouping
#' (all lengths collapsed) or `neighbor_c80_fine` for length-sensitive
#' analysis (variants distinguished).
#'
#' @export
compute_c80_variants <- function(gene_neighbors) {
  c80_variants_mapping <- gene_neighbors %>% 
    select(neighbor_c80_coarse, neighbor_gene_length) %>% unique() %>%
    group_by(neighbor_c80_coarse) %>%
    mutate(
      neighbor_gene_length_count = n_distinct(neighbor_gene_length),
      rank = if_else(neighbor_gene_length_count > 1, rank(neighbor_gene_length, ties.method = "first"), NA_real_),
      neighbor_c80_fine = if_else(
        !is.na(rank),
        paste0(neighbor_c80_coarse, "_", rank),
        neighbor_c80_coarse
      )
    ) %>%
    select(-neighbor_gene_length_count, -rank) %>%
    ungroup() %>%
    select(neighbor_c80_fine, neighbor_gene_length, everything())
  return(c80_variants_mapping)
}
