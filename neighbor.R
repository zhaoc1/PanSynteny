# ------------------------------------------------------------------------------
# neighbor.R
#
# Gene-neighborhood analysis utilities: loading, positional filtering, and
# operon-pattern extraction.
#
# Author:   Chunyu Zhao <chunyu.zhao@gladstone.ucsf.edu>
# Created:  2025-07-14
# Updated:  2026-03-18
# ------------------------------------------------------------------------------

library(tidyverse)
library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(readr)
library(purrr)
library(stringr)
library(tools)
library(parallel)


#' Load gene neighbor annotations
#'
#' Read a tab-delimited gene-neighbor file and join centroid-80 annotation
#' metadata from `gene_to_c80`.
#'
#' The input file is expected to have seven tab-separated columns in this order:
#' `gene_member`, `neighbor_gene_id`, `neighbor_contig_id`,
#' `neighbor_gene_start`, `neighbor_gene_end`, `neighbor_gene_strand`,
#' and `neighbor_gene_type`.
#'
#' If `gene_length` is present in `gene_to_c80`, missing values are filled using
#' `neighbor_gene_end - neighbor_gene_start`. The resulting column is renamed
#' to `neighbor_gene_length`.
#'
#' @export
load_gene_neighbors <- function(nfp, gene_to_c80) {
  cols_neighbors <- c("gene_member", "neighbor_gene_id", "neighbor_contig_id", 
                      "neighbor_gene_start", "neighbor_gene_end", 
                      "neighbor_gene_strand", "neighbor_gene_type") 
  
  gene_neighbors <- read_delim(nfp, delim = "\t", col_names = cols_neighbors, show_col_types = F) 
  
  if (nrow(gene_neighbors) < 2) 
    return(data.frame())
  
  gene_neighbors <- gene_neighbors %>%
    left_join(gene_to_c80, by=c("neighbor_gene_id" = "gene_id")) %>%
    dplyr::rename(neighbor_c80_coarse = centroid_80, neighbor_c80_genome_prevalence = genome_prevalence)
  
  gene_neighbors <- gene_neighbors %>% 
    mutate(gene_length = ifelse(is.na(gene_length), neighbor_gene_end - neighbor_gene_start, gene_length)) %>% 
    dplyr::rename(neighbor_gene_length = gene_length)
  
  return(gene_neighbors)
}


#' Add relative gene positions within neighborhoods
#'
#' Parse numeric suffixes from `neighbor_gene_id`, identify the focal gene
#' for each `gene_member`, and compute each neighbor's position relative to
#' the focal gene. The result is then filtered to a symmetric positional
#' window and to neighborhoods with at least `min_operon_size` genes.
#'
#' This function assumes `neighbor_gene_id` values end in an underscore
#' followed by an integer suffix reflecting contig-adjacency order, for
#' example `"PROKKA_contig1_00123"` (suffix `123`).
#'
#' @export
compute_relative_positions <- function(gene_neighbors, upper_bound = 10, min_operon_size = 5) {
  gn <- gene_neighbors %>%
    group_by(gene_member) %>%
    mutate(neighbor_gene_suffix = as.integer(sub(".*_(\\d+)$", "\\1", neighbor_gene_id))) %>%
    arrange(neighbor_gene_suffix, .by_group = TRUE) %>%
    # focal suffix is the one where neighbor_gene_id == gene_member
    mutate(
      focal_suffix = neighbor_gene_suffix[neighbor_gene_id == gene_member][1],
      relative_position = neighbor_gene_suffix - focal_suffix
    ) %>%
    ungroup() %>% 
    select(all_of(c(colnames(gene_neighbors), "relative_position")))
  
  gn <- gn %>% filter(abs(relative_position) <= upper_bound)
  
  operon_size_tbl <- gn %>%
    count(gene_member, name = "operon_size") %>%
    filter(operon_size >= min_operon_size)
  gn <- gn %>% filter(gene_member %in% operon_size_tbl$gene_member)
  
  return(gn)
}


#' Compute operon size per gene member
#'
#' Counts the number of neighboring genes associated with each `gene_member`
#' and returns those that meet a minimum operon size threshold.
#' Rows are ordered by decreasing `operon_size`.
#'
#' @export
compute_operon_size <- function(gene_neighbors, min_operon_size) {
  operon_size_tbl <- gene_neighbors %>%
    count(gene_member, name = "operon_size") %>%
    filter(operon_size >= min_operon_size) %>%
    arrange(desc(operon_size))
  return(operon_size_tbl)
}


#' Filter focal-gene neighborhoods by flanking-neighbor support
#'
#' Filters a `gene_neighbors` table to retain only `gene_member`s with sufficient
#' neighboring-gene support and sufficient operon size.
#'
#' By default (`use_strict = NULL`), the function first evaluates whether
#' enough `gene_member`s satisfy a strict flank requirement of at least
#' `min_left_neighbors` genes on the left and at least `min_right_neighbors`
#' genes on the right. If at least `min_strict_members` `gene_member`s satisfy
#' this strict criterion, only those strict cases are retained. Otherwise, a
#' relaxed criterion is used, requiring sufficient neighbors on either the
#' left or the right side. Passing `use_strict = TRUE` or `FALSE` overrides
#' the auto-detection and forces the chosen mode; `min_strict_members` is
#' then ignored.
#'
#' After flank-based filtering, the function applies a second filter using
#' [compute_operon_size()] and retains only `gene_member`s whose operon size meets
#' `min_operon_size`.
#'
#' @export
filter_by_flanking_coverage <- function(gene_neighbors, min_strict_members,
                                        min_left_neighbors = 2L, min_right_neighbors = 2L,
                                        min_operon_size = 5L, use_strict = NULL) {

  if (nrow(gene_neighbors) == 0L) return(NULL)
  stopifnot(is.null(use_strict) || (is.logical(use_strict) && length(use_strict) == 1L && !is.na(use_strict)))

  flank_counts <- gene_neighbors %>%
    group_by(gene_member) %>%
    summarise(
      n_left  = sum(relative_position < 0L),
      n_right = sum(relative_position > 0L),
      .groups = "drop"
    ) %>%
    mutate(
      meets_strict  = n_left >= min_left_neighbors & n_right >= min_right_neighbors,
      meets_relaxed = n_left >= min_left_neighbors | n_right >= min_right_neighbors
    )

  if (is.null(use_strict)) {
    use_strict <- sum(flank_counts$meets_strict, na.rm = TRUE) >= min_strict_members
  }
  keep_members <- flank_counts %>%
    filter(if (use_strict) meets_strict else meets_relaxed) %>%
    pull(gene_member)
  
  filtered <- gene_neighbors %>% filter(gene_member %in% keep_members)
  if (nrow(filtered) == 0L) return(NULL)
  
  operon_size_tbl <- compute_operon_size(filtered, min_operon_size = min_operon_size)
  
  filtered %>%
    semi_join(operon_size_tbl, by = "gene_member") %>%
    mutate(gene_member = factor(gene_member, levels = operon_size_tbl$gene_member))
}

#' Remove sparse position columns
#'
#' Keep only position columns with sufficient non-missing coverage and return
#' the continuous span between the minimum and maximum valid positions.
#'
#' Position columns must be named like `position_-3`, `position_0`, `position_4`.
#' Non-position columns are always retained.
#'
#' @export
remove_na_positions <- function(df, focal_min_genomes, min_positions) {
  # Remove sparse neighbor positions and keeps only a dense, continuous region of positions with enough coverage across genes.
  position_cols <- grep("^position_", names(df), value = TRUE)
  pos_mat <- df[, position_cols, drop = FALSE]
  # Keep only positions with enough genomes
  non_na_counts <- colSums(!is.na(pos_mat))
  valid_position_cols <- names(non_na_counts)[non_na_counts >= focal_min_genomes]
  
  if (length(valid_position_cols) < min_positions) {
    return(NULL)
  }
  # Find continuous region
  valid_positions <- as.integer(sub("^position_", "", valid_position_cols))
  min_pos <- min(valid_positions)
  max_pos <- max(valid_positions)
  
  all_positions <- as.integer(sub("^position_", "", position_cols))
  final_position_cols <- position_cols[all_positions >= min_pos & all_positions <= max_pos]
  # Keep metadata columns too
  non_position_cols <- setdiff(colnames(df), position_cols)
  final_cols <- c(non_position_cols, final_position_cols)
  
  df[, final_cols, drop = FALSE]
}


#' Orient focal-gene-centered neighborhoods
#'
#' Orient each focal-gene-centered path independently using the nearest non-missing
#' flanking `neighbor_c80_coarse` labels on both sides of the focal gene, and retain
#' only oriented path groups supported by at least `focal_min_genomes` focal genes.
#'
#' For each `gene_member`, the function identifies the nearest non-`NA` neighbor
#' to the left (`relative_position < 0`) and the nearest non-`NA` neighbor to the
#' right (`relative_position > 0`). Only paths with valid flanking neighbors on
#' both sides are retained. Orientation is then determined from the immediate
#' flanking pair:
#'
#' - if `left_anchor <= right_anchor`, the path is kept as-is (`orientation = 1L`)
#' - otherwise, the path is reversed (`orientation = -1L`)
#'
#' After orientation, the full oriented path is represented by the ordered
#' sequence of non-`NA` `neighbor_c80_coarse` labels. Paths are then grouped by this
#' full oriented representation, and only groups supported by at least
#' `focal_min_genomes` distinct `gene_member`s are retained.
#'
#' @details
#' Paths lacking a valid non-`NA` flanking neighbor on either side are discarded.
#' A message is emitted reporting how many `gene_member`s were dropped for this
#' reason.
#'
#' **Known limitation: palindromic flanks.** When the immediate-flank pair is
#' identical on both sides (`left_anchor == right_anchor`), the `<=` test
#' returns `TRUE` in both observed orientations, so both collapse to
#' `orientation = 1L` and keep their own `asc_str`. Two observations of the
#' same palindromic-core operon therefore do not share a `canonical_path` and
#' will land in different groups. If this matters for your data, fall back to
#' comparing `asc_str` vs `desc_str` on ties.
#'
#' @export
orient_focal_gene_neighbors <- function(gene_neighbors, focal_min_genomes) {
  # Orientation is decided by the immediate non-NA flanking genes (walking
  # outward until a non-NA c80 is found). More robust than full-path
  # lexicographical order, which is sensitive to long/short neighbors.
  all_members <- n_distinct(gene_neighbors$gene_member)
  
  flanking_orientation <- gene_neighbors %>%
    filter(!is.na(neighbor_c80_coarse)) %>%
    group_by(gene_member) %>%
    arrange(relative_position, .by_group = TRUE) %>%
    summarise(
      asc_str = paste(neighbor_c80_coarse, collapse = "|"), desc_str = paste(rev(neighbor_c80_coarse), collapse = "|"),
      # closest non-NA gene to the left of the focal (largest negative position)
      left_anchor  = {
        neg <- relative_position[relative_position < 0]
        if (length(neg) == 0) NA_character_ else neighbor_c80_coarse[relative_position == max(neg)][1]
      },
      # closest non-NA gene to the right of the focal (smallest positive position)
      right_anchor = {
        pos <- relative_position[relative_position > 0]
        if (length(pos) == 0) NA_character_ else neighbor_c80_coarse[relative_position == min(pos)][1]
      },
      .groups = "drop"
    ) %>%
    # require valid non-NA flanking neighbors on both sides
    filter(!is.na(left_anchor), !is.na(right_anchor)) %>%
    mutate(
      orientation    = if_else(left_anchor <= right_anchor, 1L, -1L),
      canonical_path = if_else(orientation == 1L, asc_str, desc_str)
    ) %>%
    select(gene_member, orientation, canonical_path)
  
  n_dropped <- all_members - nrow(flanking_orientation)
  if (n_dropped > 0)
    message(n_dropped, " gene_member(s) dropped: no non-NA flanking gene on one or both sides")
  
  out <- gene_neighbors %>%
    inner_join(flanking_orientation, by = "gene_member") %>%
    mutate(relative_position = relative_position * orientation) %>%
    arrange(gene_member, relative_position)
  
  out <- out %>% mutate(group_excludeNA = as.integer(factor(canonical_path))) #<- could rename to oriented_path_id
  
  out <- out %>% group_by(group_excludeNA) %>% filter(n_distinct(gene_member) >= focal_min_genomes) %>% ungroup()
  out
}


#' Extract conserved gene neighborhood patterns
#'
#' Identify conserved gene neighborhood patterns from focal-gene-centered paths.
#' This function groups identical neighborhood configurations across relative
#' positions, retains the dominant pattern within each broader path group, and
#' returns a plotting-ready long-format representation.
#'
#' @details
#' This function reads the following values from `job_config` via [cfg_get()]
#' (set by [load_job_config()]):
#' \itemize{
#'   \item `focal_min_genomes`
#'   \item `min_positions`
#'   \item `gene_padding_bp`
#' }
#'
#' Missing `neighbor_c80_coarse` values are encoded as `"NA|length"` during processing
#' to preserve gene-length information.
#'
#' **Key invariants:**
#' \itemize{
#'   \item **Two resolution levels** exist throughout:
#'     \itemize{
#'       \item `group_excludeNA` — the *broad* orientation group (set upstream
#'         by [orient_focal_gene_neighbors()]), NA-agnostic.
#'       \item `group_includeNA` — the *detailed* pattern within a broad group,
#'         distinguishing unannotated genes by length.
#'     }
#'     Only the dominant detailed pattern within each broad group survives.
#'   \item **Column-name overload:** after the long-pivot stage, the
#'     `gene_member` column no longer means "focal-gene path ID" (its meaning
#'     at input). It instead holds the c80 label of each neighbor position.
#'   \item **`gene_start` / `gene_end` are synthetic plot coordinates**, not
#'     real genomic coordinates. They are laid out with a `gene_padding_bp`
#'     gap (from `job_config`) between adjacent genes for `gggenes`-style rendering.
#' }
#'
#' @export
extract_gene_neighbor_patterns <- function(gene_neighbors, gene_to_c80){
  focal_min_genomes <- cfg_get(job_config, "focal_min_genomes")
  min_positions <- cfg_get(job_config, "min_positions")
  gene_padding_bp <- cfg_get(job_config, "gene_padding_bp")

  # identifying conserved gene neighborhood patterns
  # reshape neighborhoods into a wide position-by-position signature
  df_wide_includeNA <- gene_neighbors %>%
    select(gene_member, relative_position, neighbor_c80_coarse, neighbor_gene_length) %>%
    mutate(value = ifelse(is.na(neighbor_c80_coarse), paste(neighbor_c80_coarse, neighbor_gene_length, sep = "|"), neighbor_c80_coarse)) %>%
    select(gene_member, relative_position, value) %>%
    pivot_wider(names_from = relative_position, values_from = value, names_prefix = "position_")
  # remove poorly supported positions
  df_wide_includeNA <- remove_na_positions(df_wide_includeNA, focal_min_genomes, min_positions)
  if (is.null(df_wide_includeNA)) return(NULL)
  
  # group identical signatures
  df_group_includeNA <- df_wide_includeNA %>% 
    group_by(across(-gene_member)) %>% 
    mutate(group_includeNA = cur_group_id()) %>% 
    ungroup()
  group_assignments <- df_group_includeNA %>% select(gene_member, group_includeNA) #<- for later use
  group_assignments <- left_join(group_assignments, unique(gene_neighbors %>% select(gene_member, group_excludeNA)), by=c("gene_member"))
  
  # For each broader path group, retain the most frequent exact pattern
  df_group_includeNA <- df_group_includeNA %>%
    left_join(group_assignments, by=c("gene_member", "group_includeNA"))
  # select the dominant detailed pattern within each broader path group
  dominant_groups <- df_group_includeNA %>%
    group_by(group_excludeNA, group_includeNA) %>%
    summarise(n_genomes = n(), .groups = "drop_last") %>%
    mutate(total_genomes = sum(n_genomes)) %>%
    slice_max(order_by = n_genomes, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # Keep the most prevalent includeNA path, especially for MGEs
  df_dominant <- df_group_includeNA %>%
    select(-gene_member) %>%
    left_join(
      dominant_groups %>%
        select(group_excludeNA, group_includeNA, total_genomes),
      by = c("group_excludeNA", "group_includeNA")
    ) %>%
    filter(!is.na(total_genomes)) %>%
    unique() %>%
    filter(total_genomes >= focal_min_genomes)
  if (nrow(df_dominant) == 0) return(NULL)
  
  df <- df_dominant %>% 
    pivot_longer(cols = starts_with("position_"), names_to = "position", values_to = "gene") %>%
    mutate(
      gene_member = ifelse(str_detect(gene, "\\|"), str_extract(gene, "^[^|]+"), gene),
      neighbor_gene_length = ifelse(str_detect(gene, "\\|"),as.integer(str_extract(gene, "(?<=\\|)\\d+$")),0L)
    ) %>%
    mutate(
      gene_member = ifelse(is.na(gene_member), "NA", gene_member), 
      neighbor_gene_length = replace_na(neighbor_gene_length, 0L)
    ) %>%
    mutate(position = as.integer(str_remove(position, "position_"))) %>%
    dplyr::rename(group = group_includeNA) %>% select(-group_excludeNA) %>%
    arrange(group, position) %>%
    select(-gene)
  
  gene_length_df <- gene_to_c80 %>% 
    filter(centroid_80 == gene_id) %>% 
    select(centroid_80, gene_length)
  
  df <- df %>% left_join(gene_length_df, by=c("gene_member" = "centroid_80"))
  df <- df %>% 
    mutate(
      neighbor_gene_length = if_else(
        neighbor_gene_length == 0 & !is.na(gene_length),
        gene_length,
        neighbor_gene_length
      )
    ) %>%
    select(-gene_length)
  
  # calculate gene start/end coordinates
  df <- df %>%
    group_by(group) %>%
    arrange(group, position) %>%
    mutate(
      gene_start = cumsum(lag(neighbor_gene_length, default = 0) + gene_padding_bp),
      gene_end = gene_start + neighbor_gene_length
    ) %>%
    ungroup()
  
  # reorder groups by total gene length
  df <- df %>%
    group_by(group) %>%
    mutate(total_length = sum(neighbor_gene_length, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(group = factor(group, levels = unique(group[order(total_length)])))
  
  # create new group labels with total_genomes
  group_labels <- df %>%
    distinct(group, total_genomes) %>%
    mutate(group_label = paste("group", group, ": ", total_genomes, " genomes", sep=""))
  
  # Join labels back to df
  df <- df %>% 
    left_join(group_labels, by = c("group", "total_genomes")) %>%
    left_join(gene_to_c80, by=c("gene_member" = "gene_id"))
  
  return(list(neighbor_groups = df, group_assignments = group_assignments %>% dplyr::rename(group = group_includeNA)))
}


#' Filter operon-size records by supported neighbor groups
#'
#' Filter an operon-size table to retain only `gene_member`s belonging to
#' neighborhood-pattern groups supported by at least `focal_min_genomes` genomes.
#'
#' This function uses the output of [extract_gene_neighbor_patterns()] and removes
#' low-support neighborhood groups before subsetting `operon_size_tbl`.
#'
#' @details
#' The function first filters `neighbor_groups` to retain only sufficiently
#' supported groups, then filters `group_assignments` to retain only
#' `gene_member`s assigned to those groups, and finally subsets
#' `operon_size_tbl` accordingly.
#'
#' @export
filter_neighbor_groups <- function(gene_neighbor_patterns, operon_size_tbl, focal_min_genomes) {
  # remove nearly singleton genome groups
  ga <- gene_neighbor_patterns[["group_assignments"]]
  gd <- gene_neighbor_patterns[["neighbor_groups"]]
  gd <- gd %>% filter(total_genomes >= focal_min_genomes)
  ga <- ga %>% filter(group %in% unique(gd$group))
  
  operon_size_tbl <- operon_size_tbl %>% 
    filter(gene_member %in% unique(ga$gene_member))
  return(operon_size_tbl)
}


#' Parse and analyze focal-gene-centered neighborhoods for one focal centroid
#'
#' Run the neighborhood-analysis workflow for a single focal `centroid_80`
#' target. The workflow loads neighbor data, computes relative positions, filters
#' neighborhoods by flanking coverage, orients focal-gene-centered paths,
#' extracts conserved neighborhood patterns, generates summary figures, filters
#' supported groups, selects prevalent operon sizes, regenerates selected
#' neighborhood patterns, and saves the final neighborhood-group table.
#'
#' @details
#' This function has side effects: it writes multiple PDF figures and one RDS
#' file to output directories obtained via [get_target()].
#'
#' @export
parse_gene_neighbor <- function(in_fp, focal_c80, gene_to_c80) {
  fig_dir <- get_target("neighbor_figures")

  focal_min_genomes <- cfg_get(job_config, "focal_min_genomes")
  focal_min_total_genomes <- cfg_get(job_config, "focal_min_total_genomes")
  min_operon_size <- cfg_get(job_config, "min_positions")
  upper_bound <- cfg_get(job_config, "upper_bound")
  min_left_neighbors <- cfg_get(job_config, "min_left_neighbors")
  min_right_neighbors <- cfg_get(job_config, "min_right_neighbors")
  use_strict <- cfg_get(job_config, "use_strict")
  min_group_proportion <- cfg_get(job_config, "min_group_proportion")
  coverage_warn_threshold <- cfg_get(job_config, "coverage_warn_threshold")

  gene_neighbors <- load_gene_neighbors(in_fp, gene_to_c80)
  gene_neighbors <- compute_relative_positions(gene_neighbors, upper_bound = upper_bound, min_operon_size = min_operon_size)
  gene_neighbors <- filter_by_flanking_coverage(gene_neighbors, min_strict_members = focal_min_genomes,
                                                min_left_neighbors, min_right_neighbors, min_operon_size, use_strict = use_strict)
  if (is.null(gene_neighbors) || nrow(gene_neighbors) == 0) return(NULL)
  operon_size_tbl <- compute_operon_size(gene_neighbors, min_operon_size)

  # First pass enumerates all patterns for diagnostic plots (fig1, fig2)
  oriented_gene_neighbors <- orient_focal_gene_neighbors(gene_neighbors, focal_min_genomes)
  gene_neighbor_patterns <- extract_gene_neighbor_patterns(oriented_gene_neighbors, gene_to_c80)
  if(is.null(gene_neighbor_patterns) || length(gene_neighbor_patterns) == 0) return(NULL)
  neighbor_groups <- gene_neighbor_patterns[["neighbor_groups"]]
  
  fig_fp <- file.path(fig_dir, paste0("fig1_operon_by_gene_", focal_c80, ".pdf"))
  plot_neighbor_groups(neighbor_groups, focal_c80, fill_by="gene_member",  fig_fp=fig_fp)
  fig_fp <- file.path(fig_dir, paste0("fig2_operon_by_c80_", focal_c80, ".pdf"))
  plot_neighbor_groups(neighbor_groups, focal_c80, fill_by="centroid_80", fig_fp=fig_fp)
  
  operon_size_tbl <- filter_neighbor_groups(gene_neighbor_patterns, operon_size_tbl, focal_min_genomes)
  if(is.null(operon_size_tbl) || nrow(operon_size_tbl) == 0) return(NULL)
  
  fig_fp <- file.path(fig_dir, paste0("fig3_operon_dist_", focal_c80, ".pdf"))
  operon_size_tbl_selected <- find_most_prevalent_operon_size(operon_size_tbl, fig_fp=fig_fp)
  gene_neighbors_selected <- gene_neighbors %>%
    semi_join(operon_size_tbl_selected, by = "gene_member") %>%
    mutate(gene_member = factor(gene_member, levels = operon_size_tbl_selected$gene_member))
  if(is.null(gene_neighbors_selected)) return(NULL)
  
  # Second orient+extract pass reruns the same pipeline on the subset of members with the most prevalent
  # operon size, producing tighter groupings for the final output figures (fig4, fig5) and the saved RDS.
  oriented_gene_neighbors <- orient_focal_gene_neighbors(gene_neighbors_selected, focal_min_genomes)
  gene_neighbor_patterns <- extract_gene_neighbor_patterns(oriented_gene_neighbors, gene_to_c80)
  if(is.null(gene_neighbor_patterns) || length(gene_neighbor_patterns) == 0) return(NULL)
  
  neighbor_groups <- gene_neighbor_patterns[["neighbor_groups"]]
  
  fig_fp <- file.path(fig_dir, paste0("fig4_selected_operon_by_gene_", focal_c80, ".pdf"))
  plot_neighbor_groups(neighbor_groups, focal_c80, fill_by="gene_member", fig_fp=fig_fp)
  fig_fp <- file.path(fig_dir, paste0("fig5_selected_operon_by_prev_", focal_c80, ".pdf"))
  plot_neighbor_groups(neighbor_groups, focal_c80, fill_by="genome_prevalence", fig_fp=fig_fp)
  
  # filter neighbor groups by per group genome counts
  total_counts <- neighbor_groups %>%
    distinct(group, total_genomes) %>%
    pull(total_genomes) %>%
    sum()
  if (total_counts <= focal_min_total_genomes) return(NULL)

  groups_to_keep <- neighbor_groups %>%
    distinct(group, total_genomes) %>%
    mutate(genome_proportion = total_genomes / total_counts) %>%
    filter(total_genomes >= focal_min_genomes | genome_proportion >= min_group_proportion) %>%
    select(group, total_genomes)
  coverage <- sum(groups_to_keep$total_genomes) / total_counts
  if (coverage < coverage_warn_threshold) warning(paste0("Filtered groups only cover ", round(coverage * 100, 1), "% of total for ", focal_c80))
  
  neighbor_groups <- neighbor_groups %>% 
    inner_join(groups_to_keep, by = c("group", "total_genomes"))
  if(is.null(neighbor_groups) || nrow(neighbor_groups) == 0) return(NULL)
  
  gene_neighbors_grouped <- gene_neighbors_selected %>% 
    left_join(gene_neighbor_patterns[["group_assignments"]], by = "gene_member") %>%
    left_join(groups_to_keep %>% mutate(group = as.numeric(as.character(group))), by=c("group")) %>%
    filter(!is.na(total_genomes))
  
  gene_neighbors_grouped <- gene_neighbors_grouped %>%
    select(-one_of(c("group", "group_excludeNA")))
  
  out_dir <- get_target("neighbor_groups_by_focal")
  saveRDS(neighbor_groups, file.path(out_dir, paste(focal_c80, "_neighbor_groups.rds", sep="")))
  
  return(gene_neighbors_grouped)
}


#' Write gene-neighbor rows split into per-genome TSVs
#'
#' Split `df` by `gene_member_genome` and write one tab-separated file per
#' genome into `<outdir>/<genome>/<focal_c80>.tsv`. The per-genome subdirectory
#' is created on demand. The `gene_member_genome` column is retained in each
#' written TSV (via `.keep = TRUE` on `group_walk`), so the output is readable
#' by [load_neighbors_across_genomes()] without further column inference.
#'
#' @export
write_slices_by_genome <- function(df, outdir, focal_c80){
  df %>%
    group_by(gene_member_genome) %>%
    group_walk(~{
      genome <- .y$gene_member_genome[[1]]
      subdir <- file.path(outdir, genome)
      dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
      fp <- file.path(subdir, paste0(focal_c80, ".tsv"))
      write.table(.x, fp, sep = "\t", quote = FALSE, row.names = FALSE)
    }, .keep = TRUE) %>%
    invisible()
}


#' Write focal-gene neighborhoods as per-genome TSVs
#'
#' Derive `gene_member_genome` from `gene_member` (by stripping the trailing
#' `_<integer>` suffix), warn if any genome contributes more than one gene
#' copy to the focal cluster (a paralog signal), and write the resulting table
#' to disk as one TSV per genome under
#' `<neighbor_groups_by_genome>/<focal_label>/<genome>/<focal_c80>.tsv`. The
#' actual sharded writes are delegated to [write_slices_by_genome()].
#'
#' @details
#' **Paralog warning.** If two or more distinct `gene_member` values map to
#' the same `gene_member_genome`, the focal cluster has multiple gene copies
#' in at least one genome (paralogs, multi-copy MGE elements, or tandem
#' repeats). A `warning()` is emitted naming the focal; downstream output is
#' still produced and includes all copies.
#'
#' **Derivation assumption.** `sub("_[^_]+$", "", gene_member)` assumes the
#' Prokka-style naming convention where the last underscore separates the
#' genome ID from the per-gene index. If `gene_member` ever uses a different
#' format, the derived `gene_member_genome` will be incorrect and per-genome
#' slicing will mis-group rows.
#'
#' @export
write_gene_neighbor <- function(gene_neighbors, focal_c80, focal_label) {
  # Slice by gene_member_genome
  genomes_contain_focal_gene <- gene_neighbors %>% 
    select(gene_member) %>% distinct() %>% mutate(gene_member_genome = sub("_[^_]+$", "", gene_member))
  gene_neighbors <- gene_neighbors %>% 
    left_join(genomes_contain_focal_gene, by=c("gene_member")) %>% select(gene_member_genome, everything())
  
  has_dup <- genomes_contain_focal_gene %>% group_by(gene_member_genome) %>% filter(n() > 1) %>% nrow() > 0
  if (has_dup) warning(paste0("Multiple genes from the same genome correspond to focal centroid: ", focal_c80))
  
  out_dir <- get_target("neighbor_groups_by_genome")
  dir.create(file.path(out_dir, focal_label), recursive = T, showWarnings=F)
  write_slices_by_genome(gene_neighbors, file.path(out_dir, focal_label), focal_c80)
}


#' Extract and write per-focal-gene neighborhoods
#'
#' Iterate over each focal centroid_80 listed in `focal_c80_df`, locate its
#' neighborhood TSV under `get_target("neighbor_list")`, run the full per-focal
#' Step 1 pipeline via [parse_gene_neighbor()], and write the resulting tables
#' to per-genome TSVs via [write_gene_neighbor()]. Errors and missing files are
#' handled per focal so that one failure does not abort the whole batch.
#'
#' @details
#' **Per-focal error handling.** Missing input files emit a `warning()` and
#' the focal is skipped. Errors inside [parse_gene_neighbor()] or
#' [write_gene_neighbor()] are caught per focal via `tryCatch` and logged with
#' `message()`, so one bad focal does not abort processing of the rest. The
#' caller therefore should not wrap this function in its own `tryCatch`.
#'
#' **Configuration dependency.** Reads `get_target("neighbor_list")` for
#' inputs and `get_target("neighbor_groups_by_genome")` for outputs, so
#' `job_config` must already be loaded via [load_job_config()] before calling.
#'
#' @export
extract_and_write_per_focal_neighbors <- function(focal_c80_df, gene_to_c80) {
  focal_c80_df %>%
    group_by(focal_c80, focal_label) %>%
    group_walk(~{
      focal_c80 <- .y$focal_c80[[1]]
      focal_label <- .y$focal_label[[1]]
      message("c80=", focal_c80, " | label=", focal_label)

      in_fp <- file.path(get_target("neighbor_list"), paste0(focal_c80, ".tsv"))
      if (!file.exists(in_fp)) {
        warning("Missing neighbors file for c80=", focal_c80, ": ", in_fp)
        return(invisible(NULL))
      }

      tryCatch({
        gene_neighbors <- parse_gene_neighbor(in_fp, focal_c80, gene_to_c80)
        if (!is.null(gene_neighbors)) {
          write_gene_neighbor(gene_neighbors, focal_c80, focal_label)
        }
      }, error = function(e) {
        message("Error in group c80=", focal_c80, " | label=", focal_label,
                ":\n   ", conditionMessage(e))
      })
    }, .keep = TRUE)
  invisible(NULL)
}


#' Load gene-neighbor tables across genomes
#'
#' Recursively load all `.tsv` gene-neighbor files from an input directory and
#' combine them into a single table. Metadata are derived from the file path:
#' the file basename defines `focal_c80`, and the grandparent directory defines
#' `path_label`.
#'
#' @import data.table
#' @import dplyr
#' @importFrom parallel mclapply detectCores
#' @importFrom tools file_path_sans_ext
#' @export
load_neighbors_across_genomes <- function(input_dir, mc_cores = NULL) {
  # each neighbor_genome can have multiple gene_members corresponding to the same focal_c80
  
  all_files <- list.files(input_dir, pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE)
  
  col_classes <- c(
    gene_member_genome  = "character",
    gene_member         = "character",
    neighbor_gene_id    = "character",
    neighbor_contig_id  = "character",
    neighbor_gene_start = "integer",
    neighbor_gene_end   = "integer",
    neighbor_gene_strand = "character",
    neighbor_gene_type   = "character",
    neighbor_c80_coarse  = "character",
    neighbor_gene_length = "integer",
    neighbor_c80_length_coarse = "integer",
    neighbor_c80_genome_prevalence = "numeric",
    relative_position    = "integer",
    total_genomes        = "integer"
  )
  
  if (is.null(mc_cores)) {
    mc_cores <- max(1L, min(parallel::detectCores(), length(all_files)))
  } else {
    if (!is.numeric(mc_cores) || length(mc_cores) != 1L || is.na(mc_cores) || mc_cores < 1) {
      stop("`mc_cores` must be a single positive number.", call. = FALSE)
    }
    mc_cores <- as.integer(mc_cores)
  }
  
  read_one_file <- function(fp) {
    dt <- tryCatch(
      data.table::fread(
        fp,
        sep = "\t",
        colClasses = col_classes,
        showProgress = FALSE
      ),
      error = function(e) NULL
    )
    
    if (is.null(dt) || nrow(dt) == 0L) {
      return(NULL)
    }
    
    dt[, `:=`(
      focal_c80 = tools::file_path_sans_ext(basename(fp)),
      path_label = basename(dirname(dirname(fp)))
    )]
    
    dt
  }
  
  parts <- parallel::mclapply(
    X = all_files,
    mc.cores = mc_cores,
    FUN = read_one_file)
  
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) == 0L) {
    return(tibble::tibble())
  }
  
  gene_neighbors <- data.table::rbindlist(parts, use.names = TRUE, fill = TRUE) %>%
    relocate(focal_c80, path_label, .before = 1)

  return(gene_neighbors)
}


#' Run Step 1 — per-focal neighborhood extraction + label attachment
#'
#' Orchestrator for Step 1. Three sub-stages, all gated on the existence
#' of the cached `neighbor_groups_rds`:
#'
#' \enumerate{
#'   \item \strong{Per-focal extraction.}
#'     [extract_and_write_per_focal_neighbors()] runs
#'     [parse_gene_neighbor()] over every `is_focal == TRUE` centroid in
#'     `focal_c80_df` (orient + group + filter), sharding surviving
#'     neighborhoods to per-focal-per-genome TSVs.
#'   \item \strong{Cross-genome assembly.}
#'     [load_neighbors_across_genomes()] concatenates every shard into
#'     one long `gene_neighbors` table.
#'   \item \strong{Label attachment.} [assign_c80_to_short_genes()]
#'     synthesises focal-scoped labels for unannotated short ORFs and
#'     emits the `short_gene_prevalence` lookup; [compute_c80_variants()]
#'     adds length-variant suffixes to multi-length clusters.
#' }
#'
#' Idempotent: if `neighbor_groups_rds` already exists on disk, the
#' build/save block is skipped and all three artefacts are read back
#' from their respective caches. To force a re-run, delete
#' `get_target("neighbor_groups_rds")`.
#'
#' @param focal_c80_df Focal centroid table with an `is_focal` boolean
#'   column. Only `is_focal == TRUE` rows drive Step 1.
#' @param gene_to_c80 Per-gene-id to coarse-cluster lookup
#'   (from [load_c80_tables()] in midas.R).
#'
#' @return A list with three named entries:
#'   * `gene_neighbors` — unified neighbor table (per focal, per
#'     genome, per neighbor position) with both `neighbor_c80_coarse`
#'     and `neighbor_c80_fine` populated.
#'   * `short_gene_prevalence` — per-synthetic-c80 prevalence map for
#'     unannotated short ORFs.
#'   * `c80_variants_mapping` — `(neighbor_c80_coarse,
#'     neighbor_gene_length) -> neighbor_c80_fine` mapping.
#'
#' @export
run_step1_neighbor_extraction <- function(focal_c80_df, gene_to_c80) {
  gene_neighbors_rds <- get_target("neighbor_groups_rds")
  if (file.exists(gene_neighbors_rds)) {
    return(list(
      gene_neighbors = readRDS(get_target("neighbor_groups_rds")),
      short_gene_prevalence = readRDS(get_target("short_gene_prevalence")),
      c80_variants_mapping = readRDS(get_target("c80_variants_mapping"))
    ))
  }

  extract_and_write_per_focal_neighbors(focal_c80_df %>% filter(is_focal), gene_to_c80)
  gene_neighbors <- load_neighbors_across_genomes(get_target("neighbor_groups_by_genome"))

  # decorate with (1) synthetic smallORFs (2) length-variants
  res <- assign_c80_to_short_genes(gene_neighbors)
  gene_neighbors <- res$gene_neighbors
  short_gene_prevalence <- res$short_gene_prevalence
  c80_variants_mapping <- compute_c80_variants(gene_neighbors)

  gene_neighbors <- gene_neighbors %>%
    left_join(c80_variants_mapping, by = c("neighbor_c80_coarse", "neighbor_gene_length")) %>%
    select(focal_c80:neighbor_c80_coarse, neighbor_c80_length_coarse, neighbor_c80_fine, neighbor_gene_length, everything())

  saveRDS(gene_neighbors, get_target("neighbor_groups_rds"))
  saveRDS(short_gene_prevalence, get_target("short_gene_prevalence"))
  saveRDS(c80_variants_mapping, get_target("c80_variants_mapping"))

  list(
    gene_neighbors = gene_neighbors,
    short_gene_prevalence = short_gene_prevalence,
    c80_variants_mapping = c80_variants_mapping
  )
}
