# ------------------------------------------------------------------------------
# plot.R
#
# Plotting functions used at two stages of the pipeline:
#   * Step 1 diagnostic plotters (`plot_neighbor_groups`,
#     `find_most_prevalent_operon_size`, `compute_selected_bins`,
#     `find_peaks`, `get_fill_scale`): called from `parse_gene_neighbor`
#     in neighbor.R; render the per-focal fig1-fig5 PDFs.
#   * Step 6 operon-visualization plotters (`plot_coarse_operons`,
#     `plot_fine_operons`, `plot_*_by_component`, `run_step6_figures`)
#     Called from pipeline.R; render the global + per-component
#     gggenes PDFs from the Step 3 canonical-paths outputs.
#
# Author:  Chunyu Zhao <chunyu.zhao@gladstone.ucsf.edu>
# Created: 2025-07-14
# Updated: 2026-04-30
# ------------------------------------------------------------------------------

library(tidyverse)
library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)


find_peaks <- function(y) {
  which(diff(sign(diff(y))) == -2) + 1  # local maxima indices
}


compute_selected_bins <- function(centroid_counts) {
  d <- density(centroid_counts$operon_size, na.rm = TRUE)
  peak_indices <- find_peaks(d$y)
  peak_x <- d$x[peak_indices]
  peak_y <- d$y[peak_indices]
  # Order by decreasing peak height
  ord <- order(peak_y, decreasing = TRUE)
  peak_x_sorted <- peak_x[ord][1]
  selected_bins <- sort(unique(c(floor(peak_x_sorted), ceiling(peak_x_sorted))))
  return(selected_bins)
}


find_most_prevalent_operon_size <- function(operon_size_tbl, fig_fp) {
  n_unique <- length(unique(operon_size_tbl$operon_size))
  if (n_unique > 1) {
    f <- ggplot(operon_size_tbl, aes(x = operon_size)) +
      geom_density(fill = "purple", alpha = 0.3) +
      theme_minimal() +
      labs(title = paste("Density of Non-Missing centroid_80 Counts"),
           x = "# Neighbors with centroid_80",
           y = "Density")
    ggsave(fig_fp, f, height = 4, width = 6)
    
    # find the most represented operon
    selected_bins <- compute_selected_bins(operon_size_tbl)
    operon_size_tbl_selected <- operon_size_tbl %>%
      filter(operon_size %in% selected_bins) %>%
      arrange(operon_size)
    return(operon_size_tbl_selected)
  }
  return(operon_size_tbl)
}


load_cog_info <- function(operon_annoted) {
  
  coginfo <- operon_annoted %>%
    ungroup()%>% 
    select(centroid_80,COG_category) %>%
    unique() %>%
    separate_wider_position(COG_category,c(COG1 = 1, COG2 = 1, COG3 = 1),too_few = "align_start") %>% 
    mutate(is1 = 1) %>% 
    pivot_longer(cols = starts_with("COG")) %>% filter(!is.na(value))# %>%  
  
  cogmeta<-fread("data/cog.tsv")
  
  coginfo <- coginfo %>% left_join(cogmeta,by=c("value"="COGCAT")) %>% 
    mutate(COGDes= if_else(is.na(COGDes),"-",COGDes), COGdes2= if_else(is.na(COGdes2), "-", COGdes2))
  
  cog_levels <- names(sort(tapply(coginfo$is1, coginfo$COGDes, sum)))
  
  coginfo$cogfactor <- factor(coginfo$value, levels=cog_levels)
  coginfo <- coginfo %>% filter(!is.na(value))
  coginfo$cogfactor <- factor(coginfo$COGDes, levels=cog_levels)
  
  coginfo <- coginfo %>% select(-is1)
  
  return(coginfo)
}


add_COGdes2 <- function(operon_annoted, coginfo) {
  collapsed_tbl <- coginfo %>%
    select(centroid_80, COGdes2) %>%
    group_by(centroid_80) %>%
    summarise(COGdes2_collapsed = paste(unique(COGdes2), collapse = "\n")) %>%
    ungroup()
  
  operon_annoted <- operon_annoted %>% left_join(collapsed_tbl, by=c("centroid_80"))
  return(operon_annoted)
}



get_fill_scale <- function(df, fill_by) {
  if (grepl("prevalence", fill_by)) {
    # continuous 0–1 palette
    return(scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(0, 1)))
  }
  
  if (fill_by == "gene_member") {
    gene_ids <- sort(unique(df$gene_member), na.last = TRUE)
    n <- length(gene_ids)
    
    # at least 1 color; reserve lightgray for last level if desired
    base_cols <- colorRampPalette(brewer.pal(12, "Paired"))(max(n - 1, 1))
    vals <- if (n >= 2) c(base_cols, "lightgray") else base_cols
    
    color_mapping <- setNames(vals, gene_ids)
    return(scale_fill_manual(values = color_mapping, na.value = "grey90"))
  }
  
  # default discrete palette for other categorical fields
  scale_fill_viridis_d(na.value = "grey90")
}


plot_neighbors <- function(gene_neighbors, operon_size_tbl, focal_c80, fill_by = "neighbor_c80_coarse", fig_fp = "fig1.pdf", height_per_row = 0.03) {
  operon_size_tbl <- operon_size_tbl %>%
    mutate(gene_member = factor(gene_member, levels = .$gene_member))


  n_categories <- length(unique(gene_neighbors$gene_member))
  gg_height <- max(3, n_categories * height_per_row)
  
  fill_scale <- if (grepl("prevalence", fill_by)) {
    scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(0, 1))
  } else {
    scale_fill_viridis_d(na.value = "grey90")
  }
  
  f <- gene_neighbors %>%
    ggplot(aes(x = relative_position, y = gene_member, fill = .data[[fill_by]] )) +
    geom_tile(color = "gray30") +
    fill_scale +
    scale_y_discrete(limits = rev(levels(operon_size_tbl$gene_member))) +  # optional: top = high count
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),  # inner panel
      plot.background = element_rect(fill = "white", color = NA),   # outer plot area
      panel.grid.major = element_blank(),  # remove major grid lines
      panel.grid.minor = element_blank(),  # remove minor grid lines
      axis.line = element_line(color = "black"),  # keep axis lines
      axis.ticks.x = element_line(color = "black")  # keep ticks
    ) +
    labs(
      title = paste("Genomic Neighborhood (",focal_c80, ")"),
      x = "Relative Position",
      y = "Focal Gene (Sorted by # Non-NA centroid_80)"
    ) +
    theme(axis.text.y = element_blank()) + 
    guides(fill = guide_legend(ncol = 1)) + 
    ggtitle(focal_c80)
  ggsave(fig_fp, f, height = gg_height, width = 10, limitsize = FALSE)
}


plot_neighbor_groups <- function(gene_pats, focal_c80, fill_by = "gene_member", fig_fp = "fig1.pdf", height_per_row = 0.25) {
  df <- gene_pats %>% 
    mutate(group_label = fct_reorder(group_label, total_genomes, .desc = FALSE))
  
  offset_df <- df %>%
    filter(position == 0) %>%
    select(group_label, zero_gene_start = gene_start) 
  
  df <- df %>%
    left_join(offset_df, by = "group_label") %>%
    mutate(gene_start = gene_start - zero_gene_start, gene_end = gene_end - zero_gene_start)
  df <- df %>% filter(neighbor_gene_length > 0)
  
  n_rows <- length(unique(df$group_label))
  gg_height <- n_rows * height_per_row
  gg_height <- max(3, gg_height)
  
  f <- ggplot(df, aes(xmin = gene_start, xmax = gene_end, y = group_label, fill = .data[[fill_by]] )) +
    geom_rect(aes(ymin = as.numeric(factor(group_label)) - 0.4,  ymax = as.numeric(factor(group_label)) + 0.4)) +
    geom_text(
      aes(
        x = (gene_start + gene_end) / 2,
        y = group_label,
        label = position
      ),
      size = 2.5,
      color = "black"
    ) +
    get_fill_scale(df, fill_by) + 
    theme_minimal() +
    labs(x = "Genomic Position (pseudo-scale)", y = "Neighborhood group (with genome counts)") + 
    guides(fill = guide_legend(ncol = 1)) + 
    ggtitle(focal_c80)
  ggsave(fig_fp, f, height = gg_height, width = 14, limitsize = FALSE)
}


saturated_rainbow <- function (n, saturation_limit=0.4) {
  #https://github.com/kylebittinger/qiimer/blob/master/R/otu_table.R
  saturated_len <- floor(n * (1 - saturation_limit))
  rainbow_colors <- rev(rainbow(n - saturated_len, start=0, end=0.6))
  last_color <- tail(rainbow_colors, n=1)
  saturated_colors <- rep(last_color, saturated_len)
  colors <- c(rainbow_colors, saturated_colors)
  colors[1] <- "#FFFFFFFF"
  colors
}

# ------------------------------------------------------------------------------
# Operon visualization (gggenes): Step 6
#
# Multi-fill global plotters and per-component plotters that consume the
# canonical operon tables emitted by Step 3 (canonical_paths_c80s,
# canonical_paths_fine_c80s) plus the Step 5 selection sets. Wired
# together by `run_step6_figures` at the bottom of this file. Data-prep
# helpers (assign_c80_label, decorate_with_updated_path_type) live in
# parse.R.
# ------------------------------------------------------------------------------

# Attach position / layout columns per operon for gggenes. `is_truncated`
# and `is_fragmented` may be all-NA (coarse) or logical (fine); the
# case_when handles both. Glyph precedence (top of case_when wins): U/D
# (focal direction by beta sign) > F (fragmented) > T (truncated) > s
# (small ORF). Focal direction wins on focal rows so the trait-
# association signal isn't masked by structural tags. F beats T because
# fragmentation is the rarer, more specific signal: a row is fragmented
# only when the same coarse cluster appears at >=2 length variants
# within one isoform, and the short copy of such a pair is almost
# always also truncated -- if T won, F would never fire on the plot.
# Plain truncations still show T.
#
# `annotate_smallORF_placeholder = TRUE` appends "**" to fill_symbol on
# rows where `neighbor_c80_coarse` has 3 underscores (the 4-token
# synthetic small-ORF placeholder shape from `assign_c80_to_short_genes`);
# kept off by default since it duplicates the "s" glyph for the typical
# input shape.
.layout_operon_tracks <- function(df_annot, group_key, label_col, annotate_smallORF_placeholder = FALSE) {
  if (!"is_truncated"  %in% names(df_annot)) df_annot$is_truncated  <- NA
  if (!"is_fragmented" %in% names(df_annot)) df_annot$is_fragmented <- NA

  gene_padding_bp <- cfg_get(job_config, "gene_padding_bp")

  out <- df_annot %>%
    group_by(across(all_of(group_key))) %>%
    mutate(
      order   = row_number(),
      padded  = neighbor_gene_length + gene_padding_bp,
      start   = lag(cumsum(padded), default = 0),
      end     = start + neighbor_gene_length,
      forward = TRUE
    ) %>%
    ungroup() %>%
    mutate(
      track_label = .data[[label_col]],
      fill_symbol = case_when(
        coalesce(is_focal, FALSE) & !is.na(beta) & beta < 0 ~ "D",
        coalesce(is_focal, FALSE) & !is.na(beta) & beta > 0 ~ "U",
        !is.na(is_fragmented) & is_fragmented               ~ "F",
        !is.na(is_truncated)  & is_truncated                ~ "T",
        is_smallORF                                         ~ "s",
        TRUE                                                ~ NA_character_
      )
    )

  if (annotate_smallORF_placeholder) {
    out <- out %>%
      mutate(
        .placeholder = !is.na(neighbor_c80_coarse) & str_count(neighbor_c80_coarse, "_") == 3,
        fill_symbol  = case_when(
          .placeholder & !is.na(fill_symbol) ~ paste0(fill_symbol, "**"),
          .placeholder                       ~ "**",
          TRUE                               ~ fill_symbol
        )
      ) %>%
      select(-.placeholder)
  }

  out
}


# Internal: build a `fill_gene` categorical column on a laid-out frame.
# Colors U/D rows (focal direction) and non-focal anchor rows; all others NA.
# Suffix encodes the row class so multiple instances of the same c80 don't
# collapse to one factor level when half the rows are anchor and half are U/D.
.add_fill_gene_col <- function(pdf) {
  pdf %>%
    mutate(fill_gene = case_when(
      fill_symbol %in% c("U", "D")                                       ~ paste0(neighbor_c80_coarse, "|", fill_symbol),
      path_type == "anchor" & !coalesce(is_focal, FALSE) & !is_smallORF  ~ paste0(neighbor_c80_coarse, "|A"),
      TRUE                                                               ~ NA_character_
    ))
}


# Internal: pick a fill scale for a fill mode. `pdf` only used for fill_gene
# (it determines the palette size).
.fill_scale_for <- function(fill_by, pdf) {
  switch(
    fill_by,
    "beta" = scale_fill_viridis_c(option = "D", na.value = "white", name = "beta"),
    "sample_prevalence" = scale_fill_viridis_c(option = "D", limits = c(0, 1), na.value = "white", name = "Sample prevalence"),
    "cor_to_b"  = scale_fill_viridis_c(option = "D", limits = c(-0.5, 0.5), na.value = "white", name = "cor to b"),
    "fill_gene" = {
      gene_ids <- sort(unique(pdf$fill_gene))
      n <- length(gene_ids)
      if (n == 0) {
        scale_fill_manual(values = character(0), na.value = "white", name = "Gene")
      } else {
        cols <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n)
        scale_fill_manual(values = setNames(cols, gene_ids), na.value = "white", name = "Gene")
      }
    },
    stop("Unknown fill_by: ", fill_by)
  )
}


# Internal: render one gggenes plot from a laid-out frame and save to disk.
# `pdf` must already carry start, end, forward, fill_symbol, and the column
# named by `fill_by` (and, for facet_var, that column too). `y_var` defaults
# to "track_label".
.draw_operons <- function(pdf, fill_by, out_fp, title, caption = NULL,
                          width = 20, height_per_row = 0.3,
                          facet_var = NULL, y_var = "track_label", y_label = "Operon") {
  n_rows <- length(unique(pdf[[y_var]]))
  gg_h   <- max(3, n_rows * height_per_row)

  p <- ggplot(pdf, aes(xmin = start, xmax = end, y = .data[[y_var]], forward = forward)) +
    gggenes::geom_gene_arrow(aes(fill = .data[[fill_by]]), 
                             arrowhead_height = unit(4, "mm"),
                             arrowhead_width  = unit(1, "mm")) +
    gggenes::geom_gene_label(aes(label = fill_symbol), size = 8, fontface = "bold") +
    .fill_scale_for(fill_by, pdf) +
    gggenes::theme_genes() +
    labs(x = "Position (bp, padded)", y = y_label, title = title, caption = caption)

  if (!is.null(facet_var)) {
    p <- p + facet_wrap(stats::as.formula(paste0("~", facet_var)), scales = "free_y")
  }

  ggsave(out_fp, p, width = width, height = gg_h, limitsize = FALSE)
  invisible(p)
}


.fill_modes_validate <- function(fill_by) {
  ok <- c("beta", "sample_prevalence", "cor_to_b", "fill_gene")
  bad <- setdiff(fill_by, ok)
  if (length(bad)) stop("Unknown fill_by mode(s): ", paste(bad, collapse = ", "),
                        ". Allowed: ", paste(ok, collapse = ", "))
  fill_by
}


.prepare_operon_pdf <- function(ann, group_key, label_col, fill_by) {
  pdf <- .layout_operon_tracks(ann, group_key = group_key, label_col = label_col)

  if ("fill_gene" %in% fill_by) pdf <- .add_fill_gene_col(pdf)
  if ("sample_prevalence" %in% fill_by && "genome_prevalence" %in% names(pdf)) {
    pdf <- pdf %>% mutate(sample_prevalence = coalesce(sample_prevalence, genome_prevalence))
  }
  pdf
}


#' Plot coarse recurring operons
#'
#' Expects `canonical_paths_c80s` to have already been decorated by
#' [decorate_c80s_w_smallORFs()] (adds `is_smallORF`, used for the "s"
#' glyph). `coarse_summary` is treated as a (typically pre-filtered)
#' selection set: only its `uid` column is read, used to `semi_join`
#' the c80s frame down to the operons to plot. Filter upstream: e.g.
#' `coarse_summary %>% semi_join(selected_fine, by = "uid")`: to keep
#' the PDF a manageable size.
#'
#' `fill_by` accepts any subset of `c("beta", "sample_prevalence",
#' "cor_to_b", "fill_gene")`; one PDF is emitted per mode, named
#' `coarse_operons_<fill_by>.pdf`. Default `c("beta")` matches the
#' previous single-output behavior.
#'
#' @export
plot_coarse_operons <- function(coarse_summary, canonical_paths_c80s, 
                                out_dir, fill_by, width = 20, height_per_row = 0.3) {
  if (!nrow(coarse_summary)) return(invisible(NULL))
  stopifnot("is_smallORF" %in% names(canonical_paths_c80s))
  fill_by <- .fill_modes_validate(fill_by)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  ann <- canonical_paths_c80s %>%
    semi_join(coarse_summary %>% select(uid), by = "uid") %>%
    mutate(track_label = uid)

  pdf <- .prepare_operon_pdf(ann, group_key = "uid", label_col = "track_label", fill_by = fill_by)

  for (fb in fill_by) {
    .draw_operons(pdf, fill_by = fb,
                  out_fp = file.path(out_dir, paste0("coarse_operons_", fb, ".pdf")),
                  title = "Coarse recurring operons",
                  caption = "U = focal up (beta > 0), D = focal down (beta < 0), s = small ORF",
                  width = width, height_per_row = height_per_row,
                  y_label = "Canonical operon")
  }
  invisible(NULL)
}


#' Plot fine-isoform operons
#'
#' Expects `canonical_paths_fine_c80s` to have already been decorated by
#' both [decorate_c80s_w_smallORFs()] and [decorate_c80s_w_truncation()]
#' (adds `is_smallORF` for the "s" glyph, `is_truncated` for "T", and
#' `is_fragmented` for "F"). `fine_summary` is treated as a (typically
#' pre-filtered) selection set: only its `uid_fine` column is read,
#' used to `semi_join` the c80s frame down to the isoforms to plot.
#' Filter upstream: e.g. `fine_summary %>% filter(n_fine_genomes >= ...)`
#': to keep the PDF a manageable size.
#'
#' `fill_by` accepts any subset of `c("beta", "sample_prevalence",
#' "cor_to_b", "fill_gene")`; one PDF is emitted per mode, named
#' `fine_operons_<fill_by>.pdf`. Default `c("beta")` matches the
#' previous single-output behavior.
#'
#' @section Caveat: why decoration reads `neighbor_c80_coarse`, not `neighbor_c80_fine`:
#' [compute_c80_variants()] appends a `_<length_rank>` suffix to
#' `neighbor_c80_fine` when a c80 cluster has multiple observed lengths.
#' Decoding from that column would yield an incorrect `smallORF_type`:
#' the `_\d+$` stripper removes only one trailing `_<digits>` chunk —
#' the length-rank suffix — leaving the inner gene-type rank attached.
#' Example: `_GUT_001-CDS_1_2` would decode to `smallORF_type = "CDS_1"`
#' instead of `"CDS"`. `decorate_c80s_w_smallORFs` hardcodes
#' `neighbor_c80_coarse` as the source to avoid this.
#'
#' @export
plot_fine_operons <- function(fine_summary, canonical_paths_fine_c80s, 
                              out_dir, fill_by, width = 20, height_per_row = 0.3) {
  if (!nrow(fine_summary)) return(invisible(NULL))
  stopifnot(all(c("is_smallORF", "is_truncated") %in%
                  names(canonical_paths_fine_c80s)))
  fill_by <- .fill_modes_validate(fill_by)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  ann <- canonical_paths_fine_c80s %>%
    semi_join(fine_summary %>% select(uid_fine), by = "uid_fine") %>%
    mutate(track_label = uid_fine)

  pdf <- .prepare_operon_pdf(ann, group_key = "uid_fine", label_col = "track_label", fill_by = fill_by)

  # Order tracks so isoforms of the same coarse uid are adjacent.
  pdf <- pdf %>%
    arrange(uid, isoform_rank) %>%
    mutate(track_label = factor(track_label, levels = unique(track_label)))

  for (fb in fill_by) {
    .draw_operons(pdf, fill_by = fb,
                  out_fp = file.path(out_dir, paste0("fine_operons_", fb, ".pdf")),
                  title = "Fine-isoform operons",
                  caption = "U = focal up (beta > 0), D = focal down (beta < 0), T = truncated, F = fragmented, s = small ORF",
                  width = width, height_per_row = height_per_row,
                  y_label = "Fine isoform")
  }
  invisible(NULL)
}


#' Plot coarse operons grouped by joint component
#'
#' One PDF per `(joint_component_id, fill_by)` combination, written to
#' `<out_dir>/02_by_component_coarse/comp_<id>_<fill_by>.pdf`. Within each
#' component, tracks are ordered by LCS-similarity hclust over their
#' c80 sequence (via [`generate_path_order`] from graph.R), and the
#' panel is faceted by `updated_path_type` (anchor_pos / anchor_neg /
#' the original path_type) with `scales = "free_y"`.
#'
#' Same scoping contract as [`plot_coarse_operons`]: pass a pre-filtered
#' `coarse_summary` (e.g. `selected_coarse`) and the function will
#' `semi_join` the c80s frame down by `uid`.
#'
#' @export
plot_coarse_operons_by_component <- function(coarse_summary, canonical_paths_c80s, 
                                             out_dir, fill_by,  
                                             width = 24, height_per_row = 0.3, min_paths = 1L) {
  if (!nrow(coarse_summary)) return(invisible(NULL))
  stopifnot("is_smallORF" %in% names(canonical_paths_c80s))
  fill_by <- .fill_modes_validate(fill_by)

  comp_dir <- file.path(out_dir, "02_by_component_coarse")
  if (!dir.exists(comp_dir)) dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)

  ann <- canonical_paths_c80s %>%
    semi_join(coarse_summary %>% select(uid), by = "uid") %>%
    mutate(track_label = uid)

  pdf <- .prepare_operon_pdf(ann, group_key = "uid",
                              label_col = "track_label", fill_by = fill_by)
  pdf <- decorate_with_updated_path_type(pdf)

  comp_ids <- pdf %>%
    distinct(joint_component_id, canonical_path_id, path_type) %>%
    count(joint_component_id, name = "n_paths") %>%
    filter(n_paths >= min_paths) %>%
    pull(joint_component_id)

  for (cid in comp_ids) {
    sub <- pdf %>% filter(joint_component_id == cid)
    ordered <- sub %>%
      transmute(path_label = track_label, gene = neighbor_c80_coarse, order = order) %>%
      generate_path_order()
    sub <- sub %>%
      mutate(track_label = factor(track_label, levels = ordered))

    for (fb in fill_by) {
      .draw_operons(sub, fill_by = fb,
                    out_fp = file.path(comp_dir, paste0("comp_", cid, "_", fb, ".pdf")),
                    title = paste("Component", cid),
                    caption = "U = focal up (beta > 0), D = focal down (beta < 0), s = small ORF",
                    width = width, height_per_row = height_per_row,
                    facet_var = "updated_path_type",
                    y_label = "Canonical operon")
    }
  }
  invisible(NULL)
}


#' Plot fine isoforms grouped by joint component
#'
#' Fine analog of [`plot_coarse_operons_by_component`]; one PDF per
#' `(joint_component_id, fill_by)` written to
#' `<out_dir>/03_by_component_fine/comp_<id>_<fill_by>.pdf`. Faceted by
#' `updated_path_type` with `scales = "free_y"`. Within each component,
#' tracks are ordered by LCS-similarity hclust.
#'
#' @export
plot_fine_operons_by_component <- function(fine_summary, canonical_paths_fine_c80s, 
                                           out_dir, fill_by, 
                                           width = 24, height_per_row = 0.3, min_paths = 1L) {
  if (!nrow(fine_summary)) return(invisible(NULL))
  stopifnot(all(c("is_smallORF", "is_truncated") %in%
                  names(canonical_paths_fine_c80s)))
  fill_by <- .fill_modes_validate(fill_by)

  comp_dir <- file.path(out_dir, "03_by_component_fine")
  if (!dir.exists(comp_dir)) dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)

  ann <- canonical_paths_fine_c80s %>%
    semi_join(fine_summary %>% select(uid_fine), by = "uid_fine") %>%
    mutate(track_label = uid_fine)

  pdf <- .prepare_operon_pdf(ann, group_key = "uid_fine",
                              label_col = "track_label", fill_by = fill_by)
  pdf <- decorate_with_updated_path_type(pdf)

  comp_ids <- pdf %>%
    distinct(joint_component_id, canonical_path_id, path_type) %>%
    count(joint_component_id, name = "n_paths") %>%
    filter(n_paths >= min_paths) %>%
    pull(joint_component_id)

  for (cid in comp_ids) {
    sub <- pdf %>% filter(joint_component_id == cid)
    ordered <- sub %>%
      transmute(path_label = track_label, gene = neighbor_c80_coarse, order = order) %>%
      generate_path_order()
    sub <- sub %>%
      mutate(track_label = factor(track_label, levels = ordered))

    for (fb in fill_by) {
      .draw_operons(sub, fill_by = fb,
                    out_fp = file.path(comp_dir, paste0("comp_", cid, "_", fb, ".pdf")),
                    title = paste("Component", cid),
                    caption = "U = focal up (beta > 0), D = focal down (beta < 0), T = truncated, F = fragmented, s = small ORF",
                    width = width, height_per_row = height_per_row,
                    facet_var = "updated_path_type",
                    y_label = "Fine isoform")
    }
  }
  invisible(NULL)
}


#' Run Step 6 — render all gggenes figures
#'
#' Orchestrator for Step 6. Renders the four `plot_*` outputs at both
#' global and per-component scope:
#'
#' \itemize{
#'   \item [plot_coarse_operons()] — one PDF per `fill_by` mode at the global level.
#'   \item [plot_fine_operons()] — fine-isoform analogue.
#'   \item [plot_coarse_operons_by_component()] — one PDF per
#'     `(joint_component_id, fill_by)` under `step6_figures/02_by_component_coarse/`.
#'   \item [plot_fine_operons_by_component()] — fine analogue under
#'     `step6_figures/03_by_component_fine/`.
#' }
#'
#' Inputs are passed in explicitly so the caller (pipeline.R) is the
#' one place that loads the four TSVs from disk. Re-running Step 6 in
#' isolation (e.g., to tweak `parse.fill_modes` and re-render without
#' touching Steps 1-5) is still cheap — load the four TSVs from disk
#' first and call this function:
#'
#' \preformatted{
#' Rscript -e 'source("pipeline/config.R"); source("pipeline/model.R");
#'             source("pipeline/parse.R");  source("pipeline/plot.R");
#'             load_job_config("example.yaml");
#'             selected_coarse <- read_delim(get_target("parse_selected_coarse"), delim = "\\t", show_col_types = FALSE);
#'             selected_fine <- read_delim(get_target("parse_selected_fine"), delim = "\\t", show_col_types = FALSE);
#'             c80s_coarse <- read_delim(get_target("canonical_paths_c80s"), delim = "\\t", show_col_types = FALSE);
#'             c80s_fine <- read_delim(get_target("canonical_paths_fine_c80s"), delim = "\\t", show_col_types = FALSE);
#'             run_step6_figures(selected_coarse, selected_fine, c80s_coarse, c80s_fine)'
#' }
#'
#' Reads `fill_modes` from `job_config` via `cfg_get` and resolves the
#' two figure-output directories via `get_target` (both currently point
#' to `step6_figures/`).
#'
#' @param selected_coarse Coarse-summary selection set (read from
#'   `parse_selected_coarse`); only `uid` is consumed by the plotters as
#'   a `semi_join` filter.
#' @param selected_fine Fine-isoform selection set (read from
#'   `parse_selected_fine`); only `uid_fine` is consumed.
#' @param c80s_coarse Decorated L1 per-gene table (read from
#'   `canonical_paths_c80s`).
#' @param c80s_fine Decorated L2 per-isoform per-gene table (read from
#'   `canonical_paths_fine_c80s`).
#'
#' @return `invisible(NULL)`. Side effect: PDFs written under
#'   `step6_figures/`.
#'
#' @export
run_step6_figures <- function(selected_coarse, selected_fine, c80s_coarse, c80s_fine) {
  fill_modes <- cfg_get(job_config, "fill_modes")
  fig_coarse <- get_target("parse_coarse_figures")
  fig_fine <- get_target("parse_fine_figures")

  plot_coarse_operons(selected_coarse, c80s_coarse, fig_coarse, fill_modes)
  plot_fine_operons(selected_fine, c80s_fine, fig_fine, fill_modes)
  plot_coarse_operons_by_component(selected_coarse, c80s_coarse, fig_coarse, fill_modes)
  plot_fine_operons_by_component(selected_fine, c80s_fine, fig_fine, fill_modes)

  invisible(NULL)
}
