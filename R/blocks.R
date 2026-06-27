# ------------------------------------------------------------------------------
# blocks.R
#
# Step 6 of the strain-aware operon pipeline: focal block analysis (gated by
# blocks.skip_block; writes under step6_blocks/).
#
# Walks from per-gene trait statistics along canonical operon paths to a
# ranked set of focal block representatives with per-genome
# attribution. Four stages:
#   keep_focal_blocks()               - per-canonical-path focal gate + block clustering.
#   aggregate_blocks()                - focal-block extraction + cross-path aggregation.
#   rank_block_representatives()      - dominant-block selection + ranking.
#   map_representatives_to_genomes()  - per-genome attribution for reps.
#   run_step6_blocks()      - Step 6 orchestrator: runs all stages, writes outputs.
#
# Author:  Chunyu Zhao <chunyu.zhao@gladstone.ucsf.edu>
# Created: 2025-10-10 (extracted from pipeline.R; renumbered to Step 6 on 2026-05-15)
# ------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(purrr)


#' Identify and cluster focal blocks per canonical path
#'
#' Within each `(joint_component_id, canonical_path_id, path_type)` group,
#' walk the rows in their existing per-path order, label each as focal
#' (`is_focal == TRUE`) or non-focal, and cluster adjacent focals into
#' blocks. A new block starts whenever more than `allow_gaps` non-focal
#' rows intervene between two focals. Non-focal rows are dropped; focal
#' rows are returned tagged with a per-path `block_num`.
#'
#' @details
#' **Gap arithmetic.** `gaps[i] = hit_idx[i] - hit_idx[i-1]` counts
#' positions between consecutive focals (inclusive of the next focal).
#' The rule `gaps > allow_gaps + 1` therefore means "new block when more
#' than `allow_gaps` non-focal rows intervene."
#'
#' **Zero-focal paths.** Canonical paths with no focal rows return zero
#' rows (via `group_modify`) and are effectively dropped from the output.
#'
#' **NA handling.** Rows joined in from outside `focal_c80_df` (short
#' ORFs, neighbors absent from `focal_meta`) carry `is_focal = NA`. NA is
#' coerced to FALSE, so they correctly count as non-focal.
#'
#' @seealso [aggregate_blocks()] for the cross-canonical aggregation
#'   that consumes this output.
#' @export
keep_focal_blocks <- function(df, allow_gaps = 2) {
  # Drop synthetic `_`-prefixed ORF rows up front so block geometry (gap
  # distances) is computed over the same rows as before they inherited a
  # component. Real length-variant backbone genes are untouched.
  df <- df %>% filter(!startsWith(neighbor_c80_coarse, "_"))

  keys <- intersect(c("joint_component_id", "canonical_path_id", "path_type"), names(df))

  # Stage 1: tag each row focal-or-not, and save its global position so the
  # original row order can be restored after group_modify shuffles things.
  # NA in is_focal (e.g., short-ORF rows) is coerced to FALSE.
  df1 <- df %>%
    mutate(
      .row_id = row_number(),
      .is_focal = coalesce(is_focal, FALSE)
    )

  # Stage 2: within each canonical path, locate the focal rows and cluster
  # them into blocks. group_modify() can return any number of rows per
  # group (including zero), which is what lets us drop non-focal rows and
  # zero-focal paths in the same pass.
  out <- df1 %>%
    group_by(across(all_of(keys))) %>%
    group_modify(~{ # group_modify() can return any number of rows (including 0) per group
      g <- .x
      hit_idx <- which(g$.is_focal)  # row positions of focals within this path

      if (length(hit_idx) == 0) return(g[0, ]) # no hits in this path
      if (length(hit_idx) == 1) return(g[hit_idx, ] %>% mutate(block_num = 1L))
      
      # cluster hits into blocks: new block when gap > allow_gaps + 1
      gaps <- c(Inf, diff(hit_idx)) # differences between consecutive hit positions
      block_id <- cumsum(gaps > (allow_gaps + 1)) # at most 1 non-hit between them, a gap > 2 starts a new block
      g[hit_idx, ] %>% mutate(block_num = block_id) # you can drop block if not needed
      
    }) %>%
    ungroup()

  out %>% arrange(.row_id) %>% select(-.row_id, -.is_focal)
}


#' Extract focal blocks and aggregate across canonical paths
#'
#' Stage 1 + 2 of Step 6. Within each canonical path, find contiguous runs
#' of focal genes (`is_focal == TRUE`) via [keep_focal_blocks()], collapse
#' each run into an ordered c80 sequence, then aggregate across canonical
#' paths in the same `(joint_component_id, path_type)` group so equivalent
#' block patterns collapse into one row.
#'
#' @export
aggregate_blocks <- function(c80s_coarse, allow_gaps = 2) {
  # 1) Per-gene focal gate + block clustering within each canonical path
  c80_blocks <- keep_focal_blocks(c80s_coarse, allow_gaps = allow_gaps)

  # 2) Collapse each block into an ordered c80 sequence.
  #    `uid` is added to the group_by because it is functionally determined by
  #    (joint_component_id, path_type, canonical_path_id, n_genomes); carrying
  #    it through lets step 6 emit a parallel `canonical_uids` list alongside
  #    `canonical_paths`.
  edge_c80_per_block <- c80_blocks %>%
    group_by(joint_component_id, canonical_path_id, path_type, block_num, n_genomes, uid) %>%
    mutate(.idx = row_number()) %>%   # preserve within-block order
    summarise(seq_raw = list(neighbor_c80_coarse[order(.idx)]), .groups = "drop") %>%
    # Dedup consecutive duplicates before computing edges/paths
    mutate(
      seq_dedup        = map(seq_raw, dedup_consecutive_vec),
      n_genes          = lengths(seq_dedup),
      left_orig        = map_chr(seq_dedup, ~ .x[1]),
      right_orig       = map_chr(seq_dedup, ~ .x[length(.x)]),
      block_c80s_path  = map_chr(seq_dedup, ~ paste(.x, collapse = " → ")),
      is_singleton     = (n_genes == 1),
      edge_genes       = ifelse(is_singleton, left_orig, paste(left_orig, right_orig, sep = " → "))
    ) %>%
    select(-seq_raw, -seq_dedup)

  # 3) edge_pair identifier (left_orig;right_orig or singleton token)
  edge_table <- edge_c80_per_block %>%
    mutate(edge_pair = ifelse(is_singleton, left_orig, paste(left_orig, right_orig, sep = ";")))

  # 4) Aggregate at block-path level within each edge_pair; per-component frequencies
  block_agg <- edge_table %>%
    group_by(joint_component_id, path_type, edge_pair, block_c80s_path) %>%
    summarise(
      block_n_paths    = n(),
      block_n_genomes  = sum(n_genomes, na.rm = TRUE),
      block_n_genes    = max(n_genes, na.rm = TRUE),
      left_orig        = first(left_orig),
      right_orig       = first(right_orig),
      canonical_paths  = paste(unique(canonical_path_id), collapse = ";"),
      # Parallel list of canonical uids. Same first-occurrence order as
      # `canonical_paths` because `uid` is a deterministic function of
      # `canonical_path_id` within this grouping.
      canonical_uids   = paste(unique(uid), collapse = ";"),
      .groups = "drop_last"
    ) %>%
    group_by(joint_component_id, path_type) %>%
    mutate(
      block_total = sum(block_n_genomes),
      block_freq  = block_n_genomes / block_total
    ) %>%
    ungroup()

  block_agg
}


#' Select a reference block per component and rank everything else
#'
#' Stage 3 + 4 of Step 6. Within each `(joint_component_id, path_type)`
#' group, pick the dominant block as a reference (`selected_tbl`), then use
#' [annotate_group()] to rank every block relative to that reference based
#' on sequence-overlap relationships. Non-redundant blocks carry
#' `rep_rank >= 1`; redundant subsequences of existing reps get
#' `rep_rank = 0` and are dropped from the returned `representatives` table.
#'
#' @details
#' **Reference-selection tie-break**: within each `(component, path_type)`,
#' sort by `block_freq` desc, then `block_n_genomes` desc, then
#' `block_n_genes` desc, then `block_n_paths` desc, then `block_c80s_path`
#' alphabetically (for determinism). Take the top row.
#'
#' @export
rank_block_representatives <- function(block_agg, min_overlap = 1) {
  # 1) Per (component, path_type), pick the dominant block as reference
  selected_tbl <- block_agg %>%
    group_by(joint_component_id, path_type) %>%
    arrange(desc(block_freq), desc(block_n_genomes), desc(block_n_genes),
            desc(block_n_paths), block_c80s_path, .by_group = TRUE) %>%
    slice(1) %>%
    transmute(
      joint_component_id, path_type,
      selected_freq          = block_freq,
      selected_edge_pair     = edge_pair,
      selected_block_path    = block_c80s_path,
      selected_n_genes       = block_n_genes,
      selected_n_genomes     = block_n_genomes,
      selected_canonical_ids = canonical_paths
    ) %>%
    ungroup()

  # 2) Rank every block relative to its component's reference
  annotated <- block_agg %>%
    dplyr::group_by(joint_component_id, path_type) %>%
    dplyr::group_modify(~ annotate_group(.x, .y, selected_tbl = selected_tbl, min_overlap  = min_overlap)) %>%
    dplyr::ungroup()

  # 3) Keep only non-redundant representatives, ordered by rank. Bake the
  #    block-level uid here so it flows into downstream outputs without
  #    colliding with the canonical `uid` in canonical_paths / _c80s.
  representatives <- annotated %>%
    filter(rep_rank >= 1) %>%
    arrange(joint_component_id, path_type, rep_rank) %>%
    mutate(block_uid = paste0("cmp", joint_component_id, "-", path_type, "-rank", rep_rank, "-nge", block_n_genes)) %>%
    select(block_uid, everything())

  list(selected_tbl = selected_tbl, annotated = annotated, representatives = representatives)
}


#' Attach per-genome attribution to each representative block
#'
#' Stage 5 of Step 6. For each representative block, walk back through the
#' Step 3 provenance chain via
#' [explode_canonical_into_collapsed_paths()] to find the genomes that
#' carry its contributing canonical paths, and construct a stable block
#' identifier `uid` for downstream use.
#'
#' @details
#' **Sanity check**: `stopifnot(path_type == path_type_per_genome)` - the block's
#' `path_type` must agree with the per-genome `path_type` joined in from
#' `path_df` for every row, confirming the provenance chain stayed
#' type-consistent. The per-genome column is dropped before returning.
#'
#' Does not dedupe to `(uid, neighbor_genome, left_orig, right_orig)`; the
#' caller can apply that projection when writing the final per-genome
#' presence table.
#'
#' @export
map_representatives_to_genomes <- function(representatives, canonical_paths, collapsed_paths, path_df) {
  # 1) canonical → per-genome provenance lookup
  can_to_col <- explode_canonical_into_collapsed_paths(canonical_paths,
                                                       collapsed_paths,
                                                       path_df)

  # 2) Attach per-genome attribution to each representative block.
  #    `block_uid` is inherited from `representatives` (constructed in
  #    rank_block_representatives). `canonical_paths` and `canonical_uids` are
  #    parallel `;`-joined lists - `separate_rows` explodes them together so
  #    each exploded row pairs one canonical_path_id with its canonical uid.
  rep_slim <- representatives %>%
    select(block_uid, joint_component_id, path_type, rep_rank,
           block_n_paths, block_n_genomes, block_n_genes,
           left_orig, right_orig, canonical_paths, canonical_uids, block_freq,
           relation_to_selected) %>%
    separate_rows(canonical_paths, canonical_uids, sep = ";") %>%
    dplyr::rename(canonical_uid = canonical_uids) %>%
    left_join(
      can_to_col %>% dplyr::rename(path_type_per_genome = path_type),
      by = c("canonical_paths" = "canonical_path_id"))

  # 3) Sanity: block path_type must match per-genome path_type
  stopifnot(nrow(rep_slim %>% filter(path_type != path_type_per_genome)) == 0)

  rep_slim %>% select(-path_type_per_genome)
}


#' Diagnose whether surviving representatives contain truly-overlapping paths
#'
#' Within each `(joint_component_id, path_type)` group, count pairs of
#' representatives that share a contiguous substring of ≥ `min_shared`
#' tokens but neither contains the other. These are the cases
#' [annotate_group()]'s subset-only redundancy check cannot collapse -
#' they survive as separate reps and potentially double-count genomes in
#' `rep_slim`.
#'
#' Called as a pre-flight diagnostic from the pipeline; the summary
#' message indicates whether an overlap-annotation step (see
#' `flag_overlapping_reps`) is worth implementing on this dataset.
diagnose_rep_overlaps <- function(representatives, min_shared = 2) {
  per_group <- representatives %>%
    group_by(joint_component_id, path_type) %>%
    summarise(
      n_reps = n(),
      n_overlap_pairs = {
        p <- block_c80s_path
        if (length(p) < 2) {
          0L
        } else {
          pairs <- combn(length(p), 2)
          sum(vapply(seq_len(ncol(pairs)), function(k) {
            i <- pairs[1L, k]; j <- pairs[2L, k]
            lccs_len(p[i], p[j]) >= min_shared &&
              !is_contig_subseq(p[i], p[j]) &&
              !is_contig_subseq(p[j], p[i])
          }, logical(1)))
        }
      },
      .groups = "drop"
    )

  summary_tbl <- tibble(
    total_components    = nrow(per_group),
    comps_with_overlap  = sum(per_group$n_overlap_pairs > 0),
    total_overlap_pairs = sum(per_group$n_overlap_pairs),
    max_pairs_in_a_comp = if (nrow(per_group) == 0) 0L else max(per_group$n_overlap_pairs)
  )

  list(per_group = per_group, summary = summary_tbl)
}


#' Rank blocks within one `(component, path_type)` group by redundancy
#'
#' For the blocks of one `(joint_component_id, path_type)` group, look up
#' the pre-computed reference block in `selected_tbl`, tag each row with
#' its relation to the reference via [get_relation()], and greedily assign
#' representative ranks so contiguous-subset-redundancies are marked out.
#' Called once per group via `dplyr::group_modify` from
#' [rank_block_representatives()].
#'
#' @details
#' **Greedy rep construction, seeded by the selected reference.** Rows are
#' sorted by `(block_n_genes desc, block_n_genomes desc, block_n_paths
#' desc, block_c80s_path asc)` so longer blocks get processed before
#' shorter ones. The selected reference is installed as rank 1. Each
#' subsequent row is either (a) a contiguous subset of an already-
#' installed rep (marked redundant, `rep_rank = 0`) or (b) a new
#' representative with the next integer rank.
#'
#' **Redundancy is subset-only, forward-direction.** Uses
#' [is_contig_subseq()] which does not test the reverse of `p`. A block
#' and its exact mirror currently survive as two separate reps. A
#' reverse-aware variant would also test `is_contig_subseq(rev(p), r)`;
#' not yet wired in (see parked/ROADMAP.md R2).
#'
#' **Superpaths are not redundant.** A block that *contains* an existing
#' rep gets its own rank - it adds information (length-variant extension)
#' relative to the shorter rep.
annotate_group <- function(df_group, keys, selected_tbl, min_overlap = 1) {
  sel <- keys %>%
    dplyr::left_join(selected_tbl, by = c("joint_component_id", "path_type")) %>%
    dplyr::slice(1)
  sel_path <- sel$selected_block_path
  
  # 1) per-row relation vs selected
  df <- df_group %>%
    mutate(relation_to_selected = map_chr(block_c80s_path, ~ get_relation(.x, sel_path, min_overlap)))
  
  # 2) minimal non-redundant representatives within this group, seeded by selected
  df_ord <- df %>% arrange(desc(block_n_genes), desc(block_n_genomes), desc(block_n_paths), block_c80s_path)
  
  reps     <- character(0)
  rep_ids  <- character(nrow(df_ord))
  rep_rank <- integer(nrow(df_ord))
  
  # keep selected as rank 1
  reps <- c(reps, sel_path)
  for (i in seq_len(nrow(df_ord))) {
    p <- df_ord$block_c80s_path[i]
    if (p == sel_path) {
      rep_ids[i]  <- p
      rep_rank[i] <- 1L
      next
    }
    covered <- FALSE
    if (length(reps)) {
      for (r in reps) {
        if (is_contig_subseq(p, r)) {   # p is subset of any existing rep → redundant
          rep_ids[i]  <- r
          rep_rank[i] <- 0L
          covered     <- TRUE
          break
        }
      }
    }
    if (!covered) {
      reps <- c(reps, p)               # add new representative
      rep_ids[i] <- p
      rep_rank[i] <- length(reps)      # 2, 3, ...
    }
  }
  # map back to original row order
  df_ord$representative_path <- rep_ids
  df_ord$is_redundant <- (rep_rank == 0L)
  df_ord$rep_rank <- rep_rank
  
  df %>%
    left_join(
      df_ord %>% select(block_c80s_path, representative_path, is_redundant, rep_rank),
      by = "block_c80s_path"
    )
}


# drop only consecutive duplicates (A A B B C -> A B C)
dedup_consecutive_vec <- function(v) {
  if (length(v) <= 1) return(v)
  v[c(TRUE, v[-1] != v[-length(v)])]
}

# TRUE iff every token of short_path appears as a contiguous run inside long_path.
is_contig_subseq <- function(short_path, long_path) {
  a <- split_path_string(short_path); b <- split_path_string(long_path)
  la <- length(a); lb <- length(b)
  if (la > lb) return(FALSE)
  if (la == 0)  return(TRUE)
  for (i in seq_len(lb - la + 1L)) if (all(b[i:(i+la-1L)] == a)) return(TRUE)
  FALSE
}

# Length of the longest common contiguous token subsequence (substring on tokens).
lccs_len <- function(path1, path2) {
  a <- split_path_string(path1); b <- split_path_string(path2)
  if (!length(a) || !length(b)) return(0L)
  m <- matrix(0L, nrow = length(a) + 1L, ncol = length(b) + 1L)
  best <- 0L
  for (i in seq_along(a)) {
    ai <- a[i]
    for (j in seq_along(b)) {
      if (ai == b[j]) {
        m[i+1L, j+1L] <- m[i, j] + 1L
        if (m[i+1L, j+1L] > best) best <- m[i+1L, j+1L]
      }
    }
  }
  best
}

# Classify a block path's relation to the selected block:
#   "selected" | "subset" | "superpath" | "overlap" | "disjoint".
get_relation <- function(path, sel_path, min_overlap = 1) {
  if (path == sel_path)                        return("selected")
  if (is_contig_subseq(path, sel_path))        return("subset")
  if (is_contig_subseq(sel_path, path))        return("superpath")
  if (lccs_len(path, sel_path) >= min_overlap) return("overlap")
  "disjoint"
}


#' Run Step 6 - focal block extraction + representative ranking
#'
#' Orchestrator for Step 6. Mines `c80s_coarse` for runs of focal genes
#' ("hit blocks"), aggregates equivalent runs across canonical paths,
#' picks one dominant block per `(joint_component_id, path_type)` as
#' reference, ranks every other block by its containment / overlap
#' relation to that reference, drops subset-redundant blocks, and walks
#' the surviving reps back to per-genome attribution via the Step 3
#' provenance walker.
#'
#' Reads `allow_gaps`, `min_overlap`, `min_shared` from the global
#' `job_config` via `cfg_get`. Writes `rep_path_df` (the
#' `representative_path.tsv`) and `uid_path_df` (the slim per-genome
#' projection) via `get_target`. If at least two non-redundant reps
#' survive and the resulting block-by-genome matrix is at least 3x3, a
#' `rep_heatmap.pdf` is rendered next to the per-genome TSV. Always
#' prints a one-line `diagnose_rep_overlaps` summary to stderr.
#'
#' @param c80s_coarse L1 per-gene table from
#'   [run_step3_consolidation()]; must carry `is_focal`, `n_genomes`,
#'   `joint_component_id`, `canonical_path_id`, `path_type`, `uid`, and
#'   `neighbor_c80_coarse`.
#' @param c_paths L1 canonical-paths table (Step 3 return).
#' @param collapsed_paths Pre-canonical collapse (Step 3 return).
#' @param path_df Per-genome maximal paths (Step 2 return).
#'
#' @return The slim per-genome attribution table `rep_slim`
#'   (`block_uid` x `canonical_uid` x `neighbor_genome` x `left_orig` x
#'   `right_orig`). The other intermediate tables (`representatives`,
#'   `rep_overlap_diag`) are written to disk / stderr only.
#'
#' @export
run_step6_blocks <- function(c80s_coarse, c_paths,
                                       collapsed_paths, path_df) {

  allow_gaps <- cfg_get(job_config, "allow_gaps")
  min_overlap <- cfg_get(job_config, "min_overlap")
  min_shared <- cfg_get(job_config, "min_shared")

  # Stages 1-3: aggregate hit blocks across canonical paths, then rank
  # the per-(component, path_type) representatives, dropping subset-
  # redundant ones.
  block_agg <- aggregate_blocks(c80s_coarse, allow_gaps = allow_gaps)
  reps_result <- rank_block_representatives(block_agg, min_overlap = min_overlap)
  representatives <- reps_result$representatives
  write.table(representatives, get_target("rep_path_df"), sep = "\t", quote = FALSE, row.names = FALSE)

  # Stage 4: per-genome attribution via the Step 3 provenance walker.
  rep_slim <- map_representatives_to_genomes(representatives, c_paths, collapsed_paths, path_df) %>%
    select(block_uid, canonical_uid, neighbor_genome, left_orig, right_orig) %>%
    unique()
  write.table(rep_slim, get_target("uid_path_df"), sep = "\t", quote = FALSE, row.names = FALSE)

  # Diagnostic: how often do surviving reps within a component share a
  # contiguous substring but neither contains the other?
  rep_overlap_diag <- diagnose_rep_overlaps(representatives, min_shared = min_shared)
  message(sprintf(
    "Rep overlap diagnostic: %d/%d (component,path_type) groups have true-overlap pairs; %d pairs total, max %d in one group.",
    rep_overlap_diag$summary$comps_with_overlap,
    rep_overlap_diag$summary$total_components,
    rep_overlap_diag$summary$total_overlap_pairs,
    rep_overlap_diag$summary$max_pairs_in_a_comp
  ))

  # Optional rep_heatmap.pdf: only rendered when at least two reps
  # survive and the block-by-genome matrix is at least 3x3 (pheatmap
  # requires multiple rows/columns to cluster).
  if (n_distinct(rep_slim$block_uid) > 1) {
    fd <- dirname(get_target("uid_path_df"))
    fp <- file.path(fd, "rep_heatmap.pdf")
    d <- rep_slim %>% select(block_uid, neighbor_genome) %>% mutate(v = 1)
    w <- d %>% spread(block_uid, v, fill = 0)
    m <- as.matrix(w %>% select(-neighbor_genome))
    if (nrow(m) > 2 && ncol(m) > 2) {
      pheatmap::pheatmap(m, color = c("white", "red"),
                         width = 7, height = 14, filename = fp)
    }
  }

  rep_slim
}
