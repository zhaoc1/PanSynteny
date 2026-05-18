# ------------------------------------------------------------------------------
# parse.R
#
# Step 4 orchestrator + reusable post-processing helpers for the Step 3
# canonical-paths outputs. `run_step4_parse` is the entry point called by
# pipeline.R; the individual helpers below are also safe to call standalone
# on Step 3 outputs (canonical_paths_c80s, canonical_paths_fine_c80s,
# canonical_paths_per_genome) for ad-hoc analysis.
#
#   decorate_c80s_w_smallORFs       -  decode synthetic small-ORF labels; add is_smallORF, centroid_80,
#                                         smallORF_type, n_smallORFs, dist_to_smallORFs.
#   decorate_c80s_w_truncation      -  flag truncated and fragmented c80s; add is_truncated, truncate_ratio,
#                                         n_truncated, is_fragmented, n_fragmented_c80s, fragmented_c80s.
#   summarize_coarse_operons        -  one row per `uid` (expects decorated coarse c80s).
#   summarize_fine_isoforms         -  one row per `uid_fine` (expects decorated fine c80s).
#   sample_genome_from_fine_paths   -  fine-only sampler (one genome per surviving isoform);
#                                       long-format per-gene table with provenance +
#                                       gene_neighbors metadata merged.
#   enrich_fine_long                -  left-join per-isoform context from c80s_fine
#                                       (trait stats, small-ORF + truncation flags,
#                                       centroid_80) onto fine_long.
#   write_blast_gene_lists          -  per-(uid_fine, neighbor_genome) gene-id TSVs for BLAST.
#   assign_c80_label                -  per-row pos/neg/neu/anchor label from a trait stat;
#                                       used by decorate_with_updated_path_type and the Step 5 plotters.
#   decorate_with_updated_path_type -  per-(component, canonical_path, path_type) summary;
#                                       adds c80_label_combo, purity_status, updated_path_type
#                                       (anchor_pos / anchor_neg / anchor_mixed). Facet variable
#                                       for the Step 5 per-component plotters.
#   run_step4_parse                 -  Step 4 orchestrator: builds summaries + selection sets,
#                                       samples one exemplar genome per surviving fine isoform,
#                                       writes five TSVs and the BLAST gene-id directory.
#
# Step 5 gggenes plotters (`plot_coarse/fine_operons`, `plot_*_by_component`,
# `run_step5_figures`) and their layout/scale helpers live in plot.R.
#
# Author:   Chunyu Zhao <chunyu.zhao@gladstone.ucsf.edu>
# Created:  2026-04-24
# ------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(purrr)
library(stringr)


#' Decode synthetic small-ORF labels and attach per-focal metrics
#'
#' Synthetic small-ORF labels (set by [compute_short_gene_prevalence()])
#' have the form `"_<focal_c80>-<gene_type>_<rank>"`. This function
#' explodes that encoding into queryable columns and attaches per-operon
#' small-ORF metrics. Six columns added; original columns and row order
#' are preserved.
#'
#' * **`is_smallORF`** - `TRUE` when the gene id starts with `_`.
#' * **`centroid_80`** - for small-ORF rows, the **focal_c80** decoded from
#'   the synthetic label (everything between the leading `_` and the last
#'   `-`). For normal rows, equal to `neighbor_c80_coarse`. The shared
#'   `centroid_80` is what makes the per-focal scoping in `dist_to_smallORFs`
#'   meaningful: a focal gene and the small ORFs derived from its
#'   neighborhood land in the same `(group_key, centroid_80)` subgroup.
#' * **`smallORF_type`** - for small-ORF rows, the `<gene_type>` token
#'   (`"CDS"`, `"tRNA"`, etc.); the trailing `_<rank>` disambiguator is
#'   stripped. NA for normals.
#' * **`n_smallORFs`** - per-operon count of small-ORF rows within
#'   `group_key`, broadcast to every row of that operon (same value on all
#'   rows of a given `uid` / `uid_fine`).
#' * **`n_focal`** - per-operon count of focal rows (`is_focal == TRUE`)
#'   within `group_key`, broadcast to every row of that operon. NA values
#'   in `is_focal` are coalesced to FALSE so non-`focal_meta` rows
#'   (short ORFs, neighbors absent from `focal_meta`) correctly count as
#'   non-focal.
#' * **`dist_to_smallORFs`** - for focal rows (`is_focal == TRUE`), the
#'   smallest positional distance within `group_key` to a small-ORF row
#'   sharing its `centroid_80` (a small ORF derived from the same focal).
#'   NA for non-focal rows, small-ORF rows themselves, and focals whose
#'   subgroup has no associated small ORFs.
#'
#' Per-focal scoping is the key biological signal: a small ORF derived
#' from focal A does not contribute to focal B's distance even when both
#' share a canonical path.
#'
#' Source column is hardcoded to `neighbor_c80_coarse`. The coarse table has
#' `neighbor_c80_coarse == c80` (set by [build_canonical_paths_c80s()]); the fine
#' table only has `neighbor_c80_coarse`. `neighbor_c80_fine` (fine only) is
#' explicitly *not* used here. Pass `group_key = "uid"` for coarse,
#' `"uid_fine"` for fine.
#'
#' The focal flag column is hardcoded to `is_focal` (added to `focal_meta`
#' by the driver). NA values are coalesced to FALSE.
#'
#' For how `neighbor_c80_coarse`, `neighbor_c80_fine`, and `centroid_80` differ
#' and when to use each, see `docs/PIPELINE.md`.
#'
#' @export
decorate_c80s_w_smallORFs <- function(df, group_key = "uid") {
  # c80_col is neighbor_c80_coarse for both level
  df <- df %>%
    mutate(
      is_smallORF = startsWith(as.character(neighbor_c80_coarse), "_"),
      .stripped = if_else(is_smallORF, str_sub(as.character(neighbor_c80_coarse), 2L), as.character(neighbor_c80_coarse)),
      centroid_80 = if_else(is_smallORF, str_replace(.stripped, "-[^-]+$", ""), .stripped),
      smallORF_type = if_else(is_smallORF, str_replace(str_extract(.stripped, "[^-]+$"), "_\\d+$", ""), NA_character_)) %>%
    select(-.stripped)

  df <- df %>%
    group_by(across(all_of(group_key))) %>%
    mutate(
      .row_in_path = row_number(),
      n_smallORFs  = sum(is_smallORF),
      n_focal      = sum(coalesce(is_focal, FALSE))
    ) %>%
    ungroup()

  df <- df %>%
    group_by(across(all_of(group_key)), centroid_80) %>%
    mutate(
      .small_orf_pos = list(.row_in_path[is_smallORF]),
      dist_to_smallORFs = if_else(
        coalesce(is_focal, FALSE),
        map2_dbl(.row_in_path, .small_orf_pos, ~{
          if (length(.y) == 0) NA_real_ else min(abs(.x - .y))
        }),
        NA_real_
      )
    ) %>%
    ungroup() %>%
    select(-.row_in_path, -.small_orf_pos)

  # Group the 6 added columns next to neighbor_c80_coarse so they sit with their
  # natural sibling identity columns. In the coarse table this places them
  # immediately before joint_component_ids; in the fine table, immediately
  # before the isoform-identity block (c80_path_fine_canonical, ...).
  df %>%
    relocate(is_smallORF, smallORF_type, n_smallORFs, n_focal, dist_to_smallORFs, centroid_80,
             .after = neighbor_c80_coarse)
}



#' Flag truncated and fragmented c80s; attach per-operon metrics (fine only)
#'
#' Fine-only by design. Truncation compares `neighbor_gene_length` against
#' `neighbor_c80_length_coarse` (the database centroid reference length) -
#' a coarse-table equivalent would be misleading because that table's
#' `neighbor_gene_length` is the max across isoforms (see
#' [build_canonical_paths_c80s()]), so "shorter than cutoff" would
#' actually mean "even the longest observed isoform is shorter".
#' Fragmentation detects whether a coarse cluster shows up at multiple
#' length variants within one isoform. `truncation_cutoff` is the only knob;
#' the rest of the inputs (`uid_fine`, `neighbor_c80_coarse`,
#' `neighbor_c80_fine`, `neighbor_gene_length`, `neighbor_c80_length_coarse`)
#' are read by hardcoded name.
#'
#' Six columns added; original columns and row order preserved.
#'
#' Truncation (length-based):
#' * **`is_truncated`** - per row: `TRUE` when `neighbor_gene_length <
#'   truncation_cutoff * neighbor_c80_length_coarse`. NA-guarded: synthetic
#'   small ORFs (`neighbor_c80_length_coarse = NA`) auto-excluded.
#' * **`truncate_ratio`** - per row: `neighbor_gene_length /
#'   neighbor_c80_length_coarse`, floor-truncated to 3 decimal places. NA for
#'   rows where either length is NA (short ORFs).
#' * **`n_truncated`** - per-isoform broadcast: `sum(is_truncated)` within
#'   `uid_fine`.
#'
#' Fragmentation (label-based - independent of truncation):
#' * **`is_fragmented`** - per row: `TRUE` when this row's `neighbor_c80_coarse`
#'   shows up under ≥2 distinct `neighbor_c80_fine` values within this
#'   `uid_fine` - the same coarse cluster observed at multiple lengths in
#'   one operon (the "split-gene" signature). Synthetic small ORFs
#'   (`neighbor_c80_coarse` starts with `_`) are excluded.
#' * **`fragmented_c80s`** - per row: this row's `neighbor_c80_coarse` when
#'   `is_fragmented = TRUE`, else `NA`. Slicing rows where this is non-NA
#'   gives the rows that participate in fragmentation, labeled by which
#'   coarse cluster each one belongs to.
#' * **`n_fragmented_c80s`** - per-isoform broadcast:
#'   `n_distinct(fragmented_c80s, na.rm = TRUE)` within `uid_fine`. Number
#'   of distinct coarse clusters in this operon that are fragmented.
#'
#' @export
decorate_c80s_w_truncation <- function(df, truncation_cutoff = 0.8) {
  # Per-row is_truncated + truncate_ratio: observed length vs centroid reference.
  # truncate_ratio is floor-truncated to 3 decimal places; NA propagates naturally
  # for short ORFs (neighbor_c80_length_coarse is NA for them).
  df <- df %>%
    mutate(
      is_truncated = !is.na(neighbor_c80_length_coarse) & !is.na(neighbor_gene_length) &
                     neighbor_gene_length < truncation_cutoff * neighbor_c80_length_coarse,
      truncate_ratio = floor(neighbor_gene_length / neighbor_c80_length_coarse * 1000) / 1000
    )
  
  # n_truncated per isoform (broadcast)
  df <- df %>%
    group_by(uid_fine) %>%
    mutate(n_truncated = sum(is_truncated, na.rm = TRUE)) %>%
    ungroup()
  
  # Per-row is_fragmented: this row's coarse neighbor_c80_coarse appears under
  # >=2 distinct neighbor_c80_fine values within this uid_fine. Excludes smallORFs.
  df <- df %>%
    group_by(uid_fine, neighbor_c80_coarse) %>%
    mutate(is_fragmented = !startsWith(as.character(neighbor_c80_coarse), "_") & n_distinct(neighbor_c80_fine) >= 2) %>%
    ungroup() %>%
    # Per-row: c80 id when fragmented, NA otherwise.
    mutate(fragmented_c80s = if_else(is_fragmented, as.character(neighbor_c80_coarse), NA_character_))
  
  # n_fragmented_c80s per isoform (broadcast)
  df <- df %>%
    group_by(uid_fine) %>%
    mutate(n_fragmented_c80s = n_distinct(fragmented_c80s, na.rm = TRUE)) %>%
    ungroup()
  
  df %>%
    relocate(is_truncated, truncate_ratio, n_truncated, is_fragmented, fragmented_c80s, n_fragmented_c80s, .after = dist_to_smallORFs) %>%
    relocate(uid, canonical_path_id, fine_canonical_id, .after = c80_path_fine_canonical)
}


#' One-row-per-operon coarse summary
#'
#' Expects `canonical_paths_c80s` to have already been decorated by
#' [decorate_c80s_w_smallORFs()] so `n_smallORFs` / `n_focal` /
#' `dist_to_smallORFs` are present. All `uid`s in the input already
#' cleared the `path_min_genomes` gate during canonical-path generation, so
#' every surviving operon is "recurring" by definition - no extra
#' filter is applied here. Truncation is not reported at this level
#' because the coarse `neighbor_gene_length` is the max across
#' isoforms (see [summarize_fine_isoforms()] for truncation-aware
#' calls).
#'
#' Output columns:
#' * Carry-through identity: `uid`, `path_type`, `n_genomes`. (Other
#'   identity columns - `joint_component_ids`, `canonical_path_id` -
#'   are encoded in `uid` and recoverable by parsing or by joining to
#'   `c_paths`.) `neighbor_genomes` is intentionally **not** carried:
#'   genome-level traceback at this level is rarely useful - drop down
#'   to `summarize_fine_isoforms()` or `c_paths_per_genome` for that.
#' * `n_genes` - number of c80 positions in the operon.
#' * `n_focal` / `n_smallORFs` - per-operon counts inherited from
#'   the small-ORF decoration.
#' * `min_dist_to_smallORFs` - smallest non-NA `dist_to_smallORFs`
#'   across focal rows in the operon, integer; NA if no focal has
#'   any associated small ORF.
#' * `smallORF_type_combo` - sorted, comma-joined unique values of
#'   `smallORF_type` across the small-ORF rows of the operon
#'   (e.g. `"CDS,tRNA"`). NA if the operon has no small ORFs.
#' * `coarse_path_string` - `neighbor_c80_coarse` tokens joined with
#'   `" → "` in canonical direction.
#'
#' Sorted by `desc(n_genomes)`.
#'
#' @export
summarize_coarse_operons <- function(canonical_paths_c80s) {

  carry_cols <- intersect(c("uid", "path_type", "n_genomes"), names(canonical_paths_c80s))

  canonical_paths_c80s %>%
    group_by(across(all_of(carry_cols))) %>%
    summarise(
      n_genes = n(),
      n_focal = first(n_focal),
      n_smallORFs = first(n_smallORFs),
      min_dist_to_smallORFs = {
        x <- dist_to_smallORFs[!is.na(dist_to_smallORFs)]
        if (length(x)) as.integer(min(x)) else NA_integer_
      },
      smallORF_type_combo = {
        x <- unique(smallORF_type[!is.na(smallORF_type)])
        if (length(x)) paste(sort(x), collapse = ",") else NA_character_
      },
      coarse_path_string = paste(neighbor_c80_coarse, collapse = " → "),
      .groups = "drop"
    ) %>%
    arrange(desc(n_genomes))
}


#' One-row-per-isoform fine summary
#'
#' Expects `canonical_paths_fine_c80s` to have already been decorated by
#' both [decorate_c80s_w_smallORFs()] and [decorate_c80s_w_truncation()],
#' so `n_smallORFs` / `n_focal` / `n_truncated` / `n_fragmented_c80s` /
#' `fragmented_c80s` are present. Rolls up per-isoform flags. No
#' minimum-support filter is applied - every surviving isoform from
#' [expand_canonical_paths_to_fine()] is included; downstream callers
#' can filter by `n_fine_genomes` as needed.
#'
#' Output columns:
#' * Carry-through identity: `uid_fine`, `uid`, `path_type`,
#'   `n_fine_genomes`, `n_genomes`, `fine_neighbor_genomes`. `uid` is
#'   carried so downstream callers (e.g. `sample_genome_from_fine_paths`)
#'   can group isoforms by their coarse parent without parsing
#'   `uid_fine` strings. Other identity columns (`joint_component_ids`,
#'   `canonical_path_id`, `fine_canonical_id`, `isoform_rank`) are
#'   encoded in `uid_fine` (`<uid>-iso<rank>-ngf<n_fine_genomes>`) and
#'   recoverable by parsing or by joining to `c_paths_fine`.
#' * `n_genes` - number of c80 positions in the isoform.
#' * `n_focal` / `n_smallORFs` / `n_truncated` / `n_fragmented_c80s`
#'   - per-isoform counts inherited from the decorations.
#' * `fragmented_c80s` - sorted comma-joined coarse cluster IDs that
#'   are fragmented in this isoform; NA if none.
#' * `fine_path_string` - `neighbor_c80_fine` tokens joined with
#'   `" → "` in canonical direction.
#'
#' Sorted by `desc(n_fine_genomes)`. `fine_neighbor_genomes` is
#' relocated to the last column so the lengthy genome list does not
#' clutter the leftmost columns.
#'
#' @export
summarize_fine_isoforms <- function(canonical_paths_fine_c80s) {

  carry_cols <- intersect(
    c("uid_fine", "uid", "path_type", "n_fine_genomes", "n_genomes", "fine_neighbor_genomes"),
    names(canonical_paths_fine_c80s)
  )

  canonical_paths_fine_c80s %>%
    group_by(across(all_of(carry_cols))) %>%
    summarise(
      n_genes = n(),
      n_focal = first(n_focal),
      n_smallORFs = first(n_smallORFs),
      n_truncated = first(n_truncated),
      n_fragmented_c80s = first(n_fragmented_c80s),
      fragmented_c80s = {
        x <- unique(fragmented_c80s[!is.na(fragmented_c80s)])
        if (length(x)) paste(sort(x), collapse = ",") else NA_character_
      },
      fine_path_string  = paste(neighbor_c80_fine, collapse = " → "),
      .groups = "drop"
    ) %>%
    arrange(desc(n_fine_genomes)) %>%
    relocate(fine_neighbor_genomes, .after = last_col())
}


#' Sample one exemplar genome per fine isoform and explode to per-gene rows
#'
#' Fine-only sampler. Draws **one** random genome per surviving fine
#' isoform; coarse-level BLAST hits are derived post-hoc by aggregating
#' on `sub("-iso.*", "", uid_fine)` (the `uid` prefix). The caller is
#' responsible for the survival filter (typically
#' `n_fine_genomes >= ceiling(path_min_genomes * fine_coverage_ratio)`).
#'
#' Output is long-format: one row per gene of the sampled genome's
#' canonical-direction path, with per-gene metadata (cluster ids,
#' contig coordinates, strand, gene type, length) merged from
#' `gene_neighbors`. Direction is canonical because the explode reads
#' `gene_path_canonical` (already direction-aligned by
#' `expand_canonical_paths_per_genome` via the L2 `needs_flip` flag);
#' `position_in_path` is the 1-indexed canonical position. The full
#' path string is **not** kept on the output - reconstruct it on
#' demand via `arrange(position_in_path) %>% summarise(paste(gene_id,
#' collapse = " → "))` per `(uid_fine, neighbor_genome)`.
#'
#' Output schema (12 columns):
#' * Sample-level: `uid_fine`, `neighbor_genome`.
#' * Position: `position_in_path` (1-indexed canonical), `gene_id`.
#' * Per-gene cluster context: `neighbor_c80_coarse`, `neighbor_c80_fine`.
#' * Per-gene coordinates: `neighbor_contig_id`, `neighbor_gene_start`,
#'   `neighbor_gene_end`, `neighbor_gene_strand`.
#' * Per-gene attributes: `neighbor_gene_type`, `neighbor_gene_length`.
#'
#' @param selected_fine Pre-filtered fine summary; required cols
#'   `uid_fine`, `uid`.
#' @param canonical_paths_per_genome Output of
#'   [`expand_canonical_paths_per_genome`]; required cols `uid_fine`,
#'   `neighbor_genome`, `gene_path`, `gene_path_canonical`.
#' @param gene_neighbors Step-1 cached neighbor table; required cols
#'   `neighbor_gene_id`, `neighbor_c80_coarse`, `neighbor_c80_fine`,
#'   `neighbor_contig_id`, `neighbor_gene_start`, `neighbor_gene_end`,
#'   `neighbor_gene_strand`, `neighbor_gene_type`, `neighbor_gene_length`.
#' @param seed Global RNG seed (default 616). Used to make the
#'   per-isoform draw reproducible across runs; not recorded in the
#'   output.
#'
#' @export
sample_genome_from_fine_paths <- function(selected_fine, canonical_paths_per_genome,
                                          gene_neighbors, seed = 616) {

  pool <- canonical_paths_per_genome %>%
    semi_join(selected_fine %>% select(uid_fine), by = "uid_fine") %>%
    select(uid_fine, neighbor_genome, gene_path, gene_path_canonical) %>%
    distinct()

  n_bad <- sum(is.na(pool$gene_path) | pool$gene_path == "")
  if (n_bad) {
    warning(sprintf("sample_genome_from_fine_paths: dropping %d (uid_fine, neighbor_genome) rows with NA/empty gene_path before sampling.", n_bad))
    pool <- pool %>% filter(!is.na(gene_path) & gene_path != "")
  }

  # One random genome per surviving fine isoform.
  set.seed(seed)
  sampled <- pool %>%
    group_by(uid_fine) %>%
    slice_sample(n = 1) %>%
    ungroup()

  # Gene-intrinsic per-gene metadata; dedup on neighbor_gene_id to
  # avoid row explosion (gene_neighbors carries the same gene multiple
  # times, once per (focal, neighbor) co-occurrence).
  neighbor_lookup <- gene_neighbors %>%
    select(neighbor_gene_id, neighbor_c80_coarse, neighbor_c80_fine,
           neighbor_contig_id:neighbor_gene_type, neighbor_gene_length) %>%
    distinct(neighbor_gene_id, .keep_all = TRUE)

  sampled %>%
    mutate(.gene_id = strsplit(gene_path_canonical, " → ")) %>%
    unnest_longer(.gene_id) %>%
    group_by(uid_fine, neighbor_genome) %>%
    mutate(position_in_path = row_number()) %>%
    ungroup() %>%
    rename(gene_id = .gene_id) %>%
    left_join(neighbor_lookup, by = c("gene_id" = "neighbor_gene_id")) %>%
    select(uid_fine, neighbor_genome,
           position_in_path, gene_id,
           neighbor_c80_coarse, neighbor_c80_fine,
           neighbor_contig_id, neighbor_gene_start, neighbor_gene_end, neighbor_gene_strand,
           neighbor_gene_type, neighbor_gene_length)
}


#' Enrich `fine_long` with per-isoform context from `c80s_fine`
#'
#' Left-joins the per-(uid_fine, c80) annotations from
#' `canonical_paths_fine_c80s` onto the per-gene-row long table
#' produced by [`sample_genome_from_fine_paths()`]. Brings in
#' per-isoform-broadcast counts (`n_focal`, `n_smallORFs`,
#' `n_truncated`, `n_fragmented_c80s`), per-c80 trait stats (`beta`,
#' `cor_to_b`, `is_focal`), per-c80 small-ORF and truncation flags,
#' and `centroid_80`. The full set of added columns is the contiguous
#' range `uid_fine:centroid_80` in `c80s_fine`, minus the three
#' exceptions noted below.
#'
#' Join key is `(uid_fine, position_in_path)`. `c80s_fine` does not
#' carry a position column, but its rows are in canonical order
#' within each `uid_fine` by construction (built by
#' [`build_canonical_paths_fine_c80s()`]); `row_number()` within
#' `uid_fine` therefore matches the position semantics that
#' `fine_long` already carries. Joining on `(uid_fine,
#' neighbor_c80_fine)` alone would silently row-multiply if a fine
#' path has the same coarse cluster at two positions with the same
#' length variant - `position_in_path` disambiguates.
#'
#' Three columns are excluded from the join payload to avoid
#' colliding with `fine_long`:
#' * `neighbor_c80_coarse` and `neighbor_c80_fine` already exist in
#'   `fine_long` with identical values per joined row; dropping is
#'   purely for hygiene.
#' * `neighbor_gene_length` is **intentionally not** brought in: in
#'   `fine_long` it is the per-gene observed length (from
#'   `gene_neighbors`), while in `c80s_fine` it is the per-isoform
#'   consensus length (from the canonical-fine builder). Same name,
#'   different semantics - keeping both would invite confusion.
#'
#' @param fine_long Output of [`sample_genome_from_fine_paths`].
#'   Required cols: `uid_fine`, `position_in_path`.
#' @param c80s_fine Decorated fine c80s frame
#'   ([`canonical_paths_fine_c80s`]); must contain `uid_fine` and
#'   the contiguous column range `uid_fine:centroid_80`.
#'
#' @export
enrich_fine_long <- function(fine_long, c80s_fine) {
  c80s_fine_keyed <- c80s_fine %>%
    group_by(uid_fine) %>%
    mutate(position_in_path = row_number()) %>%
    ungroup() %>%
    select(position_in_path, uid_fine:centroid_80) %>%
    select(-any_of(c("neighbor_c80_coarse", "neighbor_c80_fine", "neighbor_gene_length")))

  out <- fine_long %>%
    left_join(c80s_fine_keyed, by = c("uid_fine", "position_in_path"))

  # Sanity: the join key is unique on the c80s_fine side by construction,
  # so row count must match. Detects any row-multiplication bug.
  stopifnot(nrow(out) == nrow(fine_long))
  out
}


# -----------------------------------------------------------------------------
# BLAST gene-list writer
# -----------------------------------------------------------------------------

#' Write one TSV per `(uid_fine, neighbor_genome)` with gene IDs, one per line
#'
#' File name pattern: `fine_{uid_fine}_{neighbor_genome}.tsv`. Consumes the
#' long-format output of [`sample_genome_from_fine_paths`] directly; no
#' re-splitting of `gene_path_canonical` is needed. Rows are written in
#' canonical order (`position_in_path` ascending).
#'
#' @export
write_blast_gene_lists <- function(sampled_long, out_dir) {
  stopifnot(all(c("uid_fine", "neighbor_genome", "position_in_path", "gene_id")
                %in% names(sampled_long)))

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  sampled_long %>%
    arrange(uid_fine, neighbor_genome, position_in_path) %>%
    group_by(uid_fine, neighbor_genome) %>%
    group_walk(function(.x, .y) {
      outfile <- file.path(
        out_dir,
        paste0("fine_", .y$uid_fine, "_", .y$neighbor_genome, ".tsv")
      )
      writeLines(.x$gene_id, con = outfile)
    })

  invisible(NULL)
}


# -----------------------------------------------------------------------------
# Plot data-prep helpers (used by the Step 5 plotters in plot.R)
# -----------------------------------------------------------------------------

#' Label c80 rows by focal-direction or anchor role
#'
#' Adds a `c80_label` column to a c80s frame for use by the Step 5
#' plotters in `plot.R`. Focal-aware: only focal rows (`is_focal == TRUE`)
#' get direction labels; non-focal rows on anchor canonical paths get the
#' `"anchor"` label so they show up as the per-path frame of reference.
#'
#' Labels:
#' * `"pos"`: focal row with `value_col > 0`.
#' * `"neg"`: focal row with `value_col < 0`.
#' * `"neu"`: focal row with `value_col == 0` (rare; a focal gene whose
#'   trait association is exactly zero).
#' * `"anchor"`: non-focal row whose canonical path is `path_type ==
#'   "anchor"`.
#' * `NA`: everything else (non-focal rows on non-anchor paths; focal
#'   rows with `NA` `value_col`).
#'
#' Step 5 glyph mapping (`.layout_operon_tracks` in plot.R): `pos` -> `"U"`,
#' `neg` -> `"D"`, `neu` -> `"N"`, `anchor` -> `"O"`. Coloring (`fill_gene`)
#' is applied to `pos`/`neg`/`anchor` rows; `neu` carries the glyph but no
#' fill.
#'
#' @param df A c80s frame; must contain `is_focal`, `path_type`, and
#'   the column named by `value_col`.
#' @param value_col Name of the trait-association column to read sign
#'   from (default `"cor_to_b"`; Step 4 passes `"beta"`).
#'
#' @export
assign_c80_label <- function(df, value_col = "cor_to_b") {
  # No trait score - use focal_label directly if available
  if (!value_col %in% names(df)) {
    if ("focal_label" %in% names(df)) {
      return(df %>% mutate(c80_label = if_else(coalesce(is_focal, FALSE), focal_label, NA_character_)))
    }
    return(df %>% mutate(c80_label = NA_character_))
  }

  stopifnot(all(c("is_focal", "path_type", value_col) %in% names(df)))

  df %>%
    mutate(
      .v        = .data[[value_col]],
      .focal    = coalesce(is_focal, FALSE),
      c80_label = case_when(
        .focal  & !is.na(.v) & .v > 0   ~ "pos",
        .focal  & !is.na(.v) & .v < 0   ~ "neg",
        .focal  & !is.na(.v) & .v == 0  ~ "neu",
        !.focal & path_type == "anchor" ~ "anchor",
        TRUE                            ~ NA_character_
      )
    ) %>%
    select(-.v, -.focal)
}

#' Decorate paths with a `updated_path_type` reflecting per-component pos/neg signal
#'
#' Walks `(joint_component_id, canonical_path_id, path_type)` groups,
#' summarizes the unique `c80_label` values, and rewrites
#' `path_type == "anchor"` to `anchor_pos` / `anchor_neg` (or
#' `anchor_mixed` when both signs co-occur). Adds three columns:
#' `c80_label_combo`, `purity_status` (`"impure"` if both pos and neg
#' present, else `"pure"`), `updated_path_type`.
#'
#' Calls [`assign_c80_label`] internally with `value_col = "beta"` if
#' `c80_label` is missing: so this helper is safe to invoke directly
#' on a c80s frame post-decorate.
#'
#' Used by the per-component plotters as the facet variable.
#'
#' @export
decorate_with_updated_path_type <- function(df) {
  if (!"c80_label" %in% names(df)) df <- assign_c80_label(df, value_col = "beta")

  per_path <- df %>%
    filter(!is.na(c80_label)) %>%
    group_by(joint_component_id, canonical_path_id, path_type) %>%
    summarise(
      c80_label_combo = paste(sort(unique(c80_label)), collapse = ","),
      .groups = "drop"
    ) %>%
    mutate(
      .n_labels = str_count(c80_label_combo, ",") + 1L,
      purity_status = if_else(.n_labels > 1L, "impure", "pure"),
      updated_path_type = case_when(
        path_type == "anchor" & purity_status == "impure" ~ paste0("anchor_mixed"),
        path_type == "anchor" & purity_status == "pure"   ~ paste0("anchor_", c80_label_combo),
        TRUE                                              ~ path_type
      )
    ) %>%
    select(-.n_labels)

  df %>%
    left_join(per_path, by = c("joint_component_id", "canonical_path_id", "path_type")) %>%
    mutate(
      purity_status     = coalesce(purity_status,     "pure"),
      updated_path_type = coalesce(updated_path_type, path_type)
    )
}


#' Run Step 4: summaries, fine-coverage selection, exemplar sampling, BLAST gene lists
#'
#' Orchestrator for Step 4. Builds per-operon and per-isoform summaries,
#' applies the fine-coverage isoform-survival filter, attaches
#' isoform-map columns to `coarse_summary`, samples one exemplar genome
#' per surviving fine isoform, enriches the long-format result with
#' per-isoform context, and writes per-`(uid_fine, neighbor_genome)`
#' gene-id TSVs for the external BLAST workflow.
#'
#' Inputs are passed in explicitly so the caller (pipeline.R) is the one
#' place that loads the three Step 3 TSVs from disk; this keeps Step 4
#' re-runnable in isolation while making the call-site signature
#' self-documenting. The plotting block lives in [run_step5_figures()]
#' (plot.R) so a re-render after editing `fill_modes` does not require
#' re-running the rest of Step 4.
#'
#' Reads `path_min_genomes`, `fine_coverage_ratio`, and `seed` from
#' `job_config` via `cfg_get`. Writes five TSVs via `get_target`:
#' `parse_coarse_summary`, `parse_fine_summary`, `parse_selected_coarse`,
#' `parse_selected_fine`, `parse_fine_long`, and a directory of per-
#' `(uid_fine, neighbor_genome)` gene-id files at `parse_genome_paths_dir`.
#'
#' @param c80s_coarse Decorated L1 per-gene table (read from
#'   `canonical_paths_c80s`).
#' @param c80s_fine Decorated L2 per-isoform per-gene table (read from
#'   `canonical_paths_fine_c80s`).
#' @param per_genome L3 per-genome master table (read from
#'   `canonical_paths_per_genome`).
#' @param gene_neighbors Step 1 output (in-memory). Used only by
#'   [sample_genome_from_fine_paths()] to attach per-gene metadata
#'   (contig coordinates, strand, observed gene length) to the chosen
#'   exemplars.
#'
#' @return The enriched per-gene long-format frame `fine_long` (one row
#'   per gene of each sampled `(uid_fine, neighbor_genome)`). The other
#'   intermediate frames (summaries, selection sets) are persisted to
#'   disk only.
#'
#' @export
run_step4_parse <- function(c80s_coarse, c80s_fine, per_genome, gene_neighbors) {
  # 1. Summaries
  coarse_summary <- summarize_coarse_operons(c80s_coarse)
  fine_summary <- summarize_fine_isoforms(c80s_fine)

  # 2. Fine-coverage isoform survival filter.
  path_min_genomes <- cfg_get(job_config, "path_min_genomes")
  fine_coverage_ratio <- cfg_get(job_config, "fine_coverage_ratio")
  selected_fine <- fine_summary %>% filter(n_fine_genomes >= ceiling(path_min_genomes * fine_coverage_ratio))
  selected_coarse <- coarse_summary %>% semi_join(selected_fine, by = "uid")

  # Attach isoform map (raw + filtered counts + surviving uid_fine list) to coarse_summary.
  coarse_summary <- coarse_summary %>%
    left_join(
      fine_summary %>% group_by(uid) %>% summarise(n_isoforms_raw = n(), .groups = "drop"),
      by = "uid"
    ) %>%
    left_join(
      selected_fine %>%
        group_by(uid) %>%
        summarise(
          n_isoforms_filtered = n(),
          n_coarse_genome_filtered = sum(n_fine_genomes),
          uid_fine_list = paste(sort(uid_fine), collapse = ";"),
          .groups = "drop"
        ),
      by = "uid"
    ) %>%
    mutate(
      n_isoforms_raw = coalesce(n_isoforms_raw, 0L),
      n_isoforms_filtered = coalesce(n_isoforms_filtered, 0L),
      n_coarse_genome_filtered = coalesce(n_coarse_genome_filtered, 0L)) %>%
    relocate(n_isoforms_raw, n_isoforms_filtered, n_coarse_genome_filtered, .after = n_smallORFs) %>%
    relocate(uid_fine_list, .after = last_col())

  write.table(coarse_summary, get_target("parse_coarse_summary"),  sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(fine_summary, get_target("parse_fine_summary"),    sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(selected_coarse, get_target("parse_selected_coarse"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(selected_fine, get_target("parse_selected_fine"),   sep = "\t", row.names = FALSE, quote = FALSE)

  # 3. Exemplar genome sampling + per-isoform context enrichment.
  fine_long <- sample_genome_from_fine_paths(
    selected_fine, per_genome, gene_neighbors,
    seed = cfg_get(job_config, "seed")
  )
  fine_long <- enrich_fine_long(fine_long, c80s_fine)
  write.table(fine_long, get_target("parse_fine_long"), sep = "\t", row.names = FALSE, quote = FALSE)

  # 4. Per-(uid_fine, neighbor_genome) gene-id TSVs for BLAST.
  write_blast_gene_lists(fine_long, get_target("parse_genome_paths_dir"))

  fine_long
}
