# ------------------------------------------------------------------------------
# path.R
#
# Canonical-path expansion: turn canonical_paths and canonical_paths_fine into
# the long-format, join-able TSVs that downstream analysis consumes.
#
#   explode_canonical_into_collapsed_paths   — shared provenance walker (canonical -> collapsed → per-genome).
#   expand_canonical_paths_to_fine           — per-isoform aggregate view.
#   expand_canonical_paths_per_genome        — per-genome master table.
#   build_canonical_paths_c80s               — coarse per-gene anchor (max length).
#   build_canonical_paths_fine_c80s          — fine per-gene anchor (per-isoform length).
#   run_step3_consolidation                  — Step 3 orchestrator: builds all L1/L2/L3 frames and writes the five TSVs.
#
# Author:  Chunyu Zhao <chunyu.zhao@gladstone.ucsf.edu>
# Created: 2025-10-10 (extracted from graph.R)
# Updated: 2026-04-28
# ------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(purrr)


#' Explode canonical -> collapsed -> per-genome provenance into long format
#'
#' Walk the Step 3 provenance chain (`canonical_path_id` ->
#' `collapsed_path_id` -> `per_genome_path_w_ids` -> `path_genome_comp`) once,
#' emitting one row per (canonical_path_id, contributing per-genome path).
#' Shared backbone for any derivation that needs per-genome provenance under
#' a canonical identity:isoform expansion, genome-lookup tables, trait-
#' block attribution, etc.
#'
#' The output includes a `needs_flip` flag comparing each collapsed row's
#' original coarse direction to the final canonical direction, so callers
#' that need direction-aware payloads (fine-path expansion, edge
#' normalization) can use it without recomputing.
#'
#' @details
#' **Coarse -> fine inheritance (why one boolean is enough).**
#' `needs_flip` is computed on the coarse strings only
#' (`c80_path_coarse != canonical_path_coarse`) before the per-genome
#' explode. After [separate_rows()][tidyr::separate_rows] expands
#' `per_genome_path_w_ids`, every per-genome row that descends from the
#' same collapsed row inherits that collapsed row's `needs_flip` value.
#' Callers that need to flip the fine string (or any other
#' direction-sensitive per-position payload) just consume the inherited
#' boolean and call `rev()` on the token vector — no fine-level direction
#' decision is required.
#'
#' This is sound because both [normalize_path()] (R2) and
#' [orient_paths_within_component()] (R3) are **orientation-only**: they pick
#' between `forward` and `reverse(forward)` of the same token vector and
#' never reorder tokens within a direction or change token content. As a
#' result, coarse and fine path strings stay **position-aligned** — the
#' fine token at position k corresponds to the coarse cluster at
#' position k — so a coarse reversal at the position level is exactly
#' the right operation at fine resolution too. One coarse boolean drives
#' all per-position flips.
#'
#' **Invariant under future edits.** If a future change introduces any
#' token-reordering operation in [normalize_path()] or
#' [orient_paths_within_component()] (sorting, dedup of non-adjacent
#' duplicates, suffix re-ranking, etc.), this inheritance breaks and
#' per-payload direction decisions would need to be reintroduced. The
#' current `clean_for_orientation` helper preserves the invariant by
#' operating on a *copy* used only for the decision; the chosen
#' direction is always applied to the original full token vector.
#'
#' @export
explode_canonical_into_collapsed_paths <- function(c_paths, collapsed_paths, path_df,
                                                   extra_path_df_cols = character(0)) {
  # 1) canonical -> collapsed (explode the comma-joined list)
  can_to_col <- c_paths %>%
    select(canonical_path_id, canonical_path_coarse = c80_path_coarse_canonical, collapsed_path_id) %>%
    separate_rows(collapsed_path_id, sep = ",")

  # 2) collapsed -> original coarse direction + per-genome provenance
  col_info <- collapsed_paths %>%
    select(collapsed_path_id, c80_path_coarse, per_genome_path_w_ids)

  # 3) flag which collapsed rows are reversed relative to the canonical direction
  can_col <- can_to_col %>%
    left_join(col_info, by = "collapsed_path_id") %>%
    mutate(needs_flip = (c80_path_coarse != canonical_path_coarse))

  # 4) explode per_genome_path_w_ids into one row per contributing per-genome path
  can_col_exp <- can_col %>%
    separate_rows(per_genome_path_w_ids, sep = ";")

  # 5) attach per-genome columns from path_df (plus any caller-requested payload).
  # Join on path_df$path_genome_comp; the left-side name per_genome_path_w_ids is preserved.
  keep_from_path_df <- unique(c("path_genome_comp", "neighbor_genome", "path_type",
                                extra_path_df_cols))
  can_col_exp %>%
    left_join(
      path_df %>% select(all_of(keep_from_path_df)),
      by = c("per_genome_path_w_ids" = "path_genome_comp")
    )
}


#' Expand surviving canonical paths into isoform-resolved variants
#'
#' For each canonical path that cleared the Step 3 `path_min_genomes` gate, walk
#' the provenance chain via [explode_canonical_into_collapsed_paths()] with
#' `c80_path_fine` as the payload, align each fine rendering to the
#' coarse canonical direction, aggregate per isoform, and attach all
#' canonical-level columns (minus the coarse `neighbor_genomes`, which the
#' per-isoform `fine_neighbor_genomes` unions back to). Produces a parallel
#' fine-grained view without touching any existing Step 3 function.
#'
#' The coarse `c80_path_coarse_canonical` column is used as the direction anchor: a
#' contributing collapsed row whose `c80_path_coarse` does not match the
#' canonical direction is considered "flipped" and its fine string is
#' token-reversed to align. Because both [normalize_path()] and
#' [orient_paths_within_component()] only reverse paths — never change token
#' content — this binary equality check is sufficient; no independent
#' normalization of the fine string is required, which sidesteps palindrome
#' and suffix-reordering ambiguity.
#'
#' @export
expand_canonical_paths_to_fine <- function(c_paths, collapsed_paths, path_df) {
  fine_rows <- explode_canonical_into_collapsed_paths(
    c_paths, collapsed_paths, path_df,
    extra_path_df_cols = "c80_path_fine"
  )

  # apply coarse-driven flip to align fine string with canonical direction
  fine_rows <- fine_rows %>%
    rowwise() %>%
    mutate(c80_path_fine_canonical = if (needs_flip) {
        collapse_path(rev(split_path_string(c80_path_fine)))
      } else {
        c80_path_fine
      }) %>%
    ungroup()

  # aggregate per (canonical_path_id, c80_path_fine_canonical) + compute isoform identity
  fine_summary <- fine_rows %>%
    group_by(canonical_path_id, c80_path_fine_canonical) %>%
    summarise(
      n_fine_genomes = n_distinct(neighbor_genome),
      fine_neighbor_genomes = paste(sort(unique(neighbor_genome)), collapse = ";"),
      .groups = "drop"
    ) %>%
    group_by(canonical_path_id) %>%
    arrange(desc(n_fine_genomes), c80_path_fine_canonical, .by_group = TRUE) %>%
    mutate(
      # Per-canonical sequential rank (1 = most-supported isoform). Ties on
      # n_fine_genomes are broken alphabetically by c80_path_fine_canonical so the
      # numbering is deterministic and isoform_rank is strictly unique within
      # canonical_path_id — required so uid_fine identifies one isoform.
      isoform_rank = row_number(),
      # Self-describing composite key: canonical_path_id + iso + rank.
      # E.g. "cp_5_iso1", "cp_5_iso2". Unique within canonical_path_id.
      fine_canonical_id = paste0(canonical_path_id, "_iso", isoform_rank)
    ) %>%
    ungroup()

  # attach canonical-level context (all columns except coarse neighbor_genomes)
  keep_from_canonical <- setdiff(colnames(c_paths), "neighbor_genomes")
  fine_summary %>%
    left_join(c_paths %>% select(all_of(keep_from_canonical)), by = "canonical_path_id") %>%
    # Hierarchical fine uid (option B): extend coarse `uid` with iso rank and
    # per-isoform genome count. Stripping `-iso\d+-ngf\d+$` recovers the coarse
    # parent uid. Requires `uid` to already exist on c_paths.
    mutate(uid_fine = paste0(uid, "-iso", isoform_rank, "-ngf", n_fine_genomes)) %>%
    select(all_of(keep_from_canonical),
           c80_path_fine_canonical, fine_canonical_id, isoform_rank, uid_fine,
           n_fine_genomes, fine_neighbor_genomes) %>%
    select(uid_fine, everything())
}


#' Build per-genome master table of canonical path contributions
#'
#' For each canonical path, emit one row per contributing per-genome
#' observation with full attribution: which genome, which collapsed group,
#' whether the observation was reverse-oriented, and raw + canonical-aligned
#' path renderings at three resolutions (gene-id, fine c80-label, coarse c80).
#' Complementary to [expand_canonical_paths_to_fine()] — the two functions
#' walk the same provenance chain but stop at different aggregation grains.
#'
#' @export
expand_canonical_paths_per_genome <- function(c_paths, collapsed_paths, path_df, c_paths_fine) {

  rows <- explode_canonical_into_collapsed_paths(c_paths, collapsed_paths, path_df,
    extra_path_df_cols = c("c80_path_fine", "path_string"))

  # Compute canonical-aligned gene-id rendering (kept in output) and a
  # fine-c80 canonical string used only as the join key into c_paths_fine
  # (dot-prefixed; dropped at the end).
  # `needs_flip` was set in the explode walker by comparing c80_path_coarse
  # vs canonical_path_coarse — one boolean drives both flips because
  # normalize_path + orient_paths_within_component are orientation-only.
  rows <- rows %>%
    rowwise() %>%
    mutate(
      gene_path_canonical = if (needs_flip) {
          collapse_path(rev(split_path_string(path_string)))
        } else {
          path_string
        },
      .c80_path_fine_canonical = if (needs_flip) {
          collapse_path(rev(split_path_string(c80_path_fine)))
        } else {
          c80_path_fine
        }
    ) %>%
    ungroup()

  # Inherit canonical-level columns from c_paths (skip ones already on the
  # table or not relevant per-row).
  keep_from_canonical <- setdiff(colnames(c_paths), 
                                 c("neighbor_genomes", "c80_path_coarse_canonical", "path_type", "collapsed_path_id"))
  rows <- rows %>%
    left_join(c_paths %>% select(all_of(keep_from_canonical)), by = "canonical_path_id")

  # Inherit isoform-level identity from c_paths_fine (uid_fine,
  # fine_canonical_id, isoform_rank, n_fine_genomes).
  rows <- rows %>%
    left_join(
      c_paths_fine %>% select(uid_fine, n_fine_genomes, canonical_path_id,
                              c80_path_fine_canonical, fine_canonical_id, isoform_rank) %>% distinct(),
      by = c("canonical_path_id", ".c80_path_fine_canonical" = "c80_path_fine_canonical"))

  # Final output (Option A — lean): keep join keys + per-genome attribution +
  # gene-id renderings + needs_flip. Drop the upstream path-string payloads
  # (each is 1:1 with an ID we kept, or recoverable from gene_path) and the
  # temp join key. Any path string can be re-attached via a join to
  # c_paths / c_paths_fine / collapsed_paths if needed downstream.
  rows %>%
    dplyr::rename(gene_path = path_string) %>%
    select(-c80_path_coarse, -c80_path_fine,
           -canonical_path_coarse, -.c80_path_fine_canonical) %>%
    relocate(
      uid_fine, uid,
      canonical_path_id,
      collapsed_path_id,
      fine_canonical_id,
      per_genome_path_w_ids,
      neighbor_genome,
      gene_path, gene_path_canonical,
      needs_flip
    )
}


#' Explode canonical paths into per-gene rows with annotations (coarse)
#'
#' Turn each canonical path into a long-format, gene-level table at coarse
#' `centroid_80` resolution and attach per-gene annotations: joint component
#' membership, path type + genome support, a representative gene length per
#' c80 cluster (max over isoforms — see Details), microslam-only gene
#' metadata, and cluster_80 metadata. Output is the main per-gene coarse
#' anchor consumed by downstream trait analysis, block aggregation, and
#' per-component plotting.
#'
#' @details
#' **`neighbor_gene_length` is max-over-variants**: each c80 cluster gets
#' its longest observed length. Use [build_canonical_paths_fine_c80s()]
#' for a parallel table that keeps per-isoform exact lengths.
#'
#' @export
build_canonical_paths_c80s <- function(c_paths, c80_variants_mapping, focal_c80_df, cluster_80, jc_map, short_gene_prevalence) {
  c80s_coarse <- c_paths %>%
    select(canonical_path_id, c80_path_coarse_canonical) %>%
    mutate(gene = strsplit(c80_path_coarse_canonical, " → ")) %>%
    unnest(gene) %>%
    left_join(jc_map, by = c("gene" = "node")) %>%
    select(-c80_path_coarse_canonical)

  # inherit canonical-path identity (uid, joint_component_ids, path_type,
  # n_genomes, neighbor_genomes)
  c80s_coarse <- c80s_coarse %>%
    select(joint_component_id, everything()) %>%
    dplyr::rename(neighbor_c80_coarse = gene) %>%
    left_join(c_paths %>% select(uid, joint_component_ids, canonical_path_id, path_type, n_genomes, neighbor_genomes) %>% unique(), by = c("canonical_path_id"))

  # Coarse gene length: collapse c80_variants_mapping to one row per c80 (max length).
  # IMPORTANT — the resulting `neighbor_gene_length` column is misnamed.
  # It is NOT a single observed gene's length. It is the MAXIMUM length
  # observed across all isoforms of this c80 cluster in the data.
  # Sometimes, the full c80 is not in the operon, => neighbor_gene_length < centroid_80_length
  gl <- c80_variants_mapping %>%
    select(-neighbor_c80_fine) %>%
    unique() %>%
    group_by(neighbor_c80_coarse) %>%
    filter(neighbor_gene_length == max(neighbor_gene_length)) %>%
    ungroup()

  c80s_coarse <- c80s_coarse %>%
    left_join(gl, by = "neighbor_c80_coarse")

  # add gene label information back
  c80s_coarse <- c80s_coarse %>%
    left_join(focal_c80_df, by = c("neighbor_c80_coarse" = "focal_c80"))

  # cluster_80 metadata (NA for small ORFs)
  c80s_coarse <- c80s_coarse %>%
    left_join(cluster_80, by = c("neighbor_c80_coarse" = "c80"))

  # Fill in genome_prevalence for short-gene rows from the short_gene_prevalence mapping.
  c80s_coarse <- c80s_coarse %>%
    left_join(
      short_gene_prevalence %>% dplyr::rename(neighbor_c80_coarse = neighbor_c80_fine, short_genome_prevalence = genome_prevalence),
      by = "neighbor_c80_coarse") %>%
    mutate(genome_prevalence = coalesce(genome_prevalence, short_genome_prevalence)) %>%
    select(-short_genome_prevalence)

  # Final column order: uid first as primary key; uid-construction inputs
  # (joint_component_ids, canonical_path_id, path_type, n_genomes) and the
  # per-c80 component id parked at the end so they're easy to spot-check
  # without cluttering the data columns. Middle columns keep their join order.
  c80s_coarse <- c80s_coarse %>%
    relocate(uid) %>%
    relocate(joint_component_ids, canonical_path_id, path_type, n_genomes, joint_component_id, .after = last_col())
  
  return(c80s_coarse)
}


#' Explode canonical paths into per-gene rows at isoform (fine) resolution
#'
#' Mirror of [build_canonical_paths_c80s()] at isoform granularity.
#' Starting from `canonical_paths_fine` (one row per
#' `(canonical_path_id, c80_path_fine_canonical)`), this function splits each
#' `c80_path_fine_canonical` path into per-position isoform tokens and attaches the
#' exact isoform-specific gene length from `c80_variants_mapping`.
#'
#' @details
#' **`isoform_rank`** is assigned via `row_number()` within each
#' `canonical_path_id`, sorted by `desc(n_fine_genomes)` then by
#' `c80_path_fine_canonical` alphabetically — rank 1 = most-supported isoform,
#' ties broken by canonical-string order so the rank is strictly unique
#' (and so `uid_fine` identifies one isoform).
#'
#' Only `neighbor_c80_fine` and `neighbor_gene_length` are isoform-specific.
#' Joint component assignments and other metadata are inherited from the
#' parent coarse `neighbor_c80_coarse`.
#'
#' @export
build_canonical_paths_fine_c80s <- function(c_paths_fine, c80_variants_mapping,
                                            focal_c80_df, cluster_80, jc_map, short_gene_prevalence) {
  # inherit identity (canonical + isoform) and tokenize c80_path_fine_canonical per position.
  # isoform_rank, fine_canonical_id, uid, and uid_fine are inherited from
  # c_paths_fine (computed once in expand_canonical_paths_to_fine).
  fine_c80s <- c_paths_fine %>%
    select(uid, joint_component_ids, canonical_path_id, path_type, n_genomes,
           c80_path_fine_canonical, fine_canonical_id, isoform_rank, uid_fine, n_fine_genomes, fine_neighbor_genomes) %>%
    mutate(neighbor_c80_fine = strsplit(c80_path_fine_canonical, " → ")) %>%
    unnest(neighbor_c80_fine)

  # Per-isoform length + coarse-c80 bridge. neighbor_gene_length here is the
  # exact per-isoform value (not max-over-isoforms; see build_canonical_paths_c80s
  # for that policy). neighbor_c80_coarse arrives from this join and serves as the key
  # for all coarse-grained annotations below.
  fine_c80s <- fine_c80s %>%
    left_join(c80_variants_mapping, by = "neighbor_c80_fine")
  # jc_map, focal_c80_df, and cluster_80 don't have isoform granularity

  # coarse-keyed annotations (component / microslam / cluster_80)
  fine_c80s <- fine_c80s %>%
    left_join(jc_map, by = c("neighbor_c80_coarse" = "node")) %>%
    left_join(focal_c80_df, by = c("neighbor_c80_coarse" = "focal_c80")) %>%
    left_join(cluster_80, by = c("neighbor_c80_coarse" = "c80"))

  # Fill in genome_prevalence for short-gene rows from the short_gene_prevalence mapping.
  fine_c80s <- fine_c80s %>%
    left_join(
      short_gene_prevalence %>% dplyr::rename(short_genome_prevalence = genome_prevalence),
      by = "neighbor_c80_fine"
    ) %>%
    mutate(genome_prevalence = coalesce(genome_prevalence, short_genome_prevalence)) %>%
    select(-short_genome_prevalence)

  # Final column order: uid_fine and uid first as primary keys; isoform-level
  # identity (c80_path_fine_canonical, fine_canonical_id, isoform_rank, n_fine_genomes,
  # fine_neighbor_genomes) and canonical-level identity (joint_component_ids,
  # canonical_path_id, path_type, n_genomes) plus the per-c80 component id
  # parked at the end so they're easy to spot-check without cluttering the data
  # columns. Middle columns keep their join order. Mirrors build_canonical_paths_c80s.
  fine_c80s <- fine_c80s %>%
    relocate(uid_fine, uid, fine_canonical_id, n_fine_genomes) %>%
    relocate(c80_path_fine_canonical, isoform_rank, fine_neighbor_genomes, 
             joint_component_ids, canonical_path_id, path_type, n_genomes, joint_component_id,
             .after = last_col())

  return(fine_c80s)
}


#' Run Step 3 — cross-genome consolidation
#'
#' Orchestrator that takes per-genome maximal paths (`path_df` from Step 2)
#' and returns the canonical-operon view at three granularity levels (L1
#' coarse, L2 per-isoform, L3 per-genome), decorated with small-ORF and
#' truncation/fragmentation flags. Writes the five canonical-paths TSVs
#' (`canonical_paths_coarse`, `canonical_paths_fine`,
#' `canonical_paths_per_genome`, `canonical_paths_c80s`,
#' `canonical_paths_fine_c80s`) to the locations resolved by
#' [`get_target()`](model.R).
#'
#' Reads `path_min_genomes` and `truncation_cutoff` from the global `job_config`
#' via `cfg_get`. Always re-runs (no on-disk cache gate, unlike Steps 1/2);
#' the five TSVs are rewritten every call.
#'
#' @param path_df Output of [stitch_paths_across_focal_genes()].
#' @param c80_variants_mapping Per-cluster length-variant labels (Step 1).
#' @param focal_c80_df Focal centroid table (`focal_c80`, `is_focal`, ...).
#' @param cluster_80 MIDAS coarse cluster metadata.
#' @param short_gene_prevalence Synthetic-c80 prevalence map for short ORFs.
#'
#' @return A list with the in-memory frames Step 4 still needs:
#'   * `c_paths` — L1 canonical-paths table
#'   * `collapsed_paths` — pre-canonicalisation collapse, used by
#'     `map_representatives_to_genomes()`
#'   * `c80s_coarse` — L1 per-gene table fed to `aggregate_blocks()`
#'
#' The other three frames (`c_paths_fine`, `c80s_fine`, `c_paths_per_genome`)
#' are written to disk only; Step 5 re-reads them from disk so this
#' function does not return them.
#'
#' @export
run_step3_consolidation <- function(path_df, c80_variants_mapping,
                                    focal_c80_df, cluster_80, short_gene_prevalence) {

  # Walks per-genome paths up to canonical operon identities, then back down to
  # per-isoform and per-genome views, producing three granularity levels of
  # operon tables (each emitted at one or two resolutions, see writes below):
  #
  #   Level 1 (coarse / canonical):
  #     - collapse_paths_across_genomes()    — one row per coarse-string + path_type
  #     - generate_canonical_path()          — unify forward/reverse, gate on path_min_genomes
  #     - compute_joint_components()         — gene-level connected components (joint, type-blind)
  #     - decorate_paths_with_components()   — attach joint_component_ids per path
  #     - orient_paths_within_component()    — align within-component direction
  #     - build_canonical_paths_c80s()       — explode to per-gene rows (coarse)
  #   Level 2 (per-isoform):
  #     - expand_canonical_paths_to_fine()       — resolve length-variant isoforms
  #     - build_canonical_paths_fine_c80s()      — explode to per-gene rows (fine)
  #   Level 3 (per-genome):
  #     - expand_canonical_paths_per_genome()    — one row per (canonical, contributing genome)
  
  path_min_genomes <- cfg_get(job_config, "path_min_genomes")
  truncation_cutoff <- cfg_get(job_config, "truncation_cutoff")

  # Collapse + canonicalise direction + joint components + within-component
  # re-orientation. This is the "consolidation" half of the step — every row
  # surviving here is a recurring operon backed by ≥ path_min_genomes strains.
  collapsed_paths <- collapse_paths_across_genomes(path_df)
  c_paths <- generate_canonical_path(collapsed_paths, path_min_genomes)
  jc_map <- compute_joint_components(c_paths, edge_types = unique(c_paths$path_type))
  c_paths <- decorate_paths_with_components(c_paths, jc_map)
  c_paths <- orient_paths_within_component(c_paths)

  # Bake the coarse uid: self-describing primary key for the L1 table.
  c_paths <- c_paths %>%
    mutate(uid = paste0("cmp", joint_component_ids, "-", path_type, "-", canonical_path_id, "-ng", n_genomes)) %>%
    select(uid, everything()) %>%
    relocate(n_genomes, .after = path_type)

  # L2 per-isoform aggregate.
  c_paths_fine <- expand_canonical_paths_to_fine(c_paths, collapsed_paths, path_df)

  # Per-gene exploders at coarse and fine resolution.
  c80s_coarse <- build_canonical_paths_c80s(c_paths, c80_variants_mapping, focal_c80_df, 
                                            cluster_80, jc_map, short_gene_prevalence) %>%
    relocate(neighbor_c80_coarse, .after = centroid_80_genome_counts)
  c80s_fine <- build_canonical_paths_fine_c80s(c_paths_fine, c80_variants_mapping, focal_c80_df, 
                                               cluster_80, jc_map, short_gene_prevalence) %>%
    relocate(neighbor_c80_coarse, .after = centroid_80_genome_counts)

  # Decorate. Truncation/fragmentation are fine-only by design — see
  # decorate_c80s_w_truncation docstring for why.
  c80s_coarse <- decorate_c80s_w_smallORFs(c80s_coarse, group_key = "uid")
  c80s_fine <- decorate_c80s_w_smallORFs(c80s_fine,   group_key = "uid_fine") %>%
    decorate_c80s_w_truncation(truncation_cutoff = truncation_cutoff)

  # L3 per-genome master table.
  c_paths_per_genome <- expand_canonical_paths_per_genome(c_paths, collapsed_paths, path_df, c_paths_fine)

  # Final column reorder + round.
  c80s_coarse <- c80s_coarse %>%
    select(uid, n_genomes, neighbor_c80_coarse, neighbor_gene_length, neighbor_c80_length_coarse,
           genome_prevalence, any_of(c("sample_prevalence", "cor_to_b", "beta")),
           neighbor_c80_coarse:dist_to_smallORFs, everything()) %>%
    mutate(genome_prevalence = round(genome_prevalence, 3))
  c80s_fine <- c80s_fine %>%
    select(uid_fine:neighbor_c80_fine, neighbor_c80_length_coarse, genome_prevalence,
           any_of(c("sample_prevalence", "cor_to_b", "beta")),
           is_focal, neighbor_c80_coarse:n_fragmented_c80s, everything()) %>%
    mutate(genome_prevalence = round(genome_prevalence, 3))

  # Persist.
  write.table(c_paths, get_target("canonical_paths"), sep = "\t", row.names = FALSE)
  write.table(c_paths_fine, get_target("canonical_paths_fine"), sep = "\t", row.names = FALSE)
  write.table(c80s_coarse, get_target("canonical_paths_c80s"), sep = "\t", row.names = FALSE)
  write.table(c80s_fine, get_target("canonical_paths_fine_c80s"),  sep = "\t", row.names = FALSE)
  write.table(c_paths_per_genome, get_target("canonical_paths_per_genome"), sep = "\t", row.names = FALSE)

  list(c_paths = c_paths, collapsed_paths = collapsed_paths, c80s_coarse = c80s_coarse)
}
