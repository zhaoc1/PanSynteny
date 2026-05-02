# Strain-Aware Operon Pipeline — User Guide

A high-level orientation to the pipeline: what it does, how to run it, and how to read the outputs. For per-step input / output / logic details, see [STEPS.md](STEPS.md). For internal architecture notes (file targets, function-level docstrings), read the source files directly — they are heavily commented.

---

## What this pipeline does

Given a set of trait-associated focal genes (centroid_80 clusters from MIDAS) for one species, the pipeline reconstructs the **operons** those genes live in, harmonizes them across the strains in the corpus, and surfaces the contiguous trait-associated blocks within each operon.

The output is three things:

1. **Canonical operons** at three granularity levels (coarse cluster path → length-variant isoforms → per-genome instances), each with attached trait statistics, small-ORF flags, and truncation/fragmentation flags.
2. **Trait-associated blocks**: the non-redundant runs of trait-correlated genes within those operons, ranked per locus.
3. **Per-genome attribution**: which strains carry which operon variant and which trait block.

The driver is [`pipeline.R`](pipeline.R). Steps 1–5 produce the analytical outputs (summaries, selection sets, BLAST gene lists). Step 6 (gggenes figures) is optional.

---

## How to run

The full workflow is two commands, both reading the same YAML:

```bash
# 1. Build focal_c80_df + enumerate any missing per-focal neighbor TSVs
Rscript prepare.R <config.yaml>
#    <if any TSVs are missing, prepare.R writes the list to gene_list.tsv;
#     materialise them externally under data_dir/neighbor_list/, then re-run
#     prepare.R until it reports "Ready to run pipeline.R">

# 2. Run the pipeline end to end
Rscript pipeline.R <config.yaml>
```

A working example config is [`example.yaml`](example.yaml).

**Step 0 — `prepare.R`.** Two responsibilities:

1. Read `corrected_genes` (RDS), filter to `species_id` and `trait`, apply `|score_col| >= inclusion_cutoff` for inclusion, mark rows with `|score_col| >= focal_cutoff` as `is_focal = TRUE`, and write the resulting `focal_c80_df` to [`gene_meta`](model.R#L34) as a TSV. `score_col` selects which trait-stat column drives the gating (typically `cor_to_b` or `beta`); both columns survive in the output regardless.
2. Walk every `is_focal == TRUE` centroid, check whether `<focal_c80>.tsv` already exists under [`neighbor_list`](model.R#L30), and write any missing centroids to [`gene_list`](model.R#L35) as a one-per-line list. If everything is present, `gene_list` is removed and a "Ready to run pipeline.R" message is printed.

Always overwrites — cheap to re-run after editing `corrected_genes`, the thresholds in YAML, or the contents of `neighbor_list/`.

**Steps 1–4 — `pipeline.R`.** Reads the TSV that `prepare.R` produced and runs the analytical pipeline. The driver consumes `focal_c80_df` as-is and does **not** apply any `|cor_to_b|` filter of its own — that decision is owned by `prepare.R`. If `gene_meta.tsv` is missing or any `is_focal` centroid still lacks a neighbor TSV, the driver aborts at startup with a pointer back to `prepare.R`.

**Re-run skipping.** On a re-run after partial completion, the driver skips work whose cached output exists: Step 1 skips re-extraction if `neighbor_groups_rds` exists, Step 2 skips if `path_df` exists. **To force a re-run of a step, delete its cache file.**

---

## Configuration (YAML)

All knobs live in one YAML file. Sections are flattened into a single `job_config` namespace at load time, so any key from any section is accessible via `cfg_get(job_config, "<key>")`. See [`config.R`](config.R) for the loader.

### `job` — required

```yaml
job:
  species_id:    "102321"             # MIDAS species id (numeric)
  trait:         "age_uni"            # column name in corrected_genes
  proj_dir:      "/path/to/results"   # output root (species_id appended)
```

### `data` — required

External data sources. Most users only ever change `data_dir`.

```yaml
data:
  midas_dir:    "..."   # MIDAS gene-by-sample matrices (per species)
  midasdb_dir:  "..."   # MIDAS database (clusters_80, genes_info)
  df_dir:       "..."   # DefenseFinder per-cluster results
  data_dir:     "..."   # neighbor_list inputs + corrected_genes
```

### `prepare` — Step 0 (`prepare.R`)

```yaml
prepare:
  score_col:        cor_to_b   # column in corrected_genes used to score genes; one of {cor_to_b, beta}
  inclusion_cutoff: 0.25       # min |score_col| for any gene to be retained in focal_c80_df
  focal_cutoff:     0.5        # min |score_col| for a gene to be marked is_focal = TRUE
```

`score_col` picks which trait-stat column the cutoffs are applied to (and which sign drives `focal_label`). `inclusion_cutoff` controls inclusion (rows below this are dropped entirely from `focal_c80_df`). `focal_cutoff` controls which rows drive Step 1 neighbor extraction. The constraint `focal_cutoff >= inclusion_cutoff` is enforced. Note that `cor_to_b ∈ [-1, 1]` while `beta` is unbounded — when switching `score_col`, scale the cutoffs accordingly.

### `neighbor` — Step 1 specifics

```yaml
neighbor:
  focal_min_genomes:       10     # min genomes a Step 1 neighborhood pattern must appear in
  focal_min_total_genomes: 30     # min total focal coverage; focals below this are dropped
  min_positions:           5      # min operon size
  upper_bound:             10     # symmetric position window around focal
  min_left_neighbors:      2      # flanking-coverage strict mode
  min_right_neighbors:     2
  use_strict:              ~      # ~ (NULL) = auto-detect; TRUE/FALSE to force
  min_group_proportion:    0.05   # min per-group genome proportion to retain a neighborhood pattern
  coverage_warn_threshold: 0.8    # warn if surviving groups cover less than this fraction of total
```

`focal_min_genomes` is the per-focal recurrence cut applied throughout Step 1 (pattern survival, flanking coverage, orientation grouping). `focal_min_total_genomes` is a separate, larger floor on each focal's total cross-genome support — focals seen in fewer total genomes are dropped before pattern extraction begins.

### `path` — Step 3 (canonical operon consolidation)

```yaml
path:
  path_min_genomes:  20    # min strains a canonical operon must be backed by
  truncation_cutoff: 0.8   # length / centroid_length below which a gene counts as truncated
```

`path_min_genomes` is the **per-canonical-operon** survival cut applied in [`generate_canonical_path()`](graph.R#L538). Conceptually distinct from Step 1's `focal_min_genomes`: Step 1 asks "is this neighborhood pattern recurring around a focal?" and Step 3 asks "does this canonical operon recur across enough strains?" — they need not have the same value, though they typically do.

### `blocks` — Step 4 block extraction (`blocks.R`)

```yaml
blocks:
  allow_gaps:  2   # max non-hit positions allowed inside one block (aggregate_blocks)
  min_overlap: 1   # min LCS length to call two blocks "overlap" (rank_block_representatives)
  min_shared:  2   # diagnostic threshold for substring-overlap pairs (diagnose_rep_overlaps)
```

`min_shared` is diagnostic-only — it does not change `representatives.tsv` or `rep.tsv`, only the `Rep overlap diagnostic` log line.

### `parse` — Step 5 + Step 6 (summaries, sampling, BLAST gene lists, gggenes figures)

Consumed by Steps 5 and 6 only. Not consumed by Steps 1–4. (`fill_modes` is read by Step 6's `run_step6_figures`; the rest by Step 5's `run_step5_parse`.)

```yaml
parse:
  fine_coverage_ratio: 0.25  # n_fine_genomes >= ceiling(path_min_genomes * ratio) to survive
  seed:                616   # RNG seed for the per-fine-isoform exemplar-genome draw
  fill_modes:                # one PDF per mode; subset of:
    - beta
    - sample_prevalence
    - cor_to_b
    - fill_gene
```

`fine_coverage_ratio` sets the fraction of the path-level genome floor each fine isoform must reach. `ceiling` semantics are preserved — only the multiplier is tunable. Note the ratio is applied to `path_min_genomes` (Step 3's gate), not to `focal_min_genomes`. Historical default 0.5 ("half-coverage rule"); the example currently uses 0.25 for permissive fine-isoform retention.

### `plot` — Step 1 + Step 6 figure layout

```yaml
plot:
  gene_padding_bp: 100   # bp gap between adjacent genes in gggenes layouts
```

Used by both the Step 1 diagnostic figures ([`extract_gene_neighbor_patterns`](neighbor.R#L323)) and the Step 6 publication figures ([`.layout_operon_tracks`](plot.R#L245)) so spacing is consistent across all gggenes outputs.

---

## Pipeline overview

| Step | What it does | Driver lines | Helper file(s) | Cache / output |
|---|---|---|---|---|
| **0** | Build `focal_c80_df` from `corrected_genes`, apply `\|score_col\|` thresholds (`inclusion_cutoff` and `focal_cutoff`), write `gene_meta.tsv`, and enumerate any missing per-focal neighbor TSVs into `gene_list.tsv`. Separate script: [`prepare.R`](prepare.R). | — | `config.R`, `model.R` | `gene_meta.tsv`, `gene_list.tsv` |
| **Setup** | Load `cluster_80`, `gene_to_c80`, and the `focal_c80_df` (with `is_focal`) produced by Step 0. Re-check that every focal has a neighbor TSV; abort if not. | 11–79 | `config.R`, `model.R` | — |
| **1** | Per-focal neighborhood extraction → cross-genome assembly → small-ORF + length-variant labels. Orchestrated by [`run_step1_neighbor_extraction`](neighbor.R). | 82–92 | `neighbor.R`, `midas.R` | `step2_neighbors/neighbor_groups.RDS` (cache) |
| **2** | Per-genome operon graphs → maximal paths. Orchestrated by [`run_step2_path_stitching`](graph.R). | 95–100 | `graph.R` | `step3_path/path_df.rds`, `step3_path/esupport_df.rds` |
| **3** | Cross-genome consolidation → three granularity levels with trait stats and structural flags. Orchestrated by [`run_step3_consolidation`](path.R). | 103–128 | `graph.R`, `path.R`, `parse.R` | `step3_path/canonical_paths*.tsv`, `step3_path/canonical_paths_*c80s.tsv` (5 TSVs) |
| **4** | Trait-associated block extraction + non-redundant representative ranking + per-genome attribution. Orchestrated by [`run_step4_block_extraction`](blocks.R). | 131–139 | `blocks.R` | `step4_block/representative_path.tsv`, `step4_block/rep.tsv`, `step4_block/rep_heatmap.pdf` |
| **5** | Summaries, fine-coverage selection, exemplar-genome sampling, BLAST gene lists. Orchestrated by [`run_step5_parse`](parse.R). | 142–153 | `parse.R` | `step5_parse/*` |
| **6** | gggenes figures: global + per-component PDFs for each fill mode. Orchestrated by [`run_step6_figures`](plot.R). | 156–163 | `plot.R` | `step6_figures/*` |

For each numbered step, [STEPS.md](STEPS.md) gives a complete input / output / logic / caveats writeup.

---

## Output reference

All outputs land under `<proj_dir>/<species_id>/`. Key files:

### Canonical operon tables (Step 3 — the analytical core)

The same operons are emitted at three granularity levels. Every level carries trait statistics, joint-component membership, and structural decorations.

- **Level 1 (coarse)** — [`step3_path/canonical_paths_coarse.tsv`](model.R#L43) (one row per operon) and [`step3_path/canonical_paths_c80s.tsv`](model.R#L46) (one row per gene-in-operon). **Primary key:** `uid = "cmp{joint_component_ids}-{path_type}-{canonical_path_id}-ng{n_genomes}"` (e.g. `cmp3-anchor-cp_42-ng18`).
- **Level 2 (per-isoform)** — [`step3_path/canonical_paths_fine.tsv`](model.R#L44) and [`step3_path/canonical_paths_fine_c80s.tsv`](model.R#L47). **Primary key:** `uid_fine = "{uid}-iso{rank}-ngf{n_fine_genomes}"`. Strip `-iso\d+-ngf\d+$` to recover the parent `uid`.
- **Level 3 (per-genome)** — [`step3_path/canonical_paths_per_genome.tsv`](model.R#L45). One row per `(canonical, contributing genome)` with raw and canonical-aligned `gene_path` strings.

**When to use which level**

| Question | Level |
|---|---|
| What recurring operons are trait-associated? | L1 coarse — group by `uid`. |
| Are there length variants worth distinguishing (truncation, fragmentation, tandem split-genes)? | L2 fine — `is_truncated`, `is_fragmented` are populated here only. |
| Which strains carry operon X? Which gene IDs do I BLAST? | L3 per-genome — `neighbor_genome` + `gene_path`. |

### Block tables (Step 4 — the trait answer)

- [`step4_block/representative_path.tsv`](model.R#L49) — non-redundant trait-associated blocks. **Primary key:** `block_uid = "cmp{component}-{type}-rank{rep_rank}-nge{block_n_genes}"`.
- [`step4_block/rep.tsv`](model.R#L50) — per-genome attribution for the reps. One row per `(block_uid, canonical_uid, neighbor_genome, left_orig, right_orig)`.
- `step4_block/rep_heatmap.pdf` — block × genome presence/absence.

---

## Three c80 columns: `c80`, `neighbor_c80_coarse`, `neighbor_c80_fine`

Easy to confuse. They differ in **scope** and **resolution**.

| Column | Resolution | Scope | Where it appears |
|---|---|---|---|
| `c80` | coarse | the operon's own focal cluster id | L1 c80s table only |
| `neighbor_c80_coarse` | coarse | a neighbor's MIDAS cluster id (a real centroid_80, or a synthetic `_<focal>-<type>_<rank>` for short ORFs) | everywhere |
| `neighbor_c80_fine` | length-variant-aware | same as `neighbor_c80_coarse` but with `_<rank>` suffix when the cluster has multiple observed lengths | L2 fine c80s table; produced by Step 1 ([`compute_c80_variants()`](midas.R#L204)) |

Two practical rules:

- **For coarse grouping** (collapse all length variants of a cluster) → use `neighbor_c80_coarse`.
- **For length-sensitive analyses** (truncation, fragmentation, isoform identity) → use `neighbor_c80_fine` and the L2 fine table.

In the L1 c80s table, `neighbor_c80_coarse == c80` after the [`build_canonical_paths_c80s()`](path.R#L264) join — they're kept as parallel columns for join-compatibility, not because they hold different information.

---

## Path-direction conventions across steps

This is the most subtle thing in the pipeline. **Each step canonicalizes direction in a different frame, and what survives downstream is the most-recent canonicalization.**

| Step | Direction frame | What it discards |
|---|---|---|
| 1 | Focal-relative (left vs right of the focal gene) | Whichever direction `right_anchor < left_anchor` triggers a flip |
| 2 | Chromosomal (per genome, sorted by `neighbor_gene_start`) | **Step 1's focal-relative orientation** — overwritten |
| 3 | Lexicographic (`normalize_path` picks min of forward / reverse) → then re-aligned within joint component to match the longest reference | Step 2's chromosomal orientation; absolute biological direction |

**Implication.** In the L1 / L2 / L3 outputs, "left-to-right" means "consistent within a joint component", not "5'→3' of any real chromosome". To recover real direction, look at L3's per-genome `gene_path` and walk back to chromosomal coordinates via `path_df` (Step 2 output) or the original neighbor TSVs.

A separate gap that remains worth tracking: Step 1's per-focal orientation is lost in Step 2; Step 3 canonicalizes path direction *before* cross-genome collapse. If you need focal-relative direction preserved across the whole pipeline, that path doesn't exist today — it would require carrying Step 1's `orientation` column through Step 2's chromosomal re-derivation.

---

## File organization

```
pipeline/
├── prepare.R                   # Step 0: build focal_c80_df from corrected_genes
├── pipeline.R           # the driver — reads it top-to-bottom
├── config.R                    # YAML loader + cfg_get
├── model.R                     # target_layout + get_target (file-path resolver)
├── neighbor.R                  # Step 1: per-focal neighborhood pipeline
├── midas.R                     # Step 1: small-ORF + length-variant labels
├── graph.R                     # Step 2 + Step 3: graph building, joint components, orientation
├── path.R                      # Step 3: canonical → fine → per-genome expansions
├── blocks.R                    # Step 4: hit blocks + reps + per-genome attribution
├── parse.R                     # Step 5 orchestrator + Step 3 c80s decorators + Step 6 plot data-prep helpers
├── plot.R                      # Step 6: global + per-component gggenes plots (and Step 1 diagnostic plots)
├── example.yaml         # example config
├── USER_GUIDE.md               # this file
├── STEPS.md                    # per-step deep dive
└── PIPELINE.md                 # earlier architecture notes (kept for reference)
```

---

## Tunables — quick reference

| Knob | Section | Default | Where it bites | Effect |
| --- | --- | --- | --- | --- |
| `score_col` | `prepare` | `cor_to_b` | Step 0 (`prepare.R`) | Trait-stat column the cutoffs apply to; one of `{cor_to_b, beta}`. Also drives the `focal_label` sign. |
| `inclusion_cutoff` | `prepare` | 0.25 | Step 0 (`prepare.R`) | Minimum `\|score_col\|` to retain a gene in `focal_c80_df` at all. |
| `focal_cutoff` | `prepare` | 0.5 | Step 0 (`prepare.R`) | Minimum `\|score_col\|` for a row to be flagged `is_focal = TRUE` and drive Step 1 neighbor extraction. |
| `focal_min_genomes` | `neighbor` | 10 | Step 1 (multiple gates: pattern survival, flanking coverage, orientation grouping) | Per-focal recurrence cut. Lower = more rare neighborhoods survive. |
| `focal_min_total_genomes` | `neighbor` | 30 | Step 1 (`parse_gene_neighbor`) | Min total focal coverage; focals with thinner support are dropped before pattern extraction. |
| `min_positions` | `neighbor` | 5 | Step 1 (`compute_relative_positions`, `parse_gene_neighbor`) | Min operon size. Lower = include shorter neighborhoods. |
| `upper_bound` | `neighbor` | 10 | Step 1 (`compute_relative_positions`) | Position window around focal (±). Larger = wider neighborhoods. |
| `min_left_neighbors` / `min_right_neighbors` | `neighbor` | 2 / 2 | Step 1 (`filter_by_flanking_coverage`) | Strict-mode flanking requirement. |
| `use_strict` | `neighbor` | `~` (auto) | Step 1 (`filter_by_flanking_coverage`) | `TRUE` forces strict; `FALSE` forces relaxed; `~` auto-selects based on `focal_min_genomes`. |
| `min_group_proportion` | `neighbor` | 0.05 | Step 1 (`parse_gene_neighbor`, second-pass filter) | A neighborhood-pattern group survives if it meets `focal_min_genomes` OR has at least this fraction of the focal's total genome support. Lower = more permissive. |
| `coverage_warn_threshold` | `neighbor` | 0.8 | Step 1 (`parse_gene_neighbor`, second-pass filter) | If surviving pattern groups cover less than this fraction of the focal's total support, emit a warning. **Diagnostic-only**; does not gate the output. |
| `path_min_genomes` | `path` | 20 | Step 3 (`generate_canonical_path` survival cut), Step 5 (driver, via `fine_coverage_ratio`) | Per-canonical-operon recurrence cut. Distinct from `focal_min_genomes`. |
| `truncation_cutoff` | `path` | 0.8 | Step 3 (`decorate_c80s_w_truncation`) | A gene shorter than this fraction of its centroid length is `is_truncated`. |
| `allow_gaps` | `blocks` | 2 | Step 4 (`aggregate_blocks`) | Max non-hit positions allowed inside one focal block. Higher = more permissive merging. |
| `min_overlap` | `blocks` | 1 | Step 4 (`rank_block_representatives`, `get_relation`) | Min LCS length to call two blocks "overlap" (vs. "disjoint") when ranking reps. |
| `min_shared` | `blocks` | 2 | Step 4 (`diagnose_rep_overlaps`) | **Diagnostic-only** threshold for substring-overlap pairs. Does not change `representatives.tsv`. |
| `fine_coverage_ratio` | `parse` | 0.25 | Step 5 (driver) | Fine isoforms survive if `n_fine_genomes >= ceiling(path_min_genomes * ratio)`. Historical default 0.5 (half-coverage rule); example currently uses 0.25 for permissive retention. |
| `seed` | `parse` | 616 | Step 5 (`sample_genome_from_fine_paths`) | RNG seed for the per-fine-isoform exemplar-genome draw. |
| `fill_modes` | `parse` | all four | Step 6 (`run_step6_figures`) | Which fill modes to render. Subset of `{beta, sample_prevalence, cor_to_b, fill_gene}`. |
| `gene_padding_bp` | `plot` | 100 | Step 1 (`extract_gene_neighbor_patterns`), Step 6 (`.layout_operon_tracks`) | bp gap between adjacent genes in gggenes layouts. Same value used by both diagnostic and publication figures. |

---

## Known issues

These are documented here so you don't trip on them.

1. **Mirror-block reps survive in Step 4.** [`is_contig_subseq`](blocks.R#L391) is forward-direction only; a block and its exact reverse end up as two separate reps. The diagnostic in `diagnose_rep_overlaps` will catch this if it happens.
2. **Step 1 orientation is not preserved into Steps 2/3.** Step 2 re-derives chromosomal order; Step 3 canonicalizes lexicographically (with synthetic small-ORF tokens stripped from the decision; see `clean_for_orientation` in `graph.R`). Within-component direction consistency is the strongest guarantee you get on the output side.

