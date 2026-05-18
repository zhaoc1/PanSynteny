# c80 column glossary

Decorated c80s tables (`canonical_paths_c80s`, `canonical_paths_fine_c80s`, after
`decorate_c80s_w_smallORFs`) carry three c80-related columns that look similar
but answer different questions about the same row.

## The three columns

### `neighbor_c80_coarse` — coarse cluster id

Length variants of the same cluster collapse together.

- For **microslam genes**: the database `centroid_80` cluster id (e.g., `"3045"`).
- For **short ORFs**: the synthetic per-focal label
  `_<focal_c80>-<gene_type>_<rank>` (e.g., `"_3045-CDS_1"`). Always starts with
  `_` — that's how `is_smallORF` is detected.

Synthetic small-ORF labels are produced by `compute_short_gene_prevalence()` in `midas.R` (Step 1); microslam c80s are inherited as-is from the input neighbor TSVs.

### `neighbor_c80_fine` — cluster + length-rank

Same coarse cluster, distinguished by length variant.

- For **microslam clusters with multiple length variants**: `<c80>_<rank>`,
  ranked first = shortest (e.g., `"3045_1"` for the 600 bp variant, `"3045_2"`
  for the 900 bp variant of the same cluster).
- For **single-length clusters and short ORFs**: equal to `neighbor_c80_coarse`.

Created by `compute_c80_variants()` in `midas.R`.

### `centroid_80` — focal context

Decoded focal scope.

- For **microslam rows**: equal to `neighbor_c80_coarse`. No new information.
- For **short-ORF rows**: the **focal_c80** decoded from the synthetic label
  (everything between the leading `_` and the last `-`). This is the focal
  cluster whose neighborhood contains this small ORF — *not* the small ORF's
  own identity.

Created by `decorate_c80s_w_smallORFs()` in `parse.R`.

## Worked examples

### Microslam cluster `3045` with two length variants

| neighbor_c80_coarse | neighbor_c80_fine | centroid_80 | neighbor_gene_length |
| --- | --- | --- | --- |
| `3045` | `3045_1` | `3045` | 600 |
| `3045` | `3045_2` | `3045` | 900 |

Three columns carry only two distinct values: the cluster id (top + bottom) and
the length variant (middle).

### Short ORF in focal `3045`'s neighborhood

| neighbor_c80_coarse | neighbor_c80_fine | centroid_80 | smallORF_type | is_smallORF |
| --- | --- | --- | --- | --- |
| `_3045-CDS_1` | `_3045-CDS_1` | `3045` | `CDS` | TRUE |

`centroid_80 ≠ neighbor_c80_coarse` here. `centroid_80 = 3045` is the **focal** that
this small ORF was derived from — the same value the focal microslam row above
carries. That shared value is what makes per-focal grouping work.

## Quick guide — which one to group on

| Question | Group on |
| --- | --- |
| All rows in this coarse cluster (length variants collapsed) | `neighbor_c80_coarse` |
| Did the same cluster show up at multiple lengths in one isoform? (split-gene check — see `decorate_c80s_w_truncation`) | `neighbor_c80_fine` |
| This focal's row plus all small ORFs derived from its neighborhood | `centroid_80` |

## Granularity hierarchy

```text
centroid_80           (focal-scoped: focal + its small ORFs share this)
   ↓
neighbor_c80_coarse   (cluster-scoped: length variants collapsed)
   ↓
neighbor_c80_fine     (variant-scoped: cluster + length-rank)
```

For **microslam genes**, the top two are identical — only `neighbor_c80_fine`
adds info (the length-rank suffix).

For **short ORFs**, the top column resolves the focal context that the
synthetic label encodes; the bottom two are equal to each other (no length
variants for short ORFs).

## Why three columns instead of one

Each is the right grouping key for a different question, and pre-computing them
all means downstream consumers don't have to parse synthetic labels. In
particular, `centroid_80` is what enables `dist_to_smallORFs` per-focal scoping
in `decorate_c80s_w_smallORFs`: a small ORF derived from focal A does not
contribute to focal B's distance even when both focals appear in the same
canonical path.

## Truncation & fragmentation columns (fine only)

Added by `decorate_c80s_w_truncation()` in `parse.R`. Two for truncation
(length-based), three for fragmentation (label-based).

| Column | Type | Grain | Meaning |
| --- | --- | --- | --- |
| `is_truncated` | logical | per row | This gene shorter than `cutoff × centroid_80_gene_length`? |
| `truncate_ratio` | double | per row | `neighbor_gene_length / centroid_80_gene_length`, floor-truncated to 3 decimals; NA for short ORFs |
| `n_truncated` | int | per isoform (broadcast) | How many genes in this operon are truncated |
| `is_fragmented` | logical | per row | This row's coarse cluster appears at ≥ 2 length variants in this isoform |
| `fragmented_c80s` | char | per row | Coarse cluster id when `is_fragmented`, else NA |
| `n_fragmented_c80s` | int | per isoform (broadcast) | How many *distinct* coarse clusters in this operon are fragmented |

Truncation and fragmentation are independent — a c80 can be neither, one, or
both. Fragmentation uses `n_distinct` (not `sum`) on `fragmented_c80s`, so a
single fragmented cluster split across multiple rows counts once per isoform.

### How they differ

Both flags target the same biological concern — this c80 may not represent a
complete copy of the underlying gene — via different evidence:

- **Truncation** is *length-based*. Compares the observed `neighbor_gene_length`
  to the database `centroid_80_gene_length`. A row is truncated if it's much
  shorter than the centroid (default cutoff: 80%).
- **Fragmentation** is *label-based*. Looks for the same coarse cluster
  appearing under multiple length-variant labels in one isoform — the
  "split-gene" signature where a database cluster is realized at multiple
  positions of the operon.

The two are independent. Truth table:

| Situation | `is_truncated` | `is_fragmented` |
| --- | --- | --- |
| Full-length, single variant | FALSE | FALSE |
| Short for its centroid, no other variants of this cluster in the operon | TRUE | FALSE |
| Multiple variants of the cluster present, this row meets the length cutoff | FALSE | TRUE |
| Multiple variants present *and* this row is short | TRUE | TRUE |

**Fine-only by design.** The coarse table's `neighbor_gene_length` is
max-over-isoforms (per `build_canonical_paths_c80s`), so a coarse-table
truncation check would mean "even the longest observed isoform is shorter
than cutoff" — not the per-isoform comparison the column name implies.
Fragmentation also doesn't apply at coarse — there's no per-isoform partition
to count length variants within.

---

## Pipeline overview — draft notes

> **Status:** holding area for future pipeline documentation. Comments below
> are extracted verbatim (with light editing) from `pipeline.R` so they
> don't get lost when the orchestration script gets refactored. Reorganize into
> proper docs later.

### Step 1 — per-focal neighborhood extraction

For each focal gene, extract its neighborhood across every genome that carries
it, orient each observed copy on a common axis, select the dominant operon-size
patterns, and write per-genome TSVs. All per-focal TSVs are then loaded back
into a single `gene_neighbors` table and enriched with both cluster-label
resolutions:

- `neighbor_c80_coarse` — coarse (cluster-level, length invariant).
- `neighbor_c80_fine` — fine (distinguishes length variants).

**Logic:** the orchestrator is `run_step1_neighbor_extraction()` in `neighbor.R`, which calls `extract_and_write_per_focal_neighbors()` then `load_neighbors_across_genomes()`.

**Cache:** delete `neighbor_groups_rds` to force regeneration.

After loading, `assign_c80_to_short_genes()` in `midas.R` resolves short-ORF
labels and produces the `short_gene_prevalence` lookup; `compute_c80_variants()`
adds the per-length-variant `neighbor_c80_fine` column.

### Step 2 — assemble cross-focal paths per genome

For each genome, build a directed graph of neighbor-gene adjacencies, aggregate
per-edge support across overlapping focal neighborhoods, and enumerate every
maximal source-to-sink path. Produces:

- `path_df` — one row per operon chain.
  - `c80_path_coarse` — coarse rendering.
  - `c80_path_fine` — length-variant rendering.
- `esupport_df` — one row per directed edge with support counts.

**Logic:** `run_step2_path_stitching()` in `graph.R` (orchestrator wrapping `stitch_paths_across_focal_genes()`).

After loading, a per-row `path_genome_comp` key is constructed as
`paste(path_genome, type, path_component_id, sep="||")` and used as the join
key into per-genome attribution downstream.

### Step 3 — three tables at increasing granularity

Step 3 produces three tables that share canonical-path identity but differ in
aggregation grain (coarsest → finest):

#### Level 1 — `canonical_paths`

One row per canonical locus. Coarse c80 resolution. This is the **identity
table**; cross-genome collapse + direction canonicalization are finished here.

Built by `collapse_paths_across_genomes()` → `generate_canonical_path()` →
`compute_joint_components()` → `decorate_paths_with_components()` →
`orient_paths_within_component()`. The coarse `uid` is then baked once:

```r
uid = paste0("cmp", joint_component_ids, "-", type, "-", canonical_path_id, "-ng", n_genomes)
```

This makes `uid` the single stable key downstream code joins/filters on,
instead of rebuilding the paste at every site.

#### Level 2 — `canonical_paths_fine`

One row per `(canonical_path_id, isoform variant)`. Derived from `path_df`'s
length-variant c80 labels; preserves the per-isoform genome counts that the
coarse collapse erased. All Level 1 columns are carried alongside for context
(except coarse `neighbor_genomes`, which the per-isoform `fine_neighbor_genomes`
unions back to).

Hierarchical fine uid: `uid_fine = paste0(uid, "-iso", isoform_rank, "-ngf", n_fine_genomes)`.
Stripping `-iso\d+-ngf\d+$` recovers the coarse parent uid.

#### Level 3 — `canonical_paths_per_genome`

One row per `(canonical_path_id, contributing per-genome observation)`. The
table answers per-genome attribution questions — *for any canonical path /
isoform, which genomes contributed and what was the genome's actual gene-id
chain at that locus?* This is the only table where **per-genome gene-id
resolution** lives.

**Output schema (lean, Option A):**

Front block (relocated):

```text
uid_fine, uid,
canonical_path_id, collapsed_path_id, fine_canonical_id, per_genome_path_w_ids,
neighbor_genome,
gene_path, gene_path_canonical,
needs_flip
```

Tail (via `everything()`): `isoform_rank`, `path_type`, `joint_component_ids`,
`n_genomes`, `n_fine_genomes`.

**Intentionally not in the output** — the four c80-resolution path strings
(`c80_path_coarse`, `c80_path_coarse_canonical`, `c80_path_fine`,
`c80_path_fine_canonical`) are dropped because each is 1:1 with an ID column
the table already carries (`collapsed_path_id`, `canonical_path_id`,
`per_genome_path_w_ids`, `uid_fine`/`fine_canonical_id`) or recoverable from
`gene_path` via the gene-id → c80-label mapping. To re-attach any of them:

```r
per_genome %>% left_join(c_paths_fine %>% select(uid_fine, c80_path_fine_canonical),
                          by = "uid_fine")
```

**`needs_flip` semantics.** Per-row boolean indicating whether this row's raw
observation was reverse-oriented relative to canonical. Computed once at the
**(canonical × collapsed_path)** grain in `explode_canonical_into_collapsed_paths`
as `(c80_path_string != canonical_path_coarse)`, then inherited by every
per-genome row under that collapsed group. Because the upstream canonicalization
(`normalize_path` + `orient_paths_within_component`) is orientation-only, a single
boolean drives the row-wise flip of `gene_path` → `gene_path_canonical`. Note
the grain: all rows sharing a `collapsed_path_id` share `needs_flip`; within
one isoform (`uid_fine`) or one canonical (`uid`), `needs_flip` can be mixed
because mirror-image collapsed groups roll up into the same isoform/canonical.

This is also why `needs_flip` cannot be promoted to `c_paths_fine` or `c_paths`
— neither has a single value at those grains.

#### Per-c80 anchor tables (also Step 3)

Two parallel "exploded" tables, one row per `(canonical path × c80 position)`:

- `canonical_paths_c80s` — coarse anchor. Uses max-over-isoforms gene length.
  Built by `build_canonical_paths_c80s()`. Consumed by `blocks.R`.
- `canonical_paths_fine_c80s` — fine anchor. Uses exact per-isoform gene length.
  Built by `build_canonical_paths_fine_c80s()`. Consumed by truncation /
  fragmentation analysis (which is fine-only by design).

Both are written to TSV at the end of Step 3.

### Step 4 — summaries, selection, exemplar sampling, BLAST gene lists

Orchestrator: `run_step4_parse()` in `parse.R`. Builds `coarse_summary` (one row per `uid`) and `fine_summary` (one row per `uid_fine`), applies the fine-coverage isoform survival filter (`n_fine_genomes >= ceiling(path_min_genomes * fine_coverage_ratio)`), samples one exemplar genome per surviving fine isoform, enriches the per-gene long-format result with per-isoform context, and writes per-`(uid_fine, neighbor_genome)` gene-id TSVs for an external BLAST workflow. See [STEPS.md §STEP 4](STEPS.md) for the full per-stage breakdown.

### Step 5 — gggenes figures (writes under `step5_figures/`)

Orchestrator: `run_step5_figures()` in `plot.R`. Renders four plotters at global and per-component scope:

1. `plot_coarse_operons()` — global coarse, one PDF per `fill_by` mode.
2. `plot_fine_operons()` — global fine.
3. `plot_coarse_operons_by_component()` — one PDF per `(joint_component_id, fill_by)` under `step5_figures/02_by_component_coarse/`.
4. `plot_fine_operons_by_component()` — fine analogue under `step5_figures/03_by_component_fine/`.

Per-c80 glyphs (precedence top→bottom): `U` / `D` (focal up/down by `beta` sign), `F` (fragmented, fine only), `T` (truncated, fine only), `s` (small ORF). `fill_modes` is read from YAML `parse.fill_modes` (subset of `{beta, sample_prevalence, cor_to_b, fill_gene}`).

### Step 6 — trait-associated block extraction (writes under `step6_blocks/`, gated by `blocks.skip`)

Orchestrator: `run_step6_blocks()` in `blocks.R`. From `canonical_paths_c80s`, identify runs of consecutive focal c80s ("blocks") and rank representative blocks per joint component:

- `aggregate_blocks(canonical_paths_c80s, allow_gaps = 2)` —
  fold per-c80 hits into runs (allowing 2-position gaps).
- `rank_block_representatives(block_agg, min_overlap = 1)` — pick representative
  blocks per `(joint_component_id, path_type)` group. Outputs `representatives`.
- `map_representatives_to_genomes()` — link each representative back to the
  genomes it appears in. Outputs `rep_slim` (block_uid × neighbor_genome).
- `diagnose_rep_overlaps()` — diagnostic: how often surviving reps within a
  component share a substring but neither contains the other.

If multiple representative blocks exist, a heatmap of block-by-genome membership is rendered to `step6_blocks/rep_heatmap.pdf` (via `pheatmap`). Skip the whole step by setting `blocks.skip: true`.

### c80s decoration (inside Step 3)

`run_step3_consolidation()` calls two annotators after building the per-gene tables:

- `decorate_c80s_w_smallORFs()` — adds `is_smallORF`, `centroid_80`,
  `smallORF_type`, `n_smallORFs`, `n_focal`, `dist_to_smallORFs`. Applied to **both**
  coarse (`group_key = "uid"`) and fine (`group_key = "uid_fine"`).
- `decorate_c80s_w_truncation()` — adds `is_truncated`, `truncate_ratio`,
  `n_truncated`, `is_fragmented`, `fragmented_c80s`, `n_fragmented_c80s`.
  Applied to **fine only** (the coarse table's `neighbor_gene_length` is
  max-over-isoforms, making a truncation check meaningless there).

Brief inline gloss: *truncation* = gene shorter than the centroid; *fragment* =
the same coarse cluster appearing under multiple length variants within one
isoform.

### Configuration knobs

| Key | Used in | Purpose |
| --- | --- | --- |
| `species_id` | Step 0a onward | Selects which MIDAS pangenome to load and which species the catalog is scoped to |
| `focal_meta` | Step 0 (`prepare.R`) | Absolute path to the user-provided focal-gene TSV |
| `midasdb_dir` | Step 0a (`build_genome_catalog`), Step 1 onward (via `load_c80_tables`) | MIDAS reference DB root |
| `data_dir` | Step 0a + Step 0 | Catalog + neighbor TSVs land under `<data_dir>/<species_id>/` |
| `n_genes` | Step 0 (`get_neighbor.sh`) | Max neighbours per side around a focal |
| `sources` (list) | Step 0a (`build_genome_catalog`) | Declares each source's membership file + `.genes` location |
| `focal_min_genomes` | Step 1 | Minimum genome support for a Step 1 neighborhood pattern |
| `path_min_genomes` | Step 3, Step 4 | Minimum genome support to retain a canonical path; also the floor that `fine_coverage_ratio` is applied to |
| `truncation_cutoff` | Step 3 (`decorate_c80s_w_truncation`) | Length-ratio threshold below which a fine c80 is flagged truncated |
| `fine_coverage_ratio` | Step 4 | Fine-isoform survival fraction of `path_min_genomes` |
| `seed` | Step 4 | RNG seed for the per-isoform exemplar-genome draw |
| `fill_modes` | Step 5 | Which gggenes fill encodings to render (subset of `{beta, sample_prevalence, cor_to_b, fill_gene}`) |
| `midas_dir` (YAML `mwas:`) | **legacy — parked** | Feeds the MWAS block in `model.R`; not read by the current pipeline |

### Inclusion thresholds

Owned by Step 0 (`prepare.R`) and YAML-controlled — not hardcoded in the driver. See [USER_GUIDE.md §`prepare`](USER_GUIDE.md) for the full schema. Two gates:

- `inclusion_cutoff` — minimum `|score_col|` for a gene to be retained in `focal_c80_df`.
- `focal_cutoff` — minimum `|score_col|` to flag `is_focal = TRUE` (drives Step 1 neighbor extraction).

`score_col` selects which trait-stat column the cutoffs apply to (default `cor_to_b`; `beta` is also valid). Direction labels: `focal_label = ifelse(score_col > 0, "pos", "neg")`, computed in `prepare.R`.
