# Per-Step Reference

Detailed input / output / logic for every stage of `pipeline.R`. For an executive overview and the run command, see [USER_GUIDE.md](USER_GUIDE.md).

All file paths below resolve through [`get_target()`](model.R#L66) against the YAML config; resolved targets are listed in [`target_layout()`](model.R#L13).

Stage order in the driver:

1. Pre-Step setup (focal selection, reference tables)
2. **Step 1** — per-focal neighborhoods → small-ORF + length-variant labels
3. **Step 2** — per-genome operon graphs → maximal paths
4. **Step 3** — cross-genome consolidation → three granularity levels
5. **Step 4** — focal block extraction + representative ranking
6. **Step 5** — summaries, fine-coverage selection, exemplar sampling, BLAST gene lists
7. **Step 6** — gggenes figures (global + per-component)

Each numbered step is a one-line orchestrator call in [`pipeline.R`](pipeline.R); the body lives in the topic-specific helper file (`run_stepN_*` in `neighbor.R` / `graph.R` / `path.R` / `blocks.R` / `parse.R` / `plot.R`).

---

## Pre-Step driver setup

Before Step 1 runs, the driver loads three reference tables and filters them down to the trait-associated focal set. These are inputs to Step 1, not part of any numbered step, but understanding them anchors everything downstream.

| Symbol | Source | Role |
|---|---|---|
| `cluster_80` | `clusters_80_updated` (TSV) | One row per coarse centroid_80 cluster: cluster id, cluster gene length, genome prevalence. |
| `gene_to_c80` | `genes_info` (TSV) joined to `cluster_80` | Maps each individual `gene_id` to its parent `centroid_80`. |
| `all_genes` | `corrected_genes` (RDS) filtered by `species_id` and `trait_model` | Per-gene trait statistics (`cor_to_b`, `beta`, `pvalue_noadj`). |
| `focal_c80_df` | Built by Step 0 ([`prepare.R`](prepare.R)) and read from [`gene_meta`](model.R#L34) | One row per candidate focal cluster, with an `is_focal` boolean column set by `prepare.R` from `\|score_col\| >= focal_cutoff` (defaults: `score_col = cor_to_b`, `focal_cutoff = 0.5`). The pipeline does not re-filter. Required columns: `focal_c80`, `sample_prevalence`, `cor_to_b`, `focal_label`, `beta`, `trait`, `genome_counts`, `is_focal`. |

The driver errors out at startup in two cases, both with a pointer to `Rscript prepare.R <config.yaml>`: (1) `gene_meta` not present on disk; (2) one or more `is_focal == TRUE` centroids still missing their per-focal neighbor TSVs under [`neighbor_list`](model.R#L30). The discovery + enumeration of missing TSVs lives in [`prepare.R`](prepare.R), which writes the list to [`gene_list`](model.R#L35) for the external materialisation job; this driver only re-checks the contract.

---

## STEP 1 — Per-focal neighborhoods → small-ORF + length-variant labels

**Driver section:** [pipeline.R:82–92](pipeline.R#L82-L92) (one-line orchestrator call to [`run_step1_neighbor_extraction`](neighbor.R#L791))
**Helper files:** [neighbor.R](neighbor.R), [midas.R](midas.R)

### Input

- `focal_c80_df` — derived from `gene_meta`. One row per focal `centroid_80`, with `c80_label` ∈ {`pos`, `neg`} from the sign of `cor_to_b`. Renames `gene_id → c80`.
- Per-focal neighbor TSVs on disk under [`neighbor_list`](model.R#L30) (one `<focal_c80>.tsv` per focal). Each TSV has 7 columns: `gene_member`, `neighbor_gene_id`, `neighbor_contig_id`, `neighbor_gene_start`, `neighbor_gene_end`, `neighbor_gene_strand`, `neighbor_gene_type`. These are produced by an external (non-R) preprocessing step; the driver only reads them.
- `gene_to_c80` — to map `neighbor_gene_id → neighbor_c80_coarse` (and pull cluster-level gene length / prevalence).

### Output

- `gene_neighbors` (in-memory, cached to RDS at [`neighbor_groups_rds`](model.R#L38)). One row per (focal, genome, neighbor position), with both `neighbor_c80_coarse` (coarse) and `neighbor_c80_fine` (length-variant-aware) columns populated.
- `short_gene_prevalence` (in-memory + cached RDS at [`short_gene_prevalence`](model.R#L36)). Per-synthetic-c80 prevalence map for unannotated short ORFs, used downstream by Step 3's c80s builders.
- `c80_variants_mapping` (in-memory + cached RDS at [`c80_variants_mapping`](model.R#L37)). Per-`(neighbor_c80_coarse, neighbor_gene_length)` mapping to the variant-resolved `neighbor_c80_fine`, used downstream by Step 3.
- Per-focal-per-genome shards on disk under [`neighbor_groups_by_genome`](model.R#L40) — one TSV per `<focal_label>/<genome>/<focal_c80>.tsv`. These are the files that `load_neighbors_across_genomes()` reassembles.
- Diagnostic figures under [`neighbor_figures`](model.R#L59): `fig1`–`fig5` per focal (operon-by-gene, operon-by-c80, operon-size distribution, selected-by-gene, selected-by-prevalence).

### Logic

Three sub-stages, all gated on the existence of the cached RDS:

**(a) Per-focal extraction** — [`extract_and_write_per_focal_neighbors()`](neighbor.R#L640)

For every focal in `focal_c80_df`, call [`parse_gene_neighbor()`](neighbor.R#L473), which walks the per-focal pipeline:

1. [`load_gene_neighbors()`](neighbor.R#L39) — read the focal's TSV, join `gene_to_c80` to attach `neighbor_c80_coarse` + cluster gene length + cluster prevalence.
2. [`compute_relative_positions()`](neighbor.R#L73) — parse the integer suffix from each `neighbor_gene_id` (Prokka convention), compute each neighbor's position relative to its focal, restrict to ±`upper_bound` and to focals with ≥ `min_positions` neighbors.
3. [`filter_by_flanking_coverage()`](neighbor.R#L133) — keep only focals with sufficient flanking-neighbor support. Auto-selects between strict (left-and-right) and relaxed (left-or-right) based on whether at least `focal_min_genomes` focals satisfy the strict criterion.
4. [`compute_operon_size()`](neighbor.R#L104) — count per-focal operon size; gate.
5. **First orient + extract pass** — [`orient_focal_gene_neighbors()`](neighbor.R#L238) decides the canonical direction per focal from its immediate non-NA flanking pair (left-anchor vs right-anchor; reverse the path if `right_anchor < left_anchor`). [`extract_gene_neighbor_patterns()`](neighbor.R#L323) then collapses identical position-by-position signatures into pattern groups (`group_excludeNA`, `group_includeNA`), keeping the dominant detailed pattern per broad group. Two diagnostic plots emit (`fig1`, `fig2`).
6. [`find_most_prevalent_operon_size()`](plot.R#L43) selects the dominant operon size; emit `fig3`.
7. **Second orient + extract pass** — re-runs orientation + pattern extraction on the subset of focals matching the dominant operon size, producing tighter groupings. Two more plots (`fig4`, `fig5`).
8. Filter neighbor groups to those with ≥ `focal_min_genomes` or ≥ `min_group_proportion` (default 0.05) of total focal coverage. Warn if surviving groups cover < `coverage_warn_threshold` (default 0.8) of total. Emit a per-focal RDS under [`neighbor_groups_by_focal`](model.R#L39).
9. [`write_gene_neighbor()`](neighbor.R#L604) shards the surviving rows into per-genome TSVs under `<neighbor_groups_by_genome>/<focal_label>/<genome>/<focal_c80>.tsv`, deriving `gene_member_genome` by stripping the trailing `_<index>` suffix from `gene_member` (Prokka convention). Warns if a genome contributes ≥ 2 focal copies (paralogs).

Errors and missing input files are caught per focal — one bad focal does not abort the batch.

**(b) Cross-genome assembly** — [`load_neighbors_across_genomes()`](neighbor.R#L680)

Recursively load every per-focal-per-genome TSV under `neighbor_groups_by_genome`, stamp each with `focal_c80` (from filename) and `path_label` (from grandparent dir = `pos`/`neg`), and `rbindlist` into one long table. Parallelized via `mclapply`. Cached to `neighbor_groups_rds`.

**(c) Label attachment** — always runs after the RDS is loaded.

[`assign_c80_to_short_genes()`](midas.R#L152) calls [`compute_short_gene_prevalence()`](midas.R#L93) to give every NA-`neighbor_c80_coarse` row a synthetic, focal-scoped label of the form `_<focal_c80>-<gene_type>_<rank>`. Per-focal scope is intentional: the same physical short gene next to two different focals receives two different synthetic labels. The function returns both the augmented `gene_neighbors` and a `short_gene_prevalence` lookup (the within-focal proportion, encoded as a negative number to distinguish from genome-wide prevalence).

[`compute_c80_variants()`](midas.R#L204) then builds a global mapping from `(neighbor_c80_coarse, neighbor_gene_length) → neighbor_c80_fine`. Clusters with one observed length keep their original label; clusters with multiple observed lengths get `<neighbor_c80_coarse>_<rank>` suffixes ordered by length. Globally scoped (no focal_c80 in the key) because cluster identities are defined once at MIDAS database build time. The mapping is left-joined back onto `gene_neighbors`.

### Tunables

From the `neighbor` and `plot` sections of the YAML:

- `focal_min_genomes` (default 10) — gate at multiple stages
- `focal_min_total_genomes` (default 30) — minimum total focal coverage
- `min_positions` (default 5) — minimum operon size
- `upper_bound` (default 10) — symmetric position window around focal
- `min_left_neighbors`, `min_right_neighbors` (default 2 each)
- `use_strict` (default `~` → auto)
- `min_group_proportion` (default 0.05) — second-pass survival threshold: a pattern group survives if it meets `focal_min_genomes` OR has at least this fraction of the focal's total genome support
- `gene_padding_bp` (default 100, from `plot` section) — bp gap between adjacent genes in the diagnostic figures `extract_gene_neighbor_patterns` lays out
- `coverage_warn_threshold` (default 0.8) — warn if surviving pattern groups together cover less than this fraction of total support; diagnostic-only

### Known caveats

- **Palindromic flanks** ([`orient_focal_gene_neighbors`](neighbor.R#L238) details). When `left_anchor == right_anchor`, both observed orientations collapse to `orientation = 1L`, so two observations of the same palindromic-core operon land in different canonical groups. Fix would compare `asc_str` vs `desc_str` on ties.
- **Paralog warning** ([`write_gene_neighbor`](neighbor.R#L604)). Multi-copy focals trigger a `warning()` but produce all copies in the output.
- **Synthetic label scope.** `_<focal>-<type>_<rank>` is per-focal — never compare these strings across focals. **Step 1 orientation is not preserved into Step 2.** Step 2 re-derives chromosomal order from `neighbor_gene_start`, intentionally discarding the per-focal orientation set here.

---

## STEP 2 — Per-genome operon graphs → maximal paths

**Driver section:** [pipeline.R:95–100](pipeline.R#L95-L100) (one-line orchestrator call to [`run_step2_path_stitching`](graph.R#L882))
**Helper file:** [graph.R](graph.R)

### Input

- `gene_neighbors` from Step 1 (one row per `(focal_c80, genome, neighbor_position)` with `neighbor_c80_coarse` and `neighbor_c80_fine`).

### Output

- `path_df` (in-memory + RDS at [`path_df`](model.R#L41)). One row per per-genome maximal path. Key columns:
  - `neighbor_genome`, `path_id`, `path_component_id` (graph-component id within that genome), `path_type` (pos/neg), `path_length`
  - `path_string` — `→`-joined `neighbor_gene_id`s
  - `c80_path_coarse` — `→`-joined coarse cluster labels
  - `c80_path_fine` — `→`-joined length-variant-aware labels
  - `c80_start`, `c80_end` — endpoints
- `esupport_df` (RDS at [`esupport_df`](model.R#L42)). Per-edge support: count of focal `gene_member`s contributing to each `(from, to, path_type)` edge per genome, with `support_genes` / `support_c80s_coarse` / `support_c80s_fine` lists.
- `path_genome_comp` — composite key added by the driver after re-loading `path_df`: `paste(path_genome, path_type, path_component_id, sep = "||")`. Step 3's join column.

### Logic

Single orchestrator: [`stitch_paths_across_focal_genes()`](graph.R#L283). Per-genome loop:

1. Re-derive an in-genome `position` from chromosomal `neighbor_gene_start`. The `anchor_index` is the row of the focal gene; `position = row_number - anchor_index`. **This overrides Step 1's focal-relative orientation.** All `gene_member`s within one genome thereby share a chromosomal frame (path direction is genome-relative, not focal-relative). Cross-genome direction unification is deferred to Step 3.
2. For each association class (pos/neg in `path_label`), call [`make_positional_edges()`](graph.R#L101) to convert positionally adjacent neighbors into directed `(from → to)` edges. Self-loops dropped; duplicate edges from different focals kept distinct.
3. Aggregate over all focals in this genome to produce `edge_support`: per `(from, to, path_type)`, count distinct supporting `gene_member`s and list them.
4. For each association class, call [`get_maximal_paths_by_type()`](graph.R#L171) to walk the support-weighted directed graph: each weakly connected component is enumerated for source-to-sink paths via DFS from each in-degree-0 node. **Errors out if the graph is not a DAG** — by construction it should be (chromosomal `position` is monotonic).
5. Translate path nodes from `neighbor_gene_id` to two parallel renderings via lookup tables: `c80_path_coarse` and `c80_path_fine`. Both are kept so downstream length-sensitive analyses don't have to re-derive either.
6. Stamp per-genome provenance, accumulate, and `bind_rows` across genomes.

After re-loading `path_df`, the driver builds `path_genome_comp` (the per-genome path key) and drops the now-redundant `path_genome` and `path_component_id`.

### Tunables

None at this stage. Inherits from Step 1.

### Known caveats

- **Both edges and paths are within-genome.** No cross-genome edges; cross-genome aggregation happens in Step 3.
- **DAG assumption.** Cyclic neighborhoods (rare; would require a tandem `A→B→A` in chromosomal order) trigger a hard error.
- **Per-edge support is not propagated into the path output.** It lives in `esupport_df` keyed by `(from, to, path_type, genome)` — re-join if needed.

---

## STEP 3 — Cross-genome consolidation → three granularity levels

**Driver section:** [pipeline.R:103–128](pipeline.R#L103-L128) (orchestrator call to [`run_step3_consolidation`](path.R#L423))
**Helper files:** [graph.R](graph.R), [path.R](path.R), [parse.R](parse.R) (decorators)

### Input

- `path_df` from Step 2 (with `path_genome_comp` key applied).
- `c80_variants_mapping`, `gene_meta`, `cluster_80`, `short_gene_prevalence` from earlier driver setup / Step 1.

### Output

Five TSVs are written to `step3_path/`. They split into three granularity **levels**:

| Level | Path-level table | Per-gene table | Purpose |
|---|---|---|---|
| L1 — coarse / canonical | [`canonical_paths`](model.R#L43) | [`canonical_paths_c80s`](model.R#L46) | One row per surviving canonical operon (coarse). |
| L2 — per-isoform | [`canonical_paths_fine`](model.R#L44) | [`canonical_paths_fine_c80s`](model.R#L47) | One row per length-variant isoform within a canonical. |
| L3 — per-genome | [`canonical_paths_per_genome`](model.R#L45) | _(none — already gene-level)_ | One row per `(canonical, contributing genome)` with full attribution. |

### Logic

Eight stages.

**1. [`collapse_paths_across_genomes()`](graph.R#L422)**
Groups `path_df` on `(c80_path_coarse, path_length, path_type)`. One row per coarse operon shape per type. Records `n_genomes`, `neighbor_genomes` (`;`-joined), and `per_genome_path_w_ids` (`;`-joined `path_genome_comp`s — the provenance pointer back to the per-genome paths). Direction is **not** canonicalized here (e.g., `A→B→C` vs `C→B→A` remain two rows).

**2. [`generate_canonical_path()`](graph.R#L537)**
Applies [`normalize_path()`](graph.R#L469) (lex-min of forward vs reverse, computed on a *cleaned* token vector — synthetic small-ORF tokens stripped, adjacent duplicates collapsed via [`clean_for_orientation()`](graph.R#L445); the chosen direction is then applied to the original full vector) to every row, assigning a surrogate `canonical_path_id`. Aggregates direction-mirror rows by summing `n_genomes` and concatenating provenance. Applies the `path_min_genomes` gate. **This is the survival cut for an operon** — anything below `path_min_genomes` is dropped from L1 onward.

**3. [`compute_joint_components()`](graph.R#L616)**
Builds an undirected, type-collapsed gene-level graph from all canonical paths and computes connected components. Returns a `(node, joint_component_id)` map — every centroid_80 that appears in any canonical path is grouped with every other c80 it's transitively adjacent to anywhere in the corpus. Type info is intentionally discarded (a hub gene that appears in both a `pos` and a `neg` path bridges them into one component).

**4. [`decorate_paths_with_components()`](graph.R#L678)**
For each canonical path, look up the `joint_component_id` of every gene token and concatenate the unique values into a string `joint_component_ids` column. In practice each canonical path lives in exactly one component; multi-component strings only arise when a gene was lost from the map.

**5. [`orient_paths_within_component()`](graph.R#L831)**
Within each joint component, pick the longest path (by *cleaned* length) as direction reference and flip every other path to maximize forward-substring overlap against the reference (via [`max_overlap()`](graph.R#L377)). Result: paths in a component all read left-to-right in a mutually consistent direction. Direction is **only** consistent within a component; the absolute reference direction is biology-agnostic (inherited from `normalize_path`'s lex rule). See the function's docstring for edge cases (longest-path tie, zero-overlap, palindromic, single-gene bridge).

**6. Bake `uid`** (inside [`run_step3_consolidation`](path.R), after the within-component re-orientation)
`uid = "cmp{joint_component_ids}-{path_type}-{canonical_path_id}-ng{n_genomes}"`. Self-describing primary key for the canonical (L1) table.

**7. Per-isoform expansion** — [`expand_canonical_paths_to_fine()`](path.R#L122)
For each surviving canonical, walk the provenance chain `canonical_path_id → collapsed_path_id → per_genome_path_w_ids → path_genome_comp` via the shared backbone [`explode_canonical_into_collapsed_paths()`](path.R#L70). The walker carries `c80_path_fine` (fine renderings); rows whose original coarse direction differs from the canonical direction get their fine string token-reversed. Aggregate per `(canonical_path_id, c80_path_fine_canonical)`, assign deterministic `isoform_rank` (sorted by `desc(n_fine_genomes), c80_path_fine_canonical`), build hierarchical `uid_fine = "{uid}-iso{isoform_rank}-ngf{n_fine_genomes}"`.

**8. Per-gene expansions + decorators**

[`build_canonical_paths_c80s()`](path.R#L264) — split each canonical path into per-gene rows at coarse resolution. Attach max-over-isoforms gene length, gene metadata, cluster_80 metadata, and small-ORF prevalence. **`neighbor_gene_length` here is the max across isoforms of the cluster** — not a single observed gene's length. Use the fine table for exact per-isoform lengths.

[`build_canonical_paths_fine_c80s()`](path.R#L342) — same but at isoform resolution. `neighbor_gene_length` is the exact per-isoform value. Joint-component, gene-meta, and cluster-80 annotations are inherited via the coarse `neighbor_c80_coarse` (no isoform granularity exists for those).

[`decorate_c80s_w_smallORFs()`](parse.R#L95) — decode synthetic small-ORF labels into queryable `is_smallORF`, `centroid_80`, `smallORF_type`, and per-operon `n_smallORFs` / `n_focal` / `dist_to_smallORFs`. Run on both coarse and fine, with `group_key = "uid"` and `"uid_fine"` respectively.

[`decorate_c80s_w_truncation()`](parse.R#L181) — fine only. Adds:
- `is_truncated` (per row) when `neighbor_gene_length < truncation_cutoff * neighbor_c80_length_coarse`
- `truncate_ratio` (per row), floor-truncated to 3 decimals
- `n_truncated` (per-isoform broadcast)
- `is_fragmented` (per row) when this row's `neighbor_c80_coarse` shows up under ≥ 2 distinct `neighbor_c80_fine` values within this `uid_fine` — same coarse cluster at multiple lengths in one operon
- `fragmented_c80s`, `n_fragmented_c80s`

**Why fine-only:** the coarse `neighbor_gene_length` is the max-over-isoforms, so "shorter than cutoff" there would mean "even the longest isoform is shorter" — not what truncation should mean.

[`expand_canonical_paths_per_genome()`](path.R#L185) — L3 master table. Re-walks the provenance chain via the shared backbone, this time carrying both `c80_path_fine` and `path_string` payloads. For each per-genome contribution, emit the raw `gene_path` and the canonical-aligned `gene_path_canonical`, plus identity columns inherited from L1 and L2 (`uid`, `uid_fine`, `canonical_path_id`, `fine_canonical_id`, `isoform_rank`, `needs_flip`).

After all expansions, the driver column-orders `c80s_coarse` / `c80s_fine` for readability and writes the five TSVs.

### Tunables

From `path`:

- `path_min_genomes` — passed to `generate_canonical_path()` as the survival cut
- `truncation_cutoff` — passed to `decorate_c80s_w_truncation()`

### Known caveats

- **Hub-gene mega-components** ([`compute_joint_components`](graph.R#L616)). A promiscuous regulator or ribosomal gene can fuse many unrelated operons into one component. No automated flag — sanity-check component sizes if results look unexpectedly merged.
- **Mixed-type warning in [`generate_canonical_path`](graph.R#L537)** is emitted but the diagnostic data is built and discarded; if the invariant is violated, the subsequent `summarise(path_type = unique(path_type))` aborts with a length-mismatch instead.
- **`n_genomes` double-counts** in `generate_canonical_path` if any genome carries both directions of the same operon (rare).
- **Mixed separators are inherited:** `collapsed_path_id` is comma-separated, `neighbor_genomes` is semicolon-separated. Not harmonized.

---

## STEP 4 — Trait-associated block extraction + representative ranking

**Driver section:** [pipeline.R:131–139](pipeline.R#L131-L139) (one-line orchestrator call to [`run_step4_block_extraction`](blocks.R#L457))
**Helper file:** [blocks.R](blocks.R)

### Input

- `c80s_coarse` from Step 3 (the L1 per-gene table). Must contain `joint_component_id`, `canonical_path_id`, `path_type`, `n_genomes`, `uid`, `is_focal`, `neighbor_c80_coarse`, and the trait stat column (default `beta`).
- `c_paths`, `collapsed_paths`, `path_df` from Step 3 — used by `map_representatives_to_genomes()` to walk back to per-genome attribution.

### Output

- `representatives` (TSV at [`rep_path_df`](model.R#L49)). One row per non-redundant trait-associated block, with `block_uid = "cmp{component}-{path_type}-rank{rep_rank}-nge{block_n_genes}"`, `relation_to_selected` ∈ {`selected`, `subset`, `superpath`, `overlap`, `disjoint`}, and `canonical_paths` / `canonical_uids` provenance lists.
- `rep_slim` (TSV at [`uid_path_df`](model.R#L50)). Per-genome attribution: one row per `(block_uid, canonical_uid, neighbor_genome, left_orig, right_orig)`, after the `select(...) %>% unique()` projection in the driver.
- Diagnostic message via [`diagnose_rep_overlaps()`](blocks.R#L264): how many `(component, path_type)` groups have substring-overlap pairs that survived (i.e., reps that share a substring but neither contains the other — the case the subset-only redundancy check can't collapse).
- Optional `rep_heatmap.pdf` (block × genome presence/absence) under the `uid_path_df` directory, only when there are ≥ 2 reps and the matrix is at least 3×3.

### Logic

Five stages, all keyed on `(joint_component_id, path_type)`:

**1. [`keep_focal_blocks()`](blocks.R#L49)** (called inside `aggregate_blocks`)
Within each `(joint_component_id, canonical_path_id, path_type)` group, walk rows in their existing per-path order. Tag each as a focal hit (`is_focal == TRUE`, gated separately on `value_col`) or non-hit. Cluster adjacent hits into blocks: a new block starts when more than `allow_gaps` non-hit rows intervene (gap rule: `gaps > allow_gaps + 1`). Drop non-hit rows. Sign (`+` / `−` `is_label`) is recorded but does **not** split blocks — a block can mix signs if they're close enough.

**2. [`aggregate_blocks()`](blocks.R#L95)**
For each block, dedup consecutive duplicate c80 tokens (`A A B B C → A B C` via [`dedup_consecutive_vec`](blocks.R#L381)) and emit `block_c80s_path` (`→`-joined), `n_genes`, `left_orig` / `right_orig` endpoints, `edge_pair`, and the singleton flag. Then aggregate at `(joint_component_id, path_type, edge_pair, block_c80s_path)`: `block_n_paths` = number of canonical paths supporting this block shape, `block_n_genomes` = sum of `n_genomes` across them, and parallel `canonical_paths` / `canonical_uids` lists. Per-component frequency: `block_freq = block_n_genomes / block_total`.

**3. Reference selection** — first half of [`rank_block_representatives()`](blocks.R#L167)
Per `(joint_component_id, path_type)`, sort by `desc(block_freq), desc(block_n_genomes), desc(block_n_genes), desc(block_n_paths), block_c80s_path` (tie-break for determinism) and take the top row as the dominant block. Stored as `selected_tbl`.

**4. [`annotate_group()`](blocks.R#L324)** (called per group via `group_modify`)
Tag every block with its [`get_relation()`](blocks.R#L416) to the selected block: `selected | subset | superpath | overlap | disjoint`. Then run the greedy rep-construction loop:

- Sort the group by `desc(block_n_genes), desc(block_n_genomes), desc(block_n_paths), block_c80s_path` (longest first).
- Install the selected reference as `rep_rank = 1`.
- For each remaining row: if its `block_c80s_path` is a contiguous subset (forward direction only — see caveat) of any already-installed rep, mark redundant (`rep_rank = 0`, `is_redundant = TRUE`). Otherwise, install it as a new rep with the next integer rank.

**5. [`map_representatives_to_genomes()`](blocks.R#L223)**
Use the shared backbone [`explode_canonical_into_collapsed_paths()`](path.R#L70) to walk each rep's contributing canonical paths back to the per-genome paths that carry them. `separate_rows` explodes the parallel `canonical_paths` / `canonical_uids` lists together so each exploded row pairs one canonical with its `uid`. **Sanity check:** `stopifnot(path_type == path_type_per_genome)` — block path_type must match per-genome type for every row.

The driver writes `representatives` and the `rep_slim` projection (`block_uid, canonical_uid, neighbor_genome, left_orig, right_orig`), then runs the diagnostic and optional heatmap.

### Tunables

From `blocks`:

- `allow_gaps` (default 2) — passed to `aggregate_blocks` (and through to `keep_focal_blocks`)
- `min_overlap` (default 1) — passed to `rank_block_representatives` (and through to `get_relation`)
- `min_shared` (default 2) — passed to `diagnose_rep_overlaps`; diagnostic-only

Still hardcoded in the driver: `value_col = "beta"` for `aggregate_blocks`.

### Known caveats

- **Redundancy is subset-only, forward-direction.** [`is_contig_subseq()`](blocks.R#L387) does not test the reverse of `p`. A block and its exact mirror survive as two separate reps. A reverse-aware variant (`is_contig_subseq(rev(p), r)`) is documented in `annotate_group` but not wired in.
- **Substring-overlap pairs survive.** The diagnostic `diagnose_rep_overlaps` exists precisely because reps that share a contiguous substring but neither contains the other can't be collapsed by the subset-only check, and may double-count genomes downstream.
- **Block-level `block_uid` collides with no other `uid`** because of the `cmp...-rank...-nge...` prefix; do not confuse with the canonical-level `uid` from Step 3 or the isoform-level `uid_fine`.

---

## STEP 5 + STEP 6 — Summaries, fine-only BLAST sampler, gggenes plots

**Driver section:** [pipeline.R:142–163](pipeline.R#L142-L163) (Step 5 / Step 6 orchestrator calls — [`run_step5_parse`](parse.R#L659) and [`run_step6_figures`](plot.R#L768))
**Helper files:** [parse.R](parse.R) (data prep), [plot.R](plot.R) (plotters)

The terminal post-processing block, split across two orchestrators: [`run_step5_parse`](parse.R) (summaries, fine-coverage selection sets, exemplar-genome sampling, per-(isoform, genome) BLAST gene-id files) and [`run_step6_figures`](plot.R) (multi-fill gggenes plots — global and per-joint-component, faceted by `updated_path_type`). Both are designed to stay re-runnable in isolation: the three Step-3 TSVs (Step 5) and the four selection-and-c80s TSVs (Step 6) are loaded at the call site in [pipeline.R](pipeline.R) and passed in as explicit arguments, so the helper functions themselves contain no `read_delim` calls.

### Input

Read from disk so this block is self-contained:
- [`canonical_paths_c80s`](model.R#L46) → `c80s_coarse`
- [`canonical_paths_fine_c80s`](model.R#L47) → `c80s_fine`
- [`canonical_paths_per_genome`](model.R#L45) → `per_genome`

In-memory from Step 1:
- `gene_neighbors` — used by the sampler to attach per-gene metadata.

### Output

- `coarse_summary` / `fine_summary` (in-memory). One row per `uid` and per `uid_fine` respectively, with per-operon counts (`n_genes`, `n_smallORFs`, `n_focal`, `n_truncated`, `n_fragmented_c80s`) and a `coarse_path_string` / `fine_path_string` rendering. No filter is applied inside the summarizer; the caller filters.
- `selected_fine` (in-memory). Fine isoforms surviving `n_fine_genomes >= ceiling(path_min_genomes * fine_coverage_ratio)`. The selection set for everything downstream in this block.
- `selected_coarse` (in-memory). `coarse_summary %>% semi_join(selected_fine, by = "uid")` — the coarse uids whose fine isoforms made it through.
- `fine_long` (in-memory). Long-format per-gene table, one row per gene of each sampled `(uid_fine, neighbor_genome)`. Built in two stages by [`sample_genome_from_fine_paths`](parse.R#L383) (12 columns) and then [`enrich_fine_long`](parse.R#L470) (joins per-isoform context from `c80s_fine`).
  - **Sampler output (12 cols):** `uid_fine`, `neighbor_genome`, `position_in_path`, `gene_id`, `neighbor_c80_coarse`, `neighbor_c80_fine`, `neighbor_contig_id`, `neighbor_gene_start`, `neighbor_gene_end`, `neighbor_gene_strand`, `neighbor_gene_type`, `neighbor_gene_length` (per-gene observed length from `gene_neighbors`).
  - **Enrichment adds** the contiguous range `uid_fine:centroid_80` from `c80s_fine`, minus three columns dropped to avoid collisions: `neighbor_c80_coarse` and `neighbor_c80_fine` (already present, identical values per joined row), and `neighbor_gene_length` (intentionally not merged — the c80s_fine version is per-isoform consensus length, not per-gene observed; same name, different semantics).
  - **Net added columns:** `n_fine_genomes`, `neighbor_c80_length_coarse`, `genome_prevalence`, `sample_prevalence`, `cor_to_b`, `focal_label`, `beta`, `trait`, `genome_counts`, `is_focal`, `is_smallORF`, `smallORF_type`, `n_smallORFs`, `n_focal`, `dist_to_smallORFs`, `is_truncated`, `truncate_ratio`, `n_truncated`, `is_fragmented`, `fragmented_c80s`, `n_fragmented_c80s`, `centroid_80_genome_counts`, `centroid_80`.
- Per-`(uid_fine, neighbor_genome)` gene-id TSVs under [`parse_genome_paths_dir`](model.R#L57), file pattern `fine_<uid_fine>_<neighbor_genome>.tsv`. Bare gene-id lists in canonical-path order; the input to an external BLAST workflow.
- gggenes plots of the selected operons / isoforms under [`parse_coarse_figures`](model.R#L60) and [`parse_fine_figures`](model.R#L61) (both resolve to `step6_figures/`). One PDF per fill mode listed in `parse.fill_modes`: `coarse_operons_<fill_by>.pdf`, `fine_operons_<fill_by>.pdf`, plus per-component variants under `02_by_component_coarse/` and `03_by_component_fine/`.

### Logic

Six stages.

**1. Summaries.** [`summarize_coarse_operons()`](parse.R#L254) and [`summarize_fine_isoforms()`](parse.R#L314) roll up the decorated c80s tables to one row per operon / isoform. Both are filter-free — the carry-through includes only the columns needed for downstream sampling and plotting (notably `uid` is carried into `fine_summary` so the post-hoc coarse-summary merge can group fine isoforms by their coarse parent without parsing `uid_fine` strings).

The driver also enriches `coarse_summary` with isoform-mapping columns post-hoc — `n_isoforms_raw` (from `fine_summary`), and `n_isoforms_filtered` / `n_coarse_genome_filtered` / `uid_fine_list` (from `selected_fine`, available after stage 2). `coarse_summary$neighbor_genomes` is intentionally not carried — genome-level traceback is a fine-isoform concern.

**2. Selection.** `path_min_genomes <- cfg_get(job_config, "path_min_genomes")` and `ratio <- cfg_get(job_config, "fine_coverage_ratio")`, then `selected_fine <- fine_summary %>% filter(n_fine_genomes >= ceiling(path_min_genomes * ratio))`. The example sets `ratio = 0.25` (permissive); historical default `0.5` reproduces the half-coverage rule (an isoform must have at least half as many supporting genomes as the canonical-path survival threshold). `ceiling` (not `floor`) keeps the bar from dropping below the configured fraction on odd `path_min_genomes`.

**3. Fine-only sampling.** [`sample_genome_from_fine_paths()`](parse.R#L383). Draws **one** random genome per surviving fine isoform (no adaptive variation by isoform count). `set.seed` is called once (with `seed` from YAML `parse`, default 616) for reproducibility, then `group_by(uid_fine) %>% slice_sample(n = 1)`. The chosen `gene_path_canonical` is exploded by `unnest_longer` to one row per gene; `position_in_path` is recorded; per-gene metadata (`neighbor_c80_coarse`, `neighbor_c80_fine`, `neighbor_contig_id`, `neighbor_gene_start`, `neighbor_gene_end`, `neighbor_gene_strand`, `neighbor_gene_type`, `neighbor_gene_length`) is merged from a deduped `gene_neighbors` lookup (deduped on `neighbor_gene_id` to prevent row explosion at join time).

Sampler output is 12 columns. Sample-level columns (`uid_fine`, `neighbor_genome`) repeat across the gene rows. Direction is canonical because the explode reads `gene_path_canonical` (already flipped by `expand_canonical_paths_per_genome` via the L2 `needs_flip` flag); `position_in_path` is therefore the 1-indexed canonical position. The full path string is not kept on the output — reconstruct it from `(gene_id, position_in_path)` per sample if needed. Per-gene coordinates (contig + start + end + strand) support nucleotide sequence extraction for BLAST workflows.

Coarse-level BLAST hits are produced **post-hoc**, not by this block: aggregate fine BLAST results on `sub("-iso.*", "", uid_fine)` to recover per-`uid` coarse hits.

**4. Per-isoform context enrichment.** [`enrich_fine_long(fine_long, c80s_fine)`](parse.R#L470). Left-joins the per-(uid_fine, c80) annotations from `c80s_fine` onto `fine_long`, bringing in trait stats (`beta`, `cor_to_b`, `is_focal`), per-isoform-broadcast counts (`n_focal`, `n_smallORFs`, `n_truncated`, `n_fragmented_c80s`), small-ORF and truncation flags, and `centroid_80`. Join key is `(uid_fine, position_in_path)` — `c80s_fine` doesn't carry a position column, so the helper derives it via `row_number()` within `uid_fine` (rows are in canonical order by construction). Three columns are dropped from the c80s_fine side before joining: `neighbor_c80_coarse` / `neighbor_c80_fine` (already in fine_long, hygiene) and `neighbor_gene_length` (semantically distinct — sampler version is per-gene observed; c80s_fine version is per-isoform consensus). A `stopifnot` asserts no row multiplication.

**5. BLAST gene lists.** [`write_blast_gene_lists()`](parse.R#L501). Long-format consumer: groups the enriched `fine_long` by `(uid_fine, neighbor_genome)`, sorts by `position_in_path`, and writes one TSV per group. The writer reads only `(uid_fine, neighbor_genome, position_in_path, gene_id)` — extra enrichment columns ride along harmlessly. This is the last stage of `run_step5_parse`; the figure-rendering stage below is owned by Step 6.

**6. Step 6 — gggenes figures.** [`run_step6_figures()`](plot.R#L768) wires four plotters at global and per-component scope: [`plot_coarse_operons()`](plot.R#L520) and [`plot_fine_operons()`](plot.R#L572) (one PDF per `fill_by` mode at the global level), [`plot_coarse_operons_by_component()`](plot.R#L617) and [`plot_fine_operons_by_component()`](plot.R#L672) (one PDF per `(joint_component_id, fill_by)` under `02_by_component_coarse/` and `03_by_component_fine/`, faceted by `updated_path_type`). Each plotter consumes its summary as a `semi_join` filter on `uid` / `uid_fine`. The shared layout helper [`.layout_operon_tracks()`](plot.R#L374) computes per-track positions and a single `fill_symbol` glyph per gene with this precedence (top wins):

| Glyph | Condition |
|---|---|
| `U` | `is_focal == TRUE & beta > 0` (focal upregulated) |
| `D` | `is_focal == TRUE & beta < 0` (focal downregulated) |
| `F` | `is_fragmented` (fine only) — beats `T` because fragmentation is the rarer, more specific signal (the short copy of a fragmented pair is almost always also truncated; if `T` won, `F` would never fire) |
| `T` | `is_truncated` (fine only — coarse `is_truncated` is forced to NA) |
| `s` | `is_smallORF` |
| _(none)_ | otherwise |

The fine plots order tracks so isoforms of the same coarse `uid` are adjacent; the coarse plots order by `uid` only. `fill_modes` is read from YAML `parse.fill_modes` (subset of `{beta, sample_prevalence, cor_to_b, fill_gene}`). All PDFs land under `step6_figures/`.

### Tunables

- `path_min_genomes` (from YAML `path`, default 20) — `fine_coverage_ratio * path_min_genomes` is the fine-isoform survival threshold.
- `fine_coverage_ratio` (from YAML `parse`, default 0.25 in the example; historical default 0.5) — fraction of `path_min_genomes` an isoform must reach to survive. 0.5 reproduces the historical half-coverage rule.
- `seed` (from YAML `parse`, default 616) — RNG seed for the per-isoform exemplar-genome draw. Reproducible across runs.
- `truncation_cutoff` (from YAML `path`, default 0.8) — already applied upstream by `decorate_c80s_w_truncation`; used here only via the inherited `is_truncated` flag.
- `gene_padding_bp` (from YAML `plot`, default 100) — bp gap between adjacent genes in the Step 6 gggenes layouts. Same value used in Step 1 diagnostic figures.

### Known caveats

- **Per-genome quality scoring is intentionally skipped.** The sampler is plain random; the function leaves a `weight_by` seam open for the day a per-genome quality column exists.
- **Synthetic small-ORF rows can surface NA** in any of the per-gene metadata columns (`neighbor_c80_coarse`, `neighbor_c80_fine`, contig coordinates, strand, gene_type, gene_length) if their PROKKA gene_id isn't in `gene_neighbors`. Spot-check if observed.
- **Plot inputs MUST be pre-filtered.** The summarizers no longer filter internally, so passing raw `fine_summary` to `plot_fine_operons` can plot hundreds of tracks. The driver always passes `selected_fine` / `selected_coarse` (the fine-coverage survival sets); standalone callers should mirror that.
- **Two `neighbor_gene_length` semantics, only one in the output.** `gene_neighbors.neighbor_gene_length` is per-gene observed; `c80s_fine.neighbor_gene_length` is per-isoform consensus. The enriched `fine_long` keeps the per-gene version (from the sampler stage); the c80s_fine version is dropped at enrichment time. Don't merge them into one column — they answer different questions.
