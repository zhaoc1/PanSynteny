# SCHEMA.md - File schemas and data formats

The single source of truth for every file the pipeline reads or writes: where it lives, what columns it has, who writes it, who reads it, and any non-obvious semantics. Every entry below is keyed by the [model.R](../R/model.R) target name (where one exists) so you can grep both sides at once.

**How this is organised.** Sections follow [model.R](../R/model.R)'s **consumer-based** Input / Output split (section 3 = `# Input` block, section 4 = `# Output` block). A file written by Step 0 (e.g., the focal_meta cache, the catalog files) lives under "Pipeline inputs" because that's how it's used downstream, even though it's produced upstream. The section 3 entries note the producer in each "Writer" row.

For pipeline behaviour and run order, see [USER_GUIDE.md](USER_GUIDE.md). For per-step input/output logic, see [STEPS.md](STEPS.md). For c80-column semantics (coarse vs fine), see [PIPELINE.md](PIPELINE.md).

---

## 1. User-provided inputs (paths in the YAML)

These are the files the user supplies; everything downstream is derived from them.

> **Path placeholders.** `load_job_config` ([config.R](../R/config.R)) expands `{proj_dir}`, `{species_id}`, `{midasdb_dir}`, and `{input_dir}` in every YAML string scalar at load time - so the paths below can be written as e.g. `"{input_dir}/focal.tsv"` or `"{midasdb_dir}/pangenomes/{species_id}/clusters_80_info_updated.tsv"`. `{proj_dir}` is the **output** root (used as-is from `job.proj_dir`); `{input_dir}` is the required **input** root holding the two user-authored files for this run (the YAML config + the focal_meta TSV), declared under `job.input_dir`. Plain absolute paths work unchanged. The same global placeholders are also recognised by `build_genome_catalog.py` for `sources:` paths; each source entry may additionally declare per-source local placeholders (any non-reserved string field, e.g. `ecor_dir: "/path"` -> `{ecor_dir}`) usable within its own `genes_info` / `genomes_dir`.

### `data.focal_meta` - focal-gene metadata

| | |
|---|---|
| **Path** | Absolute path, declared in YAML. |
| **Header** | Yes (tab-delimited). |
| **Writer** | User (upstream of this pipeline). |
| **Reader** | [`prepare.R`](../prepare.R) reads + processes; cached to `focal_meta` target. |

**Required columns**

| Column | Type | Purpose |
|---|---|---|
| `focal_c80` | str | The centroid_80 id. Used as the per-focal-TSV filename and as the join key across the pipeline. |
| `focal_label` | str | `pos` / `neg` - direction tag for Step 1 sharding and Step 5 fill mode. |
| `is_focal` | bool | Step 1 gate. Only `TRUE` rows drive neighbor extraction. (Can be derived if `prepare.score_col` is set; see below.) |
| `gene_label` | str | User-defined annotation/category for the focal gene. Used by Step 5 `gene_label` fill mode. Distinct from `.genes`-file `gene_type` (which is the GFF feature type). |

**Optional columns** (carried through to the cache; consumed when present)

| Column | Consumed by |
|---|---|
| `cor_to_b` | `prepare.score_col` threshold path; Step 5 `cor_to_b` fill mode. |
| `beta` | Step 5 `beta` fill mode. |
| `sample_prevalence` | Step 5 `sample_prevalence` fill mode. |
| `trait`, `genome_counts` | Carried for provenance; no live consumer. |

**Notes**
- Extra columns beyond this list pass through to the cache unchanged.
- `is_focal` in the input is **overwritten** when `prepare.score_col` is set (with a `warning()`). To preserve a hand-curated `is_focal` (e.g. for context rows with `is_focal = FALSE`), set `prepare.score_col: ""`.

---

### `data.clusters_80_updated` - per-c80 metadata (midasdb)

| | |
|---|---|
| **Path** | Absolute path, declared in YAML. Typically `<midasdb>/pangenomes/<species_id>/clusters_80_info_updated.tsv`. |
| **Header** | Yes (tab-delimited). |
| **Writer** | midasdb maintainer (upstream of this pipeline). |
| **Reader** | [`prepare.R`](../prepare.R) seeds a local copy (once, `overwrite = FALSE`). |

**Columns** (those `load_c80_tables` reads after the rename - first four are required; later columns are optional, pass through unused)

| Column (raw) | Renamed to (in R) | Required? |
|---|---|---|
| `c80` | `c80` | **yes** |
| `centroid_80_gene_length` | `neighbor_c80_length_coarse` | **yes** |
| `centroid_80_genome_prevalence` | `genome_prevalence` | **yes** |
| `centroid_80_genome_counts` | (unchanged) | **yes** (terminator for `select(c80:centroid_80_genome_counts)`) |
| `COG_category`, `COG_description`, `Description`, `PFAMs`, `add_annot`, `EC`, `COG_category_old`, `KEGG_Pathway` | - | no - present in UHGG but unused |

**Seed-once semantics.** `prepare.R` copies this file from the YAML path to `<data_dir>/<species_id>/clusters_80_info_updated.tsv` only if the local copy doesn't exist. Subsequent runs keep the local copy so hand edits (e.g. appending ECOR-derived c80 rows) persist. To refresh from the source, `rm` the local file and re-run `prepare.R`.

---

### `sources` (YAML list) - declares genome sources for `build_genome_catalog`

Each entry declares one source of genomes to merge into the unified catalog. Schema **per entry**:

| Field | Required? | Notes |
|---|---|---|
| `name` | yes | Free label for logging. |
| `type` | yes | `midas` or `prokka`. |
| `genes_info` | yes | Absolute path to the membership file (`{midasdb_dir}` / `{species_id}` / `{proj_dir}` and per-source locals expanded). |
| `c80_col` | yes | 1-based column index of `centroid_80` in `genes_info`. |
| `length_col` | midas: optional (default 8); **prokka: required** | 1-based column index of `gene_length`. |
| `genomes_dir` | yes | Directory of `<genome_id>/<genome_id>.genes` (midas) or `<genome_id>/<genome_id>.gff` (prokka). |
| *(any other string field)* | no | Per-source local placeholder. E.g. `ecor_dir: "/abs/path"` becomes `{ecor_dir}` usable in this entry's `genes_info` / `genomes_dir`. Scoped to its own source entry. Globals are expanded inside the local value too. |

#### `sources[type=midas].genes_info` (UHGG `genes_info.tsv`)

| Col | Name |
|---|---|
| 1 | `gene_id` |
| 2 | `centroid_99` |
| 3 | `centroid_95` |
| 4 | `centroid_90` |
| 5 | `centroid_85` |
| 6 | `centroid_80` <- `c80_col: 6` |
| 7 | `centroid_75` |
| 8 | `gene_length` <- `length_col: 8` (default) |
| 9 | `marker_id` |

#### `sources[type=prokka].genes_info` (user-curated, e.g. `ecor_gene_centroid80.tsv`)

Required columns (positions configurable via `c80_col`/`length_col`):

| Col (suggested) | Name | Required? |
|---|---|---|
| 1 | `gene_id` | **yes** - must match the prokka `locus_tag` |
| `c80_col` | `centroid_80` | **yes** - the UHGG centroid_80 this ECOR gene maps to (e.g. from BLAST) |
| `length_col` | `gene_length` | **yes** - typically `end - start + 1` from the GFF; easiest to merge in from prokka's `.tsv` `length_bp` |
| any | any extra columns | pass through unused |

Header row required.

#### `sources.genomes_dir` directory layout

```
<genomes_dir>/
├── <genome_id_1>/
│   ├── <genome_id_1>.genes   (midas: ships with midasdb; prokka: derived by gff_to_genes from <g>.gff)
│   └── <genome_id_1>.gff     (prokka only)
├── <genome_id_2>/
│   ├── ...
```

The genome_id derivation from gene_id is shared with `focal_neighbor_list.sh`: strip the trailing `_NNNNN` field. `GUT_GENOME000040_00388` -> `GUT_GENOME000040`; `GCF_900448275.1_00001` -> `GCF_900448275.1`.

---

## 2. Per-genome annotation files

### `.genes` (per-genome gene annotation TSV)

| | |
|---|---|
| **Path** | `<genomes_dir>/<genome_id>/<genome_id>.genes` |
| **Header** | Yes. |
| **Writer** | midasdb (for type=midas); [`gff_to_genes.py`](../scripts/gff_to_genes.py) in-place from `.gff` (for type=prokka). |
| **Reader** | [`get_neighbor.sh`](../scripts/get_neighbor.sh) for +/-n_genes flank extraction. The 7-column per-focal output (see below) inherits all six columns from this file. |

| Col | Name | Notes |
|---|---|---|
| 1 | `gene_id` | Matches the prokka `locus_tag`. |
| 2 | `contig_id` | Prokka-style: `gnl|Prokka|<genome_id>_<contig_number>`. |
| 3 | `start` | 1-based inclusive. |
| 4 | `end` | 1-based inclusive. |
| 5 | `strand` | `+` / `-`. |
| 6 | `gene_type` | `CDS` / `rRNA` / `tRNA` / `tmRNA` (from the GFF feature type). |

### `.gff` (Prokka GFF3 source, prokka sources only)

| | |
|---|---|
| **Path** | `<genomes_dir>/<genome_id>/<genome_id>.gff` |
| **Header** | GFF3 spec (`##gff-version 3` + comment lines). |
| **Reader** | [`gff_to_genes.py`](../scripts/gff_to_genes.py) (via `gffutils`). Extracts: `seqid` -> contig_id, `start`, `end`, `strand`, `featuretype` -> gene_type, `attributes.locus_tag` -> gene_id. Skips features without an `ID` attribute (e.g., CRISPR repeats) and contig-level `prokka`-source rows. |

---

## 3. Pipeline inputs - Step 0 caches and seeded copies

These are the model.R `# Input` targets: paths the analytical pipeline (and the Step 0 bash chain) **reads** during a run. They're produced upstream by `prepare.R`, `build_genome_catalog.py`, or `build_neighbor_lists.sh` - but from the consumer's perspective they are inputs to the work that follows. Paths reference [model.R](../R/model.R) target keys.

### `run_config` - config snapshot

| | |
|---|---|
| **Path** | `{proj_dir}/step1_setup/run_config.yaml` |
| **Format** | Verbatim copy of the CLI-arg `<config.yaml>`. |
| **Writer** | [`prepare.R`](../prepare.R), overwrites every run. |
| **Reader** | None (provenance / audit). |

### `focal_meta` (cache) - processed focal table

| | |
|---|---|
| **Path** | `{proj_dir}/step1_setup/gene_meta_full.tsv` (filename kept for back-compat) |
| **Header** | Yes (tab-delimited). |
| **Writer** | [`prepare.R`](../prepare.R), overwrites every run. |
| **Reader** | [`pipeline.R`](../pipeline.R) at startup; passed to Step 1 as `focal_c80_df`. |

Same columns as the input `data.focal_meta` (see section 1), with `is_focal` possibly overwritten by `prepare.score_col` thresholds.

### `gene_list` - missing-neighbor focals

| | |
|---|---|
| **Path** | `{proj_dir}/step1_setup/gene_list.tsv` |
| **Header** | None. One `focal_c80` per line. |
| **Writer** | [`prepare.R`](../prepare.R), written only when some `is_focal == TRUE` focals lack a per-focal neighbor TSV. Removed when all are present. |
| **Reader** | [`build_neighbor_lists.sh`](../build_neighbor_lists.sh). |

### `clusters_80_updated` (local copy)

| | |
|---|---|
| **Path** | `{data_dir}/{species_id}/clusters_80_info_updated.tsv` |
| **Schema** | Same as `data.clusters_80_updated` input (see section 1). |
| **Writer** | [`prepare.R`](../prepare.R), seeded once (`overwrite = FALSE`). |
| **Reader** | [`pipeline.R`](../pipeline.R) via [`load_c80_tables`](../R/midas.R) - provides per-c80 `genome_prevalence` + `neighbor_c80_length_coarse` on `gene_to_c80`. |

### `catalog_genes_info` - unified gene -> c80 map

| | |
|---|---|
| **Path** | `{proj_dir}/step1_setup/catalog_genes_info.tsv` |
| **Header** | Yes. |
| **Writer** | [`build_genome_catalog.py`](../build_genome_catalog.py), rebuilt every run. |
| **Reader** | [`focal_neighbor_list.sh`](../scripts/focal_neighbor_list.sh) (per-focal gene-member lookup); [`load_c80_tables`](../R/midas.R) (drives `gene_to_c80`). |

| Col | Name |
|---|---|
| 1 | `gene_id` |
| 2 | `centroid_80` |
| 3 | `gene_length` |

Union across every entry in `sources:`. No deduplication (cross-source `genome_id` collisions are caught by the dup check in `build_genome_catalog`).

### `catalog_genome_toc` - genome -> `.genes` path map

| | |
|---|---|
| **Path** | `{proj_dir}/step1_setup/catalog_genome_toc.tsv` |
| **Header** | Yes. |
| **Writer** | [`build_genome_catalog.py`](../build_genome_catalog.py), rebuilt every run. |
| **Reader** | [`focal_neighbor_list.sh`](../scripts/focal_neighbor_list.sh) (`.genes` path lookup per genome_id). |

| Col | Name |
|---|---|
| 1 | `genome_id` |
| 2 | `genes_file_path` (absolute) |

One row per unique `genome_id` across all sources. Duplicate `genome_id` across sources is a build-stop error.

### `neighbor_list/<focal_c80>.tsv` - per-focal neighbor TSV

| | |
|---|---|
| **Path** | `{data_dir}/{species_id}/list_of_neighbors/<focal_c80>.tsv` |
| **Header** | **No** (deliberately, for awk consumers; `cols_neighbors` in [neighbor.R](../R/neighbor.R) assigns names on read). |
| **Writer** | [`get_neighbor.sh`](../scripts/get_neighbor.sh) via [`focal_neighbor_list.sh`](../scripts/focal_neighbor_list.sh) ([`build_neighbor_lists.sh`](../build_neighbor_lists.sh) fans the parallel xargs). |
| **Reader** | [`load_gene_neighbors`](../R/neighbor.R) (Step 1). |
| **Idempotency** | `-s "$outfile"` skip - a non-empty existing TSV is left alone. |

| Col | Name | Source |
|---|---|---|
| 1 | `gene_member` | The gene_id (member of the focal centroid_80) the flank was extracted around. |
| 2 | `neighbor_gene_id` | From the `.genes` file (col 1). |
| 3 | `neighbor_contig_id` | From `.genes` col 2. |
| 4 | `neighbor_gene_start` | From `.genes` col 3. |
| 5 | `neighbor_gene_end` | From `.genes` col 4. |
| 6 | `neighbor_gene_strand` | From `.genes` col 5. |
| 7 | `neighbor_gene_type` | From `.genes` col 6 (`CDS`/`rRNA`/`tRNA`/`tmRNA`). |

For each gene_member, this is its row + up to `+/-n_genes` flanking rows from the same contig (`grep -C n_genes`). One block of <= `2*n_genes + 1` rows per gene_member. Multiple gene_members per file (one per genome carrying the focal).

---

## 4. Pipeline outputs - produced by Step 1 onward

These are the model.R `# Output` targets: every file written by Steps 1-6 of `pipeline.R` (caches and final outputs alike). Detailed input/output specs live in **[STEPS.md](STEPS.md)** - one section per step. Quick reference:

| Step | Output target key | On-disk path (relative to `{proj_dir}/`) | Header? | Brief |
|---|---|---|---|---|
| 1 | `neighbor_groups_rds` | `step2_neighbors/neighbor_groups.RDS` | RDS | Cached `gene_neighbors` - one row per (focal, genome, neighbor position). |
| 1 | `short_gene_prevalence` | `step1_setup/short_gene_prevalence.tsv` | RDS | Per-synthetic-c80 prevalence for unannotated short ORFs. |
| 1 | `c80_variants_mapping` | `step1_setup/c80_variants_mapping.tsv` | RDS | `(neighbor_c80_coarse, length) -> neighbor_c80_fine` mapping. |
| 1 | `neighbor_groups_by_genome` | `step2_neighbors/01_neighbor_by_genome/` | per-file TSV | Per-focal-per-genome shards. |
| 1 | `neighbor_groups_by_focal` | `step2_neighbors/02_neighbor_by_focal/` | per-file RDS | Per-focal RDS dumps post-filter. |
| 1 | `neighbor_figures` | `step5_figures/01_neighbor_by_focal/` | PDFs | `fig1`-`fig5` per focal. |
| 2 | `path_df` | `step3_path/path_df.rds` | RDS | One row per per-genome maximal path. |
| 2 | `esupport_df` | `step3_path/esupport_df.rds` | RDS | Per-edge support counts. |
| 3 | `canonical_paths` | `step3_path/canonical_paths_coarse.tsv` | yes | L1 - one row per canonical operon. Primary key `uid`. |
| 3 | `canonical_paths_fine` | `step3_path/canonical_paths_fine.tsv` | yes | L2 - one row per length-variant isoform. Primary key `uid_fine`. |
| 3 | `canonical_paths_per_genome` | `step3_path/canonical_paths_per_genome.tsv` | yes | L3 - one row per (canonical, contributing genome). |
| 3 | `canonical_paths_c80s` | `step3_path/canonical_paths_c80s.tsv` | yes | L1 per-gene rows; consumed by Steps 4-6. |
| 3 | `canonical_paths_fine_c80s` | `step3_path/canonical_paths_fine_c80s.tsv` | yes | L2 per-gene rows with truncation / fragmentation flags. |
| 6 | `rep_path_df` | `step6_blocks/representative_path.tsv` | yes | Non-redundant trait-associated blocks (`block_uid`). |
| 6 | `uid_path_df` | `step6_blocks/rep.tsv` | yes | Per-genome attribution for reps. |
| 6 | - | `step6_blocks/rep_heatmap.pdf` | PDF | Block x genome presence/absence (when >= 2 reps and matrix >= 3x3). |
| 4 | `parse_coarse_summary` | `step4_parse/coarse_recurring_operons.tsv` | yes | One row per `uid`. |
| 4 | `parse_fine_summary` | `step4_parse/fine_isoform_priorities.tsv` | yes | One row per `uid_fine`. |
| 4 | `parse_selected_coarse` | `step4_parse/selected_coarse.tsv` | yes | Coarse uids whose fine isoforms passed survival. |
| 4 | `parse_selected_fine` | `step4_parse/selected_fine.tsv` | yes | Fine isoforms surviving `n_fine_genomes >= ceil(path_min_genomes * fine_coverage_ratio)`. |
| 4 | `parse_fine_long` | `step4_parse/fine_long.tsv` | yes | Long-format per-gene table for the sampled isoforms. |
| 4 | `parse_genome_paths_dir` | `step4_parse/genome_paths/fine_<uid_fine>_<genome>.tsv` | per-file | Bare gene-id lists per sampled (isoform, genome). |
| 5 | `parse_coarse_figures` | `step5_figures/coarse_operons_<fill_by>.pdf` + `02_by_component_coarse/` | PDFs | Global + per-component coarse plots, one per `fill_mode`. |
| 5 | `parse_fine_figures` | `step5_figures/fine_operons_<fill_by>.pdf` + `03_by_component_fine/` | PDFs | Global + per-component fine plots, one per `fill_mode`. |

For column lists on each Step 1-6 output, see [STEPS.md](STEPS.md) and the roxygen docstrings on the corresponding `run_stepN_*` functions.

---

## 5. MWAS (parked) targets

Declared in [model.R](../R/model.R)'s `# MWAS (parked)` block. **Not read by the current pipeline** - placeholders for future MWAS re-integration. Listed here for completeness.

| Target | Path | Format |
|---|---|---|
| `genes_info` | `{midasdb}/pangenomes/{species_id}/genes_annotated.tsv` | UHGG genes_annotated (gene_id, centroid_99..75, gene_length, marker_id, gene_is_phage/plasmid/amr/me) |
| `clusters_80` | `{midasdb}/pangenomes/{species_id}/clusters_80_info.tsv` | Older c80 metadata sibling of `clusters_80_updated` |
| `gene_by_sample_matrix` | `{midas_dir}/{species_id}/gene_by_sample_matrix.rds` | MIDAS merge output |
| `genes_to_heatmap` | `{midas_dir}/{species_id}/genes_to_heatmap.rds` | MIDAS merge output |
| `GRM_pop` | `{midas_dir}/{species_id}/GRM.rds` | MWAS GRM |
| `pca_pop` | `{midas_dir}/{species_id}/pca_df.rds` | MWAS PCA |
| `corrected_genes` | `{data_dir}/corrected_genes_0.01.RDS` | Legacy MWAS scoring input |

---

## Where each path comes from

Three layers of indirection feed `get_target("<key>")`:

| Layer | Source | Examples |
|---|---|---|
| **YAML** (`cfg_get`) | User-provided in `<config.yaml>` | `data.focal_meta`, `data.clusters_80_updated`, `data.midasdb_dir`, `data.data_dir`, `job.proj_dir`, `job.species_id` |
| **model.R `target_layout()`** (`get_target`) | Pipeline-internal path scheme (relative to `proj_dir` or absolute) | every key in this doc |
| **Script-side derivation** | When neither cfg_get nor get_target applies, derived in code (e.g. genome_id from gene_id) | see CLAUDE.md keystone #4 for the genome_id contract |

When updating a path, the rule of thumb is: **change `model.R` for path conventions, change the YAML for external file locations.** Never hardcode a path under `step{N}_*/` in helper code - use `get_target()` and add a key to `target_layout()` if needed.
