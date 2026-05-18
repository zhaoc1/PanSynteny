# Strain-Aware Operon Pipeline — User Guide

A high-level orientation to the pipeline: what it does, how to run it, and how to read the outputs. For per-step input / output / logic details, see [STEPS.md](STEPS.md). For column-level schemas of every file the pipeline reads or writes, see [SCHEMA.md](SCHEMA.md). For internal architecture notes (file targets, function-level docstrings), read the source files directly — they are heavily commented.

---

## What this pipeline does

Given a user-curated focal-gene table (centroid_80 clusters from MIDAS) for one species, the pipeline reconstructs **de novo recurring operons** those genes live in and harmonizes them across the strains in the species-level pangenome. Within each operon, it also extracts the contiguous focal-gene blocks (grouped by `focal_label` direction).

The output is three things:

1. **Canonical operons** at three granularity levels (coarse cluster path → length-variant isoforms → per-genome instances), each with attached trait statistics, small-ORF flags, and truncation/fragmentation flags.
2. **Trait-associated blocks**: the non-redundant runs of trait-correlated genes within those operons, ranked per locus.
3. **Per-genome attribution**: which strains carry which operon variant and which trait block.

The driver is [`pipeline.R`](../pipeline.R). Steps 1–4 produce the analytical outputs (summaries, selection sets, BLAST gene lists). Step 5 renders gggenes figures (`run_step5_figures`). Step 6 extracts trait-associated blocks (`run_step6_blocks`, gated by `blocks.skip`).

---

## How to run

The full workflow is **four ordered commands**, all reading the same `<config.yaml>`. Everything runs under the `strain-aware-operon` conda env (see [SETUP.md](SETUP.md)).

### Pre-flight checklist

Sanity-check the run inputs before launching:

| Check | Why |
|---|---|
| `conda activate strain-aware-operon` | every script (R + Python) uses this env's interpreters |
| `data.focal_meta` resolves to an existing file | prepare.R aborts at startup otherwise |
| focal_meta header contains `focal_c80, focal_label, is_focal, gene_label` | minimum required schema (see "data" / "focal_meta column requirements") |
| `data.clusters_80_updated` resolves to an existing file (or local copy already exists) | prepare.R seeds the local copy from this on first run; required if the local copy isn't there yet |
| Each `sources:` entry's `genes_info` and `genomes_dir` resolve to existing paths | build_genome_catalog aborts on a missing dir/file |
| `job.proj_dir` and `data.data_dir` exist or can be created (writable) | snapshot, cache, catalog, and neighbor TSVs land under these |

#### Recommended `input_dir/` layout

Treat `input_dir/` as the self-contained bundle of *everything you authored for this run* — typically two files:

```text
input_dir/
├── my_run.yaml          # the run config (the CLI arg)
└── focal_table.tsv      # focal_meta TSV; reference it via focal_meta: "{input_dir}/focal_table.tsv"
```

Co-locating them means a single `cp -r input_dir/ archived/` saves the full inputs; together with the per-run `step1_setup/run_config.yaml` snapshot, the run is reproducible. The YAML itself can technically live anywhere (you pass its path as the CLI arg), but keeping it next to `focal_meta` is the cleanest convention.

### What's cached + what triggers a rebuild

Most reruns are cheap because the expensive artifacts are `-s`-guarded. Knowing which caches gate which step:

| Cache | Owner | Triggers rebuild |
|---|---|---|
| `{proj_dir}/step1_setup/catalog_{genes_info,genome_toc}.tsv` | `build_genome_catalog` | Always rewritten on each run (cheap apart from the UHGG `genes_info.tsv` scan). Per-run (under proj_dir), so changing `sources:` between proj_dirs doesn't clobber. If you suspect schema drift from an older version (e.g. 2-col, no headers), `rm` both files first. |
| `<genomes_dir>/<g>/<g>.genes` (prokka) | `build_genome_catalog` via `gff_to_genes` | Only re-derived if the `.genes` file is missing or empty. To force re-conversion: `rm <genomes_dir>/<g>/<g>.genes`. |
| `{proj_dir}/step1_setup/run_config.yaml` | `prepare.R` | Always overwritten (config snapshot). |
| `{proj_dir}/step1_setup/gene_meta_full.tsv` | `prepare.R` (the `focal_meta` target) | Always rewritten — re-derived from `data.focal_meta` + `prepare.*` thresholds. |
| `{data_dir}/{species_id}/clusters_80_info_updated.tsv` | `prepare.R` (the `clusters_80_updated` target) | **Seeded once** from `data.clusters_80_updated` (`overwrite = FALSE`). Subsequent runs preserve any hand edits to the local copy. To refresh from the source: `rm` the local file then re-run prepare.R. |
| `{proj_dir}/step1_setup/gene_list.tsv` | `prepare.R` | Written only if some focals still lack neighbor TSVs; removed once all are present. |
| `{data_dir}/{species_id}/list_of_neighbors/<focal_c80>.tsv` | `run_species.sh` | Only materialised if missing/empty. Force re-extraction by deleting individual files (or the whole directory). |
| `step2_neighbors/neighbor_groups.RDS` | pipeline.R Step 1 | **Cache gate for Step 1** — delete to force Step 1 to re-run. |
| `step3_path/path_df.rds` | pipeline.R Step 2 | **Cache gate for Step 2** — delete to force Step 2 to re-run. |
| `step3_path/canonical_paths*.tsv`, `step4_parse/*`, `step5_figures/*`, `step6_blocks/*` | pipeline.R Steps 3-6 | Always re-run; nothing to delete. |

### The four commands

```bash
# Step 0a — build the unified genome catalog
python  build_genome_catalog.py <config.yaml>

# Step 0  — snapshot the YAML; process focal_meta into the step1 cache;
#           enumerate any missing per-focal neighbor TSVs
Rscript prepare.R               <config.yaml>

# Step 0  — materialise the missing per-focal neighbor TSVs
bash    run_species.sh          <config.yaml>

# Steps 1-6 — analytical pipeline
Rscript pipeline.R              <config.yaml>
```

Working example config: [`example.yaml`](../example.yaml) (template). A real worked-example input bundle (config + focal_meta TSV) lives under [`examples/`](../examples/).

### Tip — pin the env's `python` / `Rscript` in shell variables

If `conda activate strain-aware-operon` doesn't behave (`which python` still pointing at the base / system install, or running from an IDE terminal that doesn't load your shell init), you can skip activation entirely by calling the env's interpreters directly:

```bash
PY=/pollard/home/czhao/miniconda3/envs/strain-aware-operon/bin/python
RSC=/pollard/home/czhao/miniconda3/envs/strain-aware-operon/bin/Rscript

$PY  build_genome_catalog.py <config.yaml>
$RSC prepare.R               <config.yaml>
bash run_species.sh          <config.yaml>
$RSC pipeline.R              <config.yaml>
```

This sidesteps `conda activate` and removes any ambiguity about which interpreter is in use. The bash scripts (`run_species.sh`, `generate_neighbor_list.sh`, `get_neighbor.sh`) call `python3` internally for their YAML helper — if `python3` doesn't resolve to the env, prepend `PATH=$(dirname $PY):$PATH` to those bash invocations, or activate the env normally.

For the R scripts: invoking the env's `Rscript` directly is enough when `LD_LIBRARY_PATH` already contains the env's lib dir (avoids the `GLIBCXX_3.4.30 not found` failure from `vroom`/`readr`). If you hit that error, prepend:

```bash
export LD_LIBRARY_PATH=/pollard/home/czhao/miniconda3/envs/strain-aware-operon/lib:$LD_LIBRARY_PATH
```

(Activation does this for you automatically — direct invocation may not.)

**Step 0a — `build_genome_catalog.py`.** Reads `sources:` from the YAML; for each declared source, normalises a membership file (`gene_id → centroid_80, gene_length`) and a genome → `.genes` path map. Output:

- `{proj_dir}/step1_setup/catalog_genes_info.tsv` — `gene_id <TAB> centroid_80 <TAB> gene_length` (union across all sources)
- `{proj_dir}/step1_setup/catalog_genome_toc.tsv` — `genome_id <TAB> genes_file_path` (union across all sources)

Catalog lives **per-run** under `proj_dir` so different `sources:` lists across iterations don't clobber a shared catalog. The expensive `.genes` files for prokka sources stay under `<genomes_dir>` and remain shared across runs via the `-s` guard — only the union/normalisation step is rerun per proj_dir.

For `type: prokka` sources it converts each `<g>.gff → <g>.genes` **in place** next to the GFF (idempotent on `-s`). For `type: midas` sources it trusts the midasdb's existing `.genes` files. **Duplicate `genome_id` across sources → warn and stop** — the genome → path map must be unambiguous.

**Step 0 — `prepare.R`.** Four responsibilities:

1. Snapshot `<config.yaml>` to [`run_config`](../model.R) (`{proj_dir}/step1_setup/run_config.yaml`) for provenance.
2. Seed the local copy of `clusters_80_updated` from `data.clusters_80_updated` into [`clusters_80_updated`](../model.R) on first run (`overwrite = FALSE`); subsequent runs preserve hand edits. To refresh from the source, `rm` the local copy.
3. Read the user-provided focal table from `data.focal_meta` (path in the YAML; `{proj_dir}` / `{species_id}` / `{midasdb_dir}` placeholders expanded), optionally apply the `|score_col|` thresholds, and cache the result to [`focal_meta`](../model.R) (`step1_setup/gene_meta_full.tsv` — filename kept for back-compat).
4. Walk every `is_focal == TRUE` centroid, check whether `<focal_c80>.tsv` already exists under [`neighbor_list`](../model.R), and write any missing centroids to [`gene_list`](../model.R) as a one-per-line list. If everything is present, `gene_list` is removed and a "Ready to run pipeline.R" message is printed.

Always overwrites — cheap to re-run.

**Step 0 — `run_species.sh`.** Consumes `gene_list.tsv` (the missing-list from `prepare.R`); fans `generate_neighbor_list.sh` over each focal in parallel. `generate_neighbor_list.sh` joins the catalog `genes_info` (gene members of the focal) to `genome_toc` (each genome's `.genes` path) and calls `get_neighbor.sh` per gene member. Output: `{data_dir}/{species_id}/list_of_neighbors/<focal_c80>.tsv` (7 cols, no header). Per-focal idempotency via `-s "$outfile"`. No-op if `gene_list.tsv` is absent.

**Steps 1–4 — `pipeline.R`.** Reads the cached focal table that `prepare.R` produced. The driver consumes `focal_c80_df` as-is and does **not** apply any threshold of its own — that decision is owned by `prepare.R`. If the focal_meta cache is missing or any `is_focal` centroid still lacks a neighbor TSV, the driver aborts at startup with a pointer back to `prepare.R`.

**Re-run skipping.** On a re-run after partial completion, the driver skips work whose cached output exists: Step 1 skips re-extraction if `neighbor_groups_rds` exists, Step 2 skips if `path_df` exists. **To force a re-run of a step, delete its cache file.** `build_genome_catalog` rebuilds the catalog files fresh every run but skips the per-genome prokka conversion via `-s`.

---

## Configuration (YAML)

All knobs live in one YAML file. Scalar sections are flattened into a single `job_config` namespace at load time, so any key from any section is accessible via `cfg_get(job_config, "<key>")`. **Exception:** `sources:` is a list (not a scalar section); `config.R` skips it, so it does not pollute `job_config`. It is consumed only by `build_genome_catalog.py`. See [`config.R`](../config.R) for the loader.

### `job` — required

```yaml
job:
  species_id:    "102321"             # MIDAS species id (numeric)
  proj_dir:      "/path/to/results"   # output root, used as-is (include species_id in the path for per-species isolation, e.g. "/path/to/results/102321")
  input_dir:     "/path/to/inputs"    # required: holds the YAML config + the focal_meta TSV for this run; usable as {input_dir} anywhere in the YAML
  parallel_jobs: 2                    # required: -P for run_species.sh's xargs fan-out across focals (typical: 2)
```

`trait` was removed in this version (unused). All four `job:` keys are **required** — `load_job_config` aborts with a clear error if `input_dir` or `parallel_jobs` is missing. The recommended convention is to keep the YAML and the `focal_meta` TSV side-by-side under `input_dir/` — see [Recommended `input_dir/` layout](#recommended-input_dir-layout) above.

`parallel_jobs` controls how many focals `run_species.sh` extracts in parallel (xargs `-P`). Cap at ~8 on shared filesystems — `get_neighbor.sh` is I/O-bound on `.genes` reads, so over-parallelising hammers storage for everyone.

### `data` — required

External data paths and the catalog flank size.

```yaml
data:
  midasdb_dir:         "/path/to/midas2db-uhgg-v2"              # MIDAS reference DB root
  data_dir:            "/path/to/data"                          # catalog + neighbor TSVs land under here
  focal_meta:          "{input_dir}/focal_table.tsv"            # user-provided focal-gene metadata; prepare.R reads from here
  clusters_80_updated: "{midasdb_dir}/pangenomes/{species_id}/clusters_80_info_updated.tsv"  # prepare.R seeds the local copy on first run
  n_genes:             20                                       # max neighbours per side around a focal (get_neighbor.sh / grep -C)
```

(`input_dir` lives in `job:` alongside `proj_dir` — see [`job`](#job--required) above.)

Both `focal_meta` and `clusters_80_updated` are paths to TSVs the user provides. `prepare.R` snapshots them into the project directory: `focal_meta` is processed into `step1_setup/gene_meta_full.tsv` on every run, while `clusters_80_updated` is seeded once into `<data_dir>/<species_id>/clusters_80_info_updated.tsv` (`overwrite = FALSE`) and then left alone — so hand edits (e.g. appending ECOR-derived c80 rows) persist across reruns. To refresh from the source, `rm` the local copy and re-run prepare.R. `n_genes` is the flank size for `get_neighbor.sh`. `midas_dir` and `df_dir` were removed from this section: `midas_dir` lives under the parked `mwas:` section (see below); `df_dir` was dropped entirely (unused).

**Path placeholders.** `{proj_dir}`, `{species_id}`, `{midasdb_dir}`, and `{input_dir}` are expanded in every YAML string scalar at load time (see `load_job_config` in [config.R](../config.R)). Use them anywhere in the YAML to avoid repeating long paths. Plain absolute paths still work unchanged.

- `{proj_dir}` — the **output** root. Used as-is; include `<species_id>` in the YAML value if you want per-species isolation.
- `{input_dir}` — **input** root for the YAML config + the focal_meta TSV (the two user-authored files for this run). Keeps inputs separate from outputs. Declared under `job.input_dir` (required).
- `{species_id}` / `{midasdb_dir}` — convenience for paths derived from the species and the MIDAS DB.

#### `focal_meta` column requirements

| Column | Required? | Used by |
|---|---|---|
| `focal_c80` | **required** | Step 1 onward — the centroid_80 id, drives per-focal extraction. |
| `focal_label` | **required** | Step 1 sharding (`pos`/`neg` dir), Step 5 (`focal_label` fill mode). |
| `is_focal` | **required** (unless `prepare.score_col` is set, in which case it can be derived) | Step 1 gate — only `TRUE` rows drive neighbor extraction. |
| `gene_label` | **required** | Step 5 (`gene_label` fill mode) — your user-defined annotation/category for the focal gene. Distinct from the `.genes`-file `gene_type` (which is the GFF feature type, CDS/rRNA/etc., for neighbors). |
| `cor_to_b` | optional | Step 5 (`cor_to_b` fill mode); also drives the threshold path when `prepare.score_col = "cor_to_b"`. |
| `beta` | optional | Step 5 (`beta` fill mode). |
| `sample_prevalence` | optional | Step 5 (`sample_prevalence` fill mode). |
| `trait`, `genome_counts` | optional | Carried through for traceability; no live consumer in the current pipeline. |

Extra columns beyond this list pass through to the cache unchanged. **`is_focal` is authoritative in the input** unless `prepare.score_col` is set — in which case the threshold derivation overwrites it (loudly, with a warning).

### `mwas` — legacy / parked

```yaml
mwas:
  midas_dir:    ""                                       # gene-by-sample matrices root (MIDAS merge output)
```

The `mwas:` section feeds the parked MWAS block in [model.R](../model.R) — `gene_by_sample_matrix`, `genes_to_heatmap`, `GRM_pop`, `pca_pop`. **None of these targets are read by the current operon pipeline.** Leave `midas_dir` empty unless you're wiring up the MWAS layer. The section exists so the plumbing is in place for re-integration without re-creating the path keys.

### `sources` — required (consumed by `build_genome_catalog`)

A list, not a section. Each entry declares one genome source:

```yaml
sources:
  - name:        uhgg                                    # free label (for logs)
    type:        midas                                   # midas | prokka
    genes_info:  "{midasdb_dir}/pangenomes/{species_id}/genes_info.tsv"  # membership file; placeholders expanded
    c80_col:     6                                       # 1-based column of centroid_80 in `genes_info`
    length_col:  8                                       # 1-based column of gene_length (default 8 for midas = UHGG; required for prokka)
    genomes_dir: "{midasdb_dir}/gene_annotations/{species_id}"           # dir of <genome>/<genome>.genes
  - name:        ecor
    type:        prokka                                  # build_genome_catalog will derive .genes from <g>.gff
    genes_info:  "/abs/path/to/ecor_to_centroid80.tsv"
    c80_col:     2
    genomes_dir: "/abs/path/to/ecor_prokka"              # dir of <genome>/<genome>.gff (+ <g>.genes after conversion)
```

`{species_id}`, `{midasdb_dir}`, and `{proj_dir}` placeholders are expanded against `job:` / `data:` at build time. `length_col` is **required for prokka sources** (their membership schema is user-controlled, so the column index must be declared explicitly); for midas it defaults to 8 (UHGG schema). UHGG-only runs drop the second entry; multi-source runs append more.

**Per-source local placeholders.** Inside one `sources:` entry, any extra string field (anything other than `name` / `type` / `genes_info` / `c80_col` / `length_col` / `genomes_dir`) is treated as a local placeholder usable inside that same entry's `genes_info` / `genomes_dir`. Useful to factor out a repeated base path:

```yaml
- name:        ecor
  type:        prokka
  ecor_dir:    "/pollard/data/.../2026-05-11-ecor72"      # local placeholder
  genes_info:  "{ecor_dir}/ecor_gene_centroid80.tsv"
  c80_col:     2
  length_col:  3
  genomes_dir: "{ecor_dir}/ecor_prokka"
```

Local placeholders are scoped to their own source entry — `{ecor_dir}` defined under `ecor` is not visible to other sources. Globals (`{species_id}`, `{midasdb_dir}`, `{proj_dir}`) are expanded inside the local value too, so `ecor_dir: "{proj_dir}/inputs"` works.

### `prepare` — Step 0 (`prepare.R`)

```yaml
prepare:
  score_col:        cor_to_b   # column in focal_meta used to score genes; one of {cor_to_b, beta}
  inclusion_cutoff: 0.25       # min |score_col| for any gene to be retained in the cache
  focal_cutoff:     0.5        # min |score_col| for a gene to be marked is_focal = TRUE
```

If `score_col` is **empty**, `prepare.R` passes `focal_meta` through unchanged — the input *must* then already contain an `is_focal` column. If `score_col` is set, prepare.R filters to rows where `|focal_meta[[score_col]]| >= inclusion_cutoff` and sets `is_focal = |focal_meta[[score_col]]| >= focal_cutoff`. The constraint `focal_cutoff >= inclusion_cutoff` is enforced by usage. Note that `cor_to_b ∈ [-1, 1]` while `beta` is unbounded — when switching `score_col`, scale the cutoffs accordingly.

**`is_focal` overwrite warning.** When `prepare.score_col` is set, the threshold derivation **overwrites** any `is_focal` column already in the input — prepare.R prints a `warning()` so it's loud. If your `focal_meta` carries hand-curated `is_focal` values (e.g. a mix of focal rows and context rows with `is_focal = FALSE` for metadata only), keep `prepare.score_col: ""` so the input wins.

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

(The Step 0 parallelism knob `parallel_jobs` lives under `job:` — see [`job` — required](#job--required) above.)

### `path` — Step 3 (canonical operon consolidation)

```yaml
path:
  path_min_genomes:  20    # min strains a canonical operon must be backed by
  truncation_cutoff: 0.8   # length / centroid_length below which a gene counts as truncated
```

`path_min_genomes` is the **per-canonical-operon** survival cut applied in [`generate_canonical_path()`](../graph.R#L538). Conceptually distinct from Step 1's `focal_min_genomes`: Step 1 asks "is this neighborhood pattern recurring around a focal?" and Step 3 asks "does this canonical operon recur across enough strains?" — they need not have the same value, though they typically do.

### `blocks` — Step 6 block extraction (`blocks.R`)

```yaml
blocks:
  skip_block:  false   # set true to skip Step 6 entirely (pipeline.R short-circuits the call)
  allow_gaps:  2       # max non-hit positions allowed inside one block (aggregate_blocks)
  min_overlap: 1       # min LCS length to call two blocks "overlap" (rank_block_representatives)
  min_shared:  2       # diagnostic threshold for substring-overlap pairs (diagnose_rep_overlaps)
```

`blocks.skip_block` is the **Step 6 toggle**. Setting it to `true` skips the entire Step 6 (block extraction) call in pipeline.R; nothing is written under `step6_blocks/`. Steps 4 (parse) and 5 (figures) do not depend on block-extraction outputs.

`min_shared` is diagnostic-only — it does not change `representatives.tsv` or `rep.tsv`, only the `Rep overlap diagnostic` log line.

### `parse` — Step 4 + Step 5 (summaries, sampling, BLAST gene lists, gggenes figures)

```yaml
parse:
  fine_coverage_ratio: 0.25  # n_fine_genomes >= ceiling(path_min_genomes * ratio) to survive
  seed:                616   # RNG seed for the per-fine-isoform exemplar-genome draw
  fill_modes:                # one PDF per mode; subset of:
    - beta
    - sample_prevalence
    - cor_to_b
    - fill_gene
    - gene_label
    - focal_type
    - focal_label
```

`fine_coverage_ratio` sets the fraction of the path-level genome floor each fine isoform must reach. `ceiling` semantics are preserved — only the multiplier is tunable. Note the ratio is applied to `path_min_genomes` (Step 3's gate), not to `focal_min_genomes`. Historical default 0.5 ("half-coverage rule"); examples currently use 0.25 for permissive fine-isoform retention.

**`fill_modes` is column-tolerant.** Any mode whose backing column is absent from the c80s tables is **skipped with a warning** (Step 5 prints `Step 5: skipping fill_mode '<m>' — no backing column found …`). For example, on a minimal focal_meta with no trait stats, listing `beta` / `cor_to_b` / `sample_prevalence` is harmless — those modes just no-op. `fill_gene` is always available (it's derived, not column-backed); `gene_label` / `focal_type` / `focal_label` need their corresponding column to be present on the c80s tables (they propagate from focal_meta via the build-canonical-paths-c80s join).

### `plot` — Step 1 + Step 5 figure layout

```yaml
plot:
  gene_padding_bp: 100   # bp gap between adjacent genes in gggenes layouts
```

Used by both the Step 1 diagnostic figures ([`extract_gene_neighbor_patterns`](../neighbor.R#L324)) and the Step 5 publication figures ([`.layout_operon_tracks`](../plot.R#L245)) so spacing is consistent across all gggenes outputs.

---

## Two common setups

The `sources:` list is the one knob that changes when you add or drop a genome source. Two recipes cover ~all real-world runs.

### Case 1 — UHGG only

Standard MIDAS-style analysis: one species, one genome catalog, no external panels.

```yaml
sources:
  - name:        uhgg
    type:        midas
    genes_info:  "{midasdb_dir}/pangenomes/{species_id}/genes_info.tsv"
    c80_col:     6
    genomes_dir: "{midasdb_dir}/gene_annotations/{species_id}"
```

**What happens:**

- `build_genome_catalog.py` reads only this entry. `{species_id}` and `{midasdb_dir}` are expanded against the `job:` / `data:` sections of your YAML.
- Every gene in the species' UHGG pangenome contributes one row to `genes_info.tsv` (gene_id, centroid_80, gene_length pulled from col 8 of the source).
- Every UHGG genome that has at least one gene in that pangenome contributes one row to `genome_toc.tsv`, pointing at the midasdb's existing `<genome>.genes` file. **No `.gff → .genes` conversion runs** (midasdb already ships `.genes`).
- Step 1 onward sees one source's worth of neighbors per focal — exactly the operons UHGG knows about for that species.

**Inputs you need on disk:**

- A working MIDAS reference DB at `data.midasdb_dir` (with `pangenomes/<species_id>/` and `gene_annotations/<species_id>/<genome>/<genome>.genes`).
- Your focal-table TSV at `data.focal_meta`.

**Run command sequence is the same as always:**

```bash
python  build_genome_catalog.py <config.yaml>
Rscript prepare.R               <config.yaml>
bash    run_species.sh          <config.yaml>
Rscript pipeline.R              <config.yaml>
```

### Case 2 — UHGG + Prokka (external panel)

You have an additional genome collection (e.g. ECOR-72) annotated with Prokka, and you want those strains to contribute to operon discovery alongside UHGG.

```yaml
sources:
  - name:        uhgg
    type:        midas
    genes_info:  "{midasdb_dir}/pangenomes/{species_id}/genes_info.tsv"
    c80_col:     6
    genomes_dir: "{midasdb_dir}/gene_annotations/{species_id}"
  - name:        ecor
    type:        prokka
    genes_info:  "/abs/path/to/ecor_gene_centroid80.tsv"  # ECOR gene -> UHGG centroid_80 + gene_length
    c80_col:     2
    length_col:  3                                       # required for prokka — gene_length column
    genomes_dir: "/abs/path/to/ecor_prokka"              # one <genome>/<genome>.gff per ECOR strain
```

**What happens:**

- `build_genome_catalog.py` walks both sources. The output `catalog_genes_info.tsv` is the **union** of UHGG + ECOR membership rows; `catalog_genome_toc.tsv` is the union of both sources' genomes.
- For each ECOR genome:
  - The script finds `<genome>/<genome>.gff` under `genomes_dir`.
  - If `<genome>/<genome>.genes` is missing/empty, [`gff_to_genes.py`](../gff_to_genes.py) is called in-process to derive it (uses `gffutils`). The file lands **next to the `.gff`** — no separate output dir.
  - Idempotent on `-s "<g>.genes"`: a non-empty file is left alone on re-run.
  - `gene_length` is **read straight from the membership file** at `length_col`. The .genes conversion no longer scans for length — the file exists only so downstream `get_neighbor.sh` can read it.
- **Duplicate `genome_id` across sources → the build stops** with an error naming which sources collided. The genome → path map must be unambiguous.
- Step 1 onward sees one merged neighbor table. ECOR neighbor genes get their c80 labels via the same `gene_to_c80` join `load_c80_tables` always uses — no special-casing.

**Inputs you need on disk:**

- Everything from Case 1, plus:
- A Prokka annotation per ECOR genome at `<ecor_genomes_dir>/<genome>/<genome>.gff` (`.fna`/`.faa` are not needed by this pipeline).
- A membership TSV mapping each ECOR gene to a UHGG `centroid_80` **with `gene_length` populated** (e.g., BLAST result post-processed; typical schema: `gene_id, centroid_80, gene_length, ...`). The easiest way to populate `gene_length` is to merge in the prokka `.tsv`'s `length_bp` column on `locus_tag` — that's already `end - start + 1` from the GFF. Column 1 = gene_id (must match the prokka `locus_tag`); centroid_80 column at `c80_col`; gene_length column at `length_col`.

**Notes specific to Prokka sources:**

- The gene_id derivation (`genome_id_from_gene_id` — strip trailing `_NNNNN`) must produce `genome_id` values that match the **directory name** under `genomes_dir`. For ECOR's `GCF_900448275.1_00001` → `GCF_900448275.1`, this matches `<ecor_genomes_dir>/GCF_900448275.1/` ✓. If you bring in a new source whose id scheme doesn't follow this, extend the derivation in both `build_genome_catalog.py` (`genome_id_from_gene_id`) and `generate_neighbor_list.sh` (the awk join).
- ECOR genes can map to centroid_80s from **any** species (not just the run's `species_id`). After `load_c80_tables` joins those to `clusters_80_updated` (which is species-scoped), off-species rows carry the correct `centroid_80` and `gene_length` but their `neighbor_c80_length_coarse` / `genome_prevalence` come through as `NA`. Downstream truncation / fragmentation flags handle this NA-tolerantly.
- The first prokka run does the conversion (~5-15 s per genome via `gffutils`). Subsequent runs skip already-converted `.genes` files; budget time accordingly on first runs.

### Switching between cases

Drop or add `sources:` entries any time. The catalog is rebuilt from scratch on every `build_genome_catalog.py` invocation, so switching from Case 1 → Case 2 (or vice versa) is just an edit + re-run. The `.genes` files produced by Case 2 stay on disk under `genomes_dir`; dropping back to Case 1 doesn't delete them.

Downstream, **delete the Step 1 cache** (`step2_neighbors/neighbor_groups.RDS`) when you switch — otherwise pipeline.R reuses the old neighbor table and the new ECOR contribution won't show up.

---

## Pipeline overview

| Step | What it does | Driver | Helper file(s) | Cache / output |
|---|---|---|---|---|
| **0a** | Build the unified genome catalog from `sources:`. Convert prokka `.gff → .genes` in place (idempotent on `-s`). Dup-check genome_id across sources. | [`build_genome_catalog.py`](../build_genome_catalog.py) | [`gff_to_genes.py`](../gff_to_genes.py) | `genome_catalog/{genes_info.tsv, genome_toc.tsv}`; per-prokka-genome `<g>.genes` |
| **0** | Snapshot the YAML. Read focal_meta from YAML, optionally apply `\|score_col\|` thresholds, cache to `focal_meta` target. Enumerate any missing per-focal neighbor TSVs into `gene_list.tsv`. | [`prepare.R`](../prepare.R) | `config.R`, `model.R` | `run_config.yaml`, `gene_meta_full.tsv` (the `focal_meta` cache), `gene_list.tsv` |
| **0** | Materialise the missing per-focal neighbor TSVs. | [`run_species.sh`](../run_species.sh) | `generate_neighbor_list.sh`, `get_neighbor.sh` | `list_of_neighbors/<focal_c80>.tsv` |
| **Setup** | Load `cluster_80`, `gene_to_c80` (from the catalog `genes_info`), and the cached focal table. Re-check every focal has a neighbor TSV; abort if not. | [`pipeline.R`](../pipeline.R) lines 11–79 | `config.R`, `model.R` | — |
| **1** | Per-focal neighborhood extraction → cross-genome assembly → small-ORF + length-variant labels. Orchestrated by [`run_step1_neighbor_extraction`](../neighbor.R). | `pipeline.R` lines 82–92 | `neighbor.R`, `midas.R` | `step2_neighbors/neighbor_groups.RDS` (cache) |
| **2** | Per-genome operon graphs → maximal paths. Orchestrated by [`run_step2_path_stitching`](../graph.R). | `pipeline.R` lines 95–100 | `graph.R` | `step3_path/path_df.rds`, `step3_path/esupport_df.rds` |
| **3** | Cross-genome consolidation → three granularity levels with trait stats and structural flags. Orchestrated by [`run_step3_consolidation`](../path.R). | `pipeline.R` lines 103–128 | `graph.R`, `path.R`, `parse.R` | `step3_path/canonical_paths*.tsv` (5 TSVs) |
| **4** | Summaries, fine-coverage selection, exemplar-genome sampling, BLAST gene lists. Orchestrated by [`run_step4_parse`](../parse.R). | `pipeline.R` | `parse.R` | `step4_parse/*` |
| **5** | gggenes figures: global + per-component PDFs for each fill mode. Orchestrated by [`run_step5_figures`](../plot.R). | `pipeline.R` | `plot.R` | `step5_figures/*` |
| **6** | Trait-associated block extraction + non-redundant representative ranking + per-genome attribution. Gated by `blocks.skip`. Orchestrated by [`run_step6_blocks`](../blocks.R). | `pipeline.R` | `blocks.R` | `step6_blocks/{representative_path.tsv, rep.tsv, rep_heatmap.pdf}` |

For each numbered step, [STEPS.md](STEPS.md) gives a complete input / output / logic / caveats writeup.

---

## Output reference

R outputs land under `<proj_dir>/`. Catalog + neighbor TSVs land under `<data_dir>/<species_id>/`.

### Canonical operon tables (Step 3 — the analytical core)

The same operons are emitted at three granularity levels. Every level carries trait statistics, joint-component membership, and structural decorations.

- **Level 1 (coarse)** — `step3_path/canonical_paths_coarse.tsv` (one row per operon) and `step3_path/canonical_paths_c80s.tsv` (one row per gene-in-operon). **Primary key:** `uid = "cmp{joint_component_ids}-{path_type}-{canonical_path_id}-ng{n_genomes}"`.
- **Level 2 (per-isoform)** — `step3_path/canonical_paths_fine.tsv` and `step3_path/canonical_paths_fine_c80s.tsv`. **Primary key:** `uid_fine = "{uid}-iso{rank}-ngf{n_fine_genomes}"`. Strip `-iso\d+-ngf\d+$` to recover the parent `uid`.
- **Level 3 (per-genome)** — `step3_path/canonical_paths_per_genome.tsv`. One row per `(canonical, contributing genome)` with raw and canonical-aligned `gene_path` strings.

**When to use which level**

| Question | Level |
|---|---|
| What recurring operons are trait-associated? | L1 coarse — group by `uid`. |
| Are there length variants worth distinguishing (truncation, fragmentation, tandem split-genes)? | L2 fine — `is_truncated`, `is_fragmented` are populated here only. |
| Which strains carry operon X? Which gene IDs do I BLAST? | L3 per-genome — `neighbor_genome` + `gene_path`. |

### Block tables (Step 6 — the trait answer)

- `step6_blocks/representative_path.tsv` — non-redundant trait-associated blocks. **Primary key:** `block_uid = "cmp{component}-{type}-rank{rep_rank}-nge{block_n_genes}"`.
- `step6_blocks/rep.tsv` — per-genome attribution for the reps. One row per `(block_uid, canonical_uid, neighbor_genome, left_orig, right_orig)`.
- `step6_blocks/rep_heatmap.pdf` — block × genome presence/absence.

---

## Three c80 columns: `c80`, `neighbor_c80_coarse`, `neighbor_c80_fine`

Easy to confuse. They differ in **scope** and **resolution**.

| Column | Resolution | Scope | Where it appears |
|---|---|---|---|
| `c80` | coarse | the operon's own focal cluster id | L1 c80s table only |
| `neighbor_c80_coarse` | coarse | a neighbor's MIDAS cluster id (a real centroid_80, or a synthetic `_<focal>-<type>_<rank>` for short ORFs) | everywhere |
| `neighbor_c80_fine` | length-variant-aware | same as `neighbor_c80_coarse` but with `_<rank>` suffix when the cluster has multiple observed lengths | L2 fine c80s table; produced by Step 1 ([`compute_c80_variants()`](../midas.R#L215)) |

Two practical rules:

- **For coarse grouping** (collapse all length variants of a cluster) → use `neighbor_c80_coarse`.
- **For length-sensitive analyses** (truncation, fragmentation, isoform identity) → use `neighbor_c80_fine` and the L2 fine table.

In the L1 c80s table, `neighbor_c80_coarse == c80` after the [`build_canonical_paths_c80s()`](../path.R#L265) join — they're kept as parallel columns for join-compatibility, not because they hold different information.

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
pangenome-operons-v2/
├── build_genome_catalog.py     # Step 0a: build genome_catalog/{genes_info,genome_toc}.tsv
├── gff_to_genes.py             # Prokka GFF3 -> .genes TSV (gffutils)
├── run_species.sh              # Step 0: materialise per-focal neighbor TSVs (entry)
├── generate_neighbor_list.sh   #   ↳ per-focal driver: joins catalog to .genes via genome_toc
├── get_neighbor.sh             #     ↳ per-genome: ±n_genes flank from one .genes file
├── prepare.R                   # Step 0: snapshot YAML, process focal_meta, list missing TSVs
├── pipeline.R                  # Steps 1-6: the driver — reads top-to-bottom
├── config.R                    # YAML loader + cfg_get
├── model.R                     # target_layout + get_target (file-path resolver)
├── neighbor.R                  # Step 1: per-focal neighborhood pipeline
├── midas.R                     # Step 1: small-ORF + length-variant labels; load_c80_tables
├── graph.R                     # Step 2 + Step 3: graph building, joint components, orientation
├── path.R                      # Step 3: canonical → fine → per-genome expansions
├── blocks.R                    # Step 6: hit blocks + reps + per-genome attribution
├── parse.R                     # Step 4 orchestrator + Step 3 c80s decorators + Step 5 plot data-prep helpers
├── plot.R                      # Step 5: global + per-component gggenes plots (and Step 1 diagnostic plots)
├── example.yaml                # template config
├── environment.yml             # conda env
├── examples/                   # worked-example input bundles (one per real run)
├── README.md                   # landing — quickstart + install
├── CLAUDE.md                   # orientation for AI agents
├── docs/
│   ├── USER_GUIDE.md           # this file
│   ├── STEPS.md                # per-step deep dive
│   ├── SCHEMA.md               # column-level schemas for every file the pipeline touches
│   ├── PIPELINE.md             # c80 glossary + truncation/fragmentation flag semantics
│   ├── diagram.md              # YAML key → consumer-step data flow
│   └── SETUP.md                # env setup + troubleshooting
└── parked/                     # supplementary docs not in the active flow
```

---

## Tunables — quick reference

| Knob | Section | Default | Where it bites | Effect |
| --- | --- | --- | --- | --- |
| `focal_meta` | `data` | — | Step 0 (`prepare.R`) | Absolute path to the user-provided focal-gene TSV. **Required.** |
| `n_genes` | `data` | 20 | Step 0 (`get_neighbor.sh`, `grep -C`) | Max neighbours extracted each side of a focal along the contig. |
| `score_col` | `prepare` | `cor_to_b` | Step 0 (`prepare.R`) | Trait-stat column the cutoffs apply to. One of `{cor_to_b, beta, ""}`. Empty = pass-through (focal_meta must already carry `is_focal`). |
| `inclusion_cutoff` | `prepare` | 0.25 | Step 0 (`prepare.R`) | Minimum `\|score_col\|` to retain a gene in `focal_c80_df` at all. |
| `focal_cutoff` | `prepare` | 0.5 | Step 0 (`prepare.R`) | Minimum `\|score_col\|` for a row to be flagged `is_focal = TRUE` and drive Step 1 neighbor extraction. |
| `focal_min_genomes` | `neighbor` | 10 | Step 1 (multiple gates) | Per-focal recurrence cut. Lower = more rare neighborhoods survive. |
| `focal_min_total_genomes` | `neighbor` | 30 | Step 1 (`parse_gene_neighbor`) | Min total focal coverage; focals with thinner support are dropped before pattern extraction. |
| `min_positions` | `neighbor` | 5 | Step 1 (`compute_relative_positions`, `parse_gene_neighbor`) | Min operon size. Lower = include shorter neighborhoods. |
| `upper_bound` | `neighbor` | 10 | Step 1 (`compute_relative_positions`) | Position window around focal (±). Larger = wider neighborhoods. |
| `min_left_neighbors` / `min_right_neighbors` | `neighbor` | 2 / 2 | Step 1 (`filter_by_flanking_coverage`) | Strict-mode flanking requirement. |
| `use_strict` | `neighbor` | `~` (auto) | Step 1 (`filter_by_flanking_coverage`) | `TRUE` forces strict; `FALSE` forces relaxed; `~` auto-selects based on `focal_min_genomes`. |
| `min_group_proportion` | `neighbor` | 0.05 | Step 1 (`parse_gene_neighbor`, second-pass filter) | A neighborhood-pattern group survives if it meets `focal_min_genomes` OR has at least this fraction of the focal's total genome support. |
| `coverage_warn_threshold` | `neighbor` | 0.8 | Step 1 (`parse_gene_neighbor`, second-pass filter) | If surviving pattern groups cover less than this fraction of the focal's total support, emit a warning. **Diagnostic-only.** |
| `path_min_genomes` | `path` | 20 | Step 3 (`generate_canonical_path`), Step 4 (driver, via `fine_coverage_ratio`) | Per-canonical-operon recurrence cut. Distinct from `focal_min_genomes`. |
| `truncation_cutoff` | `path` | 0.8 | Step 3 (`decorate_c80s_w_truncation`) | A gene shorter than this fraction of its centroid length is `is_truncated`. |
| `skip_block` | `blocks` | `false` | Step 6 (driver gate) | If `true`, pipeline.R skips the entire Step 6 call. |
| `allow_gaps` | `blocks` | 2 | Step 6 (`aggregate_blocks`) | Max non-hit positions allowed inside one focal block. |
| `min_overlap` | `blocks` | 1 | Step 6 (`rank_block_representatives`, `get_relation`) | Min LCS length to call two blocks "overlap" when ranking reps. |
| `min_shared` | `blocks` | 2 | Step 6 (`diagnose_rep_overlaps`) | **Diagnostic-only.** Substring-overlap threshold for the diagnostic log. |
| `fine_coverage_ratio` | `parse` | 0.25 | Step 4 (driver) | Fine isoforms survive if `n_fine_genomes >= ceiling(path_min_genomes * ratio)`. Historical default 0.5. |
| `seed` | `parse` | 616 | Step 4 (`sample_genome_from_fine_paths`) | RNG seed for the per-fine-isoform exemplar-genome draw. |
| `fill_modes` | `parse` | (user-supplied list) | Step 5 (`run_step5_figures`) | Which fill modes to render. Subset of `{beta, sample_prevalence, cor_to_b, fill_gene, gene_label, focal_type, focal_label}`. Column-tolerant — modes whose backing column is missing are skipped with a warning. |
| `gene_padding_bp` | `plot` | 100 | Step 1 (`extract_gene_neighbor_patterns`), Step 5 (`.layout_operon_tracks`) | bp gap between adjacent genes in gggenes layouts. |

---

## Known issues

These are documented here so you don't trip on them. Tracked items live in [parked/ROADMAP.md](../parked/ROADMAP.md).

1. **Mirror-block reps survive in Step 6.** [`is_contig_subseq`](../blocks.R#L388) is forward-direction only; a block and its exact reverse end up as two separate reps. The diagnostic in `diagnose_rep_overlaps` will catch this if it happens.
2. **Step 1 orientation is not preserved into Steps 2/3.** Step 2 re-derives chromosomal order; Step 3 canonicalizes lexicographically (with synthetic small-ORF tokens stripped from the decision; see `clean_for_orientation` in `graph.R`). Within-component direction consistency is the strongest guarantee you get on the output side.
3. **Off-species c80s have NA c80 metadata.** ECOR genes can map to centroid_80s from any species. Their `neighbor_c80_length_coarse` / `genome_prevalence` come through as NA after the `load_c80_tables` join; downstream NA-tolerant.

---

## Where to look next

- **Run end-to-end:** `python build_genome_catalog.py my_run.yaml; Rscript prepare.R my_run.yaml; bash run_species.sh my_run.yaml; Rscript pipeline.R my_run.yaml`. Outputs land in `proj_dir/{step1_setup,step2_neighbors,step3_path,step4_parse,step5_figures,step6_blocks}/` (where `proj_dir` is what you set in the YAML — include `<species_id>` if you want per-species isolation) and `data_dir/<species_id>/{list_of_neighbors,clusters_80_info_updated.tsv}`. (For a real input bundle to point your config at, see [`examples/`](../examples/).)
- **Read a step in detail:** [STEPS.md](STEPS.md) §STEP N has the full input / output / logic.
- **Read a function in detail:** every helper has a roxygen-style docstring covering arguments, behavior, and known caveats. Start from the function's call site in `pipeline.R` and follow the link.
- **Trace a column back through the pipeline:** the order of derivation is roughly Step 1 (neighbor table) → Step 2 (path strings) → Step 3 (canonical/fine/per-genome) → Step 6 (block reps). The c80 column table above explains the three c80 flavors.
