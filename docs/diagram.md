# YAML → pipeline data-flow diagram

This file maps every YAML key in [example.yaml](../example.yaml) to the script that consumes it. The R pipeline reads via `cfg_get(job_config, "<key>")`; the bash + Python step-0 chain reads via its own `yaml_get` helper (or `yaml.safe_load` in Python). `sources:` is a list, not a scalar section — it bypasses `config.R`'s flatten loop and is consumed only by `build_genome_catalog.py`.

The only YAML key consumed by more than one step is **`path.path_min_genomes`** (Steps 3 and 4). Every other key has a single owner.

---

## ASCII data-flow

```text
job.species_id ──────┬────────────────────────────────────────►  every step
data.midasdb_dir ────┤
data.data_dir ───────┤
                     │
sources: (list) ─────│──────►┌─────────────────────────┐
data.n_genes ────────│──┐    │ build_genome_catalog.py │
                     │  │    │             (Step 0a)   │
                     │  │    │  genes_info.tsv +       │
                     │  │    │  genome_toc.tsv +       │
                     │  │    │  <g>.genes in place     │
                     │  │    └────────────┬────────────┘
                     │  │                 │ catalog (Step 0a out)
                     │  │                 ▼
                     ▼  │    ┌─────────────────────────┐
data.focal_meta ────────│───►│   prepare.R    (Step 0) │
prepare.score_col ──────│───►│                         │
prepare.inclusion_… ────│───►│  snapshot run_config;   │
prepare.focal_cutoff ───│───►│  process focal_meta;    │
                        │    │  list missing TSVs      │
                        │    └────────────┬────────────┘
                        │                 │ focal_meta cache + gene_list.tsv
                        │                 ▼
                        └───►┌─────────────────────────┐
                             │ run_species.sh (Step 0) │
                             │  → generate_neighbor_   │
                             │      list.sh            │
                             │  → get_neighbor.sh      │
                             │                         │
                             │  per-focal TSVs under   │
                             │  list_of_neighbors/     │
                             └────────────┬────────────┘
                                          ▼
                              ┌─────────────────────────┐
neighbor.focal_min_genomes ──►│  neighbor.R  (Step 1)   │
neighbor.focal_min_total_… ──►│  midas.R                │
neighbor.min_positions ──────►│                         │
neighbor.upper_bound ────────►│  per-focal extraction   │
neighbor.min_left_neighbors ─►│  + small-ORF labels     │
neighbor.min_right_neighbors ►│  + length variants      │
neighbor.use_strict ─────────►│                         │
neighbor.min_group_proportion►│                         │
neighbor.coverage_warn_thres ►│                         │  ◄── diagnostic-only
plot.gene_padding_bp ────────►│  (diagnostic figures)   │
                              └────────────┬────────────┘
                                           │ neighbor_groups.RDS  ◄── CACHE GATE
                                           ▼
                              ┌─────────────────────────┐
              (no YAML keys) ►│  graph.R     (Step 2)   │
                              │  per-genome graphs ►    │
                              │  maximal paths          │
                              └────────────┬────────────┘
                                           │ path_df.rds  ◄── CACHE GATE
                                           ▼
                              ┌─────────────────────────┐
path.path_min_genomes ───────►│  graph.R / path.R       │
path.truncation_cutoff ──────►│              (Step 3)   │
                              │  L1/L2/L3 canonical     │
                              └────────────┬────────────┘
                                           ▼
                              ┌─────────────────────────┐
path.path_min_genomes ───────►│  parse.R     (Step 4)   │
parse.fine_coverage_ratio ───►│  selected_coarse / fine │
parse.seed ──────────────────►│  exemplar genome RNG    │
                              └────────────┬────────────┘
                                           ▼
                              ┌─────────────────────────┐
parse.fill_modes ────────────►│  plot.R      (Step 5)   │
plot.gene_padding_bp ────────►│  one PDF per fill mode  │
                              └────────────┬────────────┘
                                           ▼
                              ┌─────────────────────────┐
blocks.skip_block ───────────►│  blocks.R    (Step 6)   │  (gated)
blocks.allow_gaps ───────────►│  trait-block extraction │
blocks.min_overlap ──────────►│                         │
blocks.min_shared ───────────►│                         │  ◄── diagnostic-only
                              └─────────────────────────┘

job.proj_dir, data.* (paths) ──── consumed by model.R::target_layout(), used at every step
                                  (path resolution; not analytic tuning)

mwas.midas_dir (parked) ───── feeds the MWAS block in model.R; no reader in the current pipeline
```

---

## Cache-invalidation table

When you change a YAML key, you must also delete the listed cache file(s) to make the change take effect; otherwise the cached output of an earlier step is reused.

| If you change a key under… | …delete this cache |
| --- | --- |
| `sources` (list), `data.midasdb_dir`, `job.proj_dir`, `length_col` | `{proj_dir}/step1_setup/{catalog_genes_info.tsv, catalog_genome_toc.tsv}` — re-run `build_genome_catalog.py` |
| (prokka source's `.gff` updated) | also `rm <genomes_dir>/<g>/<g>.genes` for the affected genome(s); the `-s` guard skips re-conversion otherwise |
| `data.focal_meta`, `prepare.*` | `step1_setup/gene_meta_full.tsv` *(the `focal_meta` cache; re-run prepare.R)* + `step2_neighbors/neighbor_groups.RDS` *(Step 1)* |
| `data.n_genes` | every `{data_dir}/{species_id}/list_of_neighbors/<focal_c80>.tsv` *(re-run run_species.sh)* + `step2_neighbors/neighbor_groups.RDS` |
| `neighbor.*`, `plot.gene_padding_bp` (Step 1 use) | `step2_neighbors/neighbor_groups.RDS` |
| `path.path_min_genomes`, `path.truncation_cutoff` | `step3_path/path_df.rds` |
| `blocks.*` | nothing (Step 6 always re-runs when not skipped) |
| `parse.fine_coverage_ratio`, `parse.seed`, `parse.fill_modes` | nothing (Steps 4/5 always re-run) |
| `plot.gene_padding_bp` (Step 5 use) | nothing — Step 5 PDFs are always re-rendered |

Step 2 has no direct YAML keys; it inherits all tuning from Step 1 via the cached `gene_neighbors`.

---

## Diagnostic-only keys (do not gate outputs)

These keys produce log/warning messages but never filter rows out of any output:

- `neighbor.coverage_warn_threshold` — warns if surviving Step 1 pattern groups cover < this fraction of total focal support.
- `blocks.min_shared` — feeds `diagnose_rep_overlaps()` only; does not change `representatives.tsv`.

---

## Multi-step coupling

Only `path.path_min_genomes` is consumed by more than one step:

- **Step 3** ([path.R:445](../path.R#L445)): the canonical-operon survival floor. Operons backed by < `path_min_genomes` strains are dropped from L1 onward.
- **Step 4** ([parse.R:658](../parse.R#L658)): the fine-isoform survival threshold is `ceiling(path_min_genomes * fine_coverage_ratio)`.

If you tighten `path_min_genomes`, both gates get stricter together; the fine-isoform cut is relative.

---

## Verification — regenerate this diagram from code

```bash
# every cfg_get key the R pipeline reads must be present in the YAML
comm -23 \
  <(grep -hoE 'cfg_get\(job_config, "[^"]+"\)' *.R | sed -E 's/.*"([^"]+)".*/\1/' | sort -u) \
  <(grep -E '^[[:space:]]+[a-z_]+:' example.yaml | sed -E 's/^[[:space:]]+([a-z_]+):.*/\1/' | sort -u)
# (should print nothing; `proj_dir` is read directly via job_config$proj_dir, not cfg_get — known)

# every cfg_get key, grouped by file
grep -hnE 'cfg_get\(job_config, "[^"]+"\)' *.R | \
  sed -E 's|^([^:]+):([0-9]+):.*"([^"]+)".*|\3\t\1:\2|' | sort

# bash/python yaml keys (Step 0a + Step 0 chain)
grep -hnE 'yaml_get [^"]*"[a-z._]+"' *.sh | sed -E 's/.*"([^"]+)".*/\1/' | sort -u
```
