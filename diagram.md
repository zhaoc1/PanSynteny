# Data-flow Diagram

This file maps every YAML key in [example.yaml](example.yaml) to the pipeline step that consumes it. It is auto-derivable from `grep -n cfg_get *.R` plus the section names in the YAML; if you suspect drift, regenerate by running the `comm` check at the bottom.

The only YAML key consumed by more than one step is **`path_min_genomes`** (Steps 3 and 5). Every other key has a single owner. This is the structural payoff of the rename: each key now has exactly one semantic frame.

---

## ASCII data-flow

```text
                              ┌─────────────────────────┐
job.species_id ──────────────►│  prepare.R   (Step 0)   │
job.trait ───────────────────►│                         │
prepare.score_col ───────────►│  build focal_c80_df     │
prepare.inclusion_cutoff ────►│  + enumerate missing    │
prepare.focal_cutoff ────────►│  neighbor TSVs          │
                              └────────────┬────────────┘
                                           │ focal_c80_df
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
blocks.allow_gaps ───────────►│  blocks.R    (Step 4)   │
blocks.min_overlap ──────────►│  trait-block extraction │
blocks.min_shared ───────────►│                         │  ◄── diagnostic-only
                              └────────────┬────────────┘
                                           ▼
                              ┌─────────────────────────┐
path.path_min_genomes ───────►│  parse.R     (Step 5)   │
parse.fine_coverage_ratio ───►│  selected_coarse / fine │
parse.seed ──────────────────►│  exemplar genome RNG    │
                              └────────────┬────────────┘
                                           ▼
                              ┌─────────────────────────┐
parse.fill_modes ────────────►│  plot.R      (Step 6)   │
plot.gene_padding_bp ────────►│  one PDF per fill mode  │
                              └─────────────────────────┘

job.proj_dir, data.* ──── consumed by model.R::target_layout(), used at every step
                          (path resolution; not analytic tuning)
```

---

## Cache-invalidation table

When you change a YAML key, you must also delete the listed cache file(s) to make the change take effect; otherwise the cached output of an earlier step is reused.

| If you change a key under… | …delete this cache |
| --- | --- |
| `prepare.*`, `neighbor.*`, `plot.gene_padding_bp` (Step 1 use) | `step1_focal_setup/gene_meta.tsv` *(prepare.R)* + `step2_neighbors/neighbor_groups.RDS` *(Step 1)* |
| `path.path_min_genomes`, `path.truncation_cutoff` | `step3_path/path_df.rds` |
| `blocks.*` | nothing (Step 4 always re-runs) |
| `parse.fine_coverage_ratio`, `parse.seed`, `parse.fill_modes` | nothing (Steps 5/6 always re-run) |
| `plot.gene_padding_bp` (Step 6 use) | nothing — Step 6 PDFs are always re-rendered |

Step 2 has no direct YAML keys; it inherits all tuning from Step 1 via the cached `gene_neighbors`.

---

## Diagnostic-only keys (do not gate outputs)

These keys produce log/warning messages but never filter rows out of any output:

- `neighbor.coverage_warn_threshold` — warns if surviving Step 1 pattern groups cover < this fraction of total focal support.
- `blocks.min_shared` — feeds `diagnose_rep_overlaps()` only; does not change `representatives.tsv`.

---

## Multi-step coupling

Only `path.path_min_genomes` is consumed by more than one step:

- **Step 3** ([path.R:445](path.R#L445)): the canonical-operon survival floor. Operons backed by < `path_min_genomes` strains are dropped from L1 onward.
- **Step 5** ([parse.R:658](parse.R#L658)): the fine-isoform survival threshold is `ceiling(path_min_genomes * fine_coverage_ratio)`.

If you tighten `path_min_genomes`, both gates get stricter together; the fine-isoform cut is relative.

