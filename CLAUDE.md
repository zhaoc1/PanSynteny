# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Read this first — keystones

Six things that aren't obvious from skimming the code:

1. **The graph is internal plumbing, not output.** Users receive operon TSVs and gggenes PDFs. The graph layer (igraph, DFS, joint components) lives entirely inside Step 2 / Step 3 and dies before Step 4. Frame any user-facing description around operons, not graphs.

2. **`proj_dir` is used as-is** — no implicit `species_id` suffix any more. If you want per-species isolation, include the species in the YAML value, e.g. `proj_dir: "/path/.../<species_id>"`. Multi-trait or multi-species runs that share the same `proj_dir` will clobber each other; treat one `proj_dir` per run as the rule. (This changed in v0.2.0 — v0.1.0 silently appended `species_id` inside `config.R`.)

3. **YAML sections are cosmetic.** `cfg_get` flattens every scalar YAML section into one namespace; section placement is for human readability only. Two keys with the same name across sections would shadow each other silently. **Exception:** `sources:` is a list (not a scalar section) — `config.R`'s flatten loop skips it entirely. It is consumed only by `build_genome_catalog.py`.

4. **Two `focal_meta` namespaces.** `cfg_get(job_config, "focal_meta")` returns the **YAML-declared input path** (the user's raw focal table, anywhere on disk). `get_target("focal_meta")` returns the **step1 cache target** (`step1_setup/gene_meta_full.tsv` — filename kept for backward compatibility). `prepare.R` reads the first, optionally applies thresholds, and writes the second. `pipeline.R` only sees the cache.

5. **Cache deletion is the re-run signal.** Step 1 skips if `step2_neighbors/neighbor_groups.RDS` exists; Step 2 skips if `step3_path/path_df.rds` exists. Steps 3–5 always re-run. `build_genome_catalog.py` rebuilds `step1_setup/catalog_{genes_info,genome_toc}.tsv` fresh every run but skips the per-genome prokka `.gff → .genes` conversion via a `-s` guard.

6. **`focal_min_genomes` ≠ `path_min_genomes`** despite sounding interchangeable. Step 1 cut on per-focal pattern recurrence; Step 3+5 cut on per-canonical-operon recurrence. Different scopes; usually but not necessarily the same value.

**Project status.** v0.2.0 testing release published to `github.com/zhaoc1/PanSynteny` (private). Schema may shift between 0.x.y versions before 1.0.

## Run commands

The full workflow is four ordered commands, all reading the same `<config.yaml>`. Everything runs under the `strain-aware-operon` conda env (both R and Python live there):

```bash
# Step 0a — build the unified genome catalog (per-source membership + .genes paths,
#           merged into two artefacts the rest of the chain consumes).
python  build_genome_catalog.py <config.yaml>

# Step 0  — process the user-provided focal_meta TSV; cache to step1; enumerate
#           any missing per-focal neighbor TSVs.
Rscript prepare.R               <config.yaml>

# Step 0  — materialise the missing per-focal neighbor TSVs (fans
#           generate_neighbor_list.sh -> get_neighbor.sh over gene_list.tsv).
bash    run_species.sh          <config.yaml>

# Steps 1-6 — the analytical pipeline.
Rscript pipeline.R              <config.yaml>
```

Working example config: [example.yaml](example.yaml). A real worked-example input bundle (config + focal_meta TSV) lives under [examples/](examples/). `prepare.R` is cheap to re-run (always overwrites the focal_meta cache, the run_config.yaml snapshot, and gene_list.tsv). `pipeline.R` aborts at startup if the focal_meta cache is missing or any `is_focal == TRUE` centroid lacks its neighbor TSV under `neighbor_list/` — both errors point back to `prepare.R`.

**Per-focal neighbor TSVs are now materialised in-repo.** `run_species.sh` consumes `gene_list.tsv` (the missing-list `prepare.R` writes) and fans `generate_neighbor_list.sh` over each focal; per-focal idempotency lives in that script. There is no longer an external preprocessing job to coordinate (v0.1.0 required one).

**Re-run skipping by cache.** Step 1 skips re-extraction if `step2_neighbors/neighbor_groups.RDS` exists; Step 2 skips if `step3_path/path_df.rds` exists. **To force a step to re-run, delete its cache file.**

There is no test suite, no Makefile, no lint config. To smoke-test that a single helper file still parses after edits:

```bash
Rscript -e 'invisible(parse(file = "graph.R")); cat("OK\n")'
```

Setup details (conda env, R + Python packages, troubleshooting) live in [SETUP.md](docs/SETUP.md).

## Architecture

Single R script per pipeline stage; `pipeline.R` is the linear driver and reads top-to-bottom. Each numbered step is a one-line orchestrator call (`run_stepN_*`) into a topic-specific helper file. The driver passes step-N outputs into step-(N+1) explicitly — there is no implicit shared state beyond `job_config`.

### R pipeline (Steps 1–6)

| File | Role |
| --- | --- |
| [pipeline.R](pipeline.R) | Driver. Sources every helper, loads config, calls `run_step{1..6}_*` in order. |
| [prepare.R](prepare.R) | Step 0. Snapshots `<config.yaml>` to the proj_dir; reads `data.focal_meta` from YAML, optionally applies `\|score_col\|` thresholds, caches to `get_target("focal_meta")`; enumerates missing per-focal neighbor TSVs. Sources only `config.R` + `model.R`. |
| [config.R](config.R) | YAML loader. `load_job_config()` flattens every scalar YAML section into a single `job_config` env; `cfg_get(job_config, "key")` is the only accessor. List sections (e.g. `sources:`) are silently skipped. |
| [model.R](model.R) | `target_layout()` — single source of truth for every input + output file path keyed by name. `get_target("key")` resolves against active config. The `# MWAS (parked)` block at the bottom is reserved for re-integration and not read by the current pipeline. |
| [neighbor.R](neighbor.R) | Step 1: per-focal neighborhood extraction. |
| [midas.R](midas.R) | Step 1: small-ORF synthetic labels + length-variant labels; `load_c80_tables()` reads `catalog_genes_info` + `clusters_80_updated`. |
| [graph.R](graph.R) | Step 2: per-genome graphs → maximal paths. Step 3: canonicalization, joint components, orientation. |
| [path.R](path.R) | Step 3: canonical → fine → per-genome expansions. |
| [parse.R](parse.R) | Step 3 c80s decorators (small-ORF, truncation/fragmentation). Step 4 orchestrator. |
| [plot.R](plot.R) | Step 5: gggenes plotters (global + per-component). Plus Step 1 diagnostic plots. |
| [blocks.R](blocks.R) | Step 6: trait-associated block extraction + representative ranking. |

Note the v0.2.0 step renumbering vs v0.1.0: parse was Step 5 → now Step 4; figures was Step 6 → now Step 5; block extraction was Step 4 → now Step 6 (and is now skippable; see gotchas).

### Step 0a — genome catalog build

| File | Role |
| --- | --- |
| [build_genome_catalog.py](build_genome_catalog.py) | Reads `sources:` from the YAML; for each source streams membership rows to a merged `catalog_genes_info.tsv` (`gene_id <TAB> centroid_80 <TAB> gene_length`), accumulates a `catalog_genome_toc.tsv` (`genome_id <TAB> .genes path`), converts prokka `.gff → .genes` in place via `gff_to_genes.py`, dup-checks genome_ids across sources. Outputs land under `{proj_dir}/step1_setup/` (per-run). |
| [gff_to_genes.py](gff_to_genes.py) | Prokka GFF3 → MIDAS `.genes` TSV converter (uses `gffutils`). Imported in-process by `build_genome_catalog.py`. |

### Step 0 — neighbor-TSV materialisation (bash chain)

| File | Role |
| --- | --- |
| [run_species.sh](run_species.sh) | Top entry point: reads `gene_list.tsv` (the missing-list from `prepare.R`), fans `generate_neighbor_list.sh` over each focal in parallel. |
| [generate_neighbor_list.sh](generate_neighbor_list.sh) | One focal centroid → its `<query>.tsv`. Joins the catalog `genes_info` (gene members of the focal) to `genome_toc` (`.genes` path per genome), then fans `get_neighbor.sh` per gene member. Idempotent on `-s "$outfile"`. |
| [get_neighbor.sh](get_neighbor.sh) | Innermost — emits one gene's ±`n_genes` flank from its `.genes` file. Has no hardcoded paths. |

### Data flow at a glance

```text
Step 0a   sources: (UHGG midasdb + ECOR prokka + ...)
            --> {proj_dir}/step1_setup/catalog_genes_info.tsv  (gene_id, centroid_80, gene_length)
            --> {proj_dir}/step1_setup/catalog_genome_toc.tsv  (genome_id, .genes path)
            --> per-prokka-genome <g>/<g>.genes  (in place next to <g>.gff)

Step 0    data.focal_meta (user TSV)
            --> {proj_dir}/step1_setup/run_config.yaml  (config snapshot)
            --> {proj_dir}/step1_setup/gene_meta_full.tsv  (focal_meta cache)
            --> {proj_dir}/step1_setup/gene_list.tsv    (only missing focals)

          run_species.sh + generate_neighbor_list.sh + get_neighbor.sh
            --> {data_dir}/{species_id}/list_of_neighbors/<focal_c80>.tsv  (7 cols, no header)

pipeline.R Step 1    per-focal neighbor TSVs + catalog c80 tables
                                                 -->  gene_neighbors (RDS cache)
                                                    + small-ORF labels + length-variant labels
           Step 2    gene_neighbors            -->  path_df (one row per per-genome maximal path, RDS cache)
           Step 3    path_df                   -->  canonical_paths (L1) / canonical_paths_fine (L2) / canonical_paths_per_genome (L3)
                                                 + canonical_paths_c80s (L1 per-gene) / canonical_paths_fine_c80s (L2 per-gene)
           Step 4    L1/L2/L3 c80s + per-gen   -->  selected_coarse / selected_fine / fine_long + per-(uid_fine,genome) BLAST gene-id TSVs
           Step 5    selected_* + c80s tables  -->  gggenes PDFs (global + per-component, one per fill_mode)
           Step 6    canonical_paths_c80s      -->  representative_path.tsv + rep.tsv (block reps × genome)
```

All R outputs land under `<proj_dir>/{step1_setup, step2_neighbors, step3_path, step4_parse, step5_figures, step6_blocks}/` — this now includes the genome catalog (`step1_setup/catalog_{genes_info,genome_toc}.tsv`), so per-run isolation is the default. Per-focal neighbor TSVs stay under `{data_dir}/{species_id}/list_of_neighbors/` (species-shared across proj_dirs). Layout is fixed by `target_layout()` in [model.R](model.R) — **always go through `get_target("key")` rather than constructing paths inline.**

### The three granularity levels (Step 3 is the analytical core)

The same operons are emitted at three levels and three stable IDs. Downstream code joins on these IDs, never on path strings:

- **L1 coarse** — `uid = "cmp{joint_component_ids}-{path_type}-{canonical_path_id}-ng{n_genomes}"`. One row per canonical operon. Length variants collapsed.
- **L2 per-isoform** — `uid_fine = "{uid}-iso{rank}-ngf{n_fine_genomes}"`. Strip `-iso\d+-ngf\d+$` to recover the parent `uid`. Truncation and fragmentation flags are populated **here only** (the L1 table's `neighbor_gene_length` is max-over-isoforms, making per-row length checks meaningless).
- **L3 per-genome** — one row per `(canonical, contributing genome)`. Only level with per-genome `gene_id` resolution. Carries `needs_flip` (boolean computed once at the canonical × collapsed-path grain) and a `gene_path_canonical` already flipped according to it.

### Three c80 columns — easy to confuse

| Column | Resolution | Scope |
| --- | --- | --- |
| `centroid_80` / `c80` | coarse | the **focal** cluster the row belongs to (focal-scoped — focal + its synthetic small ORFs share this) |
| `neighbor_c80_coarse` | coarse | a neighbor's MIDAS cluster id, OR a synthetic `_<focal>-<type>_<rank>` for short ORFs (always starts with `_`) |
| `neighbor_c80_fine` | length-variant-aware | same as coarse + `_<rank>` suffix when the cluster has multiple observed lengths |

Synthetic labels for short ORFs are **per-focal** — the same physical short gene next to two different focals receives two different synthetic labels. Never compare these strings across focals. The leading `_` is how `is_smallORF` is detected. See [PIPELINE.md](docs/PIPELINE.md) for the full glossary including truncation/fragmentation flag semantics.

### Path direction: each step canonicalizes in a different frame

This is the most subtle invariant in the pipeline.

| Step | Direction frame | What it discards |
| --- | --- | --- |
| 1 | Focal-relative (left vs right of focal) | Ambiguity from `right_anchor < left_anchor` |
| 2 | Chromosomal (sorted by `neighbor_gene_start`) | **Step 1's focal-relative orientation is overwritten here** |
| 3 | Lexicographic (`normalize_path` = lex-min of forward vs reverse), then re-aligned within joint component to match the longest reference | Step 2's chromosomal orientation; absolute biological direction |

In the L1/L2/L3 outputs, "left-to-right" means *consistent within a joint component*, not 5'→3' of any chromosome. To recover real direction, walk back via L3's per-genome `gene_path` to chromosomal coordinates in `path_df` or the original neighbor TSVs.

**Step 3 orientation-only contract.** `normalize_path` and `orient_paths_within_component` decide direction on a *cleaned* token vector (synthetic `_`-prefixed small-ORF tokens stripped via `clean_for_orientation` in `graph.R`), then apply the chosen forward/reverse to the **full** original token list. The orientation step never reorders within a direction and never changes content — this is what lets a single boolean `needs_flip` at the canonical × collapsed-path grain drive both fine-resolution and gene-id-resolution flips downstream. **If you add a new token type that shouldn't influence direction, extend `clean_for_orientation` rather than special-casing in each caller.**

### `needs_flip` grain (per-genome table)

Computed once at the **(canonical × collapsed_path)** grain in `explode_canonical_into_collapsed_paths` as `(c80_path_string != canonical_path_coarse)`, then inherited by every per-genome row under that collapsed group. All rows sharing a `collapsed_path_id` share `needs_flip`; within one isoform (`uid_fine`) or one canonical (`uid`), `needs_flip` can be mixed because mirror-image collapsed groups roll up. This is also why `needs_flip` cannot be promoted to the L1 or L2 tables.

### Shared genome_id contract

`build_genome_catalog` and `generate_neighbor_list.sh` both derive `genome_id` from `gene_id` by stripping the trailing `_NNNNN` field. This works for both `GUT_GENOME000040_00388` and `GCF_900448275.1_00001`. The catalog's `genome_toc` is keyed by the same value. **If you add a new source whose gene_ids don't follow this convention, the TOC join in `generate_neighbor_list.sh` will silently miss those genes** (it warns to stderr on each miss but doesn't abort). Extend the derivation in both places (Python `genome_id_from_gene_id` + the awk in `generate_neighbor_list.sh`).

## Documentation map

When the user asks about pipeline behavior, prefer the source docs over re-deriving from code:

- [USER_GUIDE.md](docs/USER_GUIDE.md) — high-level orientation, run command, YAML schema, output reference, tunables quick-reference.
- [STEPS.md](docs/STEPS.md) — per-step input / output / logic / known caveats. The deep reference.
- [SCHEMA.md](docs/SCHEMA.md) — every file the pipeline reads or writes, by columns. Single source of truth for data formats.
- [PIPELINE.md](docs/PIPELINE.md) — c80 column glossary + truncation/fragmentation flag semantics.
- [SETUP.md](docs/SETUP.md) — R + Python package install, conda env, troubleshooting.
- `parked/` — supplementary docs not in the active flow: ROADMAP.md, CRITIQUE.md, VALIDATION.md, FLOWCHART.md.

The MD docs reference function definitions by line number (e.g. `[run_step1_neighbor_extraction](neighbor.R#L791)`). When you edit those functions, **update the line numbers** — there is a recent commit (6f7ba3a) titled "Remove dead functions flagged in previous sweep" and one (f140b94) titled "docs: refresh driver line numbers", so the user actively maintains this and notices when it drifts.

## Gotchas

Landmines that have caught people. Read before changing anything in these areas.

- **Step 1 orientation is overwritten by Step 2.** `orient_focal_gene_neighbors` decides direction per focal in Step 1; Step 2's chromosomal re-derivation discards it. If a task asks for focal-relative direction preserved end-to-end, that path doesn't exist today and would require carrying the Step 1 orientation column through Step 2.

- **`neighbor_gene_length` has two semantics, only one in the output.** In `gene_neighbors` (Step 1) it's per-gene observed; in `c80s_fine` (Step 3) it's per-isoform consensus. Same column name, different meaning. Don't merge them — `enrich_fine_long` deliberately drops the c80s_fine version at join time.

- **Coarse-table `is_truncated` is forced to NA.** Truncation only makes sense at fine resolution because the coarse `neighbor_gene_length` is max-over-isoforms. Don't compute truncation off the L1 table.

- **Synthetic small-ORF labels are per-focal-scoped.** A label like `_<focal>-<type>_<rank>` is only meaningful within its focal context. Never compare these strings across focals. The leading `_` is also how `is_smallORF` is detected.

- **`needs_flip` lives at the canonical × collapsed-path grain.** Per-genome rows inherit it; L1 / L2 cannot promote it because mirror-image collapsed groups roll up. Don't try to add `needs_flip` to the L1 or L2 tables.

- **`generate_canonical_path` requires `path_min_genomes` explicitly.** A previous version defaulted to 10 (shadowed by the `cfg_get` value); that default was removed in the v0.1.0 refactor. A future caller that forgets to thread the YAML value gets an explicit error rather than silent fallback.

- **Off-species c80s have NA c80 metadata.** ECOR genes can map to centroid_80s from any species (not just the run's `species_id`). After the `load_c80_tables` join, those rows carry the correct `centroid_80` and `gene_length` but `neighbor_c80_length_coarse` / `genome_prevalence` are NA. The downstream truncation/fragmentation flags fall through NA-tolerantly.

- **`is_focal` in the input can be overwritten.** When `prepare.score_col` is set, prepare.R derives `is_focal` from the `|score_col| >= focal_cutoff` threshold and **overwrites any existing `is_focal` column in the input** (with a `warning()`). Users who curate `is_focal` by hand (e.g. mixing focal rows with `is_focal = FALSE` context rows for plot metadata) should set `prepare.score_col: ""` to preserve the input.

- **`blocks.skip_block` short-circuits Step 6 (block extraction).** Setting `blocks.skip_block: true` makes pipeline.R skip the `run_step6_blocks` call entirely. Steps 4 (parse) and 5 (figures) do not depend on block-extraction outputs and proceed unchanged. Nothing is written under `step6_blocks/` in that case.

- **Step 5 `fill_modes` is column-tolerant.** `run_step5_figures` filters `fill_modes` to drop any mode whose backing column is missing from the c80s tables (warning printed). Listing `beta` / `cor_to_b` / `sample_prevalence` on a minimal focal_meta is harmless — they just no-op. The minimal set that always works: `fill_gene` (derived) and any column that's in the c80s tables.

- **Catalog `genes_info.tsv` has a header.** Awk consumers must `FNR==1 { next }`. `generate_neighbor_list.sh` already does; if you add a new awk reader, replicate it.

- **Don't `source()` from `archive/`.** Files there are intentionally stale historical versions kept for archaeology only.

## Conventions worth knowing

- **Every helper has a roxygen-style docstring** covering arguments, behavior, and known caveats. Read it before changing the function body.
- **`get_target("key")` everywhere.** Never hardcode a path under `step{N}_*/`; add a key to `target_layout()` instead.
- **`cfg_get(job_config, "key")` for every tunable.** No `cfg_get`-bypassing constants in helper files. The full set of YAML keys + defaults is enumerated in [USER_GUIDE.md §Tunables](docs/USER_GUIDE.md).
- **Caches are RDS or TSV; deletion is the re-run signal.** Document in any new step which file gates re-execution.
- **`archive/` holds historical versions of every script** (timestamped subdirs) plus completed plan MDs. Useful for "what did this used to do" archaeology, but **never `source()` from `archive/`** — those files are intentionally stale.
- **The `manuscript/` directory is gitignored.** Don't commit it without asking.
