# Setup Instructions for Strain-Aware Operon Analysis

This document covers the environment needed to run the pipeline end-to-end. See [USER_GUIDE.md](USER_GUIDE.md) for the run command and YAML schema, and [STEPS.md](STEPS.md) for per-step details.

## Files

| File | Purpose |
| --- | --- |
| `environment.yml` | Conda environment spec - pins R + the Python helpers (`pyyaml`, `gffutils`). |
| `install_packages.R` | R fallback installer for packages not in conda. |
| `SETUP.md` | This file. |

## Setup

### Option 1: Conda (recommended)

1. **Create the env (pulls R packages + `pyyaml` + `gffutils` in one shot):**
   ```bash
   conda env create -f environment.yml
   ```

2. **Activate it:**
   ```bash
   conda activate pansynteny
   ```
   Everything below - `Rscript`, `python`, the bash chain - uses this env's interpreters.

3. **Optional: verify everything imports**
   ```bash
   Rscript install_packages.R
   python -c "import yaml, gffutils; print(yaml.__version__, gffutils.__version__)"
   ```

### Option 2: R only (no Python tooling)

If you skip the Python deps, `build_genome_catalog.py` will refuse to run with a clear error. You'd need to materialise `{proj_dir}/step1_setup/{catalog_genes_info.tsv, catalog_genome_toc.tsv}` and the per-prokka-genome `.genes` files yourself. Not recommended.

```r
source("install_packages.R")
```

## PDF rendering

The pipeline writes figures directly via `ggsave` - no LaTeX needed. LaTeX is only required if you also intend to render the optional `*.Rmd` companion documents.

- **TinyTeX (easiest, R-based):** `tinytex::install_tinytex()`.
- **Full distribution:** MacTeX / `texlive-full` / MiKTeX.

## Dependencies summary

### R packages (24 total)
- **tidyverse:** dplyr, tidyr, purrr, stringr, ggplot2, readr, tibble
- **Visualization:** gggenes, ggraph, gridExtra, viridis, RColorBrewer, scales, pheatmap
- **Graph analysis:** igraph
- **Data manipulation:** data.table
- **Utilities:** glue, fs
- **Reporting:** pander, knitr, rmarkdown
- **Base R:** tools, parallel (no install needed)

### Python packages (Step 0a)
- **`pyyaml`** - config parsing by every script that reads YAML directly.
- **`gffutils`** - Prokka GFF3 -> `.genes` conversion in `gff_to_genes.py`.

Both are pinned in `environment.yml`. `build_genome_catalog.py` hard-errors at startup if either is missing and tells you to activate the env.

## Verification

```r
required_packages <- c("dplyr", "tidyr", "purrr", "stringr", "ggplot2",
                       "gggenes", "tidyverse", "ggraph", "tibble",
                       "gridExtra", "viridis", "igraph", "pander",
                       "RColorBrewer", "scales", "glue",
                       "fs", "data.table", "readr", "knitr", "rmarkdown",
                       "pheatmap")
all_installed <- sapply(required_packages, require, character.only = TRUE, quietly = TRUE)
if (all(all_installed)) cat("All packages successfully installed!\n") else
  cat("Missing packages:\n"); print(names(all_installed)[!all_installed])
```

```bash
python -c "import yaml, gffutils; print('python deps OK')"
```

## Running the analysis

Once the env is set up, the workflow is four ordered commands, all reading the same `<config.yaml>`:

```bash
# Step 0a - build the unified genome catalog (genes_info.tsv + genome_toc.tsv;
#           derives <g>.genes in place for prokka sources via gff_to_genes.py).
python build_genome_catalog.py <config.yaml>

# Step 0b - snapshot the YAML, process focal_meta into the step1 cache,
#           list any missing per-focal neighbor TSVs.
Rscript prepare.R <config.yaml>

# Step 0c - materialise the missing per-focal neighbor TSVs.
bash build_neighbor_lists.sh <config.yaml>

# Steps 1-6 - the analytical pipeline.
Rscript pipeline.R <config.yaml>
```

Working example config: `example.yaml` (template).

## Troubleshooting

### Common issues

1. **`build_genome_catalog.py` says "pyyaml not importable" or "cannot import gff_to_genes".**
   The env isn't active. Run `conda activate pansynteny` (or set `PYTHON=/path/to/env/bin/python build_genome_catalog.py ...`).

2. **R fails with `GLIBCXX_3.4.30 not found` when loading `vroom`.**
   The system `Rscript` is being used instead of the env's. Either:
   ```bash
   conda activate pansynteny
   ```
   ...or invoke the env's binary directly with the env's `lib/` on `LD_LIBRARY_PATH`:
   ```bash
   LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH \
     $CONDA_PREFIX/bin/Rscript pipeline.R config.yaml
   ```

3. **`igraph` installation fails.**
   - macOS: `brew install glpk gmp`
   - Linux: `sudo apt-get install libglpk-dev libgmp-dev`

4. **PDF rendering fails.** Check `tinytex::is_tinytex()` in R, or install a system LaTeX.

5. **Memory issues with large datasets.**
   For UHGG species 102506 the catalog scan reads a 5.6 GB `genes_info.tsv` and the `gene_to_c80` join in pipeline.R materialises a ~34 M-row table. Run on a host with >=16 GB RAM; consider an HPC node for production.

## File dependencies

The pipeline driver (`pipeline.R`) sources these R scripts from the `R/` subdirectory:

- `config.R` - YAML loader + `cfg_get` accessor
- `model.R` - `target_layout` + `get_target` (file-path resolver)
- `graph.R` - Step 2 path stitching + Step 3 canonicalization helpers
- `path.R` - Step 3 canonical -> fine -> per-genome expansions
- `neighbor.R` - Step 1 per-focal neighborhood pipeline
- `midas.R` - Step 1 small-ORF + length-variant labels; `load_c80_tables`
- `blocks.R` - Step 6 focal-block extraction + per-genome attribution
- `plot.R` - Step 5 gggenes plotters + Step 1 diagnostic plotters
- `parse.R` - Step 4 orchestrator + Step 3 c80s decorators

`prepare.R` is a separate driver that sources only `config.R` + `model.R`.

The Step 0a / 0b / 0c entry-point scripts live at the repo root, alongside the R drivers (`prepare.R`, `pipeline.R`):

- `build_genome_catalog.py` - imports `scripts/gff_to_genes.py`
- `build_neighbor_lists.sh` - calls `scripts/focal_neighbor_list.sh` -> `scripts/get_neighbor.sh`

Their helpers live under `scripts/` (each resolves siblings via its own path, so the chain works regardless of CWD):

- `scripts/gff_to_genes.py`
- `scripts/focal_neighbor_list.sh`
- `scripts/get_neighbor.sh`

## Data requirements

The pipeline expects:

- **MIDAS reference DB** (`data.midasdb_dir`): `pangenomes/<species_id>/genes_info.tsv` + `pangenomes/<species_id>/clusters_80_info_updated.tsv` + `gene_annotations/<species_id>/<genome>/<genome>.genes`.
- **User-provided focal table** (`data.focal_meta`): absolute path to a TSV.
  - **Minimum required columns** (all four): `focal_c80`, `focal_label`, `is_focal`, `gene_label`.
  - **Optional (consumed only when present):** `cor_to_b`, `beta`, `sample_prevalence`, `trait`, `genome_counts`. Step 5 `fill_modes` whose backing column is absent are skipped with a warning.
  - **`is_focal` derivation:** if `prepare.score_col` is set, prepare.R derives `is_focal` from `|score_col| >= focal_cutoff` and **overwrites** any input `is_focal` (loud `warning()`). Set `prepare.score_col: ""` to preserve a hand-curated `is_focal`.
- **Project output root** (`job.proj_dir`): writable; used as-is (include the species_id in the value if you want per-species isolation).
- **Data root** (`data.data_dir`): writable; the per-focal neighbor TSVs land under `<data_dir>/<species_id>/list_of_neighbors/`. (Shared across runs.)
- **Project root** (`job.proj_dir`): writable; the genome catalog + Step 0b cache (`step1_setup/`) and every Step 1-6 output land under `<proj_dir>/`. (Per-run.)
- **Optional prokka sources** (declared in `sources:`): one directory per genome with `<genome>/<genome>.gff` (the catalog build derives `.genes` from it).

See [USER_GUIDE.md section Configuration](USER_GUIDE.md#configuration-yaml) for the full YAML schema.
