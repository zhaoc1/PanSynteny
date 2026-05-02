# Setup Instructions for Strain-Aware Operon Analysis

This document provides instructions for setting up the R environment needed to run the strain-aware operon pipeline (`prepare.R` then `pipeline.R`, both YAML-driven). See [USER_GUIDE.md](USER_GUIDE.md) for the run command and YAML schema.

## Files Created

1. **LIBRARY_LIST.md** - Complete list of all R package dependencies with descriptions
2. **environment.yml** - Conda environment specification
3. **install_packages.R** - R script to install packages directly in R
4. **README_SETUP.md** - This file

## Setup Options

### Option 1: Using Conda (Recommended)

1. **Create the conda environment:**
   ```bash
   conda env create -f environment.yml
   ```

2. **Activate the environment:**
   ```bash
   conda activate strain-aware-operon
   ```

3. **Install remaining packages not available in conda:**
   ```bash
   R -e "install.packages(c('randomcoloR', 'gggenes'), repos='https://cloud.r-project.org/')"
   ```

4. **Optional: Run the installation script to verify all packages:**
   ```bash
   Rscript install_packages.R
   ```

### Option 2: Using R Directly

1. **Open R or RStudio**

2. **Run the installation script:**
   ```r
   source("install_packages.R")
   ```

   Or install packages manually:
   ```r
   options(repos = c(CRAN = "https://cloud.r-project.org/"))

   packages <- c("tidyverse", "gggenes", "ggraph", "gridExtra",
                 "viridis", "igraph", "pander", "RColorBrewer",
                 "scales", "randomcoloR", "glue", "fs",
                 "data.table", "knitr", "rmarkdown")

   install.packages(packages)
   ```

## PDF Rendering Support

The pipeline writes its figures directly to PDF via `ggsave` (no LaTeX needed). LaTeX is only required if you also intend to render the optional `*.html` / `*.Rmd` companion documents to PDF. Choose one:

### Option A: TinyTeX (Easiest, R-based)
```r
install.packages("tinytex")
tinytex::install_tinytex()
```

### Option B: Full LaTeX Distribution
- **macOS**: Install [MacTeX](https://www.tug.org/mactex/) or BasicTeX
- **Linux**: `sudo apt-get install texlive-full` (Debian/Ubuntu)
- **Windows**: Install [MiKTeX](https://miktex.org/)

## Dependencies Summary

### Core R Packages (24 total)
- **tidyverse ecosystem**: dplyr, tidyr, purrr, stringr, ggplot2, readr, tibble
- **Visualization**: gggenes, ggraph, gridExtra, viridis, RColorBrewer, scales, randomcoloR
- **Graph analysis**: igraph
- **Data manipulation**: data.table
- **Utilities**: glue, fs
- **Reporting**: pander, knitr, rmarkdown

### Base R (included with R)
- tools, parallel

## Verification

After installation, verify all packages are available:

```r
# Check if all required packages can be loaded
required_packages <- c("dplyr", "tidyr", "purrr", "stringr", "ggplot2",
                       "gggenes", "tidyverse", "ggraph", "tibble",
                       "gridExtra", "viridis", "igraph", "pander",
                       "RColorBrewer", "scales", "randomcoloR", "glue",
                       "fs", "data.table", "readr", "knitr", "rmarkdown")

all_installed <- sapply(required_packages, require, character.only = TRUE, quietly = TRUE)
if (all(all_installed)) {
  cat("All packages successfully installed!\n")
} else {
  cat("Missing packages:\n")
  print(names(all_installed)[!all_installed])
}
```

## Running the Analysis

Once setup is complete, run the two driver scripts in order:

```bash
# Step 0: build focal_c80_df + enumerate any missing per-focal neighbor TSVs
Rscript prepare.R <config.yaml>

# Steps 1-6: the analytical pipeline
Rscript pipeline.R <config.yaml>
```

A working example config is `example.yaml`. See [USER_GUIDE.md](USER_GUIDE.md) for the YAML schema and [STEPS.md](STEPS.md) for per-step input / output / logic.

## Troubleshooting

### Common Issues

1. **igraph installation fails**
   - Install system dependencies first:
     - macOS: `brew install glpk gmp`
     - Linux: `sudo apt-get install libglpk-dev libgmp-dev`

2. **PDF rendering fails**
   - Ensure LaTeX is installed (see PDF Rendering Support above)
   - Check: `tinytex::is_tinytex()` in R

3. **gggenes not in conda**
   - Install directly from CRAN in R:
     ```r
     install.packages("gggenes")
     ```

4. **Memory issues with large datasets**
   - Increase R memory limit:
     ```r
     memory.limit(size = 16000)  # Windows
     ```
   - Use data.table for large files
   - Consider running on HPC cluster

## File Dependencies

The pipeline driver (`pipeline.R`) sources these R scripts (must be in the same directory):

- `config.R` — YAML loader + `cfg_get` accessor
- `model.R` — `target_layout` + `get_target` (file-path resolver)
- `graph.R` — Step 2 path stitching + Step 3 canonicalization helpers
- `path.R` — Step 3 canonical → fine → per-genome expansions
- `neighbor.R` — Step 1 per-focal neighborhood pipeline
- `midas.R` — Step 1 small-ORF + length-variant labels
- `blocks.R` — Step 4 focal-block extraction + per-genome attribution
- `plot.R` — Step 6 gggenes plotters + Step 1 diagnostic plotters
- `parse.R` — Step 5 orchestrator + Step 3 c80s decorators

`prepare.R` is a separate driver that sources only `config.R` + `model.R`.

## Data Requirements

The pipeline expects data in these locations (set via the `paths` section of the YAML config; see [USER_GUIDE.md](USER_GUIDE.md)):

- MIDAS gene-by-sample matrices: `paths.midas_dir`
- MIDASDB reference (clusters_80, genes_info): `paths.midasdb_dir`
- DefenseFinder per-cluster results: `paths.df_dir`
- Per-focal neighbor TSV inputs + `corrected_genes.RDS`: `paths.data_dir`
