#!/usr/bin/env Rscript
# R package installation script for strain-aware-operon analysis
# Run this after creating the conda environment

cat("Installing R packages for strain-aware-operon analysis...\n\n")

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# List of required packages
packages <- c(
  # Core tidyverse (in case not installed via conda)
  "tidyverse",
  "dplyr",
  "tidyr",
  "purrr",
  "stringr",
  "ggplot2",
  "readr",
  "tibble",

  # Visualization
  "viridis",
  "RColorBrewer",
  "scales",
  "gridExtra",
  "ggraph",
  "gggenes",      # Gene arrow diagrams

  # Graph analysis
  "igraph",

  # Data manipulation
  "data.table",

  # Utilities
  "glue",
  "fs",

  # Reporting
  "pander",
  "knitr",
  "rmarkdown"
)

# Function to install packages if not already installed
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste0("Installing ", pkg, "...\n"))
    install.packages(pkg, dependencies = TRUE)
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat(paste0("WARNING: Failed to install ", pkg, "\n"))
      return(FALSE)
    } else {
      cat(paste0("Successfully installed ", pkg, "\n"))
      return(TRUE)
    }
  } else {
    cat(paste0(pkg, " is already installed\n"))
    return(TRUE)
  }
}

# Install all packages
results <- sapply(packages, install_if_missing)

# Summary
cat("\n=== Installation Summary ===\n")
cat(paste0("Total packages: ", length(packages), "\n"))
cat(paste0("Successfully available: ", sum(results), "\n"))
cat(paste0("Failed: ", sum(!results), "\n"))

if (sum(!results) > 0) {
  cat("\nFailed packages:\n")
  print(names(results)[!results])
}

# Optional: Install TinyTeX for PDF rendering if not already installed
cat("\n=== Checking PDF rendering capability ===\n")
if (!require("tinytex", quietly = TRUE)) {
  cat("Installing tinytex package...\n")
  install.packages("tinytex")
}

if (require("tinytex", quietly = TRUE)) {
  if (!tinytex::is_tinytex()) {
    cat("TinyTeX is not installed. Installing TinyTeX for PDF rendering...\n")
    cat("This may take several minutes...\n")
    tinytex::install_tinytex()
    cat("TinyTeX installation complete!\n")
  } else {
    cat("TinyTeX is already installed\n")
  }
}

cat("\n=== Package installation complete! ===\n")
cat("You can now run the pipeline: `Rscript prepare.R <config.yaml>` then `Rscript pipeline.R <config.yaml>`.\n")
