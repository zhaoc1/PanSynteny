# ------------------------------------------------------------------------------
# config.R
#
# Job configuration utilities: YAML loader and accessors for the `job_config`
# environment consumed by the rest of the pipeline.
#
# Author:   Chunyu Zhao <chunyu.zhao@gladstone.ucsf.edu>
# Created:  2025-07-14
# ------------------------------------------------------------------------------

load_job_config <- function(yaml_path) {
  stopifnot(file.exists(yaml_path))
  
  if (!exists("job_config", envir = .GlobalEnv)) {
    assign("job_config", new.env(), envir = .GlobalEnv)
  }
  raw <- yaml::read_yaml(yaml_path)

  # species-scope proj_dir and ensure it exists
  raw$job$proj_dir <- file.path(raw$job$proj_dir, raw$job$species_id)
  dir.create(raw$job$proj_dir, showWarnings = FALSE, recursive = TRUE)

  # Flatten sections into job_config (single flat namespace)
  for (section in names(raw)) {
    section_val <- raw[[section]]
    if (!is.list(section_val)) next
    for (key in names(section_val)) {
      assign(key, section_val[[key]], envir = job_config)
    }
  }
  cat(glue::glue("Loaded {yaml_path}: species={raw$job$species_id} trait={raw$job$trait}"), "\n")
}


cfg_get <- function(cfg, name) {
  # First look in cfg itself (inherits = FALSE); if not found,
  # walk up the parent chain (inherits = TRUE).
  if (exists(name, envir = cfg, inherits = FALSE)) {
    get(name, envir = cfg, inherits = FALSE)
  } else {
    get(name, envir = cfg, inherits = TRUE)
  }
}


print_cfg <- function(cfg, include_parent = TRUE, max_chars = 80) {
  stopifnot(is.environment(cfg))
  
  locals <- ls(envir = cfg, all.names = TRUE)
  
  parents <- character()
  if (include_parent) {
    direct_parent <- parent.env(cfg)
    if (!identical(direct_parent, emptyenv())) {
      parents <- setdiff(ls(envir = direct_parent, all.names = TRUE), locals)
    }
  }
  
  preview <- function(name) {
    val <- get(name, envir = cfg, inherits = TRUE)
    out <- paste(capture.output(str(val, max.level = 0)), collapse = " ")
    substr(out, 1, max_chars)
  }
  
  tab <- data.frame(
    name     = c(locals, parents),
    location = c(rep("local", length(locals)),
                 rep("parent", length(parents))),
    value    = vapply(c(locals, parents), preview, character(1)),
    stringsAsFactors = FALSE
  )
  
  tab <- tab[order(tab$location, tab$name), ]
  print(tab, row.names = FALSE, right = FALSE)
  invisible(tab)
}
