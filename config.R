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

  # Expand `{proj_dir}`, `{species_id}`, `{midasdb_dir}`, `{input_dir}`
  # placeholders in every string scalar in the YAML (recursively). Lets users
  # shorten repetitive paths in the config —
  # `data.focal_meta: "{input_dir}/focal_table.tsv"` instead of the full path.
  # `{input_dir}` keeps user-provided inputs separate from `{proj_dir}`
  # (the output root). Mirrors the placeholder system `build_genome_catalog`
  # uses for `sources:`.
  # input_dir + parallel_jobs are required under job: (no backward-compat fallbacks).
  if (is.null(raw$job$input_dir) || !is.character(raw$job$input_dir) || !nzchar(raw$job$input_dir)) {
    stop("job.input_dir is required in <config.yaml>. Add the absolute path to your user-provided inputs (focal_meta TSV, etc.) under job:.")
  }
  if (is.null(raw$job$parallel_jobs)) {
    stop("job.parallel_jobs is required in <config.yaml>. Add an integer (typical: 2) under job: to control run_species.sh's xargs -P fan-out.")
  }

  replacements <- c(
    "{proj_dir}"    = raw$job$proj_dir,
    "{species_id}"  = raw$job$species_id,
    "{midasdb_dir}" = if (!is.null(raw$data$midasdb_dir)) raw$data$midasdb_dir else "",
    "{input_dir}"   = raw$job$input_dir
  )
  expand_placeholders <- function(x) {
    if (is.character(x) && length(x) == 1) {
      for (ph in names(replacements)) x <- gsub(ph, replacements[[ph]], x, fixed = TRUE)
      x
    } else if (is.list(x)) {
      lapply(x, expand_placeholders)
    } else {
      x
    }
  }
  raw <- expand_placeholders(raw)

  # Ensure proj_dir exists. `proj_dir` is used as-is (no implicit species_id
  # suffix); if you want per-species isolation, include the species in the
  # YAML, e.g. `proj_dir: "/path/to/results/102506"`.
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
