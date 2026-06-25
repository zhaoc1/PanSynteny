# ------------------------------------------------------------------------------
# model.R
#
# Pipeline target layout and path resolution. Defines the mapping from file-key
# names (e.g., "path_df", "canonical_paths_c80s") to project-relative paths,
# plus a resolver that expands a key against the active `run_config`.
#
# For every column of every file (schema details, writers, readers): see
# docs/SCHEMA.md. Inline comments below note non-obvious semantics only.
#
# Author:  Chunyu Zhao <chunyu.zhao@gladstone.ucsf.edu>
# Created: 2025-07-14
# Updated: 2026-05-15
# ------------------------------------------------------------------------------


target_layout <- function(species_id, asso="uni") {
  DATA_DIR <- cfg_get(job_config, "data_dir")
  MIDASDB_DIR <- cfg_get(job_config, "midasdb_dir")
  MIDAS_DIR <- cfg_get(job_config, "midas_dir")   # legacy: declared under YAML `mwas:`; used only by the MWAS block below
  
  list(
    # Input: read by prepare.R / pipeline.R / focal_neighbor_list.sh.
    neighbor_list               = glue::glue("{DATA_DIR}/{species_id}/list_of_neighbors"),
    clusters_80_updated         = glue::glue("{DATA_DIR}/{species_id}/clusters_80_info_updated.tsv"),
    
    run_config                  = glue::glue("step1_setup/run_config.yaml"),
    focal_meta                  = glue::glue("step1_setup/gene_meta_full.tsv"),
    gene_list                   = glue::glue("step1_setup/gene_list.tsv"),
    catalog_genes_info          = glue::glue("step1_setup/catalog_genes_info.tsv"),
    catalog_genome_toc          = glue::glue("step1_setup/catalog_genome_toc.tsv"),

    # Output
    short_gene_prevalence       = glue::glue("step1_setup/short_gene_prevalence.tsv"),
    c80_variants_mapping        = glue::glue("step1_setup/c80_variants_mapping.tsv"),

    neighbor_groups_rds         = glue::glue("step2_neighbors/neighbor_groups.RDS"),
    neighbor_groups_by_genome   = glue::glue("step2_neighbors/01_neighbor_by_genome"),
    neighbor_groups_by_focal    = glue::glue("step2_neighbors/02_neighbor_by_focal"),

    path_df                     = glue::glue("step3_path/path_df.rds"),
    esupport_df                 = glue::glue("step3_path/esupport_df.rds"),
    canonical_paths             = glue::glue("step3_path/canonical_paths_coarse.tsv"),
    canonical_paths_fine        = glue::glue("step3_path/canonical_paths_fine.tsv"),
    canonical_paths_per_genome  = glue::glue("step3_path/canonical_paths_per_genome.tsv"),
    canonical_paths_c80s        = glue::glue("step3_path/canonical_paths_c80s.tsv"),
    canonical_paths_fine_c80s   = glue::glue("step3_path/canonical_paths_fine_c80s.tsv"),

    parse_coarse_summary        = glue::glue("step4_parse/coarse_recurring_operons.tsv"),
    parse_fine_summary          = glue::glue("step4_parse/fine_isoform_priorities.tsv"),
    parse_selected_coarse       = glue::glue("step4_parse/selected_coarse.tsv"),
    parse_selected_fine         = glue::glue("step4_parse/selected_fine.tsv"),
    parse_fine_long             = glue::glue("step4_parse/fine_long.tsv"),
    parse_genome_paths_dir      = glue::glue("step4_parse/genome_paths"),

    neighbor_figures            = glue::glue("step5_figures/01_neighbor_by_focal"),
    parse_coarse_figures        = glue::glue("step5_figures"),
    parse_fine_figures          = glue::glue("step5_figures"),

    rep_path_df                 = glue::glue("step6_blocks/representative_path.tsv"),
    uid_path_df                 = glue::glue("step6_blocks/rep.tsv"),

    # MWAS (parked). None of these are read by the current operon pipeline.
    genes_info                  = glue::glue("{MIDASDB_DIR}/pangenomes/{species_id}/genes_annotated.tsv"),
    clusters_80                 = glue::glue("{MIDASDB_DIR}/pangenomes/{species_id}/clusters_80_info.tsv"),
    gene_by_sample_matrix       = glue::glue("{MIDAS_DIR}/{species_id}/gene_by_sample_matrix.rds"),
    genes_to_heatmap            = glue::glue("{MIDAS_DIR}/{species_id}/genes_to_heatmap.rds"),
    GRM_pop                     = glue::glue("{MIDAS_DIR}/{species_id}/GRM.rds"),
    pca_pop                     = glue::glue("{MIDAS_DIR}/{species_id}/pca_df.rds"),
    corrected_genes             = glue::glue("{DATA_DIR}/corrected_genes_0.01.RDS")
  )
}


get_target <- function(file_key, species_id=NULL, proj_dir=NULL, asso="") {
  if (!exists("job_config", envir = .GlobalEnv)) {
    stop("`job_config` not found. Please run `load_job_config()` first.")
  }
  
  species_id <- species_id %||% job_config$species_id
  proj_dir <- proj_dir %||% job_config$proj_dir
  
  rel_path <- target_layout(species_id, asso="")[[file_key]]
  if (is.null(rel_path)) {
    stop(glue::glue("file_key '{file_key}' not found in layout"))
  }
  
  is_absolute <- fs::is_absolute_path(rel_path)
  abs_path <- if (is_absolute) {
    rel_path
  } else {
    file.path(proj_dir, rel_path)
  }
  
  has_suffix <- grepl("\\.", rel_path)
  if (has_suffix) {
    dir_path <- dirname(abs_path)
    if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  } else {
    if (!dir.exists(abs_path)) dir.create(abs_path, recursive = TRUE, showWarnings = FALSE)
  }
  return(abs_path)
}
