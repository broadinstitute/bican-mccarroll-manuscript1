# options(
#     bican.mccarroll.figures.data_root_dir =
#         "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis",
#
#     bican.mccarroll.figures.out_dir =
#         "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/figure_repository",
#
#     bican.mccarroll.figures.cache_dir =
#         "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/figure_repository/data_cache"
# )
#
# run_private_eqtl_tests_bican()
# summarize_private_eqtl_tests_bican()
#
# run_private_eqtl_tests_downsampled_bican()
# summarize_private_eqtl_tests_downsampled_bican()

#' Run cell-type-specific private eQTL tests (main data set)
#'
#' Wires \code{bican.mccarroll.eqtl::run_cell_type_specific_eqtl_tests_from_manifest()}
#' into the package's options-driven path resolution, using the primary
#' (non-downsampled) cluster assignments and eQTL results. Per-analysis
#' outputs (gathered data, covariate mapping, report, summary, linear-model
#' results, and PDF) are written under a \code{private_eqtl_analysis}
#' subdirectory of the configured cache directory (see
#' \code{\link{get_cache_dir}}), since these are intermediate analysis files
#' rather than final figures.
#'
#' @param covariates Character vector of donor-level covariates to include in
#'   the gathered analysis dataset.
#' @param force_factor_covariates Character vector of covariates that should
#'   be encoded as factors during linear-model testing.
#' @param data_cache_dir Cache directory under which per-analysis outputs are
#'   written. If \code{NULL}, resolved via configured cache directory options.
#'
#' @return Invisibly returns the list produced by
#'   \code{bican.mccarroll.eqtl::run_cell_type_specific_eqtl_tests_from_manifest()}.
#'
#' @seealso
#'   \code{\link[bican.mccarroll.eqtl]{run_cell_type_specific_eqtl_tests_from_manifest}}
#'   \code{\link{summarize_private_eqtl_tests_bican}}
#'   \code{\link{run_private_eqtl_tests_downsampled_bican}}
#'
#' @export
run_private_eqtl_tests_bican <- function(covariates = c(
                                           "biobank", "age", "imputed_sex",
                                           "PC1", "PC2", "PC3", "PC4", "PC5",
                                           "pmi_hr"
                                         ),
                                         force_factor_covariates = c("biobank", "imputed_sex"),
                                         data_cache_dir = NULL) {
  .run_private_eqtl_tests(
    manifest_file = "eqtls/metadata/private_eqtl_analysis_manifest.txt",
    cell_type_file = "eqtls/metadata/region_cell_type.tsv",
    cluster_assignment_file = "eqtls/eqtl_analysis_main_figures/cluster_assignments_qval_0.01_k13.tsv",
    eqtl_root = "eqtls/results/LEVEL_6",
    cache_subdir = "private_eqtl_analysis",
    covariates = covariates,
    force_factor_covariates = force_factor_covariates,
    data_cache_dir = data_cache_dir
  )
}


#' Run cell-type-specific private eQTL tests (downsampled data set)
#'
#' Wires \code{bican.mccarroll.eqtl::run_cell_type_specific_eqtl_tests_from_manifest()}
#' into the package's options-driven path resolution, using the downsampled
#' cluster assignments and eQTL results. Per-analysis outputs are written
#' under a \code{private_eqtl_analysis_downsampled} subdirectory of the
#' configured cache directory (see \code{\link{get_cache_dir}}).
#'
#' @inheritParams run_private_eqtl_tests_bican
#'
#' @return Invisibly returns the list produced by
#'   \code{bican.mccarroll.eqtl::run_cell_type_specific_eqtl_tests_from_manifest()}.
#'
#' @seealso
#'   \code{\link[bican.mccarroll.eqtl]{run_cell_type_specific_eqtl_tests_from_manifest}}
#'   \code{\link{summarize_private_eqtl_tests_downsampled_bican}}
#'   \code{\link{run_private_eqtl_tests_bican}}
#'
#' @export
run_private_eqtl_tests_downsampled_bican <- function(covariates = c(
                                                       "biobank", "age", "imputed_sex",
                                                       "PC1", "PC2", "PC3", "PC4", "PC5",
                                                       "pmi_hr"
                                                     ),
                                                     force_factor_covariates = c("biobank", "imputed_sex"),
                                                     data_cache_dir = NULL) {
  .run_private_eqtl_tests(
    manifest_file = "eqtls/metadata/private_eqtl_analysis_manifest.txt",
    cell_type_file = "eqtls/metadata/region_cell_type.tsv",
    cluster_assignment_file = "eqtls/eqtl_analysis_downsampled_figures/cluster_assignments_qval_0.01_k13.tsv",
    eqtl_root = "eqtls/results/LEVEL_6_downsampled",
    cache_subdir = "private_eqtl_analysis_downsampled",
    covariates = covariates,
    force_factor_covariates = force_factor_covariates,
    data_cache_dir = data_cache_dir
  )
}


.run_private_eqtl_tests <- function(manifest_file,
                                    cell_type_file,
                                    cluster_assignment_file,
                                    eqtl_root,
                                    cache_subdir,
                                    covariates,
                                    force_factor_covariates,
                                    data_cache_dir = NULL) {
  paths <- .resolve_private_eqtl_paths(
    manifest_file = manifest_file,
    cell_type_file = cell_type_file,
    cluster_assignment_file = cluster_assignment_file,
    eqtl_root = eqtl_root,
    cache_subdir = cache_subdir,
    data_cache_dir = data_cache_dir
  )

  result <- bican.mccarroll.eqtl::run_cell_type_specific_eqtl_tests_from_manifest(
    manifest_file = paths$manifest_file,
    output_dir = paths$output_dir,
    cell_type_file = paths$cell_type_file,
    cluster_assignment_file = paths$cluster_assignment_file,
    covariates = covariates,
    eqtl_root = paths$eqtl_root,
    force_factor_covariates = force_factor_covariates
  )

  logger::log_info("DONE running private eQTL tests to {paths$output_dir}")

  invisible(result)
}


#' Summarize private eQTL cell-type-specific test results (main data set)
#'
#' Combines the per-analysis \code{.summary.txt} files produced by
#' \code{\link{run_private_eqtl_tests_bican}} into a single plain-text
#' manuscript table. For each cluster, the CaH-region result is used, except
#' for cluster 9 (glutamatergic layers), which only has a DFC result. Rows
#' are ordered to match the manuscript table.
#'
#' For each of the 7 analyses included in the table, this also regenerates
#' the first two summary pages of that analysis's PDF report (fraction of
#' significant genes by out-group cell-type count, and in-group versus
#' pooled out-group slopes) directly from the cached
#' \code{.linear_model_test.txt} results and saves them as individual SVG
#' files, via \code{bican.mccarroll.eqtl::plot_fdr_fraction_by_outgroup_count()}
#' and \code{bican.mccarroll.eqtl::plot_in_vs_out_group_slopes()}. This is
#' cheap (no re-fitting of models) since it reads the results table already
#' written by \code{\link{run_private_eqtl_tests_bican}}.
#'
#' All outputs (the table and the SVGs) are written to a dedicated
#' \code{private_eqtl_analysis} subdirectory of the configured figure output
#' directory (see \code{\link{get_out_dir}}): \code{private_eqtl_summary.txt},
#' plus, per analysis, \code{<analysis_id>.fdr_fraction_by_outgroup.svg} and
#' \code{<analysis_id>.in_vs_out_group_slopes.svg}.
#'
#' @param data_cache_dir Cache directory containing the per-analysis outputs
#'   written by \code{\link{run_private_eqtl_tests_bican}}. If \code{NULL},
#'   resolved via configured cache directory options.
#' @param outDir Output directory under which the dedicated summary
#'   subdirectory is created. If \code{NULL}, resolved via configured output
#'   directory options.
#'
#' @return Invisibly returns the combined summary data.frame.
#'
#' @seealso
#'   \code{\link{run_private_eqtl_tests_bican}}
#'   \code{\link{summarize_private_eqtl_tests_downsampled_bican}}
#'
#' @export
summarize_private_eqtl_tests_bican <- function(data_cache_dir = NULL, outDir = NULL) {
  .summarize_private_eqtl_tests(
    cache_subdir = "private_eqtl_analysis",
    output_filename = "private_eqtl_summary.txt",
    data_cache_dir = data_cache_dir,
    outDir = outDir
  )
}


#' Summarize private eQTL cell-type-specific test results (downsampled data set)
#'
#' Combines the per-analysis \code{.summary.txt} files produced by
#' \code{\link{run_private_eqtl_tests_downsampled_bican}} into a single
#' plain-text manuscript table, and regenerates the per-analysis summary SVGs,
#' exactly as \code{\link{summarize_private_eqtl_tests_bican}} does for the
#' main data set. Outputs are written to a dedicated
#' \code{private_eqtl_analysis_downsampled} subdirectory of the configured
#' figure output directory (see \code{\link{get_out_dir}}), as
#' \code{private_eqtl_summary_downsampled.txt} plus the per-analysis SVGs.
#'
#' @inheritParams summarize_private_eqtl_tests_bican
#'
#' @return Invisibly returns the combined summary data.frame.
#'
#' @seealso
#'   \code{\link{run_private_eqtl_tests_downsampled_bican}}
#'   \code{\link{summarize_private_eqtl_tests_bican}}
#'
#' @export
summarize_private_eqtl_tests_downsampled_bican <- function(data_cache_dir = NULL, outDir = NULL) {
  .summarize_private_eqtl_tests(
    cache_subdir = "private_eqtl_analysis_downsampled",
    output_filename = "private_eqtl_summary_downsampled.txt",
    data_cache_dir = data_cache_dir,
    outDir = outDir
  )
}


.summarize_private_eqtl_tests <- function(cache_subdir, output_filename, data_cache_dir = NULL, outDir = NULL) {
  paths <- .resolve_private_eqtl_paths(
    cache_subdir = cache_subdir,
    data_cache_dir = data_cache_dir,
    outDir = outDir
  )

  caH_summary_files <- list.files(
    paths$output_dir,
    pattern = "CaH\\.summary\\.txt$",
    full.names = TRUE
  )

  glut_summary_file <- file.path(
    paths$output_dir,
    "cluster9_glut_L23_IT-glut_L5_IT__DFC.summary.txt"
  )
  if (!file.exists(glut_summary_file)) {
    stop(
      "Expected glutamatergic (DFC) summary file not found: ",
      glut_summary_file
    )
  }

  summary_files <- c(caH_summary_files, glut_summary_file)

  summary_list <- lapply(
    summary_files,
    utils::read.table,
    header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE
  )
  combined <- do.call(rbind, summary_list)
  analysis_ids <- sub("\\.summary\\.txt$", "", basename(summary_files))

  # Cluster/region combinations included in the summary table, in manuscript
  # row order. For each cluster, the CaH-region result is used, except for
  # cluster 9 (glutamatergic layers), which only has a DFC result.
  cluster_order <- c(10, 13, 12, 11, 9, 6, 7)

  order_index <- match(cluster_order, combined$cluster)
  combined <- combined[order_index, ]
  analysis_ids <- analysis_ids[order_index]

  missing_clusters <- cluster_order[is.na(combined$cluster)]
  if (length(missing_clusters)) {
    stop(
      "Missing expected summary result(s) for cluster(s): ",
      paste(missing_clusters, collapse = ", ")
    )
  }

  for (analysis_id in analysis_ids) {
    .write_private_eqtl_summary_svgs(
      analysis_id = analysis_id,
      results_dir = paths$output_dir,
      figure_dir = paths$figure_dir
    )
  }

  output_file <- file.path(paths$figure_dir, output_filename)
  utils::write.table(
    combined, output_file,
    row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"
  )

  logger::log_info("Saved private eQTL summary table to {output_file}")

  invisible(combined)
}


.write_private_eqtl_summary_svgs <- function(analysis_id, results_dir, figure_dir) {
  results_file <- file.path(results_dir, paste0(analysis_id, ".linear_model_test.txt"))
  if (!file.exists(results_file)) {
    stop("Missing linear model results file for SVG generation: ", results_file)
  }

  results <- data.table::fread(results_file, sep = "\t")

  p_fdr <- bican.mccarroll.eqtl::plot_fdr_fraction_by_outgroup_count(results = results)
  if (!is.null(p_fdr)) {
    ggplot2::ggsave(
      file.path(figure_dir, paste0(analysis_id, ".fdr_fraction_by_outgroup.svg")),
      plot = p_fdr, device = svglite::svglite, width = 10, height = 6
    )
  }

  p_slopes <- bican.mccarroll.eqtl::plot_in_vs_out_group_slopes(results = results)
  if (!is.null(p_slopes)) {
    ggplot2::ggsave(
      file.path(figure_dir, paste0(analysis_id, ".in_vs_out_group_slopes.svg")),
      plot = p_slopes, device = svglite::svglite, width = 7, height = 7
    )
  }

  logger::log_info("Saved private eQTL summary SVGs for {analysis_id}")

  invisible(NULL)
}


.resolve_private_eqtl_paths <- function(manifest_file = NULL,
                                        cell_type_file = NULL,
                                        cluster_assignment_file = NULL,
                                        eqtl_root = NULL,
                                        cache_subdir,
                                        data_cache_dir = NULL,
                                        outDir = NULL) {
  root <- .resolve_data_root_dir(NULL)

  pick_in <- function(x) {
    if (is.null(x)) {
      return(NULL)
    }
    .resolve_under_root(root, x)
  }

  cache <- .resolve_cache_dir(data_cache_dir)
  output_dir <- file.path(cache, cache_subdir)
  .ensure_dir(output_dir)

  out <- .resolve_out_dir(outDir)
  figure_dir <- file.path(out, cache_subdir)
  .ensure_dir(figure_dir)

  list(
    data_root_dir = root,
    manifest_file = pick_in(manifest_file),
    cell_type_file = pick_in(cell_type_file),
    cluster_assignment_file = pick_in(cluster_assignment_file),
    eqtl_root = pick_in(eqtl_root),
    output_dir = output_dir,
    data_cache_dir = cache,
    out_dir = out,
    figure_dir = figure_dir
  )
}
