# source("R/paths.R")
#
# options(
#     bican.mccarroll.figures.data_root_dir =
#         "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis",
#
#     bican.mccarroll.figures.out_dir =
#         "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/figure_repository"
# )
#
# plot_sex_de_by_chromosome_bican()
# plot_de_summary_barplot_bican_sea_ad_mtg_mmc()
# plot_de_summary_barplot_pmid_39402379()

#' Plot BICAN sex-effect differential expression by chromosome group
#'
#' Wires the BICAN LEVEL_6 sex/age region-interaction differential expression
#' results into
#' \code{bican.mccarroll.differentialexpression::plot_sex_de_by_chromosome()}.
#' Writes a single summary count-plot SVG and one signed-effect density SVG
#' per region (all cell types faceted within a single plot) to a
#' \code{sex_de_bican} subdirectory of the configured figure output
#' directory (see \code{\link{get_out_dir}}).
#'
#' @param outDir Output directory for the generated SVGs. If \code{NULL},
#'   resolved via configured output directory options.
#' @param alpha Adjusted P-value threshold used to select significant genes.
#'
#' @return Invisibly returns the list produced by
#'   \code{bican.mccarroll.differentialexpression::plot_sex_de_by_chromosome()}.
#'
#' @export
plot_sex_de_by_chromosome_bican <- function(outDir = NULL, alpha = 0.05) {
  .plot_sex_de_by_chromosome_figure(
    in_dir = paste0(
      "differential_expression/results/LEVEL_6/sex_age/",
      "cell_type_region_interaction_absolute_effects"
    ),
    file_pattern = "female_vs_male",
    dataset_id = "bican",
    outDir = outDir,
    alpha = alpha
  )
}


.plot_sex_de_by_chromosome_figure <- function(in_dir,
                                              file_pattern,
                                              dataset_id,
                                              contig_yaml_file = NULL,
                                              reduced_gtf_file = NULL,
                                              cellTypeListFile = NULL,
                                              outDir = NULL,
                                              alpha = 0.05) {
  paths <- resolve_de_summary_barplot_paths(
    in_dir = in_dir,
    contig_yaml_file = contig_yaml_file,
    reduced_gtf_file = reduced_gtf_file,
    cellTypeListFile = cellTypeListFile,
    outDir = outDir
  )

  dataset_out_dir <- file.path(
    paths$outDir,
    paste0("sex_de_", dataset_id)
  )
  .ensure_dir(dataset_out_dir)

  result <- bican.mccarroll.differentialexpression::plot_sex_de_by_chromosome(
    in_dir = paths$in_dir,
    file_pattern = file_pattern,
    contig_yaml_file = paths$contig_yaml_file,
    reduced_gtf_file = paths$reduced_gtf_file,
    cellTypeListFile = paths$cellTypeListFile,
    alpha = alpha,
    svg_barplot_file = file.path(dataset_out_dir, "sex_de_barplot.svg"),
    svg_density_dir = dataset_out_dir,
    svg_density_prefix = "sex_de_effect_size_"
  )

  logger::log_info("DONE plotting sex DE by chromosome for {dataset_id}")

  invisible(result)
}


#' Plot BICAN sea_ad_mtg_mmc differential expression effect counts
#'
#' Wires the BICAN LEVEL_6_sea_ad_mtg_mmc age differential expression results
#' into \code{bican.mccarroll.differentialexpression::barplot_de_results()},
#' restricted to the SEA-AD supertype cell types for which BICAN DE results
#' were generated. Writes a single age-effect SVG to the configured figure
#' output directory (see \code{\link{get_out_dir}}).
#'
#' Sex-effect counts by chromosome group already have dedicated handling in
#' \code{\link{plot_sex_de_by_chromosome_bican}}, so this function only
#' covers the \code{age} contrast.
#'
#' @param outDir Output directory for the generated SVG. If \code{NULL},
#'   resolved via configured output directory options.
#'
#' @return Invisibly returns the `ggplot` object produced by
#'   \code{bican.mccarroll.differentialexpression::barplot_de_results()}.
#'
#' @export
plot_de_summary_barplot_bican_sea_ad_mtg_mmc <- function(outDir = NULL) {
  .plot_de_summary_barplot_figure(
    in_dir = "differential_expression/results/LEVEL_6_sea_ad_mtg_mmc/sex_age/cell_type",
    file_pattern = "__age_DE_results\\.txt$",
    cellTypeListFile = "differential_expression/metadata/cell_types_for_sea_ad_mtg_mmc_plots.txt",
    out_file = "de_summary_barplot_bican_sea_ad_mtg_mmc_age.svg",
    outDir = outDir
  )
}


#' Plot PMID_39402379 (Gabitto et al. 2024) differential expression effect counts
#'
#' Wires the PMID_39402379 voom-like differential expression results into
#' \code{bican.mccarroll.differentialexpression::barplot_de_results()},
#' restricted to the SEA-AD supertype cell types for which BICAN DE results
#' were generated. Writes one SVG per contrast to the configured figure
#' output directory (see \code{\link{get_out_dir}}).
#'
#' @param contrast Character vector of one or more of \code{"ad_cps"},
#'   \code{"early_ad_cps"}, \code{"late_ad_cps"}, \code{"versus_all"}.
#' @param outDir Output directory for the generated SVGs. If \code{NULL},
#'   resolved via configured output directory options.
#'
#' @return Invisibly returns a named list (by contrast) of the `ggplot`
#'   objects produced by
#'   \code{bican.mccarroll.differentialexpression::barplot_de_results()}.
#'
#' @export
plot_de_summary_barplot_pmid_39402379 <- function(
  contrast = c("ad_cps", "early_ad_cps", "late_ad_cps", "versus_all"),
  outDir = NULL
) {
  results <- lapply(contrast, function(ct) {
    .plot_de_summary_barplot_figure(
      in_dir = "differential_expression/external_comparison_PMID_39402379/voom-like",
      file_pattern = sprintf("__MTG__%s_DE_results\\.txt$", ct),
      cellTypeListFile = "differential_expression/metadata/cell_types_for_sea_ad_mtg_mmc_plots.txt",
      out_file = sprintf("de_summary_barplot_PMID_39402379_%s.svg", ct),
      outDir = outDir
    )
  })
  names(results) <- contrast

  invisible(results)
}


.plot_de_summary_barplot_figure <- function(in_dir,
                                            file_pattern,
                                            cellTypeListFile = NULL,
                                            out_file,
                                            outDir = NULL) {
  paths <- resolve_de_summary_barplot_paths(
    in_dir = in_dir,
    cellTypeListFile = cellTypeListFile,
    outDir = outDir
  )

  result <- bican.mccarroll.differentialexpression::barplot_de_results(
    in_dir = paths$in_dir,
    file_pattern = file_pattern,
    cellTypeListFile = paths$cellTypeListFile,
    svg_output_file = file.path(paths$outDir, out_file)
  )

  logger::log_info("DONE plotting DE summary barplot: {out_file}")

  invisible(result)
}


resolve_de_summary_barplot_paths <- function(in_dir,
                                             contig_yaml_file = NULL,
                                             reduced_gtf_file = NULL,
                                             cellTypeListFile = NULL,
                                             outDir = NULL) {
  root <- .resolve_data_root_dir(NULL)

  rel <- list(
    contig_yaml_file =
      "metadata/GRCh38_ensembl_v43.contig_groups.yaml",
    reduced_gtf_file =
      "metadata/GRCh38_ensembl_v43.reduced.gtf.gz",
    cellTypeListFile =
      "differential_expression/metadata/cell_types_for_de_filtering_plot.txt"
  )

  pick_in <- function(x, key) {
    if (is.null(x)) {
      return(file.path(root, rel[[key]]))
    }
    .resolve_under_root(root, x)
  }

  out <- .resolve_out_dir(outDir)
  .ensure_dir(out)

  list(
    data_root_dir    = root,
    in_dir           = .resolve_under_root(root, in_dir),
    contig_yaml_file = pick_in(contig_yaml_file, "contig_yaml_file"),
    reduced_gtf_file = pick_in(reduced_gtf_file, "reduced_gtf_file"),
    cellTypeListFile = pick_in(cellTypeListFile, "cellTypeListFile"),
    outDir           = out
  )
}
