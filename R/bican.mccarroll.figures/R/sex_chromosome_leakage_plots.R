# source("R/paths.R")
#
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
# plot_sex_chromosome_leakage()

#' Plot BICAN sex-chromosome expression leakage
#'
#' Wires the BICAN LEVEL_6 sex-chromosome "leakage" analysis into
#' \code{bican.mccarroll.differentialexpression::compute_sex_chromosome_leakage_from_files()},
#' \code{bican.mccarroll.differentialexpression::plot_sex_chromosome_leakage_boxplot()}, and
#' \code{bican.mccarroll.differentialexpression::plot_leakage_numerator_denominator()}. Writes
#' three SVGs (the leakage-ratio boxplot and the two numerator/denominator
#' comparison plots) to a \code{sex_chromosome_leakage} subdirectory of the
#' configured figure output directory (see \code{\link{get_out_dir}}).
#'
#' The village-level and donor-level leakage tables are cached as plain-text
#' TSV files to speed up subsequent runs and to provide reviewer-inspectable
#' intermediate data; they are not primary outputs.
#'
#' @param force_recompute Logical scalar. If \code{TRUE}, ignore existing
#'   cache files, recompute the leakage tables, and overwrite the cache.
#'   Defaults to \code{FALSE}.
#'
#' @return Invisibly returns a list with elements \code{leakage},
#'   \code{donor_leakage}, and \code{plots} (a list of the three ggplot
#'   objects saved to SVG).
#'
#' @export
plot_sex_chromosome_leakage <- function(force_recompute = FALSE) {
  paths <- resolve_sex_leakage_paths()
  data_name <- "donor_rxn_DGEList"

  cache_file_average <- file.path(paths$data_cache_dir, "sex_chromosome_leakage_average.tsv")
  cache_file_donor <- file.path(paths$data_cache_dir, "sex_chromosome_leakage_donor.tsv")

  if (!force_recompute && file.exists(cache_file_average) && file.exists(cache_file_donor)) {
    leakage_df <- data.table::fread(cache_file_average)
    donor_leakage_df <- data.table::fread(cache_file_donor)
  } else {
    computed <- bican.mccarroll.differentialexpression::compute_sex_chromosome_leakage_from_files(
      data_dir = paths$data_dir,
      data_name = data_name,
      cellTypeListFile = paths$cellTypeListFile
    )
    leakage_df <- computed$leakage
    donor_leakage_df <- computed$donor_leakage

    utils::write.table(
      leakage_df,
      file = cache_file_average,
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE
    )
    utils::write.table(
      donor_leakage_df,
      file = cache_file_donor,
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE
    )
  }

  dataset_out_dir <- file.path(paths$outDir, "sex_chromosome_leakage")
  .ensure_dir(dataset_out_dir)

  p_box <- bican.mccarroll.differentialexpression::plot_sex_chromosome_leakage_boxplot(
    leakage_df,
    outSVG = NULL
  )
  nd_plots <- bican.mccarroll.differentialexpression::plot_leakage_numerator_denominator(leakage_df)

  save_plot_svg(
    p_box,
    out_file = "sex_chromosome_leakage_boxplot.svg",
    out_dir = dataset_out_dir,
    width = 11,
    height = 8
  )
  save_plot_svg(
    nd_plots$expression,
    out_file = "sex_chromosome_leakage_expression.svg",
    out_dir = dataset_out_dir,
    width = 11,
    height = 8
  )
  save_plot_svg(
    nd_plots$expression_by_cell_type,
    out_file = "sex_chromosome_leakage_expression_by_cell_type.svg",
    out_dir = dataset_out_dir,
    width = 26,
    height = 8
  )

  logger::log_info("DONE plotting sex chromosome leakage")

  invisible(list(
    leakage = leakage_df,
    donor_leakage = donor_leakage_df,
    plots = list(
      boxplot = p_box,
      expression = nd_plots$expression,
      expression_by_cell_type = nd_plots$expression_by_cell_type
    )
  ))
}

resolve_sex_leakage_paths <- function() {
  root <- .resolve_data_root_dir(NULL)
  out <- .resolve_out_dir(NULL)
  cache <- file.path(.resolve_cache_dir(NULL), "differential_expression")

  .ensure_dir(out)
  .ensure_dir(cache)

  list(
    data_root_dir = root,
    data_dir = file.path(root, "metacells/LEVEL_6"),
    cellTypeListFile = file.path(
      root,
      "differential_expression/metadata/cell_types_for_de_filtering_plot.txt"
    ),
    outDir = out,
    data_cache_dir = cache
  )
}
