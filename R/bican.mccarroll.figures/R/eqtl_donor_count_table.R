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
# write_eqtl_donor_counts_bican()


#' Write a table of donor counts per eQTL cell type / region (BICAN)
#'
#' Wires \code{bican.mccarroll.eqtl::count_eqtl_donors()} into the configured
#' eQTL results directory for \code{level}, and writes the resulting table
#' (columns \code{cell_type}, \code{region}, \code{n_donor}) as
#' \code{eqtl_donor_counts_<LEVEL>.txt} directly in the configured figure
#' output directory (see \code{\link{get_out_dir}}).
#'
#' @param level Character scalar naming the eQTL filtering level to scan
#'   (e.g. \code{"LEVEL 6"}), resolved the same way as \code{baseline_name} in
#'   \code{\link{plot_eqtl_filtering_examples}}.
#' @param outDir Output directory for the TSV. If \code{NULL}, resolved via
#'   configured output directory options.
#'
#' @return Invisibly returns the data.frame produced by
#'   \code{bican.mccarroll.eqtl::count_eqtl_donors()}.
#'
#' @seealso \code{\link[bican.mccarroll.eqtl]{count_eqtl_donors}}
#' @export
write_eqtl_donor_counts_bican <- function(level = "LEVEL 6", outDir = NULL) {
  paths <- .resolve_eqtl_paths(baseline_name = level, outDir = outDir)

  df <- bican.mccarroll.eqtl::count_eqtl_donors(paths$baseline_eqtl_data_dir)

  output_file <- file.path(
    paths$outDir,
    sprintf("eqtl_donor_counts_%s.txt", gsub(" ", "_", level, fixed = TRUE))
  )

  utils::write.table(
    df,
    file = output_file,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  logger::log_info("Saved eQTL donor counts to {output_file}")

  invisible(df)
}
