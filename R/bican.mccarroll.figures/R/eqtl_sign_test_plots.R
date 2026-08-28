# The .sign_test.txt files consumed by plot_eqtl_sign_test_snap200_bican()
# are produced upstream, per cell type, by running:
#
# cd /broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/eqtls/sign_test_snap200_controls/scripts
# ./run_sign_tests.sh ../metadata/sign_test_manifest.tsv
#
# which invokes the private SignTest tool once per cell type, e.g.:
#
# /broad/mccarroll/software/dropseq/priv/SignTest \
#     --INPUT /broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/eqtls/results/LEVEL_6/astrocyte__DFC/astrocyte__DFC.cis_qtl.txt.gz \
#     --UNFILTERED_EQTL_FILE /broad/mccarroll/dropulation/analysis/SNAP200_controls/eqtls/results/astrocyte__BA46/astrocyte__BA46.cis_qtl_pairs.txt.gz \
#     --FLIP_ALLELES \
#     --OUTPUT /broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/eqtls/sign_test_snap200_controls/astrocyte__DFC.sign_test.txt \
#     --UNFILTERED_GENE_COLUMN_LABEL phenotype_id \
#     --UNFILTERED_EFFECT_COLUMN_LABEL slope \
#     --UNFILTERED_VARIANT_COLUMN_LABEL variant_id

# plot_eqtl_sign_test_snap200_bican()


#' Parse SignTest output into an eQTL sign-concordance summary table
#'
#' Combines the per-cell-type rows written by the private \code{SignTest}
#' tool (see \code{eqtls/sign_test_snap200_controls/scripts/run_sign_tests.sh}
#' under the configured data root) into a summary table suitable for
#' \code{\link{.plot_eqtl_sign_concordance}}.
#'
#' @param r A data.frame/data.table combining all \code{.sign_test.txt}
#'   rows, with columns \code{FILE_PERMUTED}, \code{NUM_EQTLS_COMPARED},
#'   \code{NUM_EQTLS_SIGNS_MATCH}, \code{FRAC_EQTLS_SIGNS_MATCH}.
#'
#' @return A data.table with columns \code{cell_type}, \code{n_compared},
#'   \code{n_match}, \code{percent_match}, and \code{count_label}, ordered by
#'   \code{percent_match} ascending, with \code{cell_type} a factor in that
#'   order.
.compute_eqtl_sign_concordance <- function(r) {
  required_columns <- c(
    "FILE_PERMUTED",
    "NUM_EQTLS_COMPARED",
    "NUM_EQTLS_SIGNS_MATCH",
    "FRAC_EQTLS_SIGNS_MATCH"
  )

  missing_columns <- setdiff(required_columns, names(r))
  if (length(missing_columns) > 0L) {
    stop(
      "Missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  # Make R CMD CHECK Happy
  FILE_PERMUTED <- NUM_EQTLS_COMPARED <- NUM_EQTLS_SIGNS_MATCH <- NULL
  FRAC_EQTLS_SIGNS_MATCH <- percent_match <- n_match <- n_compared <- NULL

  plot_data <- data.table::as.data.table(r)[
    ,
    .(
      cell_type = sub(
        "__DFC.*$",
        "",
        basename(dirname(FILE_PERMUTED))
      ),
      n_compared = NUM_EQTLS_COMPARED,
      n_match = NUM_EQTLS_SIGNS_MATCH,
      percent_match = 100 * FRAC_EQTLS_SIGNS_MATCH
    )
  ]

  plot_data[, count_label := paste0(n_match, "/", n_compared)]

  plot_data <- plot_data[order(percent_match)]
  plot_data[, cell_type := factor(cell_type, levels = cell_type)]

  plot_data
}


#' Plot an eQTL sign-concordance summary barplot
#'
#' @param plot_data Output of \code{\link{.compute_eqtl_sign_concordance}}.
#' @param primary_label,secondary_label Display labels for the two datasets.
#' @param effect_name Label used in the title (e.g. \code{"eQTLs"}).
#'
#' @return A ggplot object.
.plot_eqtl_sign_concordance <- function(plot_data,
                                        primary_label = "Burger et al.",
                                        secondary_label = "Ling et al.",
                                        effect_name = "eQTLs") {
  # Make R CMD CHECK Happy
  percent_match <- cell_type <- count_label <- NULL

  title <- sprintf(
    "%s-significant %s: sign concordance in %s",
    primary_label,
    effect_name,
    secondary_label
  )

  # Matches the styling of plot_sign_test_summary() in
  # bican.mccarroll.differentialexpression (theme_bw, default geom_col fill,
  # coord_flip) so DE and eQTL sign-test barplots look identical.
  ylim_upper <- max(100, max(plot_data$percent_match, na.rm = TRUE) + 5)

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = cell_type, y = percent_match)
  ) +
    ggplot2::geom_col(color = "black") +
    ggplot2::geom_text(
      ggplot2::aes(label = count_label),
      hjust = -0.1,
      size = 4
    ) +
    ggplot2::coord_flip(ylim = c(0, ylim_upper)) +
    ggplot2::labs(
      x = NULL,
      y = "% same sign",
      title = title
    ) +
    ggplot2::theme_bw(base_size = 12)
}


#' Plot Burger et al. vs. Ling et al. (SNAP200) eQTL sign-test concordance
#'
#' Reads the per-cell-type \code{.sign_test.txt} files written by the
#' private \code{SignTest} tool (see
#' \code{eqtls/sign_test_snap200_controls/scripts/run_sign_tests.sh} under
#' the configured data root) comparing this work's eQTLs against the Ling et
#' al. SNAP200 eQTL set, and plots a sign-concordance summary barplot using
#' the same dataset display names as
#' \code{\link{plot_de_sign_test_snap200_age_bican}}. Writes the SVG and a
#' companion summary TSV to an \code{eqtl_sign_test_snap200} subdirectory of
#' the configured figure output directory (see \code{\link{get_out_dir}}).
#'
#' @param results_dir Directory containing the \code{*.sign_test.txt} files.
#'   If \code{NULL}, resolved as
#'   \code{eqtls/sign_test_snap200_controls/results} under the configured
#'   data root directory (see \code{\link{get_data_root_dir}}).
#' @param outDir Output directory for the generated SVG/TSV. If \code{NULL},
#'   resolved via configured output directory options.
#'
#' @return Invisibly returns the summary data.table used to build the plot.
#'
#' @export
plot_eqtl_sign_test_snap200_bican <- function(results_dir = NULL, outDir = NULL) {
  root <- .resolve_data_root_dir(NULL)
  results_dir <- if (is.null(results_dir)) {
    .resolve_under_root(root, "eqtls/sign_test_snap200_controls/results")
  } else {
    results_dir
  }

  sign_test_files <- list.files(
    results_dir,
    pattern = "\\.sign_test\\.txt$",
    full.names = TRUE
  )

  if (length(sign_test_files) == 0L) {
    stop("No .sign_test.txt files found in ", results_dir)
  }

  r <- do.call(
    rbind,
    lapply(sign_test_files, function(f) {
      utils::read.table(f, header = TRUE, stringsAsFactors = FALSE)
    })
  )

  plot_data <- .compute_eqtl_sign_concordance(r)

  p <- .plot_eqtl_sign_concordance(
    plot_data,
    primary_label = "Burger et al.",
    secondary_label = "Ling et al.",
    effect_name = "eQTLs"
  )

  out <- .resolve_out_dir(outDir)
  dataset_out_dir <- file.path(out, "eqtl_sign_test_snap200")
  .ensure_dir(dataset_out_dir)

  ggplot2::ggsave(
    filename = file.path(dataset_out_dir, "eqtl_sign_concordance.svg"),
    plot = p,
    device = "svg",
    width = 10,
    height = 7
  )

  utils::write.table(
    plot_data,
    file = file.path(dataset_out_dir, "eqtl_sign_concordance_summary.tsv"),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  logger::log_info("DONE plotting eQTL sign test comparison for snap200")

  invisible(plot_data)
}
