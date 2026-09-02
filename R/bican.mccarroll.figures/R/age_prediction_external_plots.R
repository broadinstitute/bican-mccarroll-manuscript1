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

# Model-transfer test: applies the BICAN DFC age-prediction model's *fitted
# weights* to the external SNAP200/BA46 data (as opposed to age_pred_figures.R's
# .age_prediction_snap200*() drivers, which train a new model on SNAP200), then
# compares those transferred predictions against a model trained independently
# on SNAP200 itself (age_pred_figures.R's .age_prediction_snap200(), which writes
# its donor predictions to data_cache/SNAP200_age_snap_DE_genes/age_prediction/).
# See bican.mccarroll.differentialexpression::predict_age_external_by_celltype()
# and ::build_age_prediction_comparison_plots().

.age_prediction_external_metacell_dir <-
  "/broad/mccarroll/dropulation/analysis/SNAP200_controls/metacells/out"

.age_prediction_external_subdir <- "SNAP200_age_bican_weights"

# SVG page size [inches] per page type. Page names from
# build_age_prediction_comparison_plots() are "metrics_heatmap" or
# "<cell_type>__<page_type>"; see .age_prediction_external_page_type().
.age_prediction_external_page_dims <- data.frame(
  page_type = c("metrics_heatmap", "age_comparison", "method_agreement"),
  width = c(8, 12, 6),
  height = c(8, 6, 6),
  stringsAsFactors = FALSE
)

#' Extract the page-type key from a build_age_prediction_comparison_plots() page name
#' @noRd
.age_prediction_external_page_type <- function(page_name) {
  sub("^.*__", "", page_name)
}

#' Compare BICAN-model-transfer age predictions to an independently-trained model
#'
#' Applies the BICAN \code{model_region} age-prediction model's fitted weights
#' to the external SNAP200/\code{external_region} data (caching the resulting
#' donor-level predictions), then compares those transferred predictions
#' against predictions from a model trained independently on SNAP200 itself.
#' Writes one SVG per comparison page (a \code{model_details} x
#' \code{cell_type} metrics heatmap, then two pages per cell type), a colored
#' Blues-gradient version of that same heatmap (in the style of
#' \code{\link[bican.mccarroll.differentialexpression]{plot_age_metric_heatmap}}),
#' a broader summary metrics TSV, and a TSV of the exact values plotted in
#' the (white-tile) metrics heatmap SVG.
#'
#' @param metacell_dir Directory containing the external SNAP200 DGEList
#'   files (see \code{loadDGEList()}). Default is the SNAP200 controls
#'   metacell directory under \code{/broad/mccarroll/dropulation}; this path
#'   is outside the configured data root and so is not resolved from it.
#' @param data_prefix Filename prefix for the external DGEList. Default
#'   \code{"donor_rxn_DGEList"}.
#' @param model_file Path to the BICAN model coefficients file (as written by
#'   \code{write_age_outputs_all()}). If \code{NULL}, resolved to
#'   \code{age_prediction/age_prediction_results_model_coefficients.txt}
#'   under the configured cache directory — the cache populated by any
#'   default-variant age-prediction wrapper, e.g.
#'   \code{\link{age_prediction_mean_residual_correlation_plots}}.
#' @param reference_file Path to the independently-trained SNAP200 model's
#'   donor-predictions file. If \code{NULL}, resolved to
#'   \code{SNAP200_age_snap_DE_genes/age_prediction/age_prediction_results_donor_predictions.txt}
#'   under the configured cache directory — the cache populated by
#'   \code{age_pred_figures.R}'s \code{.age_prediction_snap200()} driver.
#' @param external_region Region label as it appears in the external DGEList
#'   samples, and the region filter applied to \code{reference_file}.
#'   Default \code{"BA46"}.
#' @param model_region BICAN region whose fitted model weights to transfer.
#'   Default \code{"DFC"}.
#' @param min_donors Minimum donors required (after library-size QC) to
#'   attempt a fit for a given cell type. Default \code{50}.
#' @param outDir Output directory for generated SVGs. If \code{NULL},
#'   resolved via the configured output directory option. SVGs are written
#'   to a \code{SNAP200_age_bican_weights} subdirectory of this directory,
#'   alongside \code{age_external_metrics_heatmap.txt} (the heatmap's plotted
#'   values); the broader summary metrics TSV is written to the cache
#'   directory alongside the cached predictions.
#'
#' @return Invisibly, the comparison metrics data.frame (as returned by
#'   \code{bican.mccarroll.differentialexpression::build_age_prediction_comparison_plots()}).
#' @export
age_prediction_external_snap200_plots <- function(metacell_dir = NULL,
                                                   data_prefix = "donor_rxn_DGEList",
                                                   model_file = NULL,
                                                   reference_file = NULL,
                                                   external_region = "BA46",
                                                   model_region = "DFC",
                                                   min_donors = 50,
                                                   outDir = NULL) {
  paths <- .resolve_age_prediction_external_paths(outDir)

  if (is.null(metacell_dir)) {
    metacell_dir <- .age_prediction_external_metacell_dir
  }
  if (is.null(model_file)) {
    model_file <- paths$model_file_default
  }
  if (is.null(reference_file)) {
    reference_file <- paths$reference_file_default
  }

  if (!file.exists(model_file)) {
    stop(
      "Model coefficients file not found: ", model_file,
      ". Run an age_prediction_* wrapper with default paths (e.g. ",
      "age_prediction_mean_residual_correlation_plots()) first to populate this cache.",
      call. = FALSE
    )
  }
  if (!file.exists(reference_file)) {
    stop(
      "Reference donor-predictions file not found: ", reference_file,
      ". Run age_pred_figures.R's .age_prediction_snap200() driver first to populate this cache.",
      call. = FALSE
    )
  }

  predictions_file <- file.path(paths$cache_dir, "age_external_predictions.txt")

  if (file.exists(predictions_file)) {
    logger::log_info("Using cached data from {predictions_file}")
  } else {
    logger::log_info(
      "No cached data from {predictions_file}; regenerating from sources. This can take a while"
    )

    region_map <- data.frame(
      external_region = external_region,
      model_region = model_region,
      stringsAsFactors = FALSE
    )

    bican.mccarroll.differentialexpression::predict_age_external_by_celltype(
      data_dir = metacell_dir,
      data_prefix = data_prefix,
      model_file = model_file,
      region_map = region_map,
      outPDFFile = NULL,
      min_donors = min_donors,
      predictions_file = predictions_file
    )
  }

  res <- bican.mccarroll.differentialexpression::build_age_prediction_comparison_plots(
    predictions_file = predictions_file,
    reference_file = reference_file,
    model_region = external_region,
    reference_group_label = "Trained on Ling et al.",
    external_group_label = "Burger et al. model weights",
    variant_order = c("independent_model", "intercept_corrected"),
    agreement_pairs = list(c("independent_model", "intercept_corrected")),
    variant_labels = list(intercept_corrected = "Burger et al. weights")
  )

  utils::write.table(
    res$metrics,
    file = file.path(paths$cache_dir, "age_external_predictions_comparison_metrics.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  # Reshape to match the heatmap's visual layout exactly: rows top-to-bottom
  # as plotted (ggplot's first discrete y level is drawn at the bottom, so
  # reverse it), columns headed by the facet strip label actually shown in
  # the SVG (model_group; model_details itself is blanked on the x-axis).
  # Metrics/metadata columns not shown in the SVG (num_donors, r, ...) are
  # dropped.
  heatmap_df <- res$metrics_heatmap_data
  cell_type_levels <- rev(levels(heatmap_df$cell_type))
  model_levels <- levels(heatmap_df$model_details)
  group_of_model <- function(m) {
    unique(as.character(heatmap_df$model_group[as.character(heatmap_df$model_details) == m]))
  }

  heatmap_wide <- data.frame(cell_type = cell_type_levels, stringsAsFactors = FALSE)
  for (model_label in model_levels) {
    sub <- heatmap_df[heatmap_df$model_details == model_label, ]
    values <- sub$value[match(cell_type_levels, as.character(sub$cell_type))]
    # Match the "%.1f" rounding .render_metrics_heatmap() uses for the tile labels.
    heatmap_wide[[group_of_model(model_label)]] <- sprintf("%.1f", values)
  }

  utils::write.table(
    heatmap_wide,
    file = file.path(paths$dataset_out_dir, "age_external_metrics_heatmap.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  for (page_name in names(res$plots)) {
    dims <- .age_prediction_external_page_dims[
      .age_prediction_external_page_dims$page_type == .age_prediction_external_page_type(page_name),
    ]
    save_plot_svg(
      res$plots[[page_name]],
      out_file = paste0("age_external_", .sanitize_age_prediction_external_name(page_name), ".svg"),
      out_dir = paths$dataset_out_dir,
      width = dims$width,
      height = dims$height
    )
  }

  # Same colored-tile heatmap style used for the BICAN cell-type x region MAE
  # heatmap (plot_age_metric_heatmap()), applied here with the comparison
  # group ("Trained on Ling et al" / "Burger et al. model weights") standing
  # in for region.
  error_heatmap <- bican.mccarroll.differentialexpression::plot_age_metric_heatmap(
    heatmap_df,
    metric = "value",
    cell_type_col = "cell_type",
    region_col = "model_group"
  ) +
    # plot_age_metric_heatmap() maps the (already gsub()-ed to character)
    # group label to a default alphabetically-ordered x scale; restore the
    # reference-then-external column order (Burger on the right).
    ggplot2::scale_x_discrete(limits = levels(heatmap_df$model_group)) +
    # It also re-sorts rows by the mean value across columns, rather than the
    # first-column-ascending order the white heatmap/TSV use; restore that
    # order so the two representations of the same data read the same way.
    ggplot2::scale_y_discrete(limits = levels(heatmap_df$cell_type)) +
    ggplot2::labs(x = NULL, fill = "Median Abs\nError (years)")

  save_plot_svg(
    error_heatmap,
    out_file = "age_external_error_heatmap.svg",
    out_dir = paths$dataset_out_dir,
    width = 4.5, height = 8
  )

  logger::log_info("DONE plotting age prediction from external model weights for SNAP200")

  invisible(res$metrics)
}

#' Turn a build_age_prediction_comparison_plots() page name into a filename fragment
#'
#' Page names are \code{"metrics_heatmap"} or \code{"<cell_type>__<page>"};
#' this collapses the \code{"__"} separator and replaces any character
#' unsafe for a filename (cell type names may contain spaces or slashes)
#' with \code{"_"}.
#'
#' @noRd
.sanitize_age_prediction_external_name <- function(x) {
  x <- gsub("__", "_", x, fixed = TRUE)
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

.resolve_age_prediction_external_paths <- function(outDir = NULL) {
  out <- .resolve_out_dir(outDir)
  .ensure_dir(out)

  dataset_out_dir <- file.path(out, .age_prediction_external_subdir)
  .ensure_dir(dataset_out_dir)

  cache <- .resolve_cache_dir(NULL)
  cache_dir <- file.path(cache, .age_prediction_external_subdir)
  .ensure_dir(cache_dir)

  list(
    outDir = out,
    dataset_out_dir = dataset_out_dir,
    cache_dir = cache_dir,
    model_file_default = file.path(cache, "age_prediction", "age_prediction_results_model_coefficients.txt"),
    reference_file_default = file.path(
      cache, "SNAP200_age_snap_DE_genes", "age_prediction", "age_prediction_results_donor_predictions.txt"
    )
  )
}
