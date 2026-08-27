## ------------------------------------------------------------------
## Apply a pre-fit donor age prediction model (see donor_age_prediction.R /
## age_model_fit_core.R) to an external expression dataset whose region
## label(s) may not match the model's training region vocabulary.
## ------------------------------------------------------------------

########################################################################
## PUBLIC API
## Top-level, file-in/file-out entry points meant to be called directly
## from an analysis script.
########################################################################

#' Predict donor ages across cell types for an external dataset and plot QC scatterplots
#'
#' Loads an external \code{DGEList} and a pre-fit age model, maps each region
#' present in the external data to the model's region vocabulary via
#' \code{region_map}, then for every (region, cell type) combination present
#' in both datasets, predicts donor ages and writes one scatter plot per
#' cell type (per requested variant) to \code{outPDFFile}.
#'
#' @param data_dir Directory containing the external DGEList files (see
#'   \code{loadDGEList()}).
#' @param data_prefix Filename prefix for the external DGEList. Default
#'   \code{"donor_rxn_DGEList"}.
#' @param model_file Path to a model coefficients file as written by
#'   \code{\link{write_age_outputs_all}} (e.g. \code{*_model_coefficients.txt}),
#'   loaded via \code{\link{load_models}}.
#' @param region_map data.frame with columns \code{external_region} and
#'   \code{model_region}, mapping each region label present in the external
#'   data to the region label used in the fitted model.
#' @param outPDFFile Path to write the PDF of per-cell-type scatter plots to.
#' @param cell_type_list Optional character vector restricting cell types to
#'   attempt. Default \code{NULL} = all cell types present in both the
#'   external data and the model (for the mapped region).
#' @param variants Character vector of prediction variants to run per cell
#'   type; any of \code{names(.age_variant_specs)}: \code{"raw"} (model's
#'   original per-gene training mean/sd, no intercept correction),
#'   \code{"intercept_corrected"} (same per-gene mean/sd, but the overall
#'   predicted-age offset corrected to the external cohort's mean
#'   chronological age), and \code{"mean_sd_intercept_corrected"} (per-gene
#'   mean/sd recentered using the external data itself, which also corrects
#'   the overall offset). Default all three.
#' @param title_prefix Character scalar prepended to each plot title (e.g. a
#'   dataset nickname). Default is the external region label.
#' @param age_col,age_divisor Column holding chronological age, and the value
#'   to divide it by to convert to decades (e.g. \code{10} if the external
#'   data stores age in years).
#' @param min_donors Minimum donors required (after library-size QC) to
#'   attempt a fit for a given cell type.
#' @param libsize_threshold_sd,libsize_bins Library-size outlier QC
#'   parameters (set \code{libsize_threshold_sd = NULL} to skip this QC step).
#' @param prior.count Prior count used when computing log-CPM.
#' @param predictions_file Optional path to write a long-format tab-delimited
#'   file of donor-level predictions across all (cell type, variant)
#'   combinations, with columns \code{cell_type}, \code{external_region},
#'   \code{model_region}, \code{variant}, \code{donor}, \code{age}, \code{pred}.
#'   Suitable as input to \code{\link{compare_age_predictions_to_reference}}.
#' @importFrom grDevices pdf dev.off
#'
#' @return Invisibly, a named list of per-(region, cell type, variant)
#'   results, each with \code{cell_type}, \code{override_mean_sd},
#'   \code{predictions}, and \code{metrics}.
#' @export
predict_age_external_by_celltype <- function(data_dir, data_prefix = "donor_rxn_DGEList",
                                             model_file, region_map,
                                             outPDFFile,
                                             cell_type_list = NULL,
                                             variants = names(.age_variant_specs),
                                             title_prefix = NULL,
                                             age_col = "age", age_divisor = 10,
                                             min_donors = 50,
                                             libsize_threshold_sd = 1.96, libsize_bins = 50,
                                             prior.count = 1,
                                             predictions_file = NULL) {
  variants <- match.arg(variants, names(.age_variant_specs), several.ok = TRUE)

  dge <- load_external_dge_for_age_prediction(data_dir, data_prefix, age_col = age_col, age_divisor = age_divisor)
  all_models <- as.data.frame(load_models(model_file), stringsAsFactors = FALSE)

  external_regions <- unique(as.character(dge$samples$region))
  results <- list()
  predictions_all <- list()

  logger::log_info(paste0("Writing external age prediction QC plots to: ", outPDFFile))
  grDevices::pdf(outPDFFile)

  for (external_region in external_regions) {
    model_region <- resolve_model_region(region_map, external_region)
    dge_region <- dge[, dge$samples$region == external_region, keep.lib.sizes = TRUE]
    model_region_df <- all_models[all_models$region == model_region, ]

    region_title_prefix <- if (is.null(title_prefix)) external_region else title_prefix

    cts <- intersect(unique(dge_region$samples$cell_type), unique(model_region_df$cell_type))
    if (!is.null(cell_type_list)) {
      cts <- intersect(cts, cell_type_list)
    }

    missing_cts <- setdiff(unique(dge_region$samples$cell_type), unique(model_region_df$cell_type))
    if (length(missing_cts) > 0) {
      logger::log_warn(paste0(
        "No model available at region [", model_region, "] for cell types: ",
        paste(missing_cts, collapse = ", "), "; skipping."
      ))
    }

    for (cell_type in cts) {
      for (variant in variants) {
        spec <- .age_variant_specs[[variant]]

        r <- predict_age_external_celltype(
          dge_region, model_region_df, cell_type,
          override_mean_sd = spec$override_mean_sd,
          min_donors = min_donors,
          libsize_threshold_sd = libsize_threshold_sd,
          libsize_bins = libsize_bins,
          prior.count = prior.count,
          correct_intercept = spec$correct_intercept
        )
        if (is.null(r)) {
          next
        }

        title_str <- paste0(region_title_prefix, " [", gsub("_", " ", variant, fixed = TRUE), "] ", cell_type, "\n")

        p <- plot_age_pred_external(r$predictions, cell_type, title_str = title_str)
        print(p)

        key <- paste(external_region, cell_type, variant, sep = "__")
        results[[key]] <- r

        if (!is.null(predictions_file)) {
          preds <- r$predictions
          preds$cell_type <- cell_type
          preds$external_region <- external_region
          preds$model_region <- model_region
          preds$variant <- variant
          predictions_all[[key]] <- preds[, c(
            "cell_type", "external_region", "model_region", "variant", "donor", "age", "pred"
          )]
        }
      }
    }
  }

  grDevices::dev.off()

  if (!is.null(predictions_file)) {
    logger::log_info(paste0("Writing external donor predictions to: ", predictions_file))
    utils::write.table(
      do.call(rbind, predictions_all),
      file = predictions_file, sep = "\t", quote = FALSE, row.names = FALSE
    )
  }

  invisible(results)
}


#' Compare external age-prediction accuracy to a reference set of donor predictions
#'
#' Reads a long-format donor predictions file written by
#' \code{predict_age_external_by_celltype(..., predictions_file = )} and a
#' reference donor predictions file (e.g. the in-sample, cross-validated
#' predictions from \code{\link{write_age_outputs_all}}), and computes
#' per-cell-type accuracy metrics from each so they can be compared directly.
#'
#' @param predictions_file Path to a file written by
#'   \code{\link{predict_age_external_by_celltype}} via its
#'   \code{predictions_file} argument (columns \code{cell_type}, \code{variant},
#'   \code{donor}, \code{age}, \code{pred}).
#' @param reference_file Path to a reference donor predictions file (e.g.
#'   \code{*_donor_predictions.txt} as written by
#'   \code{\link{write_age_outputs_all}}), with columns \code{cell_type},
#'   \code{region}, \code{donor}, \code{age}, and \code{reference_pred_col}.
#' @param model_region Region value to restrict \code{reference_file} to
#'   (e.g. \code{"DFC"}) before computing its metrics.
#' @param reference_pred_col Column in \code{reference_file} holding predicted
#'   age. Default \code{"pred_mean"}.
#' @param age_col Column name holding chronological age in both files.
#'   Default \code{"age"}.
#' @param plot_pdf_file Optional path. If provided, also writes a PDF whose
#'   first page is a \code{model_details} x \code{cell_type} heatmap overview
#'   (see \code{overview_metric}), followed by two pages per cell type: (1)
#'   predicted vs chronological age, one panel per variant in
#'   \code{variant_order} (default: the reference plus \code{raw},
#'   \code{intercept_corrected}, and \code{mean_sd_intercept_corrected} —
#'   exactly filling a 2x2 grid), and (2) the reference
#'   (independently-trained model, x axis) vs. each variant in
#'   \code{agreement_pairs} (y axis; no chronological age) — see
#'   \code{\link{plot_age_predictions_by_method}}. All pages are laid out as
#'   a 2x2 grid with unused panels left blank.
#' @param cell_type_list Optional character vector restricting which cell
#'   types are included in \code{plot_pdf_file}. Default \code{NULL} = all
#'   cell types that have at least one non-reference (external) prediction
#'   — i.e. excludes cell types present only in \code{reference_file} (e.g.
#'   BICAN-only subtypes), which would otherwise produce all-blank pages.
#' @param variant_order Character vector giving the variants to include on
#'   \code{plot_pdf_file}'s first per-cell-type page (and as heatmap
#'   columns), in order (up to 4 for the 2x2 page). Default
#'   \code{c("independent_model", "raw", "intercept_corrected",
#'   "mean_sd_intercept_corrected")}.
#' @param agreement_pairs List of one or more length-2 character vectors
#'   \code{c(variant_x, variant_y)} to compare on \code{plot_pdf_file}'s
#'   third page. Default plots the reference (\code{"independent_model"})
#'   on the x axis against each of \code{"raw"}, \code{"intercept_corrected"},
#'   and \code{"mean_sd_intercept_corrected"} on the y axis (3 panels).
#' @param overview_metric Column of the returned metrics to show in the
#'   heatmap overview page: \code{"median_abs_error"} (default) or
#'   \code{"mean_abs_error"}.
#' @param age_scale Multiplier applied throughout \code{plot_pdf_file} for
#'   display — the heatmap's \code{overview_metric}, and every scatter
#'   panel's ages/predictions and annotated error metrics (e.g. \code{10} to
#'   convert decades to years). Does not affect the returned/written
#'   \code{metrics_file} data.frame, which stays in the input files' native
#'   units (typically decades). Default \code{10}.
#' @param metrics_file Optional path to write the returned comparison
#'   data.frame to, tab-delimited.
#'
#' @return A long-format data.frame with columns \code{cell_type},
#'   \code{model_details}, \code{num_donors}, \code{r}, \code{median_abs_error},
#'   \code{mean_abs_error} (same age units as the input files). The
#'   reference rows are always labeled \code{"independent_model"}.
#' @export
compare_age_predictions_to_reference <- function(predictions_file, reference_file, model_region,
                                                 reference_pred_col = "pred_mean",
                                                 age_col = "age",
                                                 plot_pdf_file = NULL,
                                                 cell_type_list = NULL,
                                                 variant_order = c(.reference_variant_label, names(.age_variant_specs)),
                                                 agreement_pairs = lapply(names(.age_variant_specs), function(v) c(.reference_variant_label, v)),
                                                 overview_metric = c("median_abs_error", "mean_abs_error"),
                                                 age_scale = 10,
                                                 metrics_file = NULL) {
  overview_metric <- match.arg(overview_metric)

  combined <- combine_age_predictions_for_comparison(
    predictions_file, reference_file, model_region,
    reference_pred_col = reference_pred_col,
    age_col = age_col
  )

  metrics <- compute_age_metrics_by_group(combined, group_cols = c("cell_type", "variant"), age_col = age_col, pred_col = "pred")
  colnames(metrics)[colnames(metrics) == "variant"] <- "model_details"
  colnames(metrics)[colnames(metrics) == "n"] <- "num_donors"

  if (!is.null(plot_pdf_file)) {
    # Restrict to cell types with actual external predictions; the reference file
    # typically covers additional cell types (e.g. BICAN-only subtypes) that have
    # nothing to compare against and would otherwise produce all-blank pages.
    plot_cell_types <- cell_type_list
    if (is.null(plot_cell_types)) {
      plot_cell_types <- unique(combined$cell_type[combined$variant != .reference_variant_label])
    }

    plot_age_predictions_by_method_pdf(
      combined, metrics, outPDFFile = plot_pdf_file,
      cell_type_list = plot_cell_types, variant_order = variant_order,
      agreement_pairs = agreement_pairs,
      overview_metric = overview_metric, age_scale = age_scale
    )
  }

  if (!is.null(metrics_file)) {
    logger::log_info(paste0("Writing prediction comparison metrics to: ", metrics_file))
    utils::write.table(metrics, file = metrics_file, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  metrics
}


########################################################################
## HELPERS
## Internal workhorses used by the public API functions above: functions
## that read a file only to feed another method, per-cell-type/per-group
## workhorses, plotting, and small utilities. Not meant to be called
## directly when starting a new analysis script.
########################################################################

# Maps each user-facing variant name to the (override_mean_sd, correct_intercept)
# flags passed to predict_age_external_celltype(). The name itself (with
# underscores rendered as spaces) is used directly as the display title in
# both the QC plots and the 2x2 comparison pages — see gsub("_", " ", variant).
.age_variant_specs <- list(
  raw = list(override_mean_sd = FALSE, correct_intercept = FALSE),
  intercept_corrected = list(override_mean_sd = FALSE, correct_intercept = TRUE),
  mean_sd_intercept_corrected = list(override_mean_sd = TRUE, correct_intercept = TRUE)
)

# The reference is always a model trained independently on the same external cohort
# (e.g. via bican.mccarroll.figures:::.age_prediction_snap200()), not BICAN's own
# in-sample predictions transferred to it — fixed rather than caller-configurable so
# the label in comparison outputs is stable and always means the same thing.
.reference_variant_label <- "independent_model"


#' Load an external DGEList for age prediction, rescaling age to decades
#'
#' Wraps \code{loadDGEList()} and rescales \code{dge$samples[[age_col]]} so
#' that it is in the same units (decades) the age prediction models in this
#' package are trained on. \code{\link{predict_age_from_dge}} expects
#' \code{dge$samples$donor} and \code{dge$samples$age} specifically, so
#' \code{age_col} must already be (or be renamed to) \code{"age"}.
#'
#' @param data_dir Directory containing \code{<data_prefix>_counts.tsv.gz} and
#'   \code{<data_prefix>_samples.tsv.gz}, as written by \code{saveDGEList()}.
#' @param data_prefix Filename prefix passed to \code{loadDGEList()}.
#'   Default \code{"donor_rxn_DGEList"}.
#' @param age_col Column in \code{dge$samples} holding chronological age.
#'   Default \code{"age"}.
#' @param age_divisor Value to divide \code{age_col} by to convert to decades
#'   (e.g. \code{10} if the external data stores age in years). Set to
#'   \code{1} if the data is already in decades.
#'
#' @return A \code{DGEList} with raw counts and \code{age_col} rescaled.
#' @keywords internal
load_external_dge_for_age_prediction <- function(data_dir, data_prefix = "donor_rxn_DGEList",
                                                  age_col = "age", age_divisor = 10) {
  dge <- bican.mccarroll.differentialexpression::loadDGEList(dir = data_dir, prefix = data_prefix)

  if (!(age_col %in% colnames(dge$samples))) {
    stop("age_col '", age_col, "' not found in dge$samples.")
  }

  dge$samples[[age_col]] <- as.numeric(dge$samples[[age_col]]) / age_divisor
  dge
}


#' Look up the model region matching an external region label
#'
#' @param region_map A data.frame with columns \code{external_region} and
#'   \code{model_region}, mapping each external region label to the region
#'   label used in the fitted model's \code{region} column.
#' @param external_region Character scalar; the external region label to
#'   look up.
#'
#' @return Character scalar; the corresponding model region label.
#' @keywords internal
resolve_model_region <- function(region_map, external_region) {
  stopifnot(
    is.data.frame(region_map),
    all(c("external_region", "model_region") %in% colnames(region_map))
  )

  idx <- match(external_region, region_map$external_region)
  if (is.na(idx)) {
    stop(
      "No model_region mapping found for external_region = '", external_region,
      "'. Add a row to region_map."
    )
  }
  as.character(region_map$model_region[idx])
}


#' Predict donor ages for one cell type from external expression data
#'
#' Subsets \code{dge} to \code{cell_type}, optionally removes library-size
#' outliers via \code{\link{filter_by_libsize}}, then scores donors with
#' \code{\link{predict_age_from_dge}} using the pre-fit coefficients in
#' \code{model_final}.
#'
#' @param dge A \code{DGEList} of raw counts already restricted to a single
#'   region, with \code{dge$samples$cell_type}, \code{dge$samples$donor}, and
#'   \code{dge$samples$age} (in decades).
#' @param model_final data.frame of model coefficients for a single region
#'   (across one or more cell types), in the format produced by
#'   \code{\link{load_models}} (columns \code{cell_type}, \code{feature},
#'   \code{coef}, \code{mean_train}, \code{sd_train}, \code{intercept}).
#' @param cell_type Character scalar; cell type to subset \code{dge} and
#'   \code{model_final} to.
#' @param override_mean_sd Logical; if \code{TRUE}, recenters the model's
#'   per-gene mean/sd using this external \code{dge} itself before predicting
#'   (\code{predict_age_from_dge(..., override_model_params_dge = dge)}).
#'   If \code{FALSE}, uses the model's original per-gene training mean/sd
#'   unchanged. Either way, the overall predicted-age offset is separately
#'   corrected; see \code{correct_intercept}.
#' @param min_donors Minimum donors required (after library-size QC) to
#'   attempt a fit. Below this, returns \code{NULL} with a warning.
#' @param libsize_threshold_sd SD threshold passed to
#'   \code{\link{filter_by_libsize}} for library-size QC. Set to \code{NULL}
#'   to skip this QC step.
#' @param libsize_bins \code{bins} passed to \code{\link{filter_by_libsize}}.
#' @param prior.count Prior count passed to \code{\link{predict_age_from_dge}}.
#' @param correct_intercept Logical; if \code{TRUE} (default), shifts all
#'   predictions for this cell type so their mean equals the external
#'   cohort's mean chronological age. \code{predict_age_from_dge()} only
#'   recenters the intercept when \code{override_mean_sd = TRUE}; without
#'   this correction, \code{override_mean_sd = FALSE} predictions carry the
#'   training population's baseline age forward unchanged, so any systematic
#'   difference between the training and external populations' baseline
#'   expression shows up as a pure additive bias rather than genuine
#'   prediction error. This is a no-op when \code{override_mean_sd = TRUE}
#'   (already exactly centered by construction).
#'
#' @return A list with \code{cell_type}, \code{override_mean_sd},
#'   \code{predictions} (data.frame with \code{donor}, \code{age}, \code{pred}),
#'   and \code{metrics} (as returned by \code{compute_age_metrics()}); or
#'   \code{NULL} if this cell type was skipped.
#' @keywords internal
predict_age_external_celltype <- function(dge, model_final, cell_type,
                                          override_mean_sd = FALSE,
                                          min_donors = 50,
                                          libsize_threshold_sd = 1.96,
                                          libsize_bins = 50,
                                          prior.count = 1,
                                          correct_intercept = TRUE) {
  stopifnot("DGEList" %in% class(dge))

  # model_final may be a data.table (e.g. straight from load_models()); coerce to a plain
  # data.frame so `model_final$cell_type == cell_type` below can't be captured by
  # data.table's NSE (which would resolve the bare `cell_type` on the RHS to the table's
  # own "cell_type" column instead of this function's argument).
  model_final <- as.data.frame(model_final, stringsAsFactors = FALSE)

  dge_cell <- dge[, dge$samples$cell_type == cell_type, keep.lib.sizes = TRUE]
  if (ncol(dge_cell) == 0) {
    logger::log_warn(paste0(cell_type, ": no samples found in external data; skipping."))
    return(NULL)
  }

  if (!is.null(libsize_threshold_sd)) {
    r <- bican.mccarroll.differentialexpression::filter_by_libsize(
      dge_cell, threshold_sd = libsize_threshold_sd, bins = libsize_bins, strTitlePrefix = cell_type
    )
    dge_cell <- r$dge
  }

  if (ncol(dge_cell) < min_donors) {
    logger::log_warn(paste0(
      cell_type, ": only ", ncol(dge_cell), " donors after QC (need >= ", min_donors, "); skipping."
    ))
    return(NULL)
  }

  model_cell <- model_final[model_final$cell_type == cell_type & model_final$coef != 0, ]
  if (nrow(model_cell) == 0) {
    logger::log_warn(paste0(cell_type, ": no non-zero model coefficients available; skipping."))
    return(NULL)
  }

  override_dge <- if (override_mean_sd) dge_cell else NULL

  predictions <- predict_age_from_dge(
    dge_cell, model_cell,
    prior.count = prior.count,
    override_model_params_dge = override_dge
  )

  if (correct_intercept) {
    predictions$pred <- predictions$pred + (mean(predictions$age, na.rm = TRUE) - mean(predictions$pred, na.rm = TRUE))
  }

  metrics <- compute_age_metrics(predictions$pred, predictions$age)

  list(
    cell_type = cell_type,
    override_mean_sd = override_mean_sd,
    predictions = predictions,
    metrics = metrics
  )
}


#' Plot external donor age predictions for one cell type with r/MedAE/MAE annotation
#'
#' Reuses \code{\link{plot_age_predictions}} and overlays the correlation and
#' error metrics (computed via \code{\link{compute_age_metrics}}) in the
#' upper-left corner of the plot.
#'
#' @param predictions data.frame with columns \code{donor}, \code{age},
#'   \code{pred}, as returned in \code{predict_age_external_celltype()$predictions}.
#' @param cell_type Character scalar; used in the default title.
#' @param title_str Plot title. Defaults to \code{cell_type}.
#'
#' @return A ggplot object.
#' @keywords internal
plot_age_pred_external <- function(predictions, cell_type, title_str = cell_type) {
  metrics <- compute_age_metrics(predictions$pred, predictions$age)
  label <- paste0(
    "Pearson r = ", sprintf("%.3f", metrics$r),
    "\nMedian Abs. Error = ", sprintf("%.3f", metrics$median_abs_error),
    "\nMean Abs. Error = ", sprintf("%.3f", metrics$mean_abs_error)
  )

  p <- plot_age_predictions(predictions, cell_type, titleStr = title_str)

  plot_limits <- layer_scales(p)
  x_pos <- plot_limits$x$range$range[1]
  y_pos <- plot_limits$y$range$range[2]

  p + annotate(
    "text",
    x = x_pos, y = y_pos, label = label,
    hjust = 0, vjust = 1, size = 6, color = "black"
  )
}


#' Combine external and reference donor-level predictions into one long table
#'
#' Reads a long-format donor predictions file written by
#' \code{predict_age_external_by_celltype(..., predictions_file = )} and a
#' reference donor predictions file, and stacks them into a single per-donor
#' table with a common \code{variant} column, suitable for either metric
#' computation (\code{\link{compare_age_predictions_to_reference}}) or plotting
#' (\code{plot_age_predictions_by_method()}).
#'
#' @inheritParams compare_age_predictions_to_reference
#' @return A data.frame with columns \code{cell_type}, \code{variant},
#'   \code{donor}, \code{age}, \code{pred}.
#' @keywords internal
combine_age_predictions_for_comparison <- function(predictions_file, reference_file, model_region,
                                                    reference_pred_col = "pred_mean",
                                                    age_col = "age") {
  ext <- utils::read.table(predictions_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  ext <- ext[, c("cell_type", "variant", "donor", age_col, "pred")]

  ref <- as.data.frame(
    data.table::fread(reference_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE),
    stringsAsFactors = FALSE
  )
  ref <- ref[ref$region == model_region, ]
  ref$variant <- .reference_variant_label
  ref <- ref[, c("cell_type", "variant", "donor", age_col, reference_pred_col)]
  colnames(ref)[colnames(ref) == reference_pred_col] <- "pred"

  rbind(ext, ref)
}


#' Compute donor age-prediction accuracy metrics per group
#'
#' @param df data.frame containing at least \code{group_cols}, \code{age_col}, and \code{pred_col}.
#' @param group_cols Character vector of columns to group by.
#' @param age_col,pred_col Column names holding chronological and predicted age.
#' @return data.frame with \code{group_cols}, \code{n}, \code{r}, \code{median_abs_error},
#'   \code{mean_abs_error}.
#' @keywords internal
compute_age_metrics_by_group <- function(df, group_cols, age_col = "age", pred_col = "pred") {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  groups <- unique(df[, group_cols, drop = FALSE])

  rows <- lapply(seq_len(nrow(groups)), function(i) {
    keep <- Reduce(`&`, lapply(group_cols, function(g) df[[g]] == groups[[g]][i]))
    sub <- df[keep, ]
    m <- compute_age_metrics(sub[[pred_col]], sub[[age_col]])
    cbind(groups[i, , drop = FALSE], n = nrow(sub), m)
  })

  do.call(rbind, rows)
}


#' Pivot a long-format age-prediction comparison table to a wide display table
#'
#' @param comparison_df As returned by \code{\link{compare_age_predictions_to_reference}}.
#' @param metric Column of \code{comparison_df} to pivot into wide format.
#'   Default \code{"mean_abs_error"}.
#' @param age_scale Multiplier applied to \code{metric} for display (e.g.
#'   \code{10} to convert decades to years). Default \code{1} (no conversion).
#' @param variant_order Optional character vector giving the column order for
#'   the pivoted variants. Default is the order variants first appear in
#'   \code{comparison_df}.
#' @param sort_by Variant name to sort rows by (ascending). Default is the
#'   first variant in \code{variant_order}.
#' @param add_mean_row Logical; append a row with the column-wise mean.
#'   Default \code{TRUE}.
#'
#' @return A wide data.frame: one row per cell type, one column per variant.
#' @keywords internal
pivot_age_prediction_comparison <- function(comparison_df, metric = "mean_abs_error", age_scale = 1,
                                            variant_order = NULL, sort_by = NULL, add_mean_row = TRUE) {
  if (is.null(variant_order)) {
    variant_order <- unique(comparison_df$variant)
  }

  per_variant <- lapply(variant_order, function(v) {
    sub <- comparison_df[comparison_df$variant == v, c("cell_type", metric)]
    sub[[metric]] <- sub[[metric]] * age_scale
    colnames(sub)[2] <- v
    sub
  })

  wide <- Reduce(function(x, y) merge(x, y, by = "cell_type", all = TRUE), per_variant)

  sort_col <- if (!is.null(sort_by)) sort_by else variant_order[1]
  wide <- wide[order(wide[[sort_col]]), ]

  if (add_mean_row) {
    mean_row <- c(cell_type = "mean", lapply(wide[, variant_order, drop = FALSE], function(x) round(mean(x, na.rm = TRUE), 1)))
    wide <- rbind(wide, mean_row)
  }

  rownames(wide) <- NULL
  wide
}


#' Build one scatter panel of predicted vs chronological age for one variant
#'
#' @param df_variant data.frame already restricted to one \code{variant},
#'   with columns \code{age}, \code{pred} (in decades).
#' @param panel_title Panel title (typically the variant name).
#' @param age_scale Multiplier applied to \code{age}, \code{pred}, and the
#'   annotated error metrics for display (e.g. \code{10} to convert decades
#'   to years). Default \code{10}.
#'
#' @return A ggplot object, or \code{NULL} if \code{df_variant} has no rows.
#' @keywords internal
.plot_pred_vs_age_panel <- function(df_variant, panel_title, age_scale = 10) {
  if (nrow(df_variant) == 0) {
    return(NULL)
  }

  m <- compute_age_metrics(df_variant$pred, df_variant$age)
  label <- paste0(
    "r = ", sprintf("%.3f", m$r),
    "\nMedAE = ", sprintf("%.1f", m$median_abs_error * age_scale),
    "\nMAE = ", sprintf("%.1f", m$mean_abs_error * age_scale)
  )

  df_variant$age <- df_variant$age * age_scale
  df_variant$pred <- df_variant$pred * age_scale

  # Shared x/y range so the dashed y=x line is a true diagonal, not skewed by
  # independently-scaled axes.
  rng <- range(c(df_variant$age, df_variant$pred), na.rm = TRUE)
  pad <- 0.05 * diff(rng)
  lims <- c(rng[1] - pad, rng[2] + pad)

  # Make R CMD CHECK happy
  age <- pred <- NULL

  ggplot(df_variant, aes(x = age, y = pred)) +
    geom_point(size = 1.8, alpha = 0.7, color = "steelblue") +
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.6, color = "red") +
    annotate(
      "text",
      x = lims[1], y = lims[2],
      label = label, hjust = 0, vjust = 1, size = 3.2
    ) +
    coord_equal(xlim = lims, ylim = lims) +
    labs(title = panel_title, x = "Chronological Age [years]", y = "Predicted Age [years]") +
    theme_classic(base_size = 10)
}


#' Build one scatter panel comparing predictions from two variants (no chronological age)
#'
#' Inner-joins two variants on \code{donor} within one cell type and plots
#' one method's predicted age against the other's, to assess whether the
#' methods' predictions (and hence their errors relative to chronological
#' age) are correlated or largely independent noise.
#'
#' @param df data.frame already restricted to one cell type, with columns
#'   \code{variant}, \code{donor}, \code{pred} (in decades).
#' @param variant_x,variant_y Variant names to compare (x axis, y axis).
#' @param age_scale Multiplier applied to \code{pred} and the annotated mean
#'   absolute difference for display (e.g. \code{10} to convert decades to
#'   years). Default \code{10}.
#'
#' @return A ggplot object, or \code{NULL} if the two variants share no donors.
#' @keywords internal
.plot_pred_vs_pred_panel <- function(df, variant_x, variant_y, age_scale = 10) {
  dx <- df[df$variant == variant_x, c("donor", "pred")]
  dy <- df[df$variant == variant_y, c("donor", "pred")]
  colnames(dx)[2] <- "pred_x"
  colnames(dy)[2] <- "pred_y"
  m <- merge(dx, dy, by = "donor")

  if (nrow(m) == 0) {
    return(NULL)
  }

  r <- stats::cor(m$pred_x, m$pred_y, use = "pairwise.complete.obs")
  mean_abs_diff <- mean(abs(m$pred_x - m$pred_y), na.rm = TRUE) * age_scale
  label <- paste0("r = ", sprintf("%.3f", r), "\nMean Abs. Diff = ", sprintf("%.1f", mean_abs_diff))

  m$pred_x <- m$pred_x * age_scale
  m$pred_y <- m$pred_y * age_scale

  display_x <- gsub("_", " ", variant_x, fixed = TRUE)
  display_y <- gsub("_", " ", variant_y, fixed = TRUE)

  # Shared x/y range so the dashed y=x line is a true diagonal, not skewed by
  # independently-scaled axes.
  rng <- range(c(m$pred_x, m$pred_y), na.rm = TRUE)
  pad <- 0.05 * diff(rng)
  lims <- c(rng[1] - pad, rng[2] + pad)

  # Make R CMD CHECK happy
  pred_x <- pred_y <- NULL

  ggplot(m, aes(x = pred_x, y = pred_y)) +
    geom_point(size = 1.8, alpha = 0.7, color = "darkorange3") +
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.6, color = "red") +
    annotate(
      "text",
      x = lims[1], y = lims[2],
      label = label, hjust = 0, vjust = 1, size = 3.2
    ) +
    coord_equal(xlim = lims, ylim = lims) +
    labs(
      title = paste(display_x, "vs", display_y),
      x = paste0(display_x, " [years]"),
      y = paste0(display_y, " [years]")
    ) +
    theme_classic(base_size = 10)
}


#' Arrange up to 4 panels into a titled 2x2 page, padding blanks as needed
#'
#' @param panels List of up to 4 ggplot objects (\code{NULL} entries become
#'   blank panels).
#' @param page_title Title drawn above the 2x2 grid.
#'
#' @return A combined plot object (from \code{cowplot::plot_grid()}).
#' @keywords internal
.assemble_2x2_page <- function(panels, page_title) {
  blank <- ggplot() + theme_void()
  panels <- lapply(panels, function(p) if (is.null(p)) blank else p)
  if (length(panels) < 4) {
    panels <- c(panels, replicate(4 - length(panels), blank, simplify = FALSE))
  }

  grid <- cowplot::plot_grid(plotlist = panels[seq_len(4)], ncol = 2, nrow = 2)
  title_gg <- cowplot::ggdraw() + cowplot::draw_label(page_title, fontface = "bold", size = 14)
  cowplot::plot_grid(title_gg, grid, ncol = 1, rel_heights = c(0.08, 0.92))
}


#' Plot per-cell-type prediction comparisons across methods, as two 2x2 pages
#'
#' Given the combined per-donor table from
#' \code{\link{combine_age_predictions_for_comparison}}, builds two pages for
#' one cell type, each a 2x2 grid (unused panels left blank):
#' \enumerate{
#'   \item Predicted vs chronological age, one panel per variant in
#'         \code{variant_order} (up to 4).
#'   \item Predicted-vs-predicted agreement for each \code{(variant_x,
#'         variant_y)} pair in \code{agreement_pairs} (no chronological age;
#'         up to 4 pairs) — high correlation means the two methods agree
#'         (and any shared deviation from chronological age is systematic,
#'         not noise); low correlation means their errors are closer to
#'         independent noise.
#' }
#'
#' @param combined_df As returned by \code{\link{combine_age_predictions_for_comparison}}.
#' @param cell_type Character scalar; cell type to plot.
#' @param variant_order Optional character vector giving the variants to
#'   include on page 1, in order (up to 4). Default is all variants present,
#'   in first-seen order.
#' @param agreement_pairs List of one or more length-2 character vectors
#'   \code{c(variant_x, variant_y)} to plot on page 2 (up to 4 pairs).
#'   Default plots the reference (\code{"independent_model"}, x axis) against
#'   each of \code{"raw"}, \code{"intercept_corrected"}, and
#'   \code{"mean_sd_intercept_corrected"} (y axis) — 3 panels.
#' @param age_scale Multiplier applied to ages/predictions and their
#'   annotated error metrics for display (e.g. \code{10} to convert decades
#'   to years). Default \code{10}.
#'
#' @return A list of two plot objects: \code{age_comparison} and
#'   \code{method_agreement}.
#' @keywords internal
plot_age_predictions_by_method <- function(combined_df, cell_type, variant_order = NULL,
                                           agreement_pairs = lapply(names(.age_variant_specs), function(v) c(.reference_variant_label, v)),
                                           age_scale = 10) {
  df <- combined_df[combined_df$cell_type == cell_type, ]

  if (is.null(variant_order)) {
    variant_order <- unique(df$variant)
  }
  df <- df[df$variant %in% union(variant_order, unlist(agreement_pairs)), ]

  age_panels <- lapply(variant_order, function(v) {
    .plot_pred_vs_age_panel(df[df$variant == v, ], gsub("_", " ", v, fixed = TRUE), age_scale = age_scale)
  })
  page1 <- .assemble_2x2_page(age_panels, paste(cell_type, "- predicted vs chronological age [years]"))

  agreement_panels <- lapply(agreement_pairs, function(pair) .plot_pred_vs_pred_panel(df, pair[1], pair[2], age_scale = age_scale))
  page2 <- .assemble_2x2_page(agreement_panels, paste(cell_type, "- prediction agreement between methods"))

  list(age_comparison = page1, method_agreement = page2)
}


#' Plot a model_details x cell_type heatmap overview of one metric
#'
#' @param metrics_df data.frame with columns \code{cell_type},
#'   \code{model_details}, and \code{metric} (e.g. as returned by
#'   \code{\link{compare_age_predictions_to_reference}}).
#' @param metric Column of \code{metrics_df} to display. Default
#'   \code{"median_abs_error"}.
#' @param age_scale Multiplier applied to \code{metric} for display (e.g.
#'   \code{10} to convert decades to years). Default \code{10}.
#' @param variant_order Optional character vector giving column order/labels
#'   to restrict/order \code{model_details}, and the row-sort key (rows are
#'   sorted ascending by the first entry's value). Default is all
#'   \code{model_details} present, in first-seen order.
#'
#' @return A ggplot object.
#' @keywords internal
.plot_metrics_heatmap <- function(metrics_df, metric = "median_abs_error", age_scale = 10, variant_order = NULL) {
  if (is.null(variant_order)) {
    variant_order <- unique(metrics_df$model_details)
  }
  df <- metrics_df[metrics_df$model_details %in% variant_order, ]
  df$value <- df[[metric]] * age_scale

  # Group columns so a header strip distinguishes the independently-trained reference
  # from variants derived by applying an external model's (BICAN's) fitted weights.
  group_of <- function(v) if (identical(v, .reference_variant_label)) "Independently trained" else "External model weights"
  group_order <- unique(vapply(variant_order, group_of, character(1)))
  df$model_group <- factor(vapply(df$model_details, group_of, character(1)), levels = group_order)

  display_order <- gsub("_", " ", variant_order, fixed = TRUE)
  df$model_details <- factor(gsub("_", " ", df$model_details, fixed = TRUE), levels = display_order)

  sort_values <- df[df$model_details == display_order[1], c("cell_type", "value")]
  cell_type_order <- sort_values$cell_type[order(sort_values$value)]
  cell_type_order <- c(cell_type_order, setdiff(unique(df$cell_type), cell_type_order))
  df$cell_type <- factor(df$cell_type, levels = rev(cell_type_order))

  metric_label <- switch(metric,
    median_abs_error = "Median Abs. Error",
    mean_abs_error = "Mean Abs. Error",
    metric
  )
  unit_label <- switch(as.character(age_scale), "10" = "years", "1" = "decades", paste0("x", age_scale))

  # Make R CMD CHECK happy
  model_details <- cell_type <- value <- model_group <- NULL

  ggplot(df, aes(x = model_details, y = cell_type)) +
    geom_tile(fill = "white", color = "grey80") +
    geom_text(aes(label = sprintf("%.1f", value)), size = 3, color = "black") +
    facet_grid(~model_group, scales = "free_x", space = "free_x") +
    labs(title = paste0(metric_label, " [", unit_label, "] by cell type and model"), x = NULL, y = NULL) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold")
    )
}


#' Plot predicted-vs-chronological-age and cross-method comparisons for all cell types to a PDF
#'
#' Writes a \code{model_details} x \code{cell_type} heatmap overview as the
#' first page, then calls \code{plot_age_predictions_by_method()} once per
#' cell type present in \code{combined_df}, writing its two pages
#' (predicted-vs-chronological-age, then cross-method agreement).
#'
#' @param combined_df As returned by \code{\link{combine_age_predictions_for_comparison}}.
#' @param metrics_df As returned by \code{\link{compare_age_predictions_to_reference}}
#'   (columns \code{cell_type}, \code{model_details}, and \code{overview_metric}),
#'   used to build the first-page heatmap.
#' @param outPDFFile Path to write the PDF to.
#' @param cell_type_list Optional character vector restricting cell types to
#'   plot. Default \code{NULL} = all cell types present in \code{combined_df}.
#' @param variant_order,agreement_pairs Passed to
#'   \code{plot_age_predictions_by_method()} and \code{\link{.plot_metrics_heatmap}}
#'   (as its \code{variant_order}).
#' @param overview_metric Passed to \code{\link{.plot_metrics_heatmap}}.
#' @param age_scale Multiplier applied to ages/predictions/metrics for
#'   display throughout (heatmap and per-cell-type pages alike); e.g.
#'   \code{10} to convert decades to years. Default \code{10}.
#' @importFrom grDevices pdf dev.off
#'
#' @return Invisibly, \code{NULL}.
#' @keywords internal
plot_age_predictions_by_method_pdf <- function(combined_df, metrics_df, outPDFFile, cell_type_list = NULL,
                                               variant_order = NULL,
                                               agreement_pairs = lapply(names(.age_variant_specs), function(v) c(.reference_variant_label, v)),
                                               overview_metric = "median_abs_error", age_scale = 10) {
  cell_types <- unique(combined_df$cell_type)
  if (!is.null(cell_type_list)) {
    cell_types <- intersect(cell_types, cell_type_list)
  }

  logger::log_info(paste0("Writing per-method age prediction comparison plots to: ", outPDFFile))
  grDevices::pdf(outPDFFile)

  overview_df <- metrics_df[metrics_df$cell_type %in% cell_types, ]
  print(.plot_metrics_heatmap(overview_df, metric = overview_metric, age_scale = age_scale, variant_order = variant_order))

  for (cell_type in cell_types) {
    pages <- plot_age_predictions_by_method(
      combined_df, cell_type, variant_order = variant_order, agreement_pairs = agreement_pairs,
      age_scale = age_scale
    )
    print(pages$age_comparison)
    print(pages$method_agreement)
  }
  grDevices::dev.off()

  invisible(NULL)
}
