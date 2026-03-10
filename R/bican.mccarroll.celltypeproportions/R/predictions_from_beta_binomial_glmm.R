
# figures_cache_dir = "/broad/mccarroll/yooolivi/projects/bican/manuscript_1_figures/data_cache"
#
# sample_ctp <- read.table(
#   file.path(figures_cache_dir, "donor_region.annotation.ctp.txt"),
#   sep="\t", header=TRUE, stringsAsFactors = FALSE
# )
#
# sample_metadata <- read.table(
#   file.path(figures_cache_dir, "donor_region.sample_metadata.txt"),
#   sep="\t", header=TRUE, stringsAsFactors=FALSE
# )
#
# sample_df <- prepare_data_for_glmm(sample_ctp, sample_metadata, cell_type_col="annotation", cell_types="OPC")
#
# opc_glmm <- fit_beta_binomial_glmm(
#   sample_df=sample_df,
#   cell_type="OPC",
#   cell_type_col="annotation",
#   fixed_effects=c("age_decades", "sex", "brain_region_abbreviation_simple", "mean_frac_contamination"), #, "single_cell_assay"),
#   random_effects=c("donor_external_id", "village")
# )
#
# opc_predictions <- predict_many_regions(opc_glmm, region_names=c("CaH", "Pu", "NAC", "ic", "DFC"))
#
# plot_predictions_over_data(opc_glmm, cell_type="OPC", regions=c("CaH", "Pu", "NAC", "ic", "DFC"))

#' Calculate the mode of a vector, ignoring NA values.
get_mode <- function(x) {
  ux <- na.omit(unique(x))
  ux[which.max(tabulate(match(x, ux)))]
}


#' Generate a prediction grid for a fitted model, varying specified
#' variables and holding others constant.
#'
#' @param model A fitted model object for COUNTS data (e.g., from glmmTMB).
#' @param vary A named list of variables to vary in the prediction grid,
#'  where each name corresponds to a variable in the model and each value is
#'   a vector of values to use for that variable.
#' @param n_points Number of points to use for numeric variables that are not
#'  specified in `vary`. Default is 100.
#'
#'  @return A data frame containing the prediction grid, with columns for each variable in the model.
make_prediction_grid <- function(model,
                                 vary = NULL,
                                 n_points = 100)
{

  data = model.frame(model)[,-c(1,2)] # drop response counts

  # 1. Get fixed-effect formula (drop random effects)
  form <- reformulas::nobars(formula(model))

  # 2. Get term labels
  terms_obj <- terms(form)
  vars <- attr(terms_obj, "term.labels")

  # Remove interactions (optional: keeps main effects only)
  vars <- unique(unlist(strsplit(vars, ":")))

  newdata_list <- list()

  for (v in vars) {

    # If user wants to vary this variable
    if (!is.null(vary) && v %in% names(vary)) {
      newdata_list[[v]] <- vary[[v]]
      next
    }

    if (is.numeric(data[[v]])) {
      newdata_list[[v]] <- median(data[[v]], na.rm = TRUE)
    } else if (is.factor(data[[v]])) {
      mode_val <- get_mode(data[[v]])
      newdata_list[[v]] <- factor(mode_val, levels = levels(data[[v]]))
    } else {
      newdata_list[[v]] <- get_mode(data[[v]])
    }
  }

  expand.grid(newdata_list, KEEP.OUT.ATTRS = FALSE)

}




#' Generate predictions from a fitted model for a specific brain region,
#' varying age and holding other variables constant.
#'
#' @param model A fitted model object for COUNTS data (e.g., from glmmTMB).
#' @param region_name Name of the brain region to generate predictions for.
#' @param region_col Name of the column in the model data that corresponds to brain region
#'
#' @return A data frame containing the prediction grid with columns for age,
#' brain region, predicted proportion, and standard errors.
predict_one_region <- function(model, region_name, region_col="brain_region_abbreviation_simple") {

  df <- model.frame(model)

  vary_list <- rlang::list2(
    age_decades = seq(min(df$age_decades), max(df$age_decades), by = 0.1),
    !!region_col := region_name
  )

  region_data_to_predict <- make_prediction_grid(
    model,
    vary = vary_list
  )

  predictions <- predict(model, newdata = region_data_to_predict, type = "response", se.fit=TRUE, re.form=NA)

  region_data_to_predict$fit <- predictions$fit
  region_data_to_predict$se.fit <- predictions$se.fit
  region_data_to_predict$fit_lower <- region_data_to_predict$fit - 1.96 * region_data_to_predict$se.fit
  region_data_to_predict$fit_upper <- region_data_to_predict$fit + 1.96 * region_data_to_predict$se.fit

  return(region_data_to_predict)
}


#' Generate predictions from a fitted model for multiple brain regions,
#' varying age and holding other variables constant.
#'
#' @param model A fitted model object for COUNTS data (e.g., from glmmTMB).
#' @param region_names A vector of brain region names to generate predictions for.
#' @param region_col Name of the column in the model data that corresponds to brain region.
#'
#' @return A data frame containing the prediction grid with columns for age,
#' brain region, predicted proportion, and confidence intervals for all specified brain regions.
#'
predict_many_regions <- function(model, region_names, region_col="brain_region_abbreviation_simple") {
  prediction_list <- lapply(region_names, function(region) {
    predict_one_region(model, region)
  })

  names(prediction_list) <- region_names

  combined_predictions <- dplyr::bind_rows(prediction_list, .id=region_col)

  return(combined_predictions)
}


#' Plot observed data and model predictions for a specific cell type across multiple brain regions.
#'
#' This function takes a fitted model, extracts the observed data, calculates the
#' fraction of nuclei for the specified cell type, generates predictions from the model
#' for the specified brain regions, and creates a faceted plot comparing observed
#' data points with model predictions and confidence intervals. Age is in years (not decades).
#'
#' @param model A fitted model object for COUNTS data (e.g., from glmmTMB).
#' @param cell_type The name of the cell type being analyzed (used for labeling the y-axis)
#' @param regions A vector of brain region names to include in the plot.
#' @param region_col The name of the column in the model data that corresponds to brain region
#'
#' @return A ggplot object showing observed data points and model predictions for the
#' specified cell type across the specified brain regions.
plot_predictions_over_data <- function(model, cell_type,
                                       regions=c("CaH", "Pu", "NAC", "ic", "DFC"),
                                       region_col="brain_region_abbreviation_simple") {

  df <- as.data.frame(model.frame(model))
  resp <- model.response(df)
  df$fraction_nuclei <- resp[, 1] / rowSums(resp)
  df[[region_col]] <- factor(df[[region_col]], levels = regions)
  df <- df[, -1]

  predictions <- predict_many_regions(model, regions, region_col=region_col)
  predictions[[region_col]] <- factor(predictions[[region_col]], levels = regions)

  ggplot2::ggplot(
    df |> na.omit(),
    ggplot2::aes(x=age_decades*10, y=fraction_nuclei)
  ) +
    ggplot2::geom_point() +
    ggplot2::geom_line(data=predictions, ggplot2::aes(x=age_decades*10, y=fit), color="blue", lwd=1) +
    ggplot2::geom_ribbon(
      data=predictions,
      ggplot2::aes(x=age_decades*10, y=fit, ymin=fit_lower, ymax=fit_upper),
      alpha=0.2
    ) +
    ggplot2::facet_wrap(~brain_region_abbreviation_simple, scales="free", nrow=1) +
    ggpubr::stat_cor(method="spearman", cor.coef.name = "rho") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x="age",
      y=paste(cell_type, "fraction")
    )

}
