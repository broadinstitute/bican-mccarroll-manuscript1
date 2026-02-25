#
# ctp_dir=file.path("/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis", "cell_type_proportions", "data", "LEVEL_2")
#
# sample_ctp <- read.table(
#   file.path(ctp_dir, "donor_region.annotation.cell_type_proportions.txt"),
#   sep="\t", header=TRUE, stringsAsFactors=FALSE
# )
#
# sample_metadata <- read.table(
#   file.path(ctp_dir, "donor_region_metadata.txt"),
#   sep="\t", header=TRUE, stringsAsFactors=FALSE
# ) |>
#   dplyr::mutate(
#     sample_id=paste(donor_external_id, brain_region_abbreviation_simple, sep="_"),
#     sex = ifelse(imputed_sex == 2, "F", "M"),
#     age_decades = age / 10,
#     z_PC1 = scale(PC1),
#     z_PC2 = scale(PC2),
#     z_PC3 = scale(PC3),
#     z_PC4 = scale(PC4),
#     z_PC5 = scale(PC5),
#     z_pmi_hr = scale(pmi_hr)
#   )
#
# sample_df <- prepare_data_for_glm(sample_ctp, sample_metadata, cell_type_col="annotation")
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

  region_data_to_predict <- make_prediction_grid(
    model,
    vary = list(
      age_decades = seq(min(df$age_decades),
                        max(df$age_decades),
                        by = 0.1),
      setNames(list(region_name), region_col)
  ))

  predictions <- predict(model, newdata = region_data_to_predict, type = "response", se.fit=TRUE, re.form=NA)

  region_data_to_predict$fit <- predictions$fit
  region_data_to_predict$se.fit <- predictions$se.fit

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

  combined_predictions <- dplyr::bind_rows(prediction_list, .id=region_col) |>
    dplyr::mutate(
      predicted_proportion = fit,
      predicted_proportion_lower = fit - 1.96 * se.fit,
      predicted_proportion_upper = fit + 1.96 * se.fit
    )

  return(combined_predictions)
}
