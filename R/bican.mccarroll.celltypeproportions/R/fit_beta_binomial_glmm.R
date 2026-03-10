#
# figures_cache_dir = "/broad/mccarroll/yooolivi/projects/bican/manuscript_1_figures/data_cache"
#
# sample_ctp <- read.table(
#   file.path(figures_cache_dir, "donor_region.annotation_sub_class_complete.ctp.txt"),
#   sep="\t", header=TRUE, stringsAsFactors = FALSE
# )
#
# sample_metadata <- read.table(
#   file.path(figures_cache_dir, "donor_region.sample_metadata.txt"),
#   sep="\t", header=TRUE, stringsAsFactors=FALSE
# )
#
# cell_type_col <- "annotation_sub_class_complete"
# fixed_effects <- c("age_decades", "sex", "brain_region_abbreviation_simple", "mean_frac_contamination")
# random_effects <- c("donor_external_id", "village")
# cell_types <- c("astrocyte", "microglia", "oligodendrocyte", "OPC")
#
# glmm_results <- fit_glmm_and_extract_fixed_effects(
#   sample_ctp = sample_ctp,
#   sample_metadata = sample_metadata,
#   cell_types = cell_types,
#   cell_type_col = cell_type_col,
#   fixed_effects = fixed_effects,
#   random_effects = random_effects
# )
#
# fixed_effect_volcano_plot(glmm_results$fixed_effects_df, "age_decades")
# fixed_effect_pvalue_barplot(glmm_results$fixed_effects_df, c("age_decades", "sexM"), cell_type_col)
#
# plot_donor_residual_correlation_scatter(
#   glmm_results$glmm_list[["astrocyte"]],
#   region1="CaH",
#   region2="Pu",
#   cell_type="astrocyte"
# )
#
# plot_donor_residual_correlation_heatmap(
#   glmm_results$glmm_list[["astrocyte"]],
#   regions=c("CaH", "Pu", "NAC", "ic", "DFC"),
#   cell_type="astrocyte"
# )

#' Combine cell type proportions and metadata dataframe.
#'
#' @param sample_ctp Dataframe containing cell type proportions for each sample.
#' @param sample_metadata Dataframe containing metadata for each sample.
#' @param cell_type_col Column name in `sample_ctp` that contains cell type labels.
#' @param cell_types Character vector of cell types to include in the combined dataframe.
#' @param min_nuclei Minimum number of nuclei for any of the specified cell types.
#' @param denominator_types If not using the total number of nuclei as the denominator,
#' specify which cell types to include in the denominator by providing a character
#' vector of cell type names. If NULL, the total number of nuclei for each sample
#' will be used as the denominator (i.e., n_other = total_nuclei - n_nuclei).
#'
#' @return Dataframe with combined cell type proportions and metadata for the specified cell types,
#' including a column for n_other to be used as the denominator in the GLMM.
prepare_data_for_glmm <- function(sample_ctp, sample_metadata, cell_type_col, cell_types, min_nuclei=0, denominator_types=NULL) {

  glmm_df <- sample_ctp |>
    dplyr::filter(.data[[cell_type_col]] %in% cell_types) |>
    dplyr::left_join(sample_metadata)

  if (min_nuclei != 0) {

    low_nuclei_samples <- glmm_df |>
      dplyr::filter(n_nuclei < min_nuclei) |>
      dplyr::select(sample_id) |>
      dplyr::distinct() |>
      dplyr::pull(sample_id)

    n_samples <- glmm_df$sample_id |> unique() |> length()
    n_low_nuclei_samples <- length(low_nuclei_samples)

    logger::log_info("Filtering out {n_low_nuclei_samples} samples with fewer than {min_nuclei} nuclei for any of the specified cell types (out of {n_samples} total samples)")

    glmm_df <- glmm_df |>
      dplyr::filter(!sample_id %in% low_nuclei_samples)
  }

  if (!is.null(denominator_types)) {

    logger::log_info("Calculating n_other for denominator types: {paste(denominator_types, collapse=', ')}")

    denominator_counts <- sample_ctp |>
      dplyr::filter(.data[[cell_type_col]] %in% denominator_types) |>
      dplyr::group_by(sample_id) |>
      dplyr::summarise(n_other = sum(n_nuclei))

    glmm_df <- glmm_df |>
      dplyr::left_join(denominator_counts, by="sample_id")

  } else {

    logger::log_info("Calculating n_other as total_nuclei - n_nuclei for each sample")

    glmm_df <- glmm_df |>
      dplyr::mutate(n_other = total_nuclei - n_nuclei)
  }

  return(glmm_df)
}


#' Fit beta-binomial GLMM to cell type counts.
#'
#' @param sample_df Dataframe containing cell type counts and metadata for each sample.
#' @param cell_type Cell type for which to fit the model.
#' @param cell_type_col Column name in `sample_df` that contains cell type labels
#' @param fixed_effects Character vector of column names in `sample_df` to include as fixed effects in the model.
#' @param random_effects Character vector of column names in `sample_df` to include as random effects.
#'
#' @return Fitted beta-binomial GLMM object for the specified cell type.
fit_beta_binomial_glmm <- function(sample_df, cell_type, cell_type_col, fixed_effects, random_effects) {

  cell_type_df <- sample_df |>
    dplyr::filter(.data[[cell_type_col]] == cell_type)

  formula <- as.formula(
    paste0("cbind(n_nuclei, n_other) ~ ",
           paste(fixed_effects, collapse = " + "), " + (1 | ",
           paste(random_effects, collapse = ") + (1 | "), ")")
    )

  # first fit a binomial
  cell_type_binomial_glmm <- glmmTMB::glmmTMB(formula, data = cell_type_df, family = stats::binomial(link = "logit"))

  # check for overdispersion in the binomial model fit
  rp <- residuals(cell_type_binomial_glmm, type = "pearson")
  phi_hat <- sum(rp^2) / df.residual(cell_type_binomial_glmm)

  # if overdispersion is present, fit a beta-binomial model instead
  if (phi_hat > 1.5) {
    logger::log_info("Overdispersion detected (phi_hat = {round(phi_hat, 2)}), fitting beta-binomial model for cell type {cell_type}")
    cell_type_glmm <- glmmTMB::glmmTMB(formula, data = cell_type_df, family = glmmTMB::betabinomial(link = "logit"))
  } else {
    # beta-binomial approximates binomial when overdispersion is low so we can
    # use the binomial model fit and avoid any convergence issues
    logger::log_info("No overdispersion detected (phi_hat = {round(phi_hat, 2)}), using binomial model for cell type {cell_type}")
    cell_type_glmm <- cell_type_binomial_glmm
  }

  return(cell_type_glmm)

}


#' Fit beta-binomial GLMMs for multiple cell types.
#'
#' @param sample_df Dataframe containing cell type counts and metadata for each sample.
#' @param cell_types Character vector of cell types for which to fit models.
#' @param cell_type_col Column name in `sample_df` that contains cell type labels.
#' @param fixed_effects Character vector of column names in `sample_df` to include as fixed effects in the models.
#' @param random_effects Character vector of column names in `sample_df` to include as random effects in the models.
#'
#' @return Named list of fitted beta-binomial GLMM objects for each specified cell type.
fit_many_beta_binomial_glmms <- function(sample_df, cell_types, cell_type_col, fixed_effects, random_effects) {

  glmm_list <- list()

  for (cell_type in cell_types) {
    glmm_list[[cell_type]] <- fit_beta_binomial_glmm(sample_df, cell_type, cell_type_col, fixed_effects, random_effects)
  }

  return(glmm_list)

}


#' Extract fixed effects from a fitted beta-binomial GLMM for a specific cell type.
#'
#' @param glmm Fitted beta-binomial GLMM object for a specific cell type.
#' @param cell_type Cell type corresponding to the fitted model.
#' @param cell_type_col Column name to use for the cell type in the output dataframe.
#'
#' @returns Dataframe containing the fixed effects estimates, confidence intervals, and cell type information for the specified model.
extract_beta_binomial_fixed_effects <- function(glmm, cell_type, cell_type_col) {

  fixed_effects_df <- broom.mixed::tidy(
      glmm,
      effects = "fixed",
      conf.int = TRUE
    ) |>
      as.data.frame()

  fixed_effects_df[[cell_type_col]] <- cell_type

  return(fixed_effects_df)

}

#' Extract fixed effects from fitted beta-binomial GLMMs for multiple cell types and combine into a single dataframe.
#'
#' @param glmm_list Named list of fitted beta-binomial GLMM objects for each cell type.
#' @param cell_type_col Column name to use for the cell type in the output dataframe
#' @param out_file Optional file path to save the combined fixed effects dataframe as a tab-delimited text file. If NULL, the dataframe will not be saved to a file.
#'
#' @return Dataframe containing the fixed effects estimates, confidence intervals, and cell type information for all specified models, combined into a single dataframe.
extract_many_beta_binomial_fixed_effects <- function(glmm_list, cell_type_col, out_file=NULL) {

  fixed_effects_dfs <- list()

  for (cell_type in names(glmm_list)) {
    fixed_effects_dfs[[cell_type]] <- extract_beta_binomial_fixed_effects(glmm_list[[cell_type]], cell_type, cell_type_col)
  }

  combined_fixed_effects_df <- dplyr::bind_rows(fixed_effects_dfs)

  if(!is.null(out_file)) {
    write.table(combined_fixed_effects_df, out_file,
                sep="\t", row.names=FALSE, quote=FALSE)
  }

  return(combined_fixed_effects_df)
}


#' Main workflow to fit beta-binomial GLMMs for multiple cell types and extract fixed effects into a combined dataframe.
#'
#' @param sample_ctp Dataframe containing cell type proportions for each sample.
#' @param sample_metadata Dataframe containing metadata for each sample.
#' @param cell_types Cell types for which to fit models and extract fixed effects.
#' @param cell_type_col Column name in `sample_df` that contains cell type labels, and to use for the cell type column in the output fixed effects dataframe.
#' @param fixed_effects Character vector of column names in `sample_df` to include as fixed effects in the models.
#' @param random_effects Character vector of column names in `sample_df` to include as random effects in the models.
#' @param out_file Optional file path to save the combined fixed effects dataframe as a tab-delimited text file.
#'
#' @return List containing the named list of fitted beta-binomial GLMM objects for each cell type and the combined fixed effects dataframe for all models.
fit_glmm_and_extract_fixed_effects <- function(sample_ctp, sample_metadata, cell_types, cell_type_col,
                                               fixed_effects, random_effects,
                                               min_nuclei=0, denominator_types=NULL, out_file=NULL) {

  sample_df <- prepare_data_for_glmm(sample_ctp, sample_metadata, cell_type_col, cell_types, min_nuclei, denominator_types)
  glmm_list <- fit_many_beta_binomial_glmms(sample_df, cell_types, cell_type_col, fixed_effects, random_effects)
  fixed_effects_df <- extract_many_beta_binomial_fixed_effects(glmm_list, cell_type_col, out_file)

  return(list(glmm_list=glmm_list, fixed_effects_df=fixed_effects_df))

}

#' Extract fitted values and residuals from a fitted beta-binomial GLMM.
#'
#' @param glmm_model Fitted beta-binomial GLMM object for a specific cell type.
#'
#' @return Dataframe containing the observed fractions, fitted fractions, and
#' residuals for each sample in the model fit, along with the original model data.
extract_beta_binomial_fitted_values <- function(glmm_model) {

  # get proportions from model fit
  model_data <- model.frame(glmm_model)
  response_matrix <- model_data[["cbind(n_nuclei, n_other)"]]
  observed_fraction <- response_matrix[, 1] / rowSums(response_matrix)
  model_data$observed_fraction <- observed_fraction

  # get marginal predictions from model fit
  # (i.e., predictions based on fixed effects only, without random effects)
  model_data$fitted_fraction <- predict(glmm_model, type = "response", re.form = NA)
  model_data$residual <- model_data$observed_fraction - model_data$fitted_fraction

  return(model_data)

}


#' Generate a volcano plot for a specific fixed effect from the combined fixed effects dataframe.
#'
#' @param fixed_effects_df Dataframe containing the fixed effects estimates, confidence intervals, and cell type information for all specified models.
#' @param fixed_effect Name of the fixed effect for which to generate the volcano plot
#'
#' @return ggplot object representing the volcano plot for the specified fixed effect,
#' showing effect size (log-odds) on the x-axis and -log10(p-value) on the y-axis,
#' with points colored by cell type and error bars representing confidence intervals.
fixed_effect_volcano_plot <- function(fixed_effects_df, fixed_effect) {

  plot_df <- fixed_effects_df |>
    dplyr::filter(term == fixed_effect)

  ggplot2::ggplot(plot_df,
         ggplot2::aes(x=estimate, y=-log10(p.value))) +
    ggplot2::geom_vline(xintercept=0, lty="dashed", col="black") +
    ggplot2::geom_point() +
    ggrepel::geom_text_repel(ggplot2::aes(label=!!dplyr::sym(cell_type_col)), size=3) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin=conf.low, xmax=conf.high)) +
    ggplot2::labs(title=paste("Volcano plot for", fixed_effect), x="Effect size (log-odds)", y="-log10(p-value)")

}

#' Generate a bar plot of -log10(p-values) for specified fixed effects across cell types.
#'
#' @param fixed_effects_df Dataframe containing the fixed effects estimates, confidence intervals, and cell type information for all specified models.
#' @param fixed_effects Character vector of fixed effect names for which to generate the bar
#' @param cell_type_col Column name in `fixed_effects_df` that contains cell type labels.
#'
#' @return ggplot object representing the bar plot of -log10(p-values) for the
#' specified fixed effects across cell types, with bars colored by fixed effect and grouped by cell type.
fixed_effect_pvalue_barplot <- function(fixed_effects_df, fixed_effects, cell_type_col) {

  plot_df <- fixed_effects_df |>
    dplyr::filter(term %in% fixed_effects)

  ggplot2::ggplot(plot_df,
         ggplot2::aes(x=!!rlang::sym(cell_type_col), y=-log10(p.value), fill=term)) +
    ggplot2::geom_bar(stat="identity", position=ggplot2::position_dodge()) +
    ggplot2::labs(
      title=paste("raw p-values for fixed effects:", paste(fixed_effects, collapse=", "))
    )

}


#' Extract donor-level residuals from a fitted beta-binomial GLMM and reshape into a wide format for correlation analysis.
#'
#' @param model Fitted beta-binomial GLMM object for a specific cell type.
#' @param brain_region_col Name of the column in the model data that corresponds to brain region.
#' @param donor_col Name of the column in the model data that corresponds to donor.
#'
#' @return Dataframe containing donor-level residuals for each brain region in wide format,
#' with one row per donor and columns for each brain region's residuals.
extract_donor_residuals <- function(model,
                                    brain_region_col="brain_region_abbreviation_simple",
                                    donor_col="donor_external_id") {

  model_fitted_values <- extract_beta_binomial_fitted_values(model)

  donor_residuals <- model_fitted_values |>
    dplyr::select(!!rlang::sym(donor_col), !!rlang::sym(brain_region_col), age_decades, residual) |>
    tidyr::pivot_wider(names_from = !!rlang::sym(brain_region_col), values_from = residual)

  return(donor_residuals)

}


#' Generate a scatter plot comparing donor-level residuals between two brain regions
#' for a specific cell type, with points colored by age and a linear regression line added.
#'
#' @param model Fitted beta-binomial GLMM object for a specific cell type.
#' @param region1 Name of the first brain region to compare (for x-axis).
#' @param region2 Name of the second brain region to compare (for y-axis).
#' @param cell_type Name of the cell type corresponding to the fitted model (for labeling)
#' @param brain_region_col Name of the column in the model data that corresponds to brain region.
#' @param donor_col Name of the column in the model data that corresponds to donor.
#'
#' @return ggplot object representing the scatter plot comparing donor-level residuals
#' between the two specified brain regions for the given cell type,
#' with points colored by age, a linear regression line added, and identity line.
plot_donor_residual_correlation_scatter <- function(model, region1, region2, cell_type,
                                                   brain_region_col="brain_region_abbreviation_simple",
                                                   donor_col="donor_external_id") {

  donor_residuals <- extract_donor_residuals(model, brain_region_col, donor_col)

  plot_df <- donor_residuals |>
    dplyr::select(!!rlang::sym(donor_col), age_decades, !!rlang::sym(region1), !!rlang::sym(region2)) |>
    na.omit()

  n_donors <- plot_df[[donor_col]] |> unique() |> length()
  min_residual <- min(plot_df |> dplyr::select(-!!rlang::sym(donor_col), -age_decades), na.rm = TRUE)
  max_residual <- max(plot_df |> dplyr::select(-!!rlang::sym(donor_col), -age_decades), na.rm = TRUE)

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = !!rlang::sym(region1), y = !!rlang::sym(region2))
  ) +
    ggplot2::theme_bw() +
    ggplot2::geom_abline(slope=1, intercept=0, lty="dashed", col="gray") +
    ggplot2::geom_point(ggplot2::aes(col=age_decades)) +
    ggplot2::geom_smooth(method="lm", col="red", se=FALSE) +
    ggpubr::stat_cor(method="spearman", cor.coef.name="rho") +
    ggplot2::labs(
      title=sprintf("%s abundance residuals", cell_type),
      subtitle=sprintf("N=%s donors", n_donors)
    ) +
    ggplot2::xlim(min_residual, max_residual) +
    ggplot2::ylim(min_residual, max_residual) +
    ggplot2::scale_color_viridis_c()

}


#' Generate a heatmap of Spearman correlations between donor-level residuals across
#'  multiple brain regions for a specific cell type.
#'
#' @param model Fitted beta-binomial GLMM object for a specific cell type.
#' @param regions Brain regions to include, ordered how they should appear in the heatmap.
#' @param cell_type Name of the cell type corresponding to the fitted model (for labeling)
#' @param brain_region_col Name of the column in the model data that corresponds to brain region.
#' @param donor_col Name of the column in the model data that corresponds to donor.
#'
#' @return ggplot object representing the heatmap of Spearman correlations between
#' donor-level residuals
plot_donor_residual_correlation_heatmap <- function(model, regions, cell_type,
                                                    brain_region_col="brain_region_abbreviation_simple",
                                                    donor_col="donor_external_id") {

  donor_residuals <- extract_donor_residuals(model, brain_region_col, donor_col)

  cor_df <- donor_residuals |>
    dplyr::select(-!!rlang::sym(donor_col), -age_decades) |>
    as.matrix() |>
    cor(method="spearman", use="pairwise.complete.obs") |>
    as.data.frame() |>
    tibble::rownames_to_column("region1") |>
    tidyr::pivot_longer(-region1, names_to="region2", values_to="spearman_rho")

  cor_df$region1 <- factor(cor_df$region1, levels=regions)
  cor_df$region2 <- factor(cor_df$region2, levels=regions)

  ggplot2::ggplot(
    cor_df,
    ggplot2::aes(x=region1, y=forcats::fct_rev(region2), fill=spearman_rho)
  ) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label=round(spearman_rho, 2))) +
    ggplot2::scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0, limits=c(-1, 1)) +
    ggplot2::labs(
      x=NULL,
      y=NULL,
      title=sprintf("Correlation of donor residuals for %s abundance", cell_type),
      fill="Spearman\ncorrelation"
    ) +
    ggplot2::theme(
      axis.ticks.x=ggplot2::element_blank(),
      axis.ticks.y=ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )

}

