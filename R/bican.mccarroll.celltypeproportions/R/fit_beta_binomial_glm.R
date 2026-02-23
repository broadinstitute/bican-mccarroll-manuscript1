# library(tidyverse)
# library(glmmTMB)
#
# sample_ctp <- read.table(
#   "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/cell_type_proportions/LEVEL_1/donor_region.annotation_most_specific.cell_type_proportions.txt",
#   sep="\t", header=TRUE, stringsAsFactors = FALSE
# )
#
# sample_metadata <- read.table(
#   "/broad/mccarroll/yooolivi/projects/bican/cell_type_proportions/CAP_freeze_3/beta_binomial_glm/outs/donor_region_metadata.txt",
#   sep="\t", header=TRUE, stringsAsFactors=FALSE
# )
#
# cell_type_col <- "annotation_most_specific"
# fixed_effects <- c("age_decades", "sex", "brain_region_abbreviation_simple", "z_PC1", "z_PC2", "z_PC3", "z_PC4", "z_PC5")
# random_effects <- c("donor_external_id", "villages")
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


#' Combine cell type proportions and metadata dataframe.
#'
#' @param sample_ctp Dataframe containing cell type proportions for each sample.
#' @param sample_metadata Dataframe containing metadata for each sample.
#' @param cell_type_col Column name in `sample_ctp` that contains cell type labels.
#'
#' @return Dataframe with combined cell type proportions and metadata, including a column for the number of nuclei of other cell types (total_nuclei - n_nuclei).
prepare_data_for_glm <- function(sample_ctp, sample_metadata, cell_type_col) {

  sample_df <- sample_ctp |>
    dplyr::left_join(sample_metadata) |>
    dplyr::mutate(n_other = total_nuclei - n_nuclei)

  return(sample_df)
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

  logger::log_info("Fitting beta-binomial GLM for {cell_type}")

  cell_type_glm <- glmmTMB::glmmTMB(formula, data = cell_type_df, family = glmmTMB::betabinomial(link = "logit"))

  return(cell_type_glm)

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
fit_glmm_and_extract_fixed_effects <- function(sample_ctp, sample_metadata, cell_types, cell_type_col, fixed_effects, random_effects, out_file=NULL) {

  sample_df <- prepare_data_for_glm(sample_ctp, sample_metadata, cell_type_col)
  glmm_list <- fit_many_beta_binomial_glmms(sample_df, cell_types, cell_type_col, fixed_effects, random_effects)
  fixed_effects_df <- extract_many_beta_binomial_fixed_effects(glmm_list, cell_type_col, out_file)

  return(list(glmm_list=glmm_list, fixed_effects_df=fixed_effects_df))

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



