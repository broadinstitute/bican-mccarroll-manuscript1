# library(tidyverse)
# library(variancePartition)
#
# sample_ctp <- read.table(
#   "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/cell_type_proportions/data/LEVEL_1/donor_region.annotation_most_specific.cell_type_proportions.txt",
#   sep="\t", header=TRUE, stringsAsFactors = FALSE
# )
#
# sample_metadata <- read.table(
#   "/broad/mccarroll/yooolivi/projects/bican/cell_type_proportions/CAP_freeze_3/beta_binomial_glm/outs/donor_region_metadata.txt",
#   sep="\t", header=TRUE, stringsAsFactors=FALSE
# )
#
# cell_types <- c("astrocyte", "microglia", "oligodendrocyte", "OPC")
#
# # no cell type specific QC covariates
# formula <- ~ age_decades + z_PC1 + z_PC2 + z_PC3 + z_PC4 + z_PC5 + z_pmi_hr + (1 | sex) + (1 | biobank) + (1 | hbcac_status) + (1 | villages) + (1 | brain_region_abbreviation_simple) + (1 | donor_external_id)
#
# vp_results <- format_and_run_variance_partition(
#   sample_ctp = sample_ctp,
#   sample_metadata = sample_metadata,
#   cell_type_col = "annotation_most_specific",
#   formula = formula,
#   cell_types = cell_types
# )
#
# variancePartition::plotVarPart(vp_results)
#
# vp_result_list <- list(vp_results, vp_results)
#
# combined_results <- combine_variance_partition_results(vp_result_list)
# plot_variance_partition_results(combined_results)
#
#
# # cell type specific QC covariates
# formula_qc <- ~ age_decades + z_PC1 + z_PC2 + z_PC3 + z_PC4 + z_PC5 + z_pmi_hr + z_mean_frac_contamination + z_mean_pct_intronic + (1 | sex) + (1 | biobank) + (1 | hbcac_status) + (1 | villages) + (1 | brain_region_abbreviation_simple) + (1 | donor_external_id)
#
# vp_results_by_cell <- run_variance_partition_per_cell_type(
#   sample_ctp = sample_ctp,
#   sample_metadata = sample_metadata,
#   formula = formula_qc,
#   cell_types = cell_types,
#   cell_type_col = "annotation_most_specific",
#   qc_covariate_cols = c("mean_frac_contamination", "mean_pct_intronic"),
#   pseudocount = 0.5
# )
#
# plot_variance_partition_results(vp_results_by_cell)


#' Compute logit transformations of cell type abundances.
#'
#' @param sample_ctp Dataframe containing cell type proportions for each sample,
#' @param cell_type_col Column name in `sample_ctp` that contains cell type labels.
#' @param pseudocount Value to add to zero counts to avoid issues with logit.
#'
#' @returns Dataframe of logit-transformed cell type proportions, with rows for each cell type and columns for each sample.
compute_logit_CTP <- function(sample_ctp, cell_type_col, pseudocount=0.5) {

  ctp_logit <- sample_ctp |>
    dplyr::select(sample_id, all_of(cell_type_col), n_nuclei) |>
    dplyr::mutate(n_nuclei = ifelse(n_nuclei == 0, pseudocount, n_nuclei)) |>
    dplyr::group_by(sample_id) |>
    dplyr::mutate(total_nuclei = sum(n_nuclei)) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      fraction_nuclei = n_nuclei / total_nuclei,
      logit_fraction_nuclei = log(fraction_nuclei / (1 - fraction_nuclei))
    ) |>
    dplyr::select(sample_id, all_of(cell_type_col), logit_fraction_nuclei) |>
    dplyr::arrange(sample_id) |>
    tidyr::pivot_wider(
      names_from = sample_id,
      values_from = logit_fraction_nuclei
    ) |>
    tibble::column_to_rownames(cell_type_col)

  return(ctp_logit)

}


#' Formats metadata for variance partition, ensuring sample IDs are rownames and ordered to match CTP matrix columns.
#'
#' @param sample_metadata_df Dataframe containing sample metadata, with a column for sample IDs.
#'
#' @return Dataframe with sample IDs as rownames, ordered to match the columns of the CTP matrix.
format_metadata_for_vp <- function(sample_metadata_df) {

  sample_metadata_ordered <- sample_metadata_df |>
    dplyr::arrange(sample_id) |>
    tibble::column_to_rownames("sample_id")

  return(sample_metadata_ordered)

}


#' Combines logit-transformed CTP and formatted metadata into a list for variance partition input.
#'
#' @param sample_ctp Dataframe containing cell type proportions for each sample.
#' @param sample_metadata Dataframe containing metadata for each sample.
#' @param cell_type_col Column name in `sample_ctp` that contains cell type labels.
#' @param pseudocount Value to add to zero counts to avoid issues with logit.
#'
#' @return List containing logit-transformed CTP matrix and formatted metadata dataframe, ready for variance partition analysis.
prepare_vp_input <- function(sample_ctp, sample_metadata, cell_type_col, pseudocount=0.5) {

  ctp_logit <- compute_logit_CTP(sample_ctp, cell_type_col, pseudocount)
  sample_metadata_ordered <- format_metadata_for_vp(sample_metadata)

  vp_input <- list(
    ctp_logit = ctp_logit,
    sample_metadata = sample_metadata_ordered
  )

  return(vp_input)

}


#' Runs variance partition analysis using the `variancePartition` package.
#'
#' @param vp_input List containing logit-transformed CTP matrix and formatted metadata dataframe.
#' @param formula Formula specifying the fixed and random effects for the variance partition model.
#' @param cell_types Optional character vector of cell types to include in the analysis.
#'
#' @return varPartResults object containing the results of the variance partition analysis,
#'  with variance explained for each variable and cell type.
run_variance_partition <- function(vp_input, formula, cell_types=NULL) {

  if(!is.null(cell_types)) {
    logger::log_info("Running variance partition for cell types: {paste(cell_types, collapse = ', ')}")
    ctp_logit <- vp_input$ctp_logit[rownames(vp_input$ctp_logit) %in% cell_types, ]
  } else {
    ctp_logit <- vp_input$ctp_logit
  }

  varPart <- variancePartition::fitExtractVarPartModel(ctp_logit, formula, vp_input$sample_metadata)

  return(varPart)

}


#' Formats input data and runs variance partition analysis.
#'
#' @param sample_ctp Dataframe containing cell type proportions for each sample.
#' @param sample_metadata Dataframe containing metadata for each sample.
#' @param cell_type_col Column name in `sample_ctp` that contains cell type labels.
#' @param formula Formula specifying the fixed and random effects for the variance partition model.
#' @param cell_types Optional character vector of cell types to include in the analysis.
#' @param pseudocount Value to add to zero counts to avoid issues with logit
#'
#' @return varPartResults object containing the results of the variance partition analysis,
format_and_run_variance_partition <- function(sample_ctp, sample_metadata, cell_type_col, formula, cell_types=NULL, pseudocount=0.5) {

  vp_input <- prepare_vp_input(sample_ctp, sample_metadata, cell_type_col, pseudocount)
  vp_results <- run_variance_partition(vp_input, formula, cell_types)

  return(vp_results)

}


#' Combines variance partition results from multiple analyses into a single dataframe for comparison.
#'
#' @param vp_results_list List of varPartResults objects from multiple variance partition analyses.
#'
#' @return Dataframe combining variance explained for each variable and cell type across all analyses,
#' with a column for cell type and columns for each variable's variance explained.
combine_variance_partition_results <- function(vp_results_list) {

  combined_results <- vp_results_list |>
    purrr::imap_dfr(function(x, i) {
      as.data.frame(x) |>
        tibble::rownames_to_column("cell_type")
    })

  return(combined_results)

}


#' Extracts cell type-specific QC covariates for a given cell type, scales them, and formats them for variance partition input.
#'
#' @param sample_ctp Dataframe containing cell type proportions and QC covariates for each sample.
#' @param cell_type Cell type for which to extract QC covariates.
#' @param cell_type_col Column name in `sample_ctp` that contains cell type labels.
#' @param qc_covariate_cols Character vector of column names in `sample_ctp` that contain the QC covariates to extract and scale.
#'
#' @return Dataframe containing sample IDs, cell type labels, and scaled QC covariate values for the specified cell type, formatted for variance partition input.
extract_cell_type_specific_qc_covariates <- function(sample_ctp, cell_type, cell_type_col, qc_covariate_cols) {

  cell_type_df <- sample_ctp |>
    dplyr::filter(.data[[cell_type_col]] == cell_type)

  cell_type_qc_covariates <- cell_type_df |>
    dplyr::select(sample_id, all_of(cell_type_col), all_of(qc_covariate_cols)) |>
    dplyr::distinct() |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(qc_covariate_cols),
        ~ as.numeric(scale(.)),
        .names = "z_{.col}"
      )
    ) |>
    dplyr::arrange(sample_id)

  return(cell_type_qc_covariates)

}


#' Prepares logit-transformed CTP and cell type-specific QC covariates for a single cell type, formatted for variance partition input.
#'
#' @param sample_ctp Dataframe containing cell type proportions and QC covariates for each sample.
#' @param sample_metadata Dataframe containing metadata for each sample.
#' @param cell_type Cell type for which to prepare the variance partition input.
#' @param cell_type_col Column name in `sample_ctp` that contains cell type labels.
#' @param qc_covariate_cols Character vector of column names in `sample_ctp` that contain the QC covariates to extract and scale for the specified cell type.
#' @param pseudocount Value to add to zero counts to avoid issues with logit.
#'
#' @return List containing the logit-transformed CTP matrix for the specified cell type
#' and the formatted metadata dataframe with cell type-specific QC covariates, ready for variance partition analysis.
prepare_vp_input_one_cell_type <- function(sample_ctp, sample_metadata, cell_type, cell_type_col, qc_covariate_cols, pseudocount=0.5) {

  ctp_logit <- compute_logit_CTP(sample_ctp, cell_type_col, pseudocount)
  metadata <- format_metadata_for_vp(sample_metadata)
  qc_covariates <- extract_cell_type_specific_qc_covariates(sample_ctp, cell_type, cell_type_col, qc_covariate_cols)

  # duplicate to get past number of rows check for VP
  ctp_logit_cell_type <- ctp_logit[c(cell_type, cell_type), ,]

  metadata_cell_type <- bind_cols(
    metadata[colnames(ctp_logit_cell_type), ],
    qc_covariates
  )

  return(list(
    ctp_logit = ctp_logit_cell_type,
    metadata = metadata_cell_type
  ))

}


#' Fits variance partition model for a single cell type, incorporating cell type-specific QC covariates.
#'
#' @param sample_ctp Dataframe containing cell type proportions and QC covariates for each sample.
#' @param sample_metadata Dataframe containing metadata for each sample.
#' @param formula Formula specifying the fixed and random effects for the variance partition model.
#' @param cell_type Cell type for which to fit the variance partition model.
#' @param cell_type_col Column name in `sample_ctp` that contains cell type labels.
#' @param qc_covariate_cols Character vector of column names in `sample_ctp` that contain the QC covariates to extract and scale for the specified cell type.
#' @param pseudocount Value to add to zero counts to avoid issues with logit.
#'
#' @return Dataframe containing the variance explained for each variable.
fit_vp_one_cell_type <- function(sample_ctp, sample_metadata, formula, cell_type, cell_type_col, qc_covariate_cols, pseudocount=0.5) {

  vp_input <- prepare_vp_input_one_cell_type(sample_ctp, sample_metadata, cell_type, cell_type_col, qc_covariate_cols, pseudocount)

  ctp_logit_cell_type <- vp_input$ctp_logit
  metadata_cell_type <- vp_input$metadata

  cell_type_varPart <- fitExtractVarPartModel(ctp_logit_cell_type, formula, metadata_cell_type)[cell_type,]

  return(cell_type_varPart)

}


#' Runs variance partition analysis for each specified cell type, incorporating cell type-specific QC covariates, and combines results into a single dataframe.
#'
#' @param sample_ctp Dataframe containing cell type proportions and QC covariates for each sample.
#' @param sample_metadata Dataframe containing metadata for each sample.
#' @param formula Formula specifying the fixed and random effects for the variance partition model.
#' @param cell_types Cell types for which to fit the variance partition model.
#' @param cell_type_col Column name in `sample_ctp` that contains cell type labels.
#' @param qc_covariate_cols Character vector of column names in `sample_ctp` that contain the QC covariates to extract and scale for the specified cell type.
#' @param pseudocount Value to add to zero counts to avoid issues with logit.
#'
#' @return Dataframe combining variance explained for each variable and cell type across all specified cell types.
run_variance_partition_per_cell_type <- function(sample_ctp, sample_metadata, formula, cell_types, cell_type_col, qc_covariate_cols, pseudocount=0.5) {

  vp_results_per_cell_type <- lapply(cell_types, function(cell_type) {
    fit_vp_one_cell_type(sample_ctp, sample_metadata, formula, cell_type, cell_type_col, qc_covariate_cols, pseudocount)
  })

  names(vp_results_per_cell_type) <- cell_types

  combined_results <- dplyr::bind_rows(vp_results_per_cell_type) |>
    tibble::rownames_to_column("cell_type")

  return(combined_results)

}


#' Plots variance explained for each variable across cell types for all the results.
#'
#' @param combined_results Dataframe combining variance explained for each variable and cell type across all analyses.
#'
#' @return ggplot of violin and boxplots of variance explained for each variable across cell types.
plot_variance_partition_results <- function(combined_results) {

  combined_results_long <- combined_results |>
    tidyr::pivot_longer(cols = -cell_type, names_to = "variable", values_to = "variance_explained") |>
    dplyr::mutate(
      variable = forcats::fct_relevel(variable, "Residuals", after = Inf)
    )

  base_colors <- scales::hue_pal()(length(levels(combined_results_long$variable)))
  names(base_colors) <- levels(combined_results_long$variable)
  base_colors["Residuals"] <- "gray70"

  ggplot2::ggplot(
    combined_results_long,
    ggplot2::aes(x = variable, y = variance_explained, fill=variable)) +
    ggplot2::geom_violin(scale="width") +
    ggplot2::geom_boxplot(width=0.1) +
    ggplot2::labs(
      y="Variance explained (%)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    ggplot2::scale_fill_manual(values = base_colors)
}

