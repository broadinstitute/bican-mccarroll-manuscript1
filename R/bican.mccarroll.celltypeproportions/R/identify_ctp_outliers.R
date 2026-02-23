#library(tidyverse)
#library(logger)

# sample_ctp <- read.table(
#   "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/cell_type_proportions/LEVEL_0/donor_region_village.annotation_most_specific.cell_type_proportions.txt",
#   sep="\t", header=TRUE, stringsAsFactors = FALSE
#)

#sample_neuron_ctp <- read.table(
#  "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/cell_type_proportions/LEVEL_0/donor_region_village.annotation_most_specific.neurons_only.all.cell_type_proportions.txt",
#  sep="\t", header=TRUE, stringsAsFactors = FALSE
#)

#unexpected_neuron_df <- read.table(
#  "/broad/mccarroll/yooolivi/projects/bican/manuscript_1/ctp_outliers/unexpected_neurons.annotation_most_specific.txt",
#  header = TRUE,
#  sep = "\t",
#  stringsAsFactors = FALSE
#)

#ctp_outliers <- identify_ctp_outliers(
#  ctp_df = sample_ctp,
#  cell_type_col = "annotation_most_specific",
#  region_col = "brain_region_abbreviation_simple",
#  unexpected_neuron_df = unexpected_neuron_df,
#  conformity_threshold = 3.75,
#  conformity_method = "pearson",
#  neuron_fraction_threshold = 0.2
#)


#' Identify samples to drop based on low nuclei counts.
#'
#' @param ctp_df Dataframe of cell type counts and proportions (one row per cell type).
#'
#' @return Dataframe with sample_id, total_nuclei, z_log10_total_nuclei, and nuclei_count_outlier (TRUE/FALSE).
identify_low_nuclei_count_outliers <- function(ctp_df) {

  sample_counts <- ctp_df |>
    dplyr::select(sample_id, total_nuclei) |>
    dplyr::distinct() |>
    dplyr::mutate(z_log10_total_nuclei = as.numeric(scale(log10(total_nuclei)))) |>
    dplyr::mutate(nuclei_count_outlier = z_log10_total_nuclei < -1.96)

  return(sample_counts)

}


#' Identify samples to drop based on unexpected neuron fraction.
#'
#' @param neuron_ctp_df Dataframe of neuron counts and proportions (out of neurons).
#' @param unexpected_neuron_types Vector of neuron types considered unexpected for the brain region.
#' @param cell_type_col Name of the column in neuron_ctp_df that contains the cell type.
#' @param threshold Fraction of unexpected neurons above which a sample is considered an outlier (default 0.2).
#'
#' @return Dataframe with sample_id, total_neurons, unexpected_neuron_fraction, and unexpected_neuron_outlier (TRUE/FALSE).
identify_unexpected_neuron_outliers <- function(neuron_ctp_df, unexpected_neuron_types, cell_type_col, threshold=0.2) {

  sample_unexpected_neurons <- neuron_ctp_df |>
    dplyr::mutate(unexpected_neuron = !!rlang::sym(cell_type_col) %in% unexpected_neuron_types) |>
    dplyr::group_by(sample_id, unexpected_neuron) |>
    dplyr::summarise(n_neurons = sum(n_nuclei), .groups = "drop") |>
    dplyr::group_by(sample_id) |>
    dplyr::mutate(
      total_neurons = sum(n_neurons),
      fraction_neurons = n_neurons / total_neurons,
      .groups = "drop"
      ) |>
    dplyr::filter(unexpected_neuron) |>
    dplyr::mutate(unexpected_neuron_outlier = fraction_neurons > threshold) |>
    dplyr::rename(unexpected_neuron_fraction=fraction_neurons) |>
    dplyr::select(sample_id, total_neurons, unexpected_neuron_fraction, unexpected_neuron_outlier)

  return(sample_unexpected_neurons)
}

#' Apply arcsin square root transformation to CTP fractions and pivot to wide format.
#'
#' @param ctp_df Dataframe of cell type counts and proportions (one row per cell type).
#' @param cell_type_col Name of the column in ctp_df that contains the cell
#'
#' @return Wide dataframe of arcsin square root transformed CTP fractions, with sample_id as rownames and cell types as columns.
arcsin_sqrt_transform_ctp <- function(ctp_df, cell_type_col) {

  arcsin_sqrt <- function(p) {
    return(asin(sqrt(p)))
  }

  transformed_df <- ctp_df |>
    dplyr::mutate(fraction_nuclei_transformed = arcsin_sqrt(fraction_nuclei))

  transformed_df_wide <- transformed_df |>
    dplyr::select(sample_id, !!rlang::sym(cell_type_col), fraction_nuclei_transformed) |>
    tidyr::pivot_wider(names_from = !!rlang::sym(cell_type_col), values_from = fraction_nuclei_transformed) |>
    as.data.frame()

  rownames(transformed_df_wide) <- transformed_df_wide$sample_id
  transformed_df_wide$sample_id <- NULL

  return(transformed_df_wide)
}


#' Compute conformity scores for each sample based on median correlation with other samples.
#'
#' @param transformed_ctp_df Wide dataframe of arcsin square root transformed CTP fractions, with sample_id as rownames and cell types as columns.
#' @param method Correlation method to use (default "pearson").
#' @param threshold Threshold for modified Z-score to identify conformity outliers (default 3.75).
#'
#' @return Dataframe with sample_id, conformity_score, modified_z_conformity_score, and conformity_outlier (TRUE/FALSE).
compute_ctp_conformity_scores <- function(transformed_ctp_df, method = "pearson", threshold=3.75) {

  samples <- rownames(transformed_ctp_df)

  mat <- transformed_ctp_df[samples, , drop = FALSE]
  cor_mat <- cor(t(mat), method = method)

  diag(cor_mat) <- NA

  sample_conformity_df <- data.frame(
    sample_id = rownames(cor_mat),
    conformity_score = apply(cor_mat, 1, median, na.rm = TRUE)
  ) |>
    dplyr::mutate(
      median = median(conformity_score, na.rm = TRUE),
      mad = mad(conformity_score, na.rm = TRUE),
      modified_z_conformity_score = (conformity_score - median) / mad,
    ) |>
    dplyr::select(-median, -mad) |>
    dplyr::mutate(conformity_outlier = abs(modified_z_conformity_score) > abs(threshold))

  return(sample_conformity_df)

}



#' Identify outliers in CTP data for one brain region based on low nuclei counts, unexpected neuron fractions, and conformity scores.
#'
#' @param ctp_region_df Dataframe of CTP data for one brain region (one row per cell type).
#' @param cell_type_col Name of the column in ctp_region_df that contains the cell type label.
#' @param conformity_threshold Threshold for modified Z-score to identify conformity outliers (default 3.75).
#' @param conformity_method Correlation method to use for conformity scores (default "pearson")
#' @param neuron_ctp_df Optional dataframe of neuron-only CTP data for the same brain region (one row per neuron cell type).
#' @param unexpected_neuron_types Optional vector of neuron types considered unexpected for the brain region.
#' @param neuron_fraction_threshold Threshold for fraction of unexpected neurons above which a sample is considered an outlier (default 0.2).
#'
#' @return Dataframe with sample_id, total_nuclei, z_log10_total_nuclei, nuclei_count_outlier, total_neurons, unexpected_neuron_fraction, unexpected_neuron_outlier, conformity_score, modified_z_conformity_score, conformity_outlier, and overall_outlier (TRUE/FALSE).
identify_ctp_outliers_one_region <- function(ctp_region_df, cell_type_col, conformity_threshold=3.75, conformity_method="pearson", neuron_ctp_df=NULL, unexpected_neuron_types=NULL, neuron_fraction_threshold=0.2) {

  # identify low nuclei count outliers
  nuclei_count_outliers_df <- identify_low_nuclei_count_outliers(ctp_region_df)
  low_nuclei_outliers <- nuclei_count_outliers_df |> dplyr::filter(nuclei_count_outlier) |> dplyr::pull(sample_id)

  # identify unexpected neuron outliers (if neuron CTP provided)
  if (!is.null(neuron_ctp_df) && !is.null(unexpected_neuron_types)) {
    logger::log_info("Unexpected neurons provided, identifying outliers based on unexpected neuron fraction with threshold {neuron_fraction_threshold}")
    unexpected_neuron_outliers_df <- identify_unexpected_neuron_outliers(neuron_ctp_df, unexpected_neuron_types, cell_type_col, neuron_fraction_threshold)
    unexpected_neuron_outliers <- unexpected_neuron_outliers_df |>
      dplyr::filter(unexpected_neuron_outlier) |>
      dplyr::pull(sample_id)

    ctp_region_df_filtered <- ctp_region_df |>
      dplyr::filter(!(sample_id %in% low_nuclei_outliers)) |>
      dplyr::filter(!(sample_id %in% unexpected_neuron_outliers))
  } else {
    logger::log_info("No unexpected neurons provided, no outliers will be identified based on unexpected neuron fraction")
    unexpected_neuron_outliers_df <- data.frame(
      sample_id = unique(ctp_region_df$sample_id),
      total_neurons = NA,
      unexpected_neuron_fraction = NA,
      unexpected_neuron_outlier = FALSE
    )

    ctp_region_df_filtered <- ctp_region_df |>
      dplyr::filter(!(sample_id %in% low_nuclei_outliers))
  }

  # transform CTP and compute conformity scores (dropping low nuclei and neuron outliers first)
  transformed_ctp_df <- arcsin_sqrt_transform_ctp(ctp_region_df, cell_type_col)
  conformity_scores_df <- compute_ctp_conformity_scores(transformed_ctp_df, method = conformity_method, threshold = conformity_threshold)

  # combine outlier information
  outlier_summary_df <- nuclei_count_outliers_df |>
    dplyr::left_join(unexpected_neuron_outliers_df, by = "sample_id") |>
    dplyr::left_join(conformity_scores_df, by = "sample_id") |>
    dplyr::mutate(
      overall_outlier = nuclei_count_outlier | unexpected_neuron_outlier | conformity_outlier
    )

  return(outlier_summary_df)

}


#' Identify outliers in CTP data for multiple brain regions based on low nuclei counts, unexpected neuron fractions, and conformity scores.
#'
#' @param ctp_df Dataframe of CTP data for all samples and brain regions (one row per cell type).
#' @param cell_type_col Name of the column in ctp_df that contains the cell type label.
#' @param region_col Name of the column in ctp_df that contains the brain region.
#' @param unexpected_neuron_df Dataframe of unexpected neuron types for each brain region, with columns for region and cell type.
#' @param conformity_threshold Threshold for modified Z-score to identify conformity outliers (default 3.75).
#' @param conformity_method Correlation method to use for conformity scores (default "pearson")
#' @param neuron_fraction_threshold Threshold for fraction of unexpected neurons above which a sample is considered an outlier (default 0.2).
#' @param out_file Optional file path to save the outlier summary dataframe as a tab-delimited text file (default NULL, no file saved).
#'
#' @return Dataframe with sample_id, brain_region, total_nuclei, z_log10_total_nuclei, nuclei_count_outlier, total_neurons, unexpected_neuron_fraction, unexpected_neuron_outlier, conformity_score, modified_z_conformity_score, conformity_outlier, and overall_outlier (TRUE/FALSE) for each sample and brain region.
identify_ctp_outliers <- function(ctp_df, cell_type_col, region_col, unexpected_neuron_df, conformity_threshold=3.75, conformity_method="pearson", neuron_fraction_threshold=0.2, out_file=NULL) {

  regions <- unique(ctp_df[[region_col]])

  outlier_summary_list <- lapply(regions, function(region) {

    logger::log_info("Identifying outliers for region {region}")

    ctp_region_df <- ctp_df |>
      dplyr::filter(.data[[region_col]] == region)

    neuron_ctp_region_df <- ctp_region_df |>
      dplyr::filter(grepl("neuron", .data[[cell_type_col]], ignore.case = TRUE))

    unexpected_neuron_types_region <- unexpected_neuron_df |>
      dplyr::filter(.data[[region_col]] == region) |>
      dplyr::pull(cell_type_col)

    if (length(unexpected_neuron_types_region) == 0) {
      unexpected_neuron_types_region <- NULL
    }

    outlier_summary_df <- identify_ctp_outliers_one_region(
      ctp_region_df = ctp_region_df,
      cell_type_col = cell_type_col,
      conformity_threshold = conformity_threshold,
      conformity_method = conformity_method,
      neuron_ctp_df = neuron_ctp_region_df,
      unexpected_neuron_types = unexpected_neuron_types_region,
      neuron_fraction_threshold = neuron_fraction_threshold
    )

    outlier_summary_df[[region_col]] <- region

    return(outlier_summary_df)

  })

  outlier_summary_all_regions <- dplyr::bind_rows(outlier_summary_list)

  if(!is.null(out_file)) {
    logger::log_info("Saving outlier summary to {out_file}")
    write.table(outlier_summary_all_regions, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  return(outlier_summary_all_regions)

}

