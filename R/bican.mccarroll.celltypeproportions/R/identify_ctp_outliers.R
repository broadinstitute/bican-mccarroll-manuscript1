# ctp_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/yooolivi/cell_type_proportions/LEVEL_0"

# sample_ctp <- read.table(
#   file.path(ctp_dir, "donor_region_village.annotation_most_specific.cell_type_proportions.txt"),
#   sep="\t", header=TRUE, stringsAsFactors = FALSE
# )

# sample_neuron_ctp <- read.table(
#  file.path(ctp_dir, "donor_region_village.annotation_most_specific.neurons_only.cell_type_proportions.txt"),
#  sep="\t", header=TRUE, stringsAsFactors = FALSE
# )

# unexpected_neuron_df <- read.table(
#  "/broad/mccarroll/yooolivi/projects/bican/manuscript_1/ctp_outliers/unexpected_neurons.annotation_most_specific.txt",
#  header = TRUE,
#  sep = "\t",
#  stringsAsFactors = FALSE
# )

# ctp_outliers <- identify_ctp_outliers(
#  ctp_df = sample_ctp,
#  cell_type_col = "annotation_most_specific",
#  region_col = "brain_region_abbreviation_simple",
#  neuron_ctp_df = sample_neuron_ctp,
#  unexpected_neuron_df = unexpected_neuron_df,
#  conformity_threshold = 3.75,
#  conformity_method = "pearson",
#  neuron_fraction_threshold = 0.2,
#   drop_low_nuclei_outliers = FALSE,
#   drop_unexpected_neuron_outliers = FALSE
# )

utils::globalVariables(c("z_log10_total_nuclei",
                         "unexpected_neuron", "unexpected_neurons", "total_neurons",
                         "unexpected_neuron_fraction",
                         "fraction_nuclei", "fraction_nuclei_transformed",
                         "conformity_score", "modified_z_conformity_score",
                         "nuclei_count_outlier", "unexpected_neuron_outlier", "conformity_outlier", "median", "mad"))


#' Identify samples to drop based on low nuclei counts.
#'
#' @param ctp_df Dataframe of cell type counts and proportions (one row per cell type).
#' @param z_threshold Z-score threshold below which a sample is considered an outlier 
#'    Default -1.96, corresponding to the lower 2.5\% of a normal distribution.
#'
#' @return Dataframe with sample_id, total_nuclei, z_log10_total_nuclei, and nuclei_count_outlier (TRUE/FALSE).
identify_low_nuclei_count_outliers <- function(ctp_df, z_threshold = -1.96) {

  # if z_log10_total_nuclei is already in the data, just use that to identify outliers
  if ("z_log10_total_nuclei" %in% colnames(ctp_df)) {
    sample_counts <- ctp_df |>
      dplyr::select(sample_id, total_nuclei, z_log10_total_nuclei) |>
      dplyr::distinct()
  }  else {
    # otherwise, compute z_log10_total_nuclei and identify outliers
    sample_counts <- ctp_df |>
      dplyr::select(sample_id, total_nuclei) |>
      dplyr::distinct() |>
      dplyr::mutate(z_log10_total_nuclei = as.numeric(scale(log10(total_nuclei))))
  } 
  
  sample_counts <- sample_counts |>
    dplyr::mutate(nuclei_count_outlier = z_log10_total_nuclei < z_threshold)

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

  neuron_ctp_df |>
    dplyr::mutate(unexpected_neuron = !!rlang::sym(cell_type_col) %in% unexpected_neuron_types) |>
    dplyr::group_by(sample_id) |>
    dplyr::summarise(
      total_neurons = sum(n_nuclei),
      unexpected_neurons = sum(n_nuclei[unexpected_neuron]),
      unexpected_neuron_fraction = unexpected_neurons / total_neurons,
      .groups = "drop"
    ) |>
    dplyr::mutate(unexpected_neuron_outlier = unexpected_neuron_fraction > threshold)

}


#' Apply arcsin square root transformation to CTP fractions and pivot to wide format.
#'
#' @param ctp_df Dataframe of cell type counts and proportions (one row per cell type).
#' @param cell_type_col Name of the column in ctp_df that contains the cell type.
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
#' @param threshold Threshold for absolute value of modified Z-score to identify conformity outliers (default 3.75).
#'
#' @return Dataframe with sample_id, conformity_score, modified_z_conformity_score, and conformity_outlier (TRUE/FALSE).
compute_ctp_conformity_scores <- function(transformed_ctp_df, method = "pearson", threshold=3.75) {

  samples <- rownames(transformed_ctp_df)
  cor_mat <- stats::cor(t(transformed_ctp_df), method = method)
  diag(cor_mat) <- NA

  sample_conformity_df <- data.frame(
    sample_id = rownames(cor_mat),
    conformity_score = apply(cor_mat, 1, stats::median, na.rm = TRUE)
  ) |>
    dplyr::mutate(
      median = stats::median(conformity_score, na.rm = TRUE),
      mad = stats::mad(conformity_score, na.rm = TRUE),
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
#' @param drop_low_nuclei_outliers Whether to drop low nuclei count outliers before computing conformity scores (default TRUE).
#' @param drop_unexpected_neuron_outliers Whether to drop unexpected neuron outliers before computing conformity scores (default TRUE).
#' 
#' @return Dataframe with sample_id, total_nuclei, z_log10_total_nuclei, nuclei_count_outlier, total_neurons, unexpected_neuron_fraction, unexpected_neuron_outlier, conformity_score, modified_z_conformity_score, conformity_outlier, and overall_outlier (TRUE/FALSE).
identify_ctp_outliers_one_region <- function(
  ctp_region_df, cell_type_col, conformity_threshold=3.75, conformity_method="pearson", 
  neuron_ctp_df=NULL, unexpected_neuron_types=NULL, neuron_fraction_threshold=0.2,
  drop_low_nuclei_outliers = TRUE, drop_unexpected_neuron_outliers = TRUE
) {

  # identify low nuclei count outliers
  nuclei_count_outliers_df <- identify_low_nuclei_count_outliers(ctp_region_df)
  low_nuclei_outliers <- nuclei_count_outliers_df |> dplyr::filter(nuclei_count_outlier) |> dplyr::pull(sample_id)

  # identify unexpected neuron outliers (if neuron CTP provided)
  if (!is.null(neuron_ctp_df) && !is.null(unexpected_neuron_types)) {
    # logger::log_info("Unexpected neurons provided, identifying outliers based on unexpected neuron fraction with threshold {neuron_fraction_threshold}")
    unexpected_neuron_outliers_df <- identify_unexpected_neuron_outliers(neuron_ctp_df, unexpected_neuron_types, cell_type_col, neuron_fraction_threshold)
    unexpected_neuron_outliers <- unexpected_neuron_outliers_df |>
      dplyr::filter(unexpected_neuron_outlier) |>
      dplyr::pull(sample_id)

  } else {
    # logger::log_info("No unexpected neurons provided, no outliers will be identified based on unexpected neuron fraction")
    unexpected_neuron_outliers_df <- data.frame(
      sample_id = unique(ctp_region_df$sample_id),
      total_neurons = NA,
      unexpected_neuron_fraction = NA,
      unexpected_neuron_outlier = FALSE
    )
  }

  # filter out low nuclei and unexpected neuron outliers before computing conformity scores

  if (drop_low_nuclei_outliers & drop_unexpected_neuron_outliers) {
    logger::log_info("Dropping low nuclei count and unexpected neuron outliers before computing conformity scores for region")
    ctp_region_df_filtered <- ctp_region_df |>
      dplyr::filter(
        !(sample_id %in% low_nuclei_outliers) &
        !(sample_id %in% unexpected_neuron_outliers)
      )
  } else if (drop_low_nuclei_outliers) {
    logger::log_info("Dropping low nuclei count outliers but not unexpected neuron outliers before computing conformity scores for region")
    ctp_region_df_filtered <- ctp_region_df |>
      dplyr::filter(
        !(sample_id %in% low_nuclei_outliers)
      )
  } else if (drop_unexpected_neuron_outliers) {
    logger::log_info("Dropping unexpected neuron outliers but not low nuclei count outliers before computing conformity scores for region")
    ctp_region_df_filtered <- ctp_region_df |>
      dplyr::filter(
        !(sample_id %in% unexpected_neuron_outliers)
      )
  } else {
    logger::log_info("Not dropping any samples before computing conformity scores for region")
    ctp_region_df_filtered <- ctp_region_df
  }

  # transform CTP and compute conformity scores (dropping low nuclei and neuron outliers first)
  transformed_ctp_df <- arcsin_sqrt_transform_ctp(ctp_region_df_filtered, cell_type_col)
  conformity_scores_df <- compute_ctp_conformity_scores(
    transformed_ctp_df, 
    method = conformity_method, 
    threshold = conformity_threshold
  )

  # combine outlier information
  if (!is.null(neuron_ctp_df) && !is.null(unexpected_neuron_types)) {
    # if both low nuclei and unexpected neuron outliers are identified, combine all outlier information
    outlier_summary_df <- nuclei_count_outliers_df |>
      dplyr::left_join(unexpected_neuron_outliers_df, by = "sample_id") |>
      dplyr::left_join(conformity_scores_df, by = "sample_id") |>
      dplyr::mutate(
        overall_outlier = nuclei_count_outlier | unexpected_neuron_outlier | conformity_outlier
      )
  } else {
    # if only low nuclei outliers are identified, combine with conformity outlier information
    outlier_summary_df <- nuclei_count_outliers_df |>
      dplyr::left_join(conformity_scores_df, by = "sample_id") |>
      dplyr::mutate(
        overall_outlier = nuclei_count_outlier | conformity_outlier
      )
  }

  return(outlier_summary_df)

}


#' Identify outliers in CTP data for multiple brain regions based on low nuclei counts, unexpected neuron fractions, and conformity scores.
#'
#' @param ctp_df Dataframe of CTP data for all samples and brain regions (one row per cell type).
#' @param cell_type_col Name of the column in ctp_df that contains the cell type label.
#' @param region_col Name of the column in ctp_df that contains the brain region.
#' @param conformity_threshold Threshold for modified Z-score to identify conformity outliers (default 3.75).
#' @param conformity_method Correlation method to use for conformity scores (default "pearson")
#' @param neuron_ctp_df Optional dataframe of neuron-only CTP data for all samples and brain regions (one row per neuron cell type).
#' @param unexpected_neuron_df Optional dataframe of unexpected neuron types for each brain region, with columns for region and cell type.
#' @param neuron_fraction_threshold Optional threshold for fraction of unexpected neurons above which a sample is considered an outlier (default 0.2).
#' @param drop_low_nuclei_outliers Whether to drop low nuclei count outliers before computing conformity scores (default TRUE).
#' @param drop_unexpected_neuron_outliers Whether to drop unexpected neuron outliers before computing conformity scores (default TRUE).
#' @param out_file Optional file path to save the outlier summary dataframe as a tab-delimited text file (default NULL, no file saved).
#'
#' @return Dataframe with sample_id, brain_region, total_nuclei, z_log10_total_nuclei, nuclei_count_outlier, total_neurons, unexpected_neuron_fraction, 
#' unexpected_neuron_outlier, conformity_score, modified_z_conformity_score, conformity_outlier, and overall_outlier (TRUE/FALSE) for each sample and brain region.
identify_ctp_outliers <- function(
  ctp_df, cell_type_col, region_col, 
  conformity_threshold=3.75, conformity_method="pearson", 
  neuron_ctp_df=NULL, unexpected_neuron_df=NULL, neuron_fraction_threshold=NULL, 
  drop_low_nuclei_outliers = TRUE, drop_unexpected_neuron_outliers = TRUE,
  out_file=NULL
) {

  # run outlier detection separately for each brain region, then combine results
  regions <- unique(ctp_df[[region_col]])

  outlier_summary_list <- lapply(regions, function(region) {

    logger::log_info("Identifying outliers for region {region}")

    ctp_region_df <- ctp_df |>
      dplyr::filter(.data[[region_col]] == region)

    # optional identification of outliers based on unexpected neuron fraction
    if (!is.null(neuron_ctp_df) & !is.null(unexpected_neuron_df) & !is.null(neuron_fraction_threshold)) {
      logger::log_info("Neuron CTP dataframe, unexpected neuron dataframe, and neuron fraction threshold provided, identifying unexpected neuron outliers for region {region}")
      
      neuron_ctp_region_df <- neuron_ctp_df |>
        dplyr::filter(.data[[region_col]] == region)

      unexpected_neuron_types_region <- unexpected_neuron_df |>
        dplyr::filter(.data[[region_col]] == region) |>
        dplyr::pull(cell_type_col)

      outlier_summary_df <- identify_ctp_outliers_one_region(
        ctp_region_df = ctp_region_df,
        cell_type_col = cell_type_col,
        conformity_threshold = conformity_threshold,
        conformity_method = conformity_method,
        neuron_ctp_df = neuron_ctp_region_df,
        unexpected_neuron_types = unexpected_neuron_types_region,
        neuron_fraction_threshold = neuron_fraction_threshold,
        drop_low_nuclei_outliers = drop_low_nuclei_outliers,
        drop_unexpected_neuron_outliers = drop_unexpected_neuron_outliers
      )

    } else {
      logger::log_info("No neuronal proportion information provided, no unexpected neuron outliers will be identified for region {region}")

      outlier_summary_df <- identify_ctp_outliers_one_region(
        ctp_region_df = ctp_region_df,
        cell_type_col = cell_type_col,
        conformity_threshold = conformity_threshold,
        conformity_method = conformity_method,
        drop_low_nuclei_outliers = drop_low_nuclei_outliers,
        drop_unexpected_neuron_outliers = drop_unexpected_neuron_outliers
      )
    }

    outlier_summary_df[[region_col]] <- region

    return(
      outlier_summary_df |>
        dplyr::select(sample_id, dplyr::all_of(region_col), dplyr::everything())
      )

  })

  outlier_summary_all_regions <- dplyr::bind_rows(outlier_summary_list)

  if(!is.null(out_file)) {
    logger::log_info("Saving outlier summary to {out_file}")
    utils::write.table(outlier_summary_all_regions, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  return(outlier_summary_all_regions)

}

