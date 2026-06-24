# Compute cell type proportions from cell metadata file

## development, comment out in package build

# input_file <- "/broad/bican_um1_mccarroll/RNAseq/analysis/cellarium_upload/CAP_freeze_3/CAP_cell_metadata.annotated.with_sub_clusters.txt.gz"
# filters <- c("brain_region_abbreviation == 'DFC'", "single_cell_assay == '10X-GEMX-3P'")
# group_cols <- c("donor_external_id", "village")
# cell_type_cols <- c("annotation", "sub_cluster")
# metric_cols <- c("pct_intronic", "frac_contamination")
#
# cell_metadata <- read.table(
#   input_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE
# )
#
# a <- compute_ctp_and_metrics(cell_metadata, group_cols, cell_type_col, metric_cols, filters)
# b <- compute_ctp_and_metrics(cell_metadata, group_cols, cell_type_col, metric_cols, filters=c(filters, "!is.na(annotation)"))
#
# sample_ctp <- compute_ctp_and_metrics(
#   cell_metadata,
#   group_cols = c("donor_external_id", "brain_region_abbreviation_simple"),
#   cell_type_col = "annotation",
#   metric_cols = NULL,
#   filters = c("is.na(ctp_donor_outlier)", "is.na(ctp_sample_outlier)")
# )
#
# ratio_df <- generate_ctp_ratio_table(
#   sample_ctp,
#   numerator=c("MSN_D1_matrix", "MSN_D1_striosome"),
#   denominator=c("MSN_D2_matrix", "MSN_D2_striosome"),
#   z_score_threshold = 5,
#   ratio_name="D1_D2_MSN"
# )


#' Filters dataframe.
#'
#' @param df Dataframe to filter.
#' @param filters Character vector of filtering expressions.
#'
#' @return Filtered dataframe.
filter_df <- function(df, filters=NULL) {

  if (is.null(filters) || length(filters) == 0) {
    logger::log_info("No filtering criteria found.")
    return(df)
  }

  if(!is.character(filters)) {
    stop("`filters` must be a character vector")
  }

  logger::log_info("Applying the following filters: {paste(filters, collapse = ', ')}")

  filtered_df <- df |>
    dplyr::filter(!!!rlang::parse_exprs(filters))

  return(filtered_df)

}


#' Computes cell type proportions for a given grouping.
#'
#' @param df Dataframe containing cell metadata.
#' @param group_cols Character vector of columns to group by (e.g., donor ID)
#' @param cell_type_cols Column name representing cell type annotations.
#'  If multiple columns are provided, then values will be concatenated.
#'
#' @return Dataframe with cell type proportions and total nuclei counts (raw & z-scored) per grouping.
compute_ctp <- function(df, group_cols, cell_type_cols) {

  # Check if specified columns exist
  missing_cols <- setdiff(c(group_cols, cell_type_cols), colnames(df))
  if (length(missing_cols) > 0) {
    stop(paste("The following columns are missing from the dataframe:", paste(missing_cols, collapse = ", ")))
  }

  logger::log_info("Computing cell type proportions grouped by: [{paste(group_cols, collapse = ', ')}], using cell type columns: [{paste(cell_type_cols, collapse =', ')}]")

  cell_type_col = paste(cell_type_cols, collapse="__")

  # Create sample ID by concatenating grouping variables
  # Concatenate cell type columns if there are multiple
  sample_df <- df |>
    tidyr::unite(
      sample_id, all_of(group_cols), sep = "__", remove = FALSE
      ) |>
    tidyr::unite(
      !!rlang::sym(cell_type_col),
      all_of(cell_type_cols), sep="__", remove=FALSE
    )

  # Extract sample information for later joining
  sample_info <- sample_df |>
    dplyr::select(sample_id, all_of(group_cols)) |>
    dplyr::distinct()

  ctp_df <- sample_df |>
    tidyr::unite(sample_id, all_of(group_cols), sep = "__", remove = FALSE) |>
    dplyr::group_by(sample_id, !!rlang::sym(cell_type_col)) |>
    dplyr::summarise(n_nuclei = dplyr::n(), .groups = 'drop') |>
    dplyr::group_by(sample_id) |>
    dplyr::mutate(
      total_nuclei = sum(n_nuclei),
      fraction_nuclei = n_nuclei / total_nuclei
    ) |>
    dplyr::ungroup() |>
    tidyr::complete(sample_id, !!rlang::sym(cell_type_col), fill=list(n_nuclei = 0, fraction_nuclei = 0))

  nuclei_counts <- ctp_df |>
    dplyr::select(sample_id, total_nuclei) |>
    dplyr::distinct() |>
    na.omit() |>
    dplyr::mutate(log10_nuclei = log10(total_nuclei)) |>
    dplyr::mutate(z_log10_nuclei = as.numeric(scale(log10_nuclei))) |>
    dplyr::select(-log10_nuclei)

  final_df <- ctp_df |>
    dplyr::select(-total_nuclei) |>
    dplyr::left_join(sample_info, by ="sample_id") |>
    dplyr::left_join(nuclei_counts, by = "sample_id") |>
    dplyr::select(sample_id, all_of(group_cols), everything())

  return(final_df)
}


#' Computes mean metrics per cell type for a given grouping.
#'
#' @param df Dataframe containing cell metadata.
#' @param group_cols Character vector of columns to group by (e.g., donor ID)
#' @param cell_type_cols Column name representing cell type annotations.
#'  If multiple columns are provided, then values will be concatenated.
#' @param metric_cols Character vector of metric columns to compute means for.
#'
#' @return Dataframe with mean metrics per cell type and grouping.
compute_mean_cell_type_metrics <- function(df, group_cols, cell_type_cols, metric_cols) {

  # Check if specified columns exist in the data frame
  missing_cols <- setdiff(c(group_cols, cell_type_cols, metric_cols), colnames(df))
  if (length(missing_cols) > 0) {
    stop(paste("The following columns are missing from the dataframe:", paste(missing_cols, collapse = ", ")))
  }

  logger::log_info("Computing mean cell type metrics grouped by: [{paste(group_cols, collapse = ', ')}], using cell type columns: [{paste(cell_type_cols, collapse =', ')}]")

  cell_type_col <- paste(cell_type_cols, collapse="__")

  # Compute mean metrics
  mean_metrics_df <- df |>
    tidyr::unite(sample_id, all_of(group_cols), sep = "__", remove = FALSE) |>
    tidyr::unite(!!rlang::sym(cell_type_col), all_of(cell_type_cols), sep="__", remove=FALSE) |>
    dplyr::group_by(sample_id, !!rlang::sym(cell_type_col)) |>
    dplyr::summarise(
      across(
        all_of(metric_cols),
        ~ mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    ) |>
    dplyr::rename_with(
      ~ paste0("mean_", .x),
      all_of(metric_cols)
    )

  return(mean_metrics_df)
}


#' Main workflow to compute cell type proportions and optional metrics.
#'
#' @param df Dataframe containing cell metadata.
#' @param group_cols Character vector of columns to group by (e.g., donor ID)
#' @param cell_type_cols Column name representing cell type annotations.
#'  If multiple columns are provided, then values will be concatenated.
#' @param metric_cols Optional; Character vector of metric columns to compute means for.
#' @param filters Optional; Character vector of filtering expressions.
#' @param out_file Optional; Output file to save the results
#'
#' @return Dataframe with cell type proportions and optional metrics.
compute_ctp_and_metrics <- function(df, group_cols, cell_type_cols, metric_cols = NULL, filters = NULL, out_file=NULL) {

  # Apply filters
  filtered_df <- filter_df(df, filters)

  # Compute cell type proportions
  ctp_df <- compute_ctp(filtered_df, group_cols, cell_type_cols)

  # Compute mean metrics if specified
  if (!is.null(metric_cols) && length(metric_cols) > 0) {
    mean_metrics_df <- compute_mean_cell_type_metrics(filtered_df, group_cols, cell_type_cols, metric_cols)
    ctp_df <- dplyr::full_join(
      ctp_df, mean_metrics_df,
      by=c("sample_id", paste(cell_type_cols, collapse="__"))
    )
  }

  # Save CTP if specified
  if (!is.null(out_file)) {
    logger::log_info("Saving cell type proportions and metrics to {out_file}")
    write.table(
      ctp_df, file = out_file, sep = "\t", row.names = FALSE, quote = FALSE
    )
  }

  return(ctp_df)

}


#' Main workflow to compute cell type proportions and optional metrics from a metadata file.
#'
#' @param cell_metadata_file Path to the cell metadata file (one row for each cell barcode),
#' with columns for grouping (e.g., donor ID, brain region), cell type, and optional metrics.
#' @param group_cols Character vector of columns to group by (e.g., donor ID)
#' @param cell_type_cols Column name representing cell type annotations.
#'  If multiple columns are provided, then values will be concatenated.
#' @param metric_cols Optional; Character vector of metric columns to compute means for.
#' @param filters Optional; Character vector of filtering expressions.
#' @param out_file Optional; Output file to save the results
#'
#' @return Dataframe with cell type proportions and optional metrics.
load_and_compute_ctp_and_metrics <- function(cell_metadata_file, group_cols, cell_type_cols, out_file=NULL, metric_cols = NULL, filters = NULL) {

  # Read input file
  logger::log_info("Loading cell metadata from {cell_metadata_file}")
  df <- read.table(cell_metadata_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

  # Compute CTP and metrics (if specified)
  ctp_df <- compute_ctp_and_metrics(df, group_cols, cell_type_cols, metric_cols, filters, out_file)

  return(ctp_df)
}


#' Generate Cell-Type Proportion (CTP) Ratio Table
#'
#' This function calculates ratios between specific cell-type populations across
#' defined brain regions. It allows for the comparison of a subset of cells
#' (numerator) against either another subset or the total nuclei count (denominator).
#'
#' @param ctp_df A data frame containing cell type counts.
#'  Typically generated by `compute_ctp`
#'  Must include columns: `sample_id`, `total_nuclei`, `n_nuclei`
#' @param numerator Character vector of cell types to include in the numerator of the ratio.
#' @param denominator Character vector of cell types to include in the denominator of the ratio.
#'  If NULL, the denominator will be the total nuclei count for the sample.
#' @param cell_type_col Column name in `ctp_df` that contains cell type annotations
#' @param group_cols Character vector of columns to group by (e.g., donor ID, brain region).
#' @param filters Optional; Character vector of filtering expressions to apply to `ctp_df` before calculating ratios.
#' @param ratio_name Name to assign to the calculated ratio
#' @param brain_region_col Optional; Column name in `ctp_df` that contains brain region annotations. 
#'   If provided, z-scores of ratios will be computed per brain region.
#' @param z_score_threshold Optional; Numeric threshold for flagging outliers based on z-scores of ratios. 
#'   If provided, a new column `outlier` will be added to the output indicating whether each sample is an outlier.
#' 
#' @return A data frame with calculated ratios and optional z-scores and outlier flags, grouped by specified columns.
generate_ctp_ratio_table <- function(
  ctp_df, 
  numerator, 
  denominator,
  cell_type_col,
  group_cols,
  filters=NULL,
  ratio_name,
  brain_region_col=NULL,
  z_score_threshold=NULL,
  out_file=NULL
) {

  # check that columns are present in dataframe 
  missing_cols <- setdiff(c(group_cols, cell_type_col, "n_nuclei", "total_nuclei"), colnames(ctp_df))
  if (length(missing_cols) > 0) {
    stop(paste("The following required columns are missing from `ctp_df`:", paste(missing_cols, collapse = ", ")))
  }

  ctp_df_filtered = filter_df(ctp_df, filters)

  ratio_df <- ctp_df_filtered |>

    # group by sample id
    dplyr::group_by(
      dplyr::across(dplyr::all_of(c("sample_id", group_cols, "total_nuclei")))
      # dplyr::all_of(c("sample_id", group_cols, "total_nuclei")) 
    ) |>
    
    # calculate sums for numerator and denominator subsets
    dplyr::summarise(
      # sum nuclei for cell types defined in 'numerator'
      numerator_sum = sum(
        n_nuclei[.data[[cell_type_col]] %in% numerator],
        na.rm = TRUE
      ),

      # if denominator is NULL, use all nuclei in the group; otherwise, sum specific subset
      denominator_sum = if (is.null(denominator)) {
        sum(n_nuclei, na.rm = TRUE)
      } else {
        sum(
          .data$n_nuclei[.data[[cell_type_col]] %in% denominator],
          na.rm = TRUE
        )
      },
      .groups = "drop"
    ) |>
    
    # compute the final ratio
    dplyr::mutate(ratio = numerator_sum / denominator_sum) |>
    
    # add the ratio name as a column
    dplyr::mutate(ratio_name = ratio_name) |>
    
    # return columns needed for plotting
    dplyr::select(dplyr::all_of(c("sample_id", group_cols, "total_nuclei", "ratio", "ratio_name")))
  
  # compute z-scores per brain region
  if (!is.null(brain_region_col)) {
    logger::log_info("Brain region column provided: {brain_region_col}. Computing z-scores of ratios per brain region.")
    ratio_df <- ratio_df |>
      dplyr::group_by(dplyr::across(all_of(brain_region_col))) |>
      dplyr::mutate(
        z_ratio = as.numeric(scale(ratio))
      ) |>
      dplyr::ungroup()
  } else {
    logger::log_info("No brain region column provided, skipping z-score calculation.")
  }

  if (!is.null(z_score_threshold)) {
    logger::log_info("Z-score threshold provided: {z_score_threshold}. Flagging outliers based on z-scores.")
    ratio_df <- ratio_df |>
      dplyr::mutate(outlier=abs(z_ratio) > z_score_threshold)
  } else {
    logger::log_info("No z-score threshold provided, skipping outlier flagging.")
  }

  if(!is.null(out_file)) {
    logger::log_info("Saving CTP ratio table to {out_file}")
    write.table(ratio_df, file = out_file, sep = "\t", row.names = FALSE, quote = FALSE)
  }
  
  return(ratio_df)
}

#' Main workflow to load CTP data and generate ratio table.
#'
#' @param ctp_file Path to the CTP data file (output from `compute_ctp`), containing columns for sample ID, total nuclei, cell type counts, and grouping variables.
#' @param numerator Character vector of cell types to include in the numerator of the ratio.
#' @param denominator Character vector of cell types to include in the denominator of the ratio.
#'  If NULL, the denominator will be the total nuclei count for the sample.
#' @param cell_type_col Column name in `ctp_file` that contains cell type annotations
#' @param group_cols Character vector of columns to group by (e.g., donor ID, brain region).
#' @param filters Optional; Character vector of filtering expressions to apply to the loaded CTP data before calculating ratios.
#' @param ratio_name Name to assign to the calculated ratio
#' @param brain_region_col Optional; Column name in `ctp_file` that contains brain region annotations. 
#'   If provided, z-scores of ratios will be computed per brain region.
#' @param z_score_threshold Optional; Numeric threshold for flagging outliers based on z-scores of ratios. 
#'   If provided, a new column `outlier` will be added to the output indicating whether each sample is an outlier.
#' @param out_file Optional; Output file to save the generated ratio table.
#' 
#' @return A data frame with calculated ratios and optional z-scores and outlier flags, grouped by specified columns.
load_and_generate_ctp_ratio_table <- function(
  ctp_file,
  numerator,
  denominator,
  cell_type_col,
  group_cols,
  filters=NULL,
  ratio_name,
  brain_region_col=NULL,
  z_score_threshold=NULL,
  out_file=NULL
) {
  
  logger::log_info("Loading CTP data from {ctp_file}")
  ctp_df <- read.table(ctp_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  ratio_df <- generate_ctp_ratio_table(
    ctp_df = ctp_df,
    numerator = numerator,
    denominator = denominator,
    cell_type_col = cell_type_col,
    group_cols = group_cols,
    filters = filters,
    ratio_name = ratio_name,
    brain_region_col = brain_region_col,
    z_score_threshold = z_score_threshold,
    out_file = out_file
  )
  
  return(ratio_df)
}


#' Extract sample-level metadata from cell-level metadata, applying optional filters and grouping.
#'
#' This function takes cell-level metadata and extracts sample-level metadata by
#' grouping based on specified columns (e.g., donor ID, brain region) and applying optional filters.
#' It can also aggregate certain metadata columns by concatenating unique values within each group.
#'
#' @param df Dataframe containing cell-level metadata, with columns for grouping (e.g., donor ID, brain region) and other metadata.
#' @param group_cols Character vector of columns to group by (e.g., donor ID, brain region).
#' @param donor_metadata_cols Character vector of donor-level metadata columns to include in the output.
#' @param metadata_cols_to_group Character vector of metadata columns to aggregate by concatenating unique values within each group.
#' @param filters Optional; Character vector of filtering expressions to apply to the cell-level metadata before extracting sample-level metadata.
#' @param out_file Optional; Output file to save the extracted sample metadata.
#'
#' @return Dataframe with one row per sample (defined by `group_cols`), including the specified donor-level metadata and aggregated metadata columns.
extract_sample_metadata <- function(df, group_cols, donor_metadata_cols, metadata_cols_to_group, filters=NULL, out_file=NULL) {

  filtered_df <- filter_df(df, filters)

  sample_metadata <- filtered_df |>
    dplyr::mutate(
      sample_id = paste(!!!rlang::syms(group_cols), sep = "_")
    ) |>
    dplyr::select(
      sample_id, all_of(group_cols), all_of(donor_metadata_cols), all_of(metadata_cols_to_group)
    ) |>
    dplyr::distinct() |>
    dplyr::group_by(sample_id) |>
    dplyr::mutate(
      dplyr::across(
        all_of(metadata_cols_to_group),
        ~ paste(sort(unique(.x)), collapse=":")
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::distinct()

  if (!is.null(out_file)) {
    logger::log_info("Saving sample metadata to {out_file}")
    write.table(sample_metadata, file = out_file, sep = "\t", row.names = FALSE, quote = FALSE)
  }

  return(sample_metadata)
}


#' Main workflow to extract sample-level metadata from a cell metadata file.
#'
#' This function reads cell-level metadata from a file, applies optional filters,
#' and extracts sample-level metadata by grouping based on specified columns.
#'  It can also aggregate certain metadata columns by concatenating unique values within each group.
#'   The resulting sample metadata can be saved to an output file if specified.
#'
#' @param cell_metadata_file Path to the cell metadata file (one row for each cell barcode),
#' with columns for grouping (e.g., donor ID, brain region) and other metadata.
#' @param group_cols Character vector of columns to group by (e.g., donor ID, brain region).
#' @param donor_metadata_cols Character vector of donor-level metadata columns to include in the output.
#' @param metadata_cols_to_group Character vector of metadata columns to aggregate by concatenating unique values within each group.
#' @param filters Optional; Character vector of filtering expressions to apply to the cell-level metadata before extracting sample-level metadata.
#' @param out_file Optional; Output file
#'
#' @return Dataframe with one row per sample (defined by `group_cols`), including the specified donor-level metadata and aggregated metadata columns.
load_and_extract_sample_metadata <- function(cell_metadata_file, group_cols, donor_metadata_cols, metadata_cols_to_group, filters=NULL, out_file=NULL) {

  # read input file
  logger::log_info("Loading cell metadata from {cell_metadata_file}")
  df <- read.table(cell_metadata_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

  # extract sample metadata
  sample_metadata <- extract_sample_metadata(df, group_cols, donor_metadata_cols, metadata_cols_to_group, filters, out_file)

  return(sample_metadata)
}
