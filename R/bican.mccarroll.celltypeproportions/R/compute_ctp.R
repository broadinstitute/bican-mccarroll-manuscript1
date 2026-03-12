# library(dplyr)
# library(rlang)
# library(logger)

# input_file <- "/broad/bican_um1_mccarroll/RNAseq/analysis/cellarium_upload/CAP_freeze_3/CAP_cell_metadata.annotated.with_sub_clusters.txt.gz"
# filters <- c("brain_region_abbreviation == 'DFC'", "single_cell_assay == '10X-GEMX-3P'")
# group_cols <- c("donor_external_id", "village")
# cell_type_col <- "annotation"
# metric_cols <- c("pct_intronic", "frac_contamination")
#
# cell_metadata <- read.table(
#   input_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE
# )
#
# a <- compute_ctp_and_metrics(cell_metadata, group_cols, cell_type_col, metric_cols, filters)
# b <- compute_ctp_and_metrics(cell_metadata, group_cols, cell_type_col, metric_cols, filters=c(filters, "!is.na(annotation)"))


#' Filters dataframe.
#'
#' @param df Dataframe to filter.
#' @param filters Character vector of filtering expressions.
#' @param group_cols Optional; Character vector of columns to group by before filtering
#' @param group_filters Optional; Character vector of filtering expressions to apply within groups.
#'
#' @return Filtered dataframe.
filter_df <- function(df, filters=NULL, group_cols=NULL, group_filters=NULL) {

  if (is.null(filters) || length(filters) == 0) {
    logger::log_info("No filtering criteria found.")
    return(df)
  }

  if(!is.character(filters)) {
    stop("`filters` must be a character vector")
  }

  logger::log_info("Applying the following filters: {paste(filters, collapse = ', ')}")

  if(!is.null(group_cols) & !is.null(group_filters)) {
    logger::log_info("Applying group-wise filters: {paste(group_filters, collapse = ', ')} within groups defined by: {paste(group_cols, collapse = ', ')}")
    filtered_df <- df |>
      dplyr::filter(!!!rlang::parse_exprs(filters)) |>
      dplyr::group_by(across(all_of(group_cols))) |>
      dplyr::filter(!!!rlang::parse_exprs(group_filters)) |>
      dplyr::ungroup()
  } else {
    filtered_df <- df |>
      dplyr::filter(!!!rlang::parse_exprs(filters))
  }

  return(filtered_df)

}


#' Computes cell type proportions for a given grouping.
#'
#' @param df Dataframe containing cell metadata.
#' @param group_cols Character vector of columns to group by (e.g., donor ID)
#' @param cell_type_col Column name representing cell type annotations.
#'
#' @return Dataframe with cell type proportions and total nuclei counts (raw & z-scored) per grouping.
compute_ctp <- function(df, group_cols, cell_type_col) {

  # Check if specified columns exist in the dataframe
  missing_cols <- setdiff(c(group_cols, cell_type_col), colnames(df))
  if (length(missing_cols) > 0) {
    stop(paste("The following columns are missing from the dataframe:", paste(missing_cols, collapse = ", ")))
  }

  # add sample label
  sample_df <- df |>
    tidyr::unite(sample_id, all_of(group_cols), sep = "_", remove = FALSE)

  # extract sample info
  sample_info <- sample_df |>
    dplyr::select(sample_id, all_of(group_cols)) |>
    dplyr::distinct()

  # Compute cell type proportions
  ctp_df <- sample_df |>
    tidyr::unite(sample_id, all_of(group_cols), sep = "_", remove = FALSE) |>
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
#' @param cell_type_col Column name representing cell type annotations.
#' @param metric_cols Character vector of metric columns to compute means for.
#'
#' @return Dataframe with mean metrics per cell type and grouping.
compute_mean_cell_type_metrics <- function(df, group_cols, cell_type_col, metric_cols) {

  # Check if specified columns exist in the data frame
  missing_cols <- setdiff(c(group_cols, metric_cols), colnames(df))
  if (length(missing_cols) > 0) {
    stop(paste("The following columns are missing from the dataframe:", paste(missing_cols, collapse = ", ")))
  }

  # Compute mean metrics
  mean_metrics_df <- df |>
    tidyr::unite(sample_id, all_of(group_cols), sep = "_", remove = FALSE) |>
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

#' Save CTP results.
#'
#' @param ctp_df Dataframe containing cell type proportions and optional metrics.
#' @param out_dir Output directory to save the results.
#' @param output_prefix Prefix for the output file name.
#'
#' @return None
save_ctp <- function(ctp_df, out_file) {
  # save file
  write.table(ctp_df, file = out_file, sep = "\t", row.names = FALSE, quote = FALSE)
}


#' Main workflow to compute cell type proportions and optional metrics.
#'
#' @param df Dataframe containing cell metadata.
#' @param group_cols Character vector of columns to group by (e.g., donor ID)
#' @param cell_type_col Column name representing cell type annotations.
#' @param metric_cols Optional; Character vector of metric columns to compute means for.
#' @param filters Optional; Character vector of filtering expressions.
#' @param group_filters Optional;  Filtering expressions to apply within groups.
#' @param out_file Optional; Output file to save the results
#'
#' @return Dataframe with cell type proportions and optional metrics.
compute_ctp_and_metrics <- function(df, group_cols, cell_type_col, metric_cols = NULL, filters = NULL, group_filters=NULL, out_file=NULL) {

  # Step 0: drop NA values
  #df_complete <- df[!is.na(df[[cell_type_col]]), ]

  # Step 1: Filter the data frame
  filtered_df <- filter_df(df, filters, group_cols, group_filters)

  # Step 2: Compute cell type proportions
  ctp_df <- compute_ctp(filtered_df, group_cols, cell_type_col)

  # Step 3: Compute mean metrics if specified
  if (!is.null(metric_cols) && length(metric_cols) > 0) {
    mean_metrics_df <- compute_mean_cell_type_metrics(filtered_df, group_cols, cell_type_col, metric_cols)

    ctp_df <- dplyr::full_join(
      ctp_df, mean_metrics_df,
      by=c("sample_id", cell_type_col)
    )

  }

  # Optional: save CTP if specified
  if (!is.null(out_file)) {
    save_ctp(ctp_df, out_file)
  }

  return(ctp_df)

}


#' Main workflow to compute cell type proportions and optional metrics from a metadata file.
#'
#' @param cell_metadata_file Path to the cell metadata file (one row for each cell barcode),
#' with columns for grouping (e.g., donor ID, brain region), cell type, and optional metrics.
#' @param group_cols Character vector of columns to group by (e.g., donor ID)
#' @param cell_type_col Column name representing cell type annotations.
#' @param metric_cols Optional; Character vector of metric columns to compute means for.
#' @param filters Optional; Character vector of filtering expressions.
#' @param group_filters Optional; Character vector of filtering expressions to apply within groups.
#' @param out_file Optional; Output file to save the results
#'
#' @return Dataframe with cell type proportions and optional metrics.
load_and_compute_ctp <- function(cell_metadata_file, group_cols, cell_type_col, out_file=NULL, metric_cols = NULL, filters = NULL, group_filters=NULL) {

  # read input file
  logger::log_info("Loading cell metadata from {cell_metadata_file}")
  df <- read.table(cell_metadata_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

  # compute CTP and metrics
  logger::log_info("Computing cell type proportions...")
  ctp_df <- compute_ctp_and_metrics(df, group_cols, cell_type_col, metric_cols, filters, group_filters)

  # save results
  if(!is.null(out_file)) {
    logger::log_info("Saving results to {out_file}")
    save_ctp(ctp_df, out_file)
  }

  return(ctp_df)
}

#' Recompute CTP from an existing CTP dataframe - either with a subset of the
#' existing grouping variables, or a new cell type label.
#'
#' This is useful for quickly recomputing CTP from an existing dataframe
#' rather than having to load the entire cell metadata file in again. The
#' tradeoff is that you cannot compute any single-cell metrics.
#'
#' @param ctp_df Dataframe of cell type proportions.
#' @param group_cols Character vector of columns to group by (e.g., donor ID)
#' @param cell_type_col Column name representing cell type annotations.
#' @param cell_type_label_map Dataframe with 2 columns: the first column should
#' match the existing cell type labels in the `cell_type_col` of the `ctp_df`,
#' and the second column should have the new cell type labels to group by.
#'
#' @return Dataframe with recomputed cell type proportions based on the new grouping.
regroup_ctp <- function(ctp_df, group_cols, cell_type_col, cell_type_label_map) {

  # check format of new annotation map
  if(ncol(annotation_map) != 2) {
    stop("`cell_type_label_map` must have exactly 2 columns: the first for the {cell_type_col} and the second for the new grouping variable.")
  }

  new_cell_type_col <- colnames(annotation_map)[2]

  # compute CTP using the counts from the original
  new_ctp_df <- ctp_df |>
    dplyr::left_join(annotation_map) |>
    dplyr::group_by(across(all_of(group_cols)), !!rlang::sym(new_cell_type_col)) |>
    dplyr::summarise(
      n_nuclei = sum(n_nuclei)
    ) |>
    dplyr::group_by(across(all_of(group_cols))) |>
    dplyr::mutate(
      total_nuclei = sum(n_nuclei),
      fraction_nuclei = n_nuclei / total_nuclei
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      sample_id = paste(!!!rlang::syms(group_cols), sep = "_")
    ) |>
    dplyr::select(sample_id, everything())

  return(new_ctp_df)

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

