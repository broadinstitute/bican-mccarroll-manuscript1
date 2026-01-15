#library(dplyr)
#library(rlang)
#library(logger)

#input_file <- "/broad/bican_um1_mccarroll/RNAseq/analysis/cellarium_upload/CAP_freeze_3/CAP_cell_metadata.annotated.txt.gz"
#df <- read.table(input_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#filters <- c("brain_region_abbreviation == 'DFC'", "single_cell_assay == '10X-GEMX-3P'")
#group_cols <- c("donor_external_id", "village")
#cell_type_col <- "annotation_most_specific"
#metric_cols <- c("pct_intronic", "frac_contamination")
#out_dir <- "/broad/mccarroll/yooolivi/test/celltypeproportions"
#prefix <- "DFC_10X-GEMX-3P"

#ctp <- compute_ctp_and_metrics(df, group_cols, cell_type_col, metric_cols, filters)
#save_ctp(ctp, out_dir, prefix)


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
#' @param cell_type_col Column name representing cell type annotations.
#'
#' @return Dataframe with cell type proportions and total nuclei counts (raw & z-scored) per grouping.
compute_ctp <- function(df, group_cols, cell_type_col) {

  # Check if specified columns exist in the dataframe
  missing_cols <- setdiff(c(group_cols, cell_type_col), colnames(df))
  if (length(missing_cols) > 0) {
    stop(paste("The following columns are missing from the dataframe:", paste(missing_cols, collapse = ", ")))
  }

  # Compute cell type proportions
  ctp_df <- df |>
    dplyr::group_by(across(all_of(group_cols)), !!rlang::sym(cell_type_col)) |>
    dplyr::summarise(n_nuclei = dplyr::n(), .groups = 'drop') |>
    dplyr::group_by(across(all_of(group_cols))) |>
    dplyr::mutate(total_nuclei = sum(n_nuclei),
           fraction_nuclei = n_nuclei / total_nuclei)

  nuclei_counts <- ctp_df |>
    dplyr::select(all_of(group_cols), total_nuclei) |>
    dplyr::distinct() |>
    dplyr::summarise(z_log10_nuclei = scale(log10(total_nuclei)))

  sample_df <- ctp_df |>
    dplyr::left_join(nuclei_counts, by = group_cols)

  return(sample_df)
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
    dplyr::group_by(across(all_of(group_cols)), !!rlang::sym(cell_type_col)) |>
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
#' @param cell_type_col Column name representing cell type annotations.
#' @param metric_cols (Optional) Character vector of metric columns to compute means for.
#' @param filters (Optional) Character vector of filtering expressions.
#'
#' @return Dataframe with cell type proportions and optional metrics.
compute_ctp_and_metrics <- function(df, group_cols, cell_type_col, metric_cols = NULL, filters = NULL) {

  # Step 1: Filter the data frame
  filtered_df <- filter_df(df, filters)

  # Step 2: Compute cell type proportions
  ctp_df <- compute_ctp(filtered_df, group_cols, cell_type_col)

  # Step 3: Compute mean metrics if specified
  if (!is.null(metric_cols) && length(metric_cols) > 0) {
    mean_metrics_df <- compute_mean_cell_type_metrics(filtered_df, group_cols, cell_type_col, metric_cols)

    combined_df <- dplyr::full_join(
      ctp_df, mean_metrics_df,
      by=c(group_cols, cell_type_col)
    )

    return(combined_df)
  } else {
    return(ctp_df)
  }
}


#' Save CTP results.
#'
#' @param ctp_df Dataframe containing cell type proportions and optional metrics.
#' @param out_dir Output directory to save the results.
#' @param output_prefix Prefix for the output file name.
#'
#' @return None
save_ctp <- function(ctp_df, out_dir, output_prefix) {

  # round floats
  ctp_df_out <- ctp_df |>
    dplyr::mutate(across(where(is.double), ~ round(.x, 6)))

  # save file
  ctp_file <- file.path(out_dir, paste0(output_prefix, ".cell_type_proportions.txt"))
  write.table(ctp_df_out, file = ctp_file, sep = "\t", row.names = FALSE, quote = FALSE)

}








