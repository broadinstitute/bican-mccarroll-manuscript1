
# sample_ctp <- read.table(
#   "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/cell_type_proportions/LEVEL_2/donor_region.annotation.cell_type_proportions.txt",
#   sep="\t", header=TRUE, stringsAsFactors = FALSE
# )
#
# ratio_df <- generate_ctp_ratio_table(
#   sample_ctp,
#   numerator=c("MSN_D1_matrix", "MSN_D1_striosome"),
#   denominator=c("MSN_D2_matrix", "MSN_D2_striosome"),
#   z_score_threshold = 5,
#   ratio_name="D1_D2_MSN"
# )


#' Generate Cell-Type Proportion (CTP) Ratio Table
#'
#' This function calculates ratios between specific cell-type populations across
#' defined brain regions. It allows for the comparison of a subset of cells
#' (numerator) against either another subset or the total nuclei count (denominator).
#'
#' @param ctp_df A data frame containing cell type counts. Must include columns:
#'   `brain_region_abbreviation_simple`, `sample_id`, `donor_external_id`,
#'   `total_nuclei`, `n_nuclei`, and the column specified in `cell_type_col`.
#' @param numerator A character vector of cell type names to be summed for the numerator.
#' @param denominator A character vector of cell type names to be summed for the denominator.
#'   If `NULL`, the total nuclei count for the sample/region is used as the denominator.
#' @param cell_type_col Character string specifying the column name in `ctp_df`
#'   that contains cell type annotations. Defaults to `"annotation"`.
#' @param region_list A character vector of brain regions to include in the analysis.
#'   Defaults to `c("CaH", "Pu", "NAC")`.
#' @param z_score_threshold Threshold for calling outliers based on the z-score of
#' the ratio for a given region.
#'
#' @return A tibble (data frame) containing `sample_id`, `brain_region_abbreviation_simple`,
#'   `donor_external_id`, `total_nuclei`, and the calculated `ratio`.
#'
#' @importFrom dplyr filter group_by summarise mutate select .data
#' @export
generate_ctp_ratio_table <- function(
    ctp_df,
    numerator,
    denominator,
    cell_type_col = "annotation",
    region_list = c("CaH", "Pu", "NAC"),
    brain_region_col="brain_region_abbreviation_simple",
    donor_col="donor_external_id",
    ratio_name = NA,
    z_score_threshold=NULL)
{

  ratio_df <- ctp_df |>

    # Filter to brain regions of interest
    dplyr::filter(.data[[region_col]]%in% region_list) |>

    # Group by sample_id
    dplyr::group_by(
      sample_id,
      brain_region_abbreviation_simple,
      donor_external_id,
      total_nuclei
    ) |>

    # Calculate sums for numerator and denominator subsets
    dplyr::summarise(
      # Sum nuclei for cell types defined in 'numerator'
      numerator_sum = sum(
        n_nuclei[.data[[cell_type_col]] %in% numerator],
        na.rm = TRUE
      ),

      # If denominator is NULL, use all nuclei in the group; otherwise, sum specific subset
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

    # Compute the final ratio
    dplyr::mutate(ratio = numerator_sum / denominator_sum) |>

    # add the ratio name
    dplyr::mutate(ratio_name = ratio_name) |>

    # Return columns needed for plotting
    dplyr::select(
      sample_id,
      all_of(c(region_col, donor_col)),
      total_nuclei,
      ratio,
      ratio_name
    )

  # compute z-scores of ratios per brain region
  ratio_df <- ratio_df |>
    dplyr::group_by(brain_region_abbreviation_simple) |>
    dplyr::mutate(
      z_ratio = as.numeric(scale(ratio))
    ) |>
    dplyr::ungroup()

  if (!is.null(z_score_threshold)) {
    ratio_df <- ratio_df |>
      dplyr::mutate(outlier=abs(z_ratio) > z_score_threshold)
  } else {
    ratio_df <- ratio_df |>
      dplyr::mutate(outlier=NA)
  }

  return(ratio_df)
}


