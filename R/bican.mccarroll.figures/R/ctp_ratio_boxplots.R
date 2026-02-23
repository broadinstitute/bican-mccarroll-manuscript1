

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
    region_list = c("CaH", "Pu", "NAC")) 
{
  
  ratio_df <- ctp_df |>
    # Filter to brain regions of interest
    dplyr::filter(.data$brain_region_abbreviation_simple %in% region_list) |>
    
    # Group by sample_id
    dplyr::group_by(
      .data$sample_id, 
      .data$brain_region_abbreviation_simple, 
      .data$donor_external_id, 
      .data$total_nuclei
    ) |>
    
    # Calculate sums for numerator and denominator subsets
    dplyr::summarise(
      # Sum nuclei for cell types defined in 'numerator'
      numerator_sum = sum(
        .data$n_nuclei[.data[[cell_type_col]] %in% numerator],
        na.rm = TRUE
      ),
      
      # If denominator is NULL, use all nuclei in the group; otherwise, sum specific subset
      denominator_sum = if (is.null(denominator)) {
        sum(.data$n_nuclei, na.rm = TRUE)
      } else {
        sum(
          .data$n_nuclei[.data[[cell_type_col]] %in% denominator],
          na.rm = TRUE
        )
      },
      .groups = "drop" 
    ) |>
    
    # Compute the final ratio
    dplyr::mutate(ratio = .data$numerator_sum / .data$denominator_sum) |>
    
    # Return columns needed for plotting
    dplyr::select(
      .data$sample_id, 
      .data$brain_region_abbreviation_simple, 
      .data$donor_external_id, 
      .data$total_nuclei, 
      .data$ratio
    )
  
  return(ratio_df)
}


#' Perform Wilcoxon Rank-Sum Test with Effect Sizes
#'
#' This function performs a Mann-Whitney U test (Wilcoxon rank-sum test) and 
#' calculates non-parametric effect sizes, including the Probability of 
#' Superiority and the Rank-Biserial Correlation.
#'
#' @param x A numeric vector of observations from the first group.
#' @param y A numeric vector of observations from the second group.
#'
#' @details 
#' The probability of superiority is calculated as:
#' \deqn{PS = \frac{U}{n_1 n_2}}
#' The rank-biserial correlation is calculated as:
#' \deqn{r_{rb} = 1 - \frac{2U}{n_1 n_2}}
#' where \eqn{U} is the Mann-Whitney U statistic.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{prob_superiority}: The probability that a randomly selected 
#'     value from \code{x} is greater than a randomly selected value from \code{y}.
#'   \item \code{rank_biserial}: The effect size ranging from -1 to 1.
#'   \item \code{p_value}: The p-value from the Wilcoxon test.
#' }
#' 
#' @importFrom stats wilcox.test
#' @export
#'
#' @examples
#' group1 <- c(10, 12, 15, 18)
#' group2 <- c(8, 9, 11, 13)
#' wilcoxon_test(group1, group2)
wilcoxon_test <- function(x, y) {
  # Clean input data by removing missing values to prevent errors in test calculation
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  
  # Execute the non-parametric Mann-Whitney U test using a normal approximation
  test_result <- stats::wilcox.test(x, y, exact = FALSE)
  p_value <- test_result$p.value
  
  # Extract the W statistic, which is functionally equivalent to the U statistic
  U <- test_result$statistic
  
  # Determine sample sizes for both groups to use as denominators for effect sizes
  n1 <- length(x)
  n2 <- length(y)
  
  # Calculate the probability of superiority
  prob_superiority <- U / (n1 * n2)
  
  # Calculate Rank-Biserial Correlation to represent the direction and strength of the difference
  rank_biserial <- 1 - (2 * U) / (n1 * n2)
  
  # Return results in a list
  return(list(
    prob_superiority = prob_superiority,
    rank_biserial = rank_biserial,
    p_value = p_value
  ))
}


#' Create Boxplots for Cell-Type Proportion Ratios
#'
#' This function generates a publication-quality boxplot comparing cell-type 
#' proportion ratios across different brain regions. It optionally includes 
#' individual data points (jitter) and pairwise statistical significance labels.
#'
#' @param ctp_df A data frame containing cell type counts and metadata.
#' @param numerator Character vector of cell types for the ratio numerator.
#' @param denominator Character vector of cell types for the ratio denominator.
#' @param title Character string for the plot title.
#' @param cell_type_col Column name containing cell type annotations. Defaults to `"annotation"`.
#' @param ylabel Character string for the y-axis label. Defaults to `"Ratio"`.
#' @param region_list Character vector of regions to include in the plot.
#' @param region_order Character vector defining the factor levels/order of the x-axis.
#' @param cmap A named character vector for manual color mapping of regions.
#' @param display_pvalues Logical; if TRUE, performs pairwise Wilcoxon tests and adds 
#'   significance bars to the plot.
#' @param jitter Logical; if TRUE, overlays individual data points on the boxplots.
#'
#' @return A \code{ggplot} object.
#' 
#' @importFrom dplyr filter mutate .data
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter scale_fill_manual labs theme_minimal theme element_text
#' @importFrom ggsignif geom_signif
#' @export
make_ctp_ratio_boxplot <- function(
    ctp_df,
    numerator,
    denominator,
    title,
    cell_type_col = "annotation",
    ylabel = "Ratio",
    region_list = c("CaH", "Pu", "NAC"),
    region_order = c("CaH", "Pu", "NAC", "ic", "DFC"),
    cmap = c(
      "CaH" = "#E69F00",
      "Pu" = "#56B4E9",
      "NAC" = "#009E73",
      "ic" = "#F0E442",
      "DFC" = "purple"
    ),
    display_pvalues = TRUE,
    jitter = TRUE
) {
  # Generate the underlying ratio table and clean up missing values for plotting
  ratio_df <- generate_ctp_ratio_table(ctp_df, numerator, denominator, cell_type_col, region_list) |>
    dplyr::filter(!is.na(.data$ratio)) |> 
    # Ensure regions appear in the specific order defined by the user for the x-axis
    dplyr::mutate(
      brain_region_abbreviation_simple = factor(
        .data$brain_region_abbreviation_simple, 
        levels = region_order
      )
    )
  
  # Initialize the ggplot object with aesthetics for region and ratio
  p <- ggplot2::ggplot(
    ratio_df, 
    ggplot2::aes(
      x = .data$brain_region_abbreviation_simple, 
      y = .data$ratio,
      fill = .data$brain_region_abbreviation_simple
    )
  ) 
  
  # Add layers conditionally based on whether the user wants to see individual distribution points
  if (jitter) {
    p <- p +
      ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.85) +
      ggplot2::geom_jitter(width = 0.2, height = 0, alpha = 0.5, seed=42)
  } else {
    p <- p +
      ggplot2::geom_boxplot(alpha = 0.7, outlier.size = 1, width = 0.85)
  }
  
  # Calculate and display p-values for all pairwise comparisons among the selected regions
  if (display_pvalues) {
    # Determine every possible unique pairing of the requested regions
    region_pairs <- combn(region_list, 2, simplify = FALSE)
    
    # Iterate through pairs to calculate p-values using the helper Wilcoxon function
    p_values <- sapply(region_pairs, function(pair) {
      group1 <- ratio_df$ratio[ratio_df$brain_region_abbreviation_simple == pair[1]]
      group2 <- ratio_df$ratio[ratio_df$brain_region_abbreviation_simple == pair[2]]
      test_result <- wilcoxon_test(group1, group2)
      return(test_result$p_value)
    })
    
    # Format p-values into reader-friendly labels (Scientific notation or "NS")
    p_labels <- sapply(p_values, function(p) {
      if (is.na(p)) "NA"
      else if (p >= 0.05) "NS"
      else sprintf("p = %.2e", p)
    })
    
    # Append the significance bars and text labels to the plot
    p <- p + ggsignif::geom_signif(
      comparisons = region_pairs,
      annotations = p_labels,
      step_increase = 0.1,
      textsize = 3,
      vjust = 0
    )
  }
  
  # Apply final aesthetic refinements, custom colors, and labels
  p <- p +
    ggplot2::scale_fill_manual(values = cmap) +
    ggplot2::labs(
      title = title,
      x = "Brain Region",
      y = ylabel
    ) +
    ggplot2::theme_minimal() +
    # Center x-axis text and remove legend as color redundant with x-axis labels
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      legend.position = "none"
    )
  
  return(p)
}


# ---------------- EXAMPLE USAGE ------------------

annotations_df <- data.table::fread(
  file.path("/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/cell_type_proportions/data/LEVEL_1/donor_region.annotation.neurons_only.cell_type_proportions.txt")
)

make_ctp_ratio_boxplot(
  annotations_df,
  numerator="MSN_D1_matrix",
  denominator=NULL, 
  title="D1 matrix proportion in different regions",
  ylabel="D1 prop",
  jitter=TRUE
)

make_ctp_ratio_boxplot(
  annotations_df,
  numerator="MSN_D1_matrix",
  denominator="MSN_D2_matrix", 
  title="D1/D2 matrix ratios",
  ylabel="D1/D2 ratio",
  jitter=TRUE
)

make_ctp_ratio_boxplot(
  annotations_df,
  numerator="MSN_D1_matrix",
  denominator=c("MSN_D1_matrix", "MSN_D2_matrix"), 
  title="D1 fraction of matrix MSNs",
  ylabel="D1 fraction",
  jitter=TRUE
)



cell_types <- unique(annotations_df$annotation)

# Define neurons that do not fit the primary Interneuron/Projection groupings
OOD_neurons <- c("VTR-HTH_glut", "SN-VTR-HTH_GABA", "other_GABA")

# Interneurons: All GABAergic types plus striatal cholinergic, excluding OOD
interneurons <- c(grep("GABA", cell_types, value = TRUE),
                  "striatal_cholinergic") |> setdiff(OOD_neurons)


MSNs <- grep("MSN", cell_types, value = TRUE)

make_ctp_ratio_boxplot(
  annotations_df,
  numerator=MSNs,
  denominator=interneurons, 
  title="MSN / Interneuron ratios across striatal gray matter",
  ylabel="MSN / Interneuron ratio",
  jitter=TRUE
)

make_ctp_ratio_boxplot(
  annotations_df,
  numerator=interneurons,,
  denominator=MSNs, 
  title="Interneuron / MSN ratios across striatal gray matter",
  ylabel="Interneuron / MSN ratio",
  jitter=TRUE
)












