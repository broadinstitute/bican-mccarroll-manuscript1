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
      ggplot2::geom_jitter(width = 0.2, height = 0, alpha = 0.5)
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
      legend.position = "none",
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank()
    )
  
  return(p)
}



#' Create Faceted Boxplots for Cell-Type Proportion Ratios
#'
#' This function generates multiple boxplots faceted by ratio types, using
#' the base `generate_ctp_ratio_table` function iteratively.
#'
#' @inheritParams make_ctp_ratio_boxplot
#' @param adjust_pvalues Logical; if TRUE, applies Bonferroni correction to p-values within each facet.
#'
#' @return A patchwork object containing multiple ggplot facets.
#' @importFrom purrr map2_dfr
#' @importFrom patchwork wrap_plots plot_layout plot_annotation
#' @export
make_ctp_ratio_boxplot_faceted <- function(
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
    adjust_pvalues = FALSE,
    jitter = TRUE
) {
  
  # --- 1. Standardize Inputs ---
  # Convert to lists if they aren't already to allow uniform iteration
  if (!is.list(numerator)) {
    numerator_list <- list("Ratio" = numerator)
    denominator_list <- list("Ratio" = denominator)
  } else {
    numerator_list <- numerator
    denominator_list <- denominator
    
    # If one denominator is provided for a list of numerators, recycle it
    if (!is.list(denominator_list)) {
      denominator_list <- rep(list(denominator_list), length(numerator_list))
      names(denominator_list) <- names(numerator_list)
    }
  }
  
  # --- 2. Generate Data using the base function ---
  full_ratio_df <- purrr::map2_dfr(
    numerator_list, 
    denominator_list, 
    function(num, den) {
      generate_ctp_ratio_table(ctp_df, num, den, cell_type_col, region_list)
    }, 
    .id = "ratio_name"
  ) |>
    dplyr::filter(!is.na(.data$ratio)) |>
    dplyr::mutate(ratio_name = factor(.data$ratio_name, levels = names(numerator_list)))
  
  ratio_names <- levels(full_ratio_df$ratio_name)
  global_y_max <- if(nrow(full_ratio_df) > 0) max(full_ratio_df$ratio, na.rm = TRUE) else 1
  
  # --- 3. Build Individual Plots ---
  plot_list <- lapply(seq_along(ratio_names), function(i) {
    r_name <- ratio_names[i]
    is_leftmost <- (i == 1)
    
    sub_df <- full_ratio_df |>
      dplyr::filter(.data$ratio_name == r_name) |>
      dplyr::mutate(
        brain_region_abbreviation_simple = factor(
          .data$brain_region_abbreviation_simple, 
          levels = region_order
        )
      )
    
    p <- ggplot2::ggplot(sub_df, ggplot2::aes(
      x = .data$brain_region_abbreviation_simple, 
      y = .data$ratio, 
      fill = .data$brain_region_abbreviation_simple
    ))
    
    # Plotting Layers
    if (jitter) {
      p <- p +
        ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.85) +
        ggplot2::geom_jitter(width = 0.2, height = 0, alpha = 0.5)
    } else {
      p <- p + ggplot2::geom_boxplot(alpha = 0.7, outlier.size = 1, width = 0.85)
    }
    
    # Scale Consistency
    p <- p + ggplot2::expand_limits(y = c(0, global_y_max))
    
    # Statistics
    if (display_pvalues && nrow(sub_df) > 0) {
      region_pairs <- combn(region_list, 2, simplify = FALSE)
      p_values <- sapply(region_pairs, function(pair) {
        g1 <- sub_df$ratio[sub_df$brain_region_abbreviation_simple == pair[1]]
        g2 <- sub_df$ratio[sub_df$brain_region_abbreviation_simple == pair[2]]
        if(length(g1) < 2 || length(g2) < 2) return(NA)
        return(wilcoxon_test(g1, g2)$p_value)
      })
      
      if (adjust_pvalues) p_values <- p.adjust(p_values, method = "bonferroni")
      
      p_labels <- sapply(p_values, function(p) {
        if (is.na(p)) "NA" else if (p >= 0.05) "NS" else sprintf("p = %.2e", p)
      })
      
      p <- p + ggsignif::geom_signif(
        comparisons = region_pairs, 
        annotations = p_labels,
        step_increase = 0.1, 
        textsize = 3, 
        vjust = 0
      )
    }
    
    # Facet Styling
    p <- p +
      ggplot2::scale_fill_manual(values = cmap, name = "Brain Region") +
      ggplot2::facet_wrap(~ratio_name) +
      ggplot2::theme_bw() + 
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        strip.background = ggplot2::element_rect(fill = "grey90", color = "black"),
        strip.text = ggplot2::element_text(face = "bold"),
        panel.grid.major.x = ggplot2::element_blank()
      )
    
    # Conditional Axis Labels
    if (is_leftmost) {
      p <- p + ggplot2::labs(x = NULL, y = ylabel)
    } else {
      p <- p + ggplot2::labs(x = NULL, y = NULL) +
        ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank())
    }
    
    return(p)
  })
  
  # --- 4. Combine via Patchwork ---
  combined_plot <- patchwork::wrap_plots(plot_list, nrow = 1) + 
    patchwork::plot_layout(guides = "collect") & 
    ggplot2::theme(legend.position = "right")
  
  combined_plot <- combined_plot + 
    patchwork::plot_annotation(title = title)
  
  return(combined_plot)
}


# ----------------------------------

generate_manuscript_ratio_boxplots <-function(ctp_table_fp, outDir){
  annotations_df_level2 <- data.table::fread(ctp_table_fp)
  
  # ---------- Single plots: ----------
  
  p <- make_ctp_ratio_boxplot(
    annotations_df_level2,
    numerator=c("MSN_D1_striosome", "MSN_D2_striosome"),
    denominator=c("MSN_D1_matrix", "MSN_D2_matrix"), 
    title="Striosome : Matrix ratios between regions",
    ylabel="Striosome : Matrix ratio",
    jitter=TRUE
  )
  
  ggplot2::ggsave(
    filename = file.path(outDir, "striosome_matrix_ratios.svg"),
    plot = p,
    width = 6,
    height = 4
  )
  
  p <- make_ctp_ratio_boxplot(
    annotations_df_level2,
    numerator=c("MSN_D1_striosome", "MSN_D1_matrix"),
    denominator=c("MSN_D2_striosome", "MSN_D2_matrix"), 
    title="D1 : D2 ratios between regions",
    ylabel="D1 : D2 ratio",
    jitter=TRUE
  )
  
  ggplot2::ggsave(
    filename = file.path(outDir, "D1_D2_ratios.svg"),
    plot = p,
    width = 6,
    height = 4
  )
  
  # ---------- Faceted plots: ----------
  numerators <- list(
    "D1 : D2 ratio in striosome" = c("MSN_D1_striosome"),
    "D1 : D2 ratio in matrix" = c("MSN_D1_matrix")
  )
  
  denominators <- list(
    "D1 : D2 ratio in striosome" = c("MSN_D2_striosome"),
    "D1 : D2 ratio in matrix" = c("MSN_D2_matrix")
  )
  
  p <- make_ctp_ratio_boxplot_faceted(
    annotations_df_level2,
    numerator = numerators,
    denominator = denominators,
    title = "D1 : D2 ratios in striosome and matrix",
    ylabel = "D1 : D2 ratio",
    jitter = FALSE,
    adjust_pvalues = TRUE
  )
  ggplot2::ggsave(
    filename = file.path(outDir, "D1_D2_ratios_faceted.svg"),
    plot = p,
    width = 8,
    height = 4
  )
  
  numerators <- list(
    "Matrix : Striosome ratio of D1 MSNs" = c("MSN_D1_matrix"),
    "Matrix : Striosome ratio of D2 MSNs" = c("MSN_D2_matrix")
  )
  denominators <- list(
    "Matrix : Striosome ratio of D1 MSNs" = c("MSN_D1_striosome"),
    "Matrix : Striosome ratio of D2 MSNs" = c("MSN_D2_striosome")
  )
  
  p <- make_ctp_ratio_boxplot_faceted(
    annotations_df_level2,
    numerator = numerators,
    denominator = denominators,
    title = "Matrix : Striosome ratios for D1 and D2 MSNs",
    ylabel = "Matrix : Striosome ratio",
    jitter = FALSE,
    adjust_pvalues = TRUE
  )
  ggplot2::ggsave(
    filename = file.path(outDir, "matrix_striosome_ratios_faceted.svg"),
    plot = p,
    width = 8,
    height = 4
  )
  
  
  
  cell_types <- unique(annotations_df_level2$annotation)
  
  OOD_neurons <- c("VTR-HTH_glut", "SN-VTR-HTH_GABA", "other_GABA")
  
  # Interneurons: All GABAergic types plus striatal cholinergic, excluding OOD
  interneurons <- c(grep("GABA", cell_types, value = TRUE),
                    "striatal_cholinergic") |> setdiff(OOD_neurons)
  
  numerators <- list(
    "GABA PTHLH-PVALB" = "striatal_GABA_MGE_PTHLH-PVALB",
    "GABA TAC3-PLPP4" = "striatal_GABA_MGE_TAC3-PLPP4",
    "GABA SST-CHODL" = "GABA_CGE_SST-CHODL",
    "GABA cholinergic" = "striatal_cholinergic"
  )
  denominators <- list(
    "GABA PTHLH-PVALB" = interneurons,
    "GABA TAC3-PLPP4" = interneurons,
    "GABA SST-CHODL" = interneurons,
    "GABA cholinergic" = interneurons
  )
  

  
  p <- make_ctp_ratio_boxplot_faceted(
    annotations_df_level2,
    numerator = numerators,
    denominator = denominators,
    title = "Interneuron subtype ratios across striatal gray matter",
    ylabel = "Proportion of interneurons",
    jitter = FALSE,
    adjust_pvalues = TRUE
  )
  
  ggplot2::ggsave(
    filename = file.path(outDir, "interneuron_ratios_faceted.svg"),
    plot = p,
    width = 12,
    height = 4
  )

}



outDir <- "/Users/emuratog/Documents/misc_outs/manuscript1_boxplots"
ctp_table_fp <- file.path("/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/cell_type_proportions/data/LEVEL_2/donor_region.annotation.neurons_only.cell_type_proportions.txt")

generate_manuscript_ratio_boxplots(ctp_table_fp, outDir)


