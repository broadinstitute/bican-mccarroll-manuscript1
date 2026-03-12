#
# root_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis"
# figures_out_dir = "/broad/mccarroll/yooolivi/projects/bican/manuscript_1_figures/figures"
# figures_cache_dir = "/broad/mccarroll/yooolivi/projects/bican/manuscript_1_figures/data_cache"
#
# sample_ctp <- read.table(
#   file.path(figures_cache_dir, "donor_region.annotation.ctp.txt"),
# )
#
# d1_d2_ratio_df <- read.table(
#   file.path(figures_cache_dir, "donor_region.D1_D2_MSN_ratio.txt"),
#   sep="\t", header=TRUE, stringsAsFactors = FALSE
# )
#
# glia_neuron_ratio_df <- read.table(
#   file.path(figures_cache_dir, "donor_region.glial_neuron_ratios.txt"),
#   sep="\t", header=TRUE, stringsAsFactors = FALSE
# )


# plot_ctp_correlation(
#   sample_ctp,
#   "CaH",
#   "annotation",
#   "astrocyte",
#   "extreme_ventral_MSN"
# )
#
# plot_ctp_region_correlation(
#   ctp_df=sample_ctp,
#   cell_type="OPC",
#   cell_type_col="annotation",
#   region1="CaH",
#   region2="Pu",
#   region_col="brain_region_abbreviation_simple",
#   donor_col="donor_external_id",
#   compute_correlation=TRUE
# )
#
# plot_ratio_region_correlation(
#   ratio_df=d1_d2_ratio_df,
#   region1="CaH",
#   region2="Pu",
#   region_col="brain_region_abbreviation_simple",
#   donor_col="donor_external_id",
#   metric_name="D1/D2 MSN ratio",
#   drop_outliers=TRUE,
#   compute_correlation=TRUE
# )
#
#
# plot_correlation_heatmap(
#   ctp_df=sample_ctp,
#   cell_types=c("astrocyte", "OPC", "oligodendrocyte", "microglia"),
#   regions=c("CaH", "Pu", "NAC")
# )


#' Make a wide-format dataframe of donor cell type proportions, with columns
#' for each cell type or brain region.
#'
#' @param ctp_df Dataframe containing cell type proportions, with columns for
#' cell type, brain region, donor ID, and the CTP metric.
#' @param filter_value Optional value to filter the dataframe by
#' (e.g., a specific cell type or brain region).
#' @param filter_col Column of dataframe to apply the filter on (e.g., cell type
#' column or brain region column).
#' @param pivot_col Column to pivot on to create wide format (e.g., cell type
#' column or brain region column).
#' @param cell_type_col Column of dataframe containing cell type labels.
#' @param region_col Column of dataframe containing brain region labels.
#' @param donor_col Column of dataframe containing donor IDs.
#' @param metric_col Column of dataframe containing the CTP metric to be plotted (e.g., fraction of nuclei).
#' @param drop_outliers Whether to filter out rows marked as outliers
#' (if an "outlier" column exists in the dataframe). Default is TRUE.
#'
#' @return A wide-format dataframe with rows for each donor and columns for each
#' cell type or brain region (depending on the pivot_col), containing the specified CTP metric.
make_donor_ctp_wide <- function(
    ctp_df,
    filter_value=NULL,
    filter_col=NULL,
    pivot_cols=NULL,
    cell_type_col="annotation",
    region_col="brain_region_abbreviation_simple",
    donor_col="donor_external_id",
    metric_col="fraction_nuclei",
    drop_outliers=TRUE
) {

  if (!is.null(filter_value) & !is.null(filter_col)) {
    ctp_df <- ctp_df |>
      dplyr::filter(.data[[filter_col]] == filter_value)
  }

  if (drop_outliers & "outlier" %in% colnames(ctp_df)) {
    ctp_df <- ctp_df |>
      dplyr::filter(!outlier | is.na(outlier))
  }

  if (!is.null(pivot_cols)) {
    wide_df <- ctp_df |>
      dplyr::select(all_of(c(donor_col, pivot_cols, metric_col))) |>
      tidyr::pivot_wider(
        names_from = all_of(pivot_cols),
        values_from = all_of(metric_col)
      )
  } else {
    wide_df <- ctp_df |>
      dplyr::select(all_of(c(donor_col, metric_col)))
  }

  return(wide_df)

}

#
# make_donor_ctp_ratio_wide <- function(
#     ratio_df,
#     region_col="brain_region_abbreviation_simple",
#     donor_col="donor_external_id",
#     drop_outliers=TRUE
# ) {
#
#   if (drop_outliers) {
#     ratio_df_filtered <- ratio_df |>
#       dplyr::filter(!outlier)
#   } else {
#     ratio_df_filtered <- ratio_df
#   }
#
#   wide_df <- ratio_df_filtered |>
#     dplyr::select(all_of(donor_col, region_col), ratio) |>
#     tidyr::pivot_wider(names_from = all_of(region_col), values_from = ratio)
#
#   return(wide_df)
#
# }

#' Generate a scatterplot comparing a specified CTP metric between two variables
#' (e.g., cell types or brain regions), with optional annotation of correlation statistics.
#'
#' @param ctp_df_wide A wide-format dataframe containing the CTP metric for each
#' donor, with columns for the two variables to compare.
#' @param var1 Name of the first variable (e.g., cell type or brain region) to compare on the x-axis.
#' @param var2 Name of the second variable (e.g., cell type or brain region) to compare on the y-axis.
#' @param metric_name Name of the CTP metric being compared (used for axis labels).
#' @param compute_correlation Whether to compute and annotate the scatterplot with Spearman correlation statistics. Default is TRUE.
#' @param correlation_rho Optional pre-computed Spearman correlation coefficient.
#' @param correlation_pvalue Optional pre-computed p-value for the Spearman correlation.
#'
#' @return A ggplot object of the scatterplot comparing the CTP metric between
#' the two variables, with optional correlation annotation.
plot_metric_correlation <- function(
    ctp_df_wide,
    var1,
    var2,
    metric_name,
    compute_correlation = TRUE,
    correlation_rho = NULL,
    correlation_pvalue = NULL
) {

  plot_df <- ctp_df_wide |>
    dplyr::select(all_of(c(var1, var2))) |>
    tidyr::drop_na()

  max_value <- max(plot_df[[var1]], plot_df[[var2]], na.rm = TRUE)
  min_value <- min(plot_df[[var1]], plot_df[[var2]], na.rm = TRUE)

  num_samples <- nrow(plot_df)
  logger::log_info("N={num_samples}.")

  ctp_scatterplot <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = !!ggplot2::sym(var1), y = !!ggplot2::sym(var2))
  ) +
    ggplot2::geom_abline(slope=1, intercept=0, lty="dashed", col="gray") +
    ggplot2::geom_point() +
    ggplot2::xlim(min_value, max_value) +
    ggplot2::ylim(min_value, max_value) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x=sprintf("%s (%s)", metric_name, var1),
      y=sprintf("%s (%s)", metric_name, var2)
    )

  if (!is.null(correlation_rho) & !is.null(correlation_pvalue)) {

    p_label <- ifelse(
      correlation_pvalue < 1e-16,
      "p < 1e-16",
      ifelse(
        correlation_pvalue < 1e-4,
        paste0("p = ",
               format(correlation_pvalue,
                      scientific = TRUE,
                      digits = 2)),
        paste0("p = ",
               format(correlation_pvalue,
                      scientific = FALSE,
                      digits = 3))
      )
    )

    logger::log_info("Annotating scatterplot with provided correlation results: rho = {round(correlation_rho, 3)}, {p_label}")

    rho_label <- paste0("\u03C1 = ", format(correlation_rho, digits=3))

    ctp_scatterplot +
      ggplot2::annotate(
        "text", x=-Inf, y=Inf,
        hjust = -0.2, vjust = 3.5,
        label = paste0(rho_label, ", ", p_label)
      )

    ctp_scatterplot <- ctp_scatterplot +
      ggplot2::annotate(
        "text", x=-Inf, y=Inf,
        hjust = -0.2, vjust = 3.5,
        label = paste0(rho_label, ", ", p_label)
      )
  } else if (compute_correlation) {
    logger::log_info("Computing Spearman correlation for scatterplot annotation.")

    ctp_scatterplot <- ctp_scatterplot +
      ggpubr::stat_cor(method="spearman", cor.coef.name="rho")
  } else if (is.null(correlation_rho) | is.null(correlation_pvalue)) {
    logger::log_info("Correlation results not provided. Scatterplot will not be annotated with correlation statistics.")
  }

  return(ctp_scatterplot)

}


#' Generate a scatterplot comparing CTP between two cell types within a specified
#' brain region, with optional annotation of correlation statistics.
#'
#' @param ctp_df Dataframe containing cell type proportions, with columns for
#' cell type, brain region, donor ID, and fraction of nuclei.
#' @param region Brain region to filter the dataframe by for the scatterplot.
#' @param cell_type_col Column of dataframe containing cell type labels.
#' @param cell_type1 First cell type to compare on the x-axis.
#' @param cell_type2 Second cell type to compare on the y-axis.
#' @param region_col Column of dataframe containing brain region labels.
#' @param donor_col Column of dataframe containing donor IDs.
#' @param drop_outliers Whether to filter out rows marked as outliers (if any).
#' @param compute_correlation Whether to compute and annotate the scatterplot with Spearman correlation.
#' @param correlation_rho Optional pre-computed Spearman correlation.
#' @param correlation_pvalue Optional pre-computed p-value for the Spearman correlation.
#'
#' @return A ggplot object of the scatterplot comparing CTP between the two cell
#'  types within the specified brain region, with optional correlation annotation.
plot_ctp_correlation <- function(
    ctp_df,
    region,
    cell_type_col,
    cell_type1,
    cell_type2,
    region_col="brain_region_abbreviation_simple",
    donor_col="donor_external_id",
    metric_col="fraction_nuclei",
    drop_outliers=TRUE,
    compute_correlation=TRUE,
    correlation_rho=NULL,
    correlation_pvalue=NULL
) {

  logger::log_info("Plotting correlation of CTP between cell types {cell_type1} and {cell_type2} in region {region}.")

  ctp_df_wide <- make_donor_ctp_wide(
    ctp_df=ctp_df,
    filter_value=region,
    filter_col=region_col,
    pivot_cols=cell_type_col,
    cell_type_col=cell_type_col,
    region_col=region_col,
    donor_col=donor_col,
    metric_col=metric_col,
    drop_outliers=drop_outliers
  )

  plot_metric_correlation(
    ctp_df_wide=ctp_df_wide,
    var1=cell_type1,
    var2=cell_type2,
    metric_name=region,
    compute_correlation=compute_correlation,
    correlation_rho=correlation_rho,
    correlation_pvalue=correlation_pvalue
  )

}


#' Generate a scatterplot comparing CTP for a specific cell type between two brain
#' regions, with optional annotation of correlation statistics.
#'
#' @param ctp_df Dataframe containing cell type proportions, with columns for
#' cell type, brain region, donor ID, and fraction of nuclei.
#' @param cell_type Cell type to filter the dataframe by for the scatterplot.
#' @param cell_type_col Column of dataframe containing cell type labels.
#' @param region1 First brain region to compare on the x-axis.
#' @param region2 Second brain region to compare on the y-axis.
#' @param region_col Column of dataframe containing brain region labels.
#' @param donor_col Column of dataframe containing donor IDs.
#' @param drop_outliers Whether to filter out rows marked as outliers (if any).
#' @param compute_correlation Whether to compute and annotate the scatterplot with Spearman correlation.
#' @param correlation_rho Optional pre-computed Spearman correlation.
#' @param correlation_pvalue Optional pre-computed p-value for the Spearman correlation.
#'
#' @return A ggplot object of the scatterplot comparing CTP for the specified cell type
#' between the two brain regions, with optional correlation annotation.
plot_ctp_region_correlation <- function(
    ctp_df,
    cell_type,
    cell_type_col,
    region1,
    region2,
    region_col="brain_region_abbreviation_simple",
    donor_col="donor_external_id",
    metric_name="fraction",
    metric_col="fraction_nuclei",
    drop_outliers=TRUE,
    compute_correlation=TRUE,
    correlation_rho=NULL,
    correlation_pvalue=NULL
) {

  logger::log_info("Plotting correlation of CTP for cell type {cell_type} between regions {region1} and {region2}.")

  ctp_df_wide <- make_donor_ctp_wide(
    ctp_df=ctp_df,
    filter_value=cell_type,
    filter_col=cell_type_col,
    pivot_cols=region_col,
    cell_type_col=cell_type_col,
    region_col=region_col,
    donor_col=donor_col,
    metric_col=metric_col,
    drop_outliers=drop_outliers
  )

  plot_metric_correlation(
    ctp_df_wide=ctp_df_wide,
    var1=region1,
    var2=region2,
    metric_name=paste(cell_type, metric_name),
    compute_correlation=compute_correlation,
    correlation_rho=correlation_rho,
    correlation_pvalue=correlation_pvalue
  )

}



# plot_ratio_region_correlation <- function(
#     ratio_df,
#     region1,
#     region2,
#     region_col="brain_region_abbreviation_simple",
#     donor_col="donor_external_id",
#     metric_name,
#     drop_outliers=TRUE,
#     compute_correlation=TRUE,
#     correlation_rho=NULL,
#     correlation_pvalue=NULL
# ) {
#
#   logger::log_info("Plotting correlation of {metric_name} between regions {region1} and {region2}.")
#
#   ratio_df_wide <- make_ctp_ratio_wide(
#     ratio_df=ratio_df,
#     region_col=region_col,
#     donor_col=donor_col,
#     drop_outliers=drop_outliers
#   )
#
#   plot_metric_correlation(
#     ctp_df_wide=ratio_df_wide,
#     var1=region1,
#     var2=region2,
#     metric_name=metric_name,
#     compute_correlation=compute_correlation,
#     correlation_rho=correlation_rho,
#     correlation_pvalue=correlation_pvalue
#   )
#
# }


#' Compute pairwise Spearman correlations of a specified CTP metric between combinations
#' of cell types and brain regions, returning a long-format dataframe of correlation results.
#'
#' @param ctp_df Dataframe containing cell type proportions, with columns for
#' cell type, brain region, donor ID, and the CTP metric.
#' @param cell_types Vector of cell types to include in the correlation analysis.
#' @param regions Vector of brain regions to include in the correlation analysis.
#' @param cell_type_col Column of dataframe containing cell type labels.
#' @param region_col Column of dataframe containing brain region labels.
#' @param donor_col Column of dataframe containing donor IDs.
#' @param metric_col Column of dataframe containing the CTP metric to be analyzed (e.g. fraction of nuclei).
#' @param drop_outliers Whether to filter out rows marked as outliers (if any) before computing correlations. Default is TRUE.
#'
#' @return A long-format dataframe containing pairwise Spearman correlation
#'  results between combinations of cell types and brain regions.
compute_ctp_correlations <- function(
    ctp_df,
    cell_types,
    regions,
    cell_type_col="annotation",
    region_col="brain_region_abbreviation_simple",
    donor_col="donor_external_id",
    metric_col="fraction_nuclei",
    drop_outliers=TRUE
) {

  # drop outliers if present
  if (drop_outliers & "outlier" %in% colnames(ctp_df)) {
    ctp_df <- ctp_df |>
      dplyr::filter(!outlier | is.na(outlier))
  }

  wide_ctp_df <- ctp_df |>
    dplyr::filter(.data[[cell_type_col]] %in% cell_types) |>
    dplyr::filter(.data[[region_col]] %in% regions) |>
    dplyr::mutate(celltype_region = paste(.data[[cell_type_col]], .data[[region_col]], sep="__")) |>
    dplyr::select(all_of(c(donor_col, metric_col)), celltype_region) |>
    tidyr::pivot_wider(names_from = celltype_region, values_from = all_of(metric_col))

  # compute correlations
  cor_res <- Hmisc::rcorr(
    as.matrix(wide_ctp_df |> dplyr::select(-all_of(donor_col))),
    type="spearman")

  # pull correlation coefficients and p-values into long format dataframes
  cor_mat_long <- cor_res$r |>
    as.data.frame() |>
    tibble::rownames_to_column(var="var1") |>
    tidyr::pivot_longer(cols=-var1, names_to="var2", values_to="rho")

  p_mat_long <- cor_res$P |>
    as.data.frame() |>
    tibble::rownames_to_column(var="var1") |>
    tidyr::pivot_longer(cols=-var1, names_to="var2", values_to="p_value")

  # combine
  res_mat_long <- cor_mat_long |>
    dplyr::left_join(p_mat_long, by=c("var1", "var2")) |>
    dplyr::filter(as.character(var1) < as.character(var2)) |>
    tidyr::separate(var1, into=c("cell_type1", "region1"), sep="__") |>
    tidyr::separate(var2, into=c("cell_type2", "region2"), sep="__")


  return(res_mat_long)

}


#' Compute pairwise Spearman correlations of a specified CTP metric between
#' combinations of brain regions for a specific cell type. Comparisons are
#' only made within the same cell type.
#'
#' @param ctp_df Dataframe containing cell type proportions, with columns for
#' cell type, brain region, donor ID, and the CTP metric.
#' @param cell_types Cell types to use.
#' @param regions Vector of brain regions to include.
#' @param cell_type_col Column of dataframe containing cell type labels.
#' @param region_col Column of dataframe containing brain region labels.
#' @param donor_col Column of dataframe containing donor IDs.
#' @param metric_col Column of dataframe containing the CTP metric to be analyzed (e.g. fraction of nuclei).
#' @param drop_outliers Whether to filter out rows marked as outliers (if any) before computing correlations. Default is TRUE.
#'
#' @return A long-format dataframe containing pairwise Spearman correlations.
compute_ctp_correlations_within_type <- function(
    ctp_df,
    cell_types,
    regions,
    cell_type_col="annotation",
    region_col="brain_region_abbreviation_simple",
    donor_col="donor_external_id",
    metric_col="fraction_nuclei",
    drop_outliers=TRUE
) {

  cell_type_correlations <- lapply(cell_types, function(cell_type) {
    logger::log_info("Computing correlations for cell type {cell_type} across regions.")
    res_mat_long <- compute_ctp_correlations(
      ctp_df=ctp_df,
      cell_types=cell_type,
      regions=regions,
      cell_type_col=cell_type_col,
      region_col=region_col,
      donor_col=donor_col,
      metric_col=metric_col,
      drop_outliers=drop_outliers
    )
    return(res_mat_long)
  }) |>
    dplyr::bind_rows()

  return(cell_type_correlations)

}


#' Compute pairwise Spearman correlations of a specified CTP metric between
#' combinations of cell types within a specific brain region. Comparisons are
#' only made within the same brain region.
#'
#' @param ctp_df Dataframe containing cell type proportions, with columns for
#' cell type, brain region, donor ID, and the CTP metric.
#' @param regions Vector of brain regions to include.
#' @param cell_types Cell types to use.
#' @param cell_type_col Column of dataframe containing cell type labels.
#' @param region_col Column of dataframe containing brain region labels.
#' @param donor_col Column of dataframe containing donor IDs.
#' @param metric_col Column of dataframe containing the CTP metric to be analyzed (e.g. fraction of nuclei).
#' @param drop_outliers Whether to filter out rows marked as outliers (if any) before computing correlations. Default is TRUE.
#'
#' @return A long-format dataframe containing pairwise Spearman correlations.
compute_ctp_correlations_within_region <- function(
    ctp_df,
    regions,
    cell_types,
    cell_type_col="annotation",
    region_col="brain_region_abbreviation_simple",
    donor_col="donor_external_id",
    metric_col="fraction_nuclei",
    drop_outliers=TRUE
) {

  region_correlations <- lapply(regions, function(region) {
    logger::log_info("Computing correlations for region {region} across cell types.")
    res_mat_long <- compute_ctp_correlations(
      ctp_df=ctp_df,
      cell_types=cell_types,
      regions=region,
      cell_type_col=cell_type_col,
      region_col=region_col,
      donor_col=donor_col,
      metric_col=metric_col,
      drop_outliers=drop_outliers
    )
    return(res_mat_long)
  }) |>
    dplyr::bind_rows()

  return(region_correlations)

}


#' Generate heatmap of Spearman corrleations.
#'
#' @param ctp_df Dataframe containing cell type proportions, with columns for
#' cell type, brain region, donor ID, and the CTP metric.
#' @param cell_types Cell types to use.
#' @param regions Vector of brain regions to include.
#' @param cell_type_col Column of dataframe containing cell type labels.
#' @param region_col Column of dataframe containing brain region labels.
#' @param donor_col Column of dataframe containing donor IDs.
#' @param metric_col Column of dataframe containing the CTP metric to be analyzed (e.g. fraction of nuclei).
#' @param drop_outliers Whether to filter out rows marked as outliers (if any) before computing correlations. Default is TRUE.
#'
#' @return A ggplot object of the correlation heatmap.
plot_correlation_heatmap <- function(
    ctp_df,
    cell_types,
    regions,
    cell_type_col="annotation",
    region_col="brain_region_abbreviation_simple",
    donor_col="donor_external_id",
    metric_col="fraction_nuclei",
    drop_outliers=TRUE
) {

  # setup for plots
  axis_order <- paste(
    rep(cell_types, each=length(regions)),
    rep(regions, times=length(cell_types)),
    sep="__"
  )

  # compute correlations
  cor_mat_long <- compute_ctp_correlations(
    ctp_df=ctp_df,
    cell_types=cell_types,
    regions=regions,
    cell_type_col=cell_type_col,
    region_col=region_col,
    donor_col=donor_col,
    metric_col=metric_col,
    drop_outliers=drop_outliers
  ) |>
    dplyr::mutate(
      var1 = paste(cell_type1, region1, sep="__"),
      var2 = paste(cell_type2, region2, sep="__")
    ) |>
    dplyr::select(var1, var2, rho)

  # complete them, so heatmap is full
  cor_mat_long <- compute_ctp_correlations(
    ctp_df, cell_types, regions
  ) |>
    dplyr::mutate(
      var1 = paste(cell_type1, region1, sep="__"),
      var2 = paste(cell_type2, region2, sep="__")
    ) |>
    dplyr::select(var1, var2, rho)

  cor_mat_long_complete <- cor_mat_long |>
    dplyr::select(var1, var2, rho)  |>

    # add reversed correlations
    dplyr::bind_rows(
      cor_mat_long |>
        dplyr::rename(var1 = var2, var2 = var1)
    ) |>

    # add diagonals
    tidyr::complete(var1 = axis_order, var2 = axis_order) |>

    # fill diagonal
    dplyr::mutate(
      rho = dplyr::case_when(
        var1 == var2 ~ 1,
        TRUE ~ rho
      )
    ) |>

    # order factors
    dplyr::mutate(
      var1 = factor(var1, levels=axis_order),
      var2 = factor(var2, levels=axis_order)
    )

  # setup for adding horizontal lines to plot
  axis_df <- data.frame(
    var = axis_order,
    pos = seq_along(axis_order)
  )

  # find where cell type changes along the axis
  boundary_positions <- axis_df |>
    dplyr::mutate(cell_type = sub("__.*$", "", var)) |>
    dplyr::group_by(cell_type) |>
    dplyr::summarise(max_pos = max(pos)) |>
    dplyr::pull(max_pos)

  # add zero at the bottom/start
  boundary_positions <- c(0, boundary_positions)
  n_tiles <- length(axis_order)

  # heatmap
  ctp_cor_heatmap <- ggplot2::ggplot(
    cor_mat_long_complete,
    ggplot2::aes(x=var1, y=forcats::fct_rev(var2))
  ) +
    ggplot2::geom_tile(ggplot2::aes(fill=rho)) +
    ggplot2::scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0, limits=c(-1, 1)) +
    # vertical separator lines
    ggplot2::geom_segment(
      data = data.frame(x = boundary_positions + 0.5),
      ggplot2::aes(x = x, xend = x, y = 0.5, yend = n_tiles + 0.5),
      color = "black"
    ) +

    # horizontal separator lines
    ggplot2::geom_segment(
      data = data.frame(y = boundary_positions + 0.5),
      ggplot2::aes(x = 0.5, xend = n_tiles + 0.5, y = y, yend = y),
      color = "black"
    ) +
    ggplot2::labs(
      x=NULL,
      y=NULL,
      fill="Spearman rho"
    ) +
    ggplot2::scale_x_discrete(
      labels = \(x) sub("__", " (", x) |> paste0(")")
    ) +
    ggplot2::scale_y_discrete(
      labels = \(x) sub("__", " (", x) |> paste0(")")
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle=90, vjust=0.5, hjust=1),
      axis.ticks.x=ggplot2::element_blank(),
      axis.ticks.y=ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank()
    )

  return(ctp_cor_heatmap)

}



