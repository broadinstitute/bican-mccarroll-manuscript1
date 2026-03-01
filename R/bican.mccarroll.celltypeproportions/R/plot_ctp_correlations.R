#library(dplyr)
#library(ggplot2)
#library(cowplot)
#library(logger)
#library(tidyr)
#
# ctp_df <- read.table("/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/cell_type_proportions/LEVEL_1/donor_region.cell_type_proportions.txt", header=TRUE, sep="\t")
#
# cell_type_col <- "annotation_most_specific"
# region_col <- "brain_region_abbreviation_simple"
# donor_col <- "donor_external_id"
#
# region1 <- "CaH"
# region2 <- "DFC"
# cell_type <- "astrocyte"
# regions <- c("CaH", "DFC", "ic", "NAC", "Pu")
#
# get_ctp_by_region(ctp_df, "astrocyte", cell_type_col, region_col, donor_col)
# plot_ctp_region_pair(ctp_df, region1, region2, cell_type, cell_type_col, region_col, donor_col)
# plot_ctp_region_pairs(ctp_df, regions, cell_type, cell_type_col, region_col, donor_col)
# plot_ctp_region_pair_by_celltype(ctp_df, region1, region2, c("astrocyte", "microglia", "oligodendrocyte", "OPC"), cell_type_col, region_col, donor_col)


#' Retrieves cell type proportions for a specific cell type across different brain regions.
#'
#' @param ctp_df Dataframe containing cell type proportions, with columns for cell type, brain region, donor ID, and fraction of nuclei.
#' @param cell_type Cell type of interest.
#' @param cell_type_col Column of dataframe containing cell type labels.
#' @param region_col Column of dataframe containing brain region labels.
#' @param donor_col Column of dataframe containing donor IDs.
#'
#' @return wide dataframe of CTP, with rows for each donor and columns for each brain region.
get_ctp_by_region <- function(ctp_df, cell_type, cell_type_col, region_col, donor_col) {

  cell_type_df <- ctp_df |>
    dplyr::filter(.data[[cell_type_col]] == cell_type) |>
    dplyr::select(!!dplyr::sym(donor_col), !!dplyr::sym(region_col), fraction_nuclei) |>
    tidyr::pivot_wider(names_from = !!dplyr::sym(region_col), values_from = fraction_nuclei)

  return(cell_type_df)
}

#' Generates a scatterplot comparing cell type proportions between two brain regions for a specific cell type.
#'
#' @param ctp_df Dataframe containing cell type proportions, with columns for cell type, brain region, donor ID, and fraction of nuclei.
#' @param region1 First brain region to compare.
#' @param region2 Second brain region to compare.
#' @param cell_type Cell type of interest.
#' @param cell_type_col Column of dataframe containing cell type labels.
#' @param region_col Column of dataframe containing brain region labels.
#' @param donor_col Column of dataframe containing donor IDs.
#'
#' @return ggplot object of scatterplot comparing cell type proportions between the two regions.
plot_ctp_region_pair <- function(ctp_df, region1, region2, cell_type, cell_type_col, region_col, donor_col) {

  cell_type_df <- get_ctp_by_region(ctp_df, cell_type, cell_type_col, region_col, donor_col)
  plot_df <- cell_type_df |>
    dplyr::select(!!dplyr::sym(donor_col), !!dplyr::sym(region1), !!dplyr::sym(region2)) |>
    tidyr::drop_na()

  max_prop <- max(unlist(cell_type_df[sapply(cell_type_df, is.numeric)]), na.rm = TRUE)
  n_donors <- nrow(plot_df)

  region_scatterplot <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x=!!dplyr::sym(region1), y=!!dplyr::sym(region2))
  ) +
    ggplot2::geom_abline(slope=1, intercept=0, linetype="dashed", color="gray") +
    ggplot2::geom_point(alpha=0.75) +
    ggpubr::stat_cor(method="spearman", cor.coef.name="rho") +
    ggplot2::labs(
      title=sprintf("%s fraction", cell_type),
      subtitle=sprintf("N=%s donors", n_donors)
    ) +
    ggplot2::xlim(0, max_prop) +
    ggplot2::ylim(0, max_prop) +
    ggplot2::theme_bw()

  return(region_scatterplot)

}


#' Generates scatterplots comparing cell type proportions between all pairs of specified brain regions for a specific cell type.
#'
#' @param ctp_df Dataframe containing cell type proportions, with columns for cell type, brain region, donor ID, and fraction of nuclei.
#' @param regions Character vector of brain regions to compare.
#' @param cell_type Cell type of interest.
#' @param cell_type_col Column of dataframe containing cell type labels.
#' @param region_col Column of dataframe containing brain region labels.
#' @param donor_col Column of dataframe containing donor IDs.
#' @param nrow Number of rows to use in combined plot layout (default is 2).
#'
#' @return cowplot object of combined scatterplots comparing cell type proportions between all pairs of specified brain regions.
plot_ctp_region_pairs <- function(ctp_df, regions, cell_type, cell_type_col, region_col, donor_col, nrow=2) {

  region_pairs <- combn(regions, 2)
  n_region_pairs <- ncol(region_pairs)

  region_plot_list <- list()

  for (i in 1:n_region_pairs) {
    region1 <- region_pairs[1, i]
    region2 <- region_pairs[2, i]

    plot <- plot_ctp_region_pair(ctp_df, region1, region2, cell_type, cell_type_col, region_col, donor_col)
    region_plot_list[[paste(region1, region2, sep = "_")]] <- plot
  }

  combined_plot <- cowplot::plot_grid(
    plotlist=region_plot_list,
    nrow=nrow
  )

  return(combined_plot)

}


#' Generates scatterplots comparing cell type proportions between two specified brain regions for multiple cell types.
#'
#' @param ctp_df Dataframe containing cell type proportions, with columns for cell type, brain region, donor ID, and fraction of nuclei.
#' @param region1 First brain region to compare.
#' @param region2 Second brain region to compare.
#' @param cell_types Character vector of cell types to compare.
#' @param cell_type_col Column of dataframe containing cell type labels.
#' @param region_col Column of dataframe containing brain region labels.
#' @param donor_col Column of dataframe containing donor IDs.
#' @param nrow Number of rows to use in combined plot layout (default is 1).
#'
#' @return cowplot object of combined scatterplots comparing cell type proportions between the two regions for multiple cell types.
plot_ctp_region_pair_by_celltype <- function(ctp_df, region1, region2, cell_types, cell_type_col, region_col, donor_col, nrow=1) {

  region_plot_list <- list()

  for (cell_type in cell_types) {
    plot <- plot_ctp_region_pair(ctp_df, region1, region2, cell_type, cell_type_col, region_col, donor_col)
    region_plot_list[[cell_type]] <- plot
  }

  combined_plot <- cowplot::plot_grid(
    plotlist=region_plot_list,
    nrow=nrow
  )

  return(combined_plot)

}
