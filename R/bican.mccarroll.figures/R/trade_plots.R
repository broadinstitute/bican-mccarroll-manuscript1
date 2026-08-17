# source("R/paths.R")
#
# options(
#     bican.mccarroll.figures.data_root_dir =
#         "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis",
#
#     bican.mccarroll.figures.out_dir =
#         "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/figure_repository",
#
#     bican.mccarroll.figures.cache_dir =
#         "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/figure_repository/data_cache"
# )

# de_dir <- de_region_subset_dir <- de_region_interaction_dir <- gene_to_chr_path <- ct_file <- data_cache_dir <- outDir <- NULL

#' Generate the TRADE manuscript figure for the BICAN analysis
#'
#' Plots the BICAN TRADE analysis of autosome age effects across cell types.
#' This adapter preserves the manuscript-specific paths, cache name, cell type
#' file, plot formatting, and output dimensions while delegating the analysis
#' and plotting to the shared private implementation.
#'
#' @param force_recompute Logical scalar. If \code{TRUE}, ignore an existing
#'   cache file, recompute TRADE results, and overwrite the cache. Defaults to
#'   \code{FALSE}.
#' @param data_cache_dir Directory used to store cached TRADE results. If
#'   \code{NULL}, the directory is resolved from
#'   \code{options("bican.mccarroll.figures.cache_dir")}.
#' @param outDir Output directory used to write the SVG file. If \code{NULL},
#'   the directory is resolved from
#'   \code{options("bican.mccarroll.figures.out_dir")}.
#'
#' @return Invisibly returns \code{NULL}. This function is called for its side
#'   effects.
#'
#' @seealso
#'   \code{bican.mccarroll.differentialexpression::load_trade_data},
#'   \code{bican.mccarroll.differentialexpression::run_trade},
#'   \code{bican.mccarroll.differentialexpression::trade_barplot}
#' @export
plot_trade_analysis_bican <- function(force_recompute = FALSE,
                                      data_cache_dir = NULL,
                                      outDir = NULL) {
  .plot_trade_analysis(
    de_dir = NULL,
    cache_name = "trade_BICAN_age_autosomes.tsv",
    dataset_id = "BICAN",
    gene_to_chr_path = NULL,
    ct_file = NULL,
    data_cache_dir = data_cache_dir,
    outDir = outDir,
    hide_cell_type_axis = TRUE,
    x_breaks = c(0, 0.001, 0.002),
    x_labels = c("0.000", "0.001", "0.002"),
    axis_title_x_size = ggplot2::rel(1.5),
    axis_text_x_size = ggplot2::rel(2.25),
    axis_tick_x_linewidth = 0.7,
    plot_margin = ggplot2::margin(5.5, 5.5, 5.5, 12),
    bar_fill = "black",
    x_label = "Transcriptome-wide impact (TRADE)",
    width = 5,
    height = 9,
    force_recompute = force_recompute
  )

  # this was purely to provide some labels to compare individual plots
  # vs other papers.
  .plot_trade_analysis(
    de_dir = NULL,
    cache_name = "trade_BICAN_age_autosomes.tsv",
    dataset_id = "BICAN_legend",
    gene_to_chr_path = NULL,
    ct_file = NULL,
    data_cache_dir = data_cache_dir,
    outDir = outDir,
    hide_cell_type_axis = FALSE,
    x_breaks = c(0, 0.001, 0.002),
    x_labels = c("0.000", "0.001", "0.002"),
    axis_title_x_size = ggplot2::rel(1.5),
    axis_text_x_size = ggplot2::rel(2.25),
    axis_tick_x_linewidth = 0.7,
    plot_margin = ggplot2::margin(5.5, 5.5, 5.5, 12),
    bar_fill = "black",
    x_label = "Transcriptome-wide impact (TRADE)",
    width = 5,
    height = 9,
    force_recompute = force_recompute
  )

  .plot_trade_analysis_regions(
    de_region_interaction_dir = NULL,
    cache_name = "trade_BICAN_age_region_interaction_autosomes.tsv",
    dataset_id = "BICAN",
    region_order = c("CaH", "Pu", "NAC", "ic", "DFC"),
    filter_ic_to_non_neurons = TRUE,
    gene_to_chr_path = NULL,
    ct_file = NULL,
    data_cache_dir = data_cache_dir,
    outDir = outDir,
    width = 5,
    height = 8,
    force_recompute = force_recompute
  )
}

plot_trade_analysis_bican_downsampled <- function(force_recompute = FALSE,
                                                  data_cache_dir = NULL,
                                                  outDir = NULL) {
  de_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/results/LEVEL_6_downsampled/sex_age/cell_type"
  de_region_interaction_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/results/LEVEL_6_downsampled/sex_age/cell_type_region_interaction_absolute_effects"

  .plot_trade_analysis(
    de_dir = de_dir,
    cache_name = "trade_BICAN_age_autosomes_downsampled.tsv",
    dataset_id = "BICAN_legend_downsampled",
    gene_to_chr_path = NULL,
    ct_file = NULL,
    data_cache_dir = data_cache_dir,
    outDir = outDir,
    hide_cell_type_axis = FALSE,
    x_breaks = c(0, 0.001, 0.002),
    x_labels = c("0.000", "0.001", "0.002"),
    axis_title_x_size = ggplot2::rel(1.5),
    axis_text_x_size = ggplot2::rel(2.25),
    axis_tick_x_linewidth = 0.7,
    plot_margin = ggplot2::margin(5.5, 5.5, 5.5, 12),
    bar_fill = "black",
    x_label = "Transcriptome-wide impact (TRADE)",
    width = 5,
    height = 9,
    force_recompute = force_recompute
  )

  .plot_trade_analysis_regions(
    de_region_interaction_dir = de_region_interaction_dir,
    cache_name = "trade_BICAN_age_region_interaction_autosomes_downsampled.tsv",
    dataset_id = "BICAN_downsampled_interaction",
    region_order = c("CaH", "Pu", "NAC", "ic", "DFC"),
    filter_ic_to_non_neurons = TRUE,
    gene_to_chr_path = NULL,
    ct_file = NULL,
    data_cache_dir = data_cache_dir,
    outDir = outDir,
    width = 5,
    height = 8,
    force_recompute = force_recompute
  )

  de_region_subset_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/results/LEVEL_6_downsampled/sex_age/cell_type_subset_region"

  .plot_trade_analysis_regions(
    de_region_interaction_dir = de_region_subset_dir,
    cache_name = "trade_BICAN_age_region_subset_autosomes_downsampled.tsv",
    dataset_id = "BICAN_downsampled_subset",
    region_order = c("CaH", "Pu", "NAC", "ic", "DFC"),
    filter_ic_to_non_neurons = TRUE,
    gene_to_chr_path = NULL,
    ct_file = NULL,
    data_cache_dir = data_cache_dir,
    outDir = outDir,
    width = 5,
    height = 8,
    force_recompute = force_recompute
  )
}

plot_trade_analysis_PMID_39227716 <- function(force_recompute = FALSE,
                                              data_cache_dir = NULL,
                                              outDir = NULL) {
  de_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/external_comparison_PMID_39227716/voom-like"

  .plot_trade_analysis(
    de_dir = de_dir,
    cache_name = "trade_PMID_39227716_age_autosomes.tsv",
    dataset_id = "PMID_39227716",
    gene_to_chr_path = NULL,
    ct_file = NULL,
    data_cache_dir = data_cache_dir,
    outDir = outDir,
    hide_cell_type_axis = FALSE,
    x_breaks = c(0, 0.00001, 0.00002),
    x_labels = c("0.000", "0.00001", "0.00002"),
    axis_title_x_size = ggplot2::rel(1.5),
    axis_text_x_size = ggplot2::rel(2.25),
    axis_tick_x_linewidth = 0.7,
    plot_margin = ggplot2::margin(5.5, 5.5, 5.5, 12),
    bar_fill = "black",
    x_label = "Transcriptome-wide impact (TRADE)",
    width = 5,
    height = 9,
    force_recompute = force_recompute
  )
}

plot_trade_analysis_Ling <- function(force_recompute = FALSE,
                                     data_cache_dir = NULL,
                                     outDir = NULL) {
  de_dir <- "/broad/mccarroll/dropulation/analysis/SNAP200_controls/differential_expression/results/sex_age/cell_type"

  .plot_trade_analysis(
    de_dir = de_dir,
    cache_name = "trade_Ling_age_autosomes.tsv",
    dataset_id = "Ling",
    gene_to_chr_path = NULL,
    ct_file = NULL,
    data_cache_dir = data_cache_dir,
    outDir = outDir,
    hide_cell_type_axis = FALSE,
    x_breaks = c(0, 0.00001, 0.00002),
    x_labels = c("0.000", "0.00001", "0.00002"),
    axis_title_x_size = ggplot2::rel(1.5),
    axis_text_x_size = ggplot2::rel(2.25),
    axis_tick_x_linewidth = 0.7,
    plot_margin = ggplot2::margin(5.5, 5.5, 5.5, 12),
    bar_fill = "black",
    x_label = "Transcriptome-wide impact (TRADE)",
    width = 5,
    height = 9,
    force_recompute = force_recompute
  )
}


#' Generate the TRADE figure for the BICAN sea_ad_mtg_mmc analysis
#'
#' Plots the BICAN TRADE analysis of autosome age effects across the SEA-AD
#' supertype cell types for which BICAN
#' \code{LEVEL_6_sea_ad_mtg_mmc} DE results were generated.
#'
#' @param force_recompute Logical scalar. If \code{TRUE}, ignore an existing
#'   cache file, recompute TRADE results, and overwrite the cache. Defaults to
#'   \code{FALSE}.
#' @param data_cache_dir Directory used to store cached TRADE results. If
#'   \code{NULL}, the directory is resolved from
#'   \code{options("bican.mccarroll.figures.cache_dir")}.
#' @param outDir Output directory used to write the SVG file. If \code{NULL},
#'   the directory is resolved from
#'   \code{options("bican.mccarroll.figures.out_dir")}.
#'
#' @return Invisibly returns \code{NULL}. This function is called for its side
#'   effects.
#'
#' @export
plot_trade_analysis_bican_sea_ad_mtg_mmc <- function(force_recompute = FALSE,
                                                     data_cache_dir = NULL,
                                                     outDir = NULL) {
  de_dir <- "differential_expression/results/LEVEL_6_sea_ad_mtg_mmc/sex_age/cell_type"
  ct_file <- "differential_expression/metadata/cell_types_for_sea_ad_mtg_mmc_plots.txt"

  .plot_trade_analysis(
    de_dir = de_dir,
    cache_name = "trade_BICAN_sea_ad_mtg_mmc_age_autosomes.tsv",
    dataset_id = "BICAN_sea_ad_mtg_mmc",
    gene_to_chr_path = NULL,
    ct_file = ct_file,
    data_cache_dir = data_cache_dir,
    outDir = outDir,
    hide_cell_type_axis = FALSE,
    bar_fill = "black",
    x_label = "Transcriptome-wide impact (TRADE)",
    width = 5,
    height = 9,
    contrast = "__age_DE_results\\.txt$",
    contrast_id = "age",
    force_recompute = force_recompute
  )
}


#' Generate the TRADE figure for the PMID_39402379 (Gabitto et al. 2024) analysis
#'
#' Plots the TRADE analysis of autosome effects across the SEA-AD supertype
#' cell types for which BICAN \code{LEVEL_6_sea_ad_mtg_mmc} DE results were
#' generated, for each of the four PMID_39402379 contrasts.
#'
#' @param contrast Character vector of one or more of \code{"ad_cps"},
#'   \code{"early_ad_cps"}, \code{"late_ad_cps"}, \code{"versus_all"}.
#' @param force_recompute Logical scalar. If \code{TRUE}, ignore an existing
#'   cache file, recompute TRADE results, and overwrite the cache. Defaults to
#'   \code{FALSE}.
#' @param data_cache_dir Directory used to store cached TRADE results. If
#'   \code{NULL}, the directory is resolved from
#'   \code{options("bican.mccarroll.figures.cache_dir")}.
#' @param outDir Output directory used to write the SVG files. If \code{NULL},
#'   the directory is resolved from
#'   \code{options("bican.mccarroll.figures.out_dir")}.
#'
#' @return Invisibly returns \code{NULL}. This function is called for its side
#'   effects.
#'
#' @export
plot_trade_analysis_PMID_39402379 <- function(
  contrast = c("ad_cps", "early_ad_cps", "late_ad_cps", "versus_all"),
  force_recompute = FALSE,
  data_cache_dir = NULL,
  outDir = NULL
) {
  de_dir <- "differential_expression/external_comparison_PMID_39402379/voom-like"
  ct_file <- "differential_expression/metadata/cell_types_for_sea_ad_mtg_mmc_plots.txt"

  for (ct in contrast) {
    .plot_trade_analysis(
      de_dir = de_dir,
      cache_name = sprintf("trade_PMID_39402379_%s_autosomes.tsv", ct),
      dataset_id = "PMID_39402379",
      gene_to_chr_path = NULL,
      ct_file = ct_file,
      data_cache_dir = data_cache_dir,
      outDir = outDir,
      hide_cell_type_axis = FALSE,
      bar_fill = "black",
      x_label = "Transcriptome-wide impact (TRADE)",
      width = 5,
      height = 9,
      contrast = sprintf("__MTG__%s_DE_results\\.txt$", ct),
      contrast_id = ct,
      force_recompute = force_recompute
    )
  }

  invisible(NULL)
}


#' Plot TRADE results as a barplot across regions (expects input already filtered)
#'
#' @param trade_results A data.table produced by run_trade(), already filtered to a single cell type.
#' @param regions_use Optional character vector specifying the order and subset of regions to display.
#' @param region_var Character scalar naming the region column.
#' @param value_var Character scalar naming the TRADE statistic to plot.
#'
#' @return A ggplot object.
#' @export
trade_barplot_regions <- function(trade_results,
                                  regions_use = NULL,
                                  region_var = "region",
                                  value_var = "trade_twi") {
  # Make R CMD CHECK Happy
  region <- value <- N <- .N <- NULL

  dt <- data.table::as.data.table(trade_results)

  if (!(region_var %in% names(dt))) {
    stop("trade_barplot_regions(): region_var not found in trade_results.")
  }
  if (!(value_var %in% names(dt))) {
    stop("trade_barplot_regions(): value_var not found in trade_results.")
  }

  # Sanity: require exactly one cell type in input
  if ("cell_type" %in% names(dt)) {
    n_ct <- data.table::uniqueN(dt[["cell_type"]])
    if (n_ct != 1L) {
      stop("trade_barplot_regions(): requires trade_results to contain exactly one cell_type. Filter before calling.")
    }
  }

  # Require one row per region
  dup_rg <- dt[, .N, by = region_var][N > 1L]
  if (nrow(dup_rg) > 0L) {
    stop("trade_barplot_regions(): requires one row per region. Filter/aggregate first.")
  }

  dt <- dt[, list(region = get(region_var), value = get(value_var))]

  if (!is.null(regions_use)) {
    dt <- dt[region %in% regions_use]
    dt[, region := factor(region, levels = rev(regions_use))]
  } else {
    dt <- dt[order(value)]
    dt[, region := factor(region, levels = region)]
  }

  ggplot2::ggplot(dt, ggplot2::aes(x = value, y = region)) +
    ggplot2::geom_col() +
    ggplot2::labs(x = paste0("TRADE ", value_var), y = NULL) +
    ggplot2::theme_classic()
}


.plot_trade_analysis <- function(de_dir,
                                 cache_name,
                                 dataset_id,
                                 gene_to_chr_path = NULL,
                                 ct_file = NULL,
                                 data_cache_dir = NULL,
                                 outDir = NULL,
                                 hide_cell_type_axis = FALSE,
                                 x_breaks = NULL,
                                 x_labels = NULL,
                                 axis_title_x_size = NULL,
                                 axis_text_x_size = NULL,
                                 axis_tick_x_linewidth = NULL,
                                 plot_margin = NULL,
                                 bar_fill = "black",
                                 x_label = "Transcriptome-wide impact (TRADE)",
                                 width = 5,
                                 height = 9,
                                 contrast = "age",
                                 contrast_id = NULL,
                                 force_recompute = FALSE) {
  if (is.null(contrast_id)) {
    contrast_id <- contrast
  }

  paths <- resolve_trade_paths(
    de_dir = de_dir,
    gene_to_chr_path = gene_to_chr_path,
    ct_file = ct_file,
    data_cache_dir = data_cache_dir,
    outDir = outDir
  )

  cell_types_use_order <- scan(paths$ct_file, what = character(), quiet = TRUE)
  cache_file <- file.path(paths$data_cache_dir, cache_name)

  if (!force_recompute && file.exists(cache_file)) {
    # The TRADE functions assume data.table rather than data.frame.
    trade_auto <- data.table::fread(cache_file)
  } else {
    de_dt <- bican.mccarroll.differentialexpression::load_trade_data(
      data_path = paths$de_dir,
      contrast = contrast,
      gene_to_chr_path = paths$gene_to_chr_path,
      cellTypeListFile = paths$ct_file,
      regions_use = NULL
    )

    trade_auto <- .run_trade_autosomes(de_dt)

    utils::write.table(
      trade_auto,
      file = cache_file,
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE
    )
  }

  p_bar_age <- bican.mccarroll.differentialexpression::trade_barplot(
    trade_auto,
    cell_types_use = cell_types_use_order,
    value_var = "trade_twi"
  )

  theme_args <- list()

  if (hide_cell_type_axis) {
    theme_args$axis.text.y <- ggplot2::element_blank()
    theme_args$axis.ticks.y <- ggplot2::element_blank()
  }
  if (!is.null(axis_title_x_size)) {
    theme_args$axis.title.x <- ggplot2::element_text(size = axis_title_x_size)
  }
  if (!is.null(axis_text_x_size)) {
    theme_args$axis.text.x <- ggplot2::element_text(size = axis_text_x_size)
  }
  if (!is.null(axis_tick_x_linewidth)) {
    theme_args$axis.ticks.x <- ggplot2::element_line(
      linewidth = axis_tick_x_linewidth
    )
  }
  if (!is.null(plot_margin)) {
    theme_args$plot.margin <- plot_margin
  }
  if (length(theme_args) > 0L) {
    p_bar_age <- p_bar_age + do.call(ggplot2::theme, theme_args)
  }

  p_bar_age <- p_bar_age +
    ggplot2::geom_col(fill = bar_fill) +
    ggplot2::xlab(x_label)

  if (!is.null(x_breaks)) {
    scale_args <- list(breaks = x_breaks)
    if (!is.null(x_labels)) {
      scale_args$labels <- x_labels
    }
    p_bar_age <- p_bar_age +
      do.call(ggplot2::scale_x_continuous, scale_args)
  }

  out_file <- paste0(
    "trade_", dataset_id, "_", contrast_id, "_autosomes_barplot.svg"
  )

  save_plot_svg(
    p_bar_age,
    out_file = out_file,
    out_dir = paths$outDir,
    width = width,
    height = height
  )

  logger::log_info("DONE plotting Trade")
  invisible(NULL)
}


.plot_trade_analysis_regions <- function(de_region_interaction_dir,
                                         cache_name,
                                         dataset_id,
                                         region_order = NULL,
                                         filter_ic_to_non_neurons = FALSE,
                                         gene_to_chr_path = NULL,
                                         ct_file = NULL,
                                         data_cache_dir = NULL,
                                         outDir = NULL,
                                         width = 5,
                                         height = 8,
                                         force_recompute = FALSE) {
  paths <- resolve_trade_paths(
    de_region_interaction_dir = de_region_interaction_dir,
    gene_to_chr_path = gene_to_chr_path,
    ct_file = ct_file,
    data_cache_dir = data_cache_dir,
    outDir = outDir
  )

  cell_types_use_order <- scan(paths$ct_file, what = character(), quiet = TRUE)
  cache_file <- file.path(paths$data_cache_dir, cache_name)

  if (!force_recompute && file.exists(cache_file)) {
    trade_auto <- data.table::fread(cache_file)
  } else {
    de_dt <- bican.mccarroll.differentialexpression::load_trade_data(
      data_path = paths$de_region_interaction_dir,
      contrast = "age",
      gene_to_chr_path = paths$gene_to_chr_path,
      cellTypeListFile = paths$ct_file,
      regions_use = NULL
    )

    if (filter_ic_to_non_neurons) {
      # Make R CMD CHECK Happy
      region <- cell_type <- NULL

      non_neuron_types <- c(
        "astrocyte", "OPC", "oligodendrocyte", "microglia"
      )

      de_dt <- de_dt[
        !(region == "ic" & !(cell_type %in% non_neuron_types))
      ]
    }

    trade_auto <- .run_trade_autosomes(de_dt)

    utils::write.table(
      trade_auto,
      file = cache_file,
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE
    )
  }

  p_heat_age <- bican.mccarroll.differentialexpression::trade_heatmap(
    trade_auto,
    cell_types_use = cell_types_use_order,
    region_order = region_order,
    value_var = "trade_twi",
    display_numbers = TRUE
  )

  out_file <- paste0(
    "trade_", dataset_id,
    "_age_region_interaction_autosomes_heatmap.svg"
  )

  save_plot_svg(
    p_heat_age,
    out_file = out_file,
    out_dir = paths$outDir,
    width = width,
    height = height
  )

  logger::log_info("DONE plotting regional Trade")
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# Legacy commented analysis blocks retained for reference.
# -----------------------------------------------------------------------------
# Dataset 2: region subset, age (AUTOSOMES ONLY)
#
# This block used de_region_subset_dir and is intentionally not re-enabled.
# The active per-region BICAN analysis now uses de_region_interaction_dir via
# .plot_trade_analysis_regions().
#
# region_order <- c("CaH", "Pu", "NAC", "ic", "DFC")
# cache_file <- file.path(paths$data_cache_dir,
#                         "trade_dataset2_age_subset_region_autosomes.tsv")
#
# if (!force_recompute && file.exists(cache_file)) {
#   trade_auto <- data.table::fread(cache_file)
# } else {
#   de_dt <- bican.mccarroll.differentialexpression::load_trade_data(
#     data_path = paths$de_region_subset_dir, contrast = "age",
#     gene_to_chr_path = paths$gene_to_chr_path,
#     cellTypeListFile = paths$ct_file, regions_use = NULL)
#
#   trade_auto <- .run_trade_autosomes(de_dt)
#
#   utils::write.table(trade_auto, file = cache_file, sep = "\t",
#                      row.names = FALSE, col.names = TRUE, quote = FALSE)
# }
#
# p_heat_age <- bican.mccarroll.differentialexpression::trade_heatmap(
#   trade_auto, cell_types_use = cell_types_use_order,
#   region_order = region_order, value_var = "trade_twi")
#
# save_plot_svg(p_heat_age,
#               out_file = "trade_dataset2_age_subset_region_autosomes_heatmap.svg",
#               out_dir = paths$outDir, width = 5, height = 8)


.run_trade_autosomes <- function(de_dt) {
  # Make R CMD CHECK Happy
  chr <- NULL

  de_auto <- de_dt[chr %in% 1:22]
  bican.mccarroll.differentialexpression::run_trade(de_auto)
}


save_plot_svg <- function(plot, out_file, out_dir = ".", width = 14, height = 7) {
  out_svg <- file.path(out_dir, out_file)

  svglite::svglite(file = out_svg, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)

  print(plot)

  invisible(out_svg)
}


resolve_trade_paths <- function(de_dir = NULL,
                                de_region_subset_dir = NULL,
                                de_region_interaction_dir = NULL,
                                gene_to_chr_path = NULL,
                                ct_file = NULL,
                                outDir = NULL,
                                data_cache_dir = NULL) {
  root <- .resolve_data_root_dir(NULL)

  rel <- list(
    de_dir =
      "differential_expression/results/LEVEL_6/sex_age/cell_type",
    de_region_subset_dir =
      "differential_expression/results/LEVEL_6/sex_age/cell_type_subset_region",
    de_region_interaction_dir =
      "differential_expression/results/LEVEL_6/sex_age/cell_type_region_interaction_absolute_effects",
    gene_to_chr_path =
      "metadata/gene_to_chromosome.txt",
    ct_file =
      "differential_expression/metadata/cell_types_for_de_filtering_plot.txt"
  )

  pick_in <- function(x, key) {
    if (is.null(x)) {
      return(file.path(root, rel[[key]]))
    }
    .resolve_under_root(root, x)
  }

  out <- .resolve_out_dir(outDir)
  cache <- .resolve_cache_dir(data_cache_dir)

  # If a cache was not set, use the differential_expression subdirectory.
  if (is.null(data_cache_dir)) {
    cache <- file.path(cache, "differential_expression")
  }

  .ensure_dir(out)
  .ensure_dir(cache)

  list(
    data_root_dir = root,
    de_dir = pick_in(de_dir, "de_dir"),
    de_region_subset_dir = pick_in(
      de_region_subset_dir, "de_region_subset_dir"
    ),
    de_region_interaction_dir = pick_in(
      de_region_interaction_dir, "de_region_interaction_dir"
    ),
    gene_to_chr_path = pick_in(gene_to_chr_path, "gene_to_chr_path"),
    ct_file = pick_in(ct_file, "ct_file"),
    outDir = out,
    data_cache_dir = cache
  )
}
