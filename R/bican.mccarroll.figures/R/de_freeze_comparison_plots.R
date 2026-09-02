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

# Compares age DE effect sizes and FDR between CAP freeze 3 and CAP freeze 3.1,
# for a given region, for three representative cell types. Replaces the ad
# hoc, base-R version in adhoc_compare_DE_by_freeze.R with a ggplot2 figure
# wrapper following the style of
# bican.mccarroll.differentialexpression::plot_effect_comparison
# (age_vs_pmi_rqs.R).

.de_freeze_comparison_old_dir <-
  "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results/LEVEL_3/sex_age/cell_type_region_interaction_absolute_effects"

.de_freeze_comparison_new_dir <-
  "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/results/LEVEL_6/sex_age/cell_type_region_interaction_absolute_effects"

#' Compare CAP freeze 3 vs CAP freeze 3.1 age DE results for one region
#'
#' For each of D1 matrix SPNs, astrocytes, and microglia, reads the age
#' differential expression results for the given region from both CAP freeze
#' 3 and CAP freeze 3.1, and writes a two-panel SVG (effect size on the left,
#' FDR on the right) comparing the two freezes gene by gene. Also writes a
#' summary TSV with correlation and slope statistics per cell type.
#'
#' Cell types are matched between freezes by name: freeze 3 called the D1
#' matrix SPN population "MSN_D1_matrix", freeze 3.1 renamed it to
#' "SPN_D1_matrix"; the astrocyte and microglia names are unchanged. All
#' three have region-interaction results for every region (CaH, DFC, ic,
#' NAC, Pu) in both freezes.
#'
#' @param region Brain region to compare, e.g. "CaH" or "DFC". Must match a
#'   region token in the region-interaction DE filenames.
#' @param outDir Output directory for generated SVG plots. If \code{NULL},
#'   resolved via the configured output directory option. SVGs and the
#'   summary TSV are written to a \code{pmi_comparison} subdirectory of this
#'   directory.
#' @return Invisibly returns a data.table with one summary row per cell type.
#' @export
de_freeze_comparison_scatter_plots <- function(region = "CaH", outDir = NULL) {
  # Make R CMD CHECK happy
  cell_type <- interaction <- gene <- NULL

  paths <- .resolve_de_freeze_comparison_paths(outDir)

  cell_type_map <- data.frame(
    old = c("MSN_D1_matrix", "astrocyte", "microglia"),
    new = c("SPN_D1_matrix", "astrocyte", "microglia"),
    stringsAsFactors = FALSE
  )

  file_pattern <- sprintf("__%s__age_DE_results\\.txt$", region)

  old_de <- data.table::as.data.table(
    bican.mccarroll.differentialexpression::parse_de_inputs(
      .de_freeze_comparison_old_dir,
      file_pattern = file_pattern
    )
  )

  new_de <- data.table::as.data.table(
    bican.mccarroll.differentialexpression::parse_de_inputs(
      .de_freeze_comparison_new_dir,
      file_pattern = file_pattern
    )
  )

  summary_list <- lapply(seq_len(nrow(cell_type_map)), function(i) {
    .plot_one_freeze_comparison(
      old_de = old_de,
      new_de = new_de,
      old_cell_type = cell_type_map$old[i],
      new_cell_type = cell_type_map$new[i],
      region = region,
      out_dir = paths$dataset_out_dir
    )
  })

  summary_dt <- data.table::rbindlist(summary_list[!vapply(summary_list, is.null, logical(1))])

  if (nrow(summary_dt) == 0L) {
    stop(
      "None of the requested cell types have age DE results in region '", region, "' in both freezes.",
      call. = FALSE
    )
  }

  utils::write.table(
    summary_dt,
    file = file.path(paths$dataset_out_dir, sprintf("de_freeze_comparison_%s_summary.tsv", region)),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  logger::log_info("DONE plotting CAP freeze 3 vs 3.1 DE comparison for region {region}")

  invisible(summary_dt)
}

#' Build and save the two-panel freeze comparison figure for one cell type
#'
#' @param old_de CAP freeze 3 DE data.table from parse_de_inputs.
#' @param new_de CAP freeze 3.1 DE data.table from parse_de_inputs.
#' @param old_cell_type Cell type label as it appears in old_de.
#' @param new_cell_type Cell type label as it appears in new_de.
#' @param region Brain region being compared, used for the plot title and
#'   output filename.
#' @param out_dir Directory to write the SVG to.
#' @return A data.table summary row for this cell type, or \code{NULL} if
#'   either freeze has no rows for this cell type in this region (e.g. a
#'   striatal cell type in a cortical region).
#' @noRd
.plot_one_freeze_comparison <- function(old_de,
                                        new_de,
                                        old_cell_type,
                                        new_cell_type,
                                        region,
                                        out_dir) {
  # Make R CMD CHECK happy
  cell_type <- gene <- logFC <- adj.P.Val <- NULL

  o <- old_de[cell_type == old_cell_type, ]
  n <- new_de[cell_type == new_cell_type, ]

  if (nrow(o) == 0L || nrow(n) == 0L) {
    logger::log_warn(
      "SKIP {old_cell_type} vs {new_cell_type} in region {region}: ",
      "no rows in {if (nrow(o) == 0L) 'CAP 3' else 'CAP 3.1'}"
    )
    return(NULL)
  }

  missing_from_new <- setdiff(o$gene, n$gene)
  missing_from_old <- setdiff(n$gene, o$gene)

  logger::log_info(
    "{old_cell_type} vs {new_cell_type}: {length(missing_from_new)} genes CAP 3 only, ",
    "{length(missing_from_old)} genes CAP 3.1 only"
  )

  merged <- merge(
    o[, list(gene, old_logFC = logFC, old_fdr = adj.P.Val)],
    n[, list(gene, new_logFC = logFC, new_fdr = adj.P.Val)],
    by = "gene"
  )

  if (nrow(merged) == 0L) {
    logger::log_warn("SKIP {old_cell_type} vs {new_cell_type} in region {region}: no shared genes")
    return(NULL)
  }

  merged[, old_neglog10_fdr := -log10(pmax(old_fdr, .Machine$double.xmin))]
  merged[, new_neglog10_fdr := -log10(pmax(new_fdr, .Machine$double.xmin))]

  display_name <- .format_cell_type_label(new_cell_type)

  p_effect <- .plot_freeze_comparison_panel(
    dt = merged,
    x_col = "old_logFC",
    y_col = "new_logFC",
    x_lab = "Age effect (log2 FC)",
    y_lab = "Age effect [with PMI] (log2 FC)",
    title = "Effect size",
    ref_line = 0
  )

  p_fdr <- .plot_freeze_comparison_panel(
    dt = merged,
    x_col = "old_neglog10_fdr",
    y_col = "new_neglog10_fdr",
    x_lab = "-log10(FDR)",
    y_lab = "-log10(FDR) [with PMI]",
    title = "FDR",
    ref_line = -log10(0.05)
  )

  plot_grid <- cowplot::plot_grid(
    p_effect,
    p_fdr,
    ncol = 2,
    labels = c("A", "B")
  )

  page_title <- paste(display_name, region, sep = " - ")

  page <- cowplot::ggdraw() +
    cowplot::draw_plot(
      plot_grid,
      x = 0,
      y = 0,
      width = 1,
      height = 0.94
    ) +
    cowplot::draw_label(
      page_title,
      x = 0.5,
      y = 0.985,
      hjust = 0.5,
      vjust = 1,
      fontface = "bold",
      size = 16
    )

  fileStr <- paste0(
    "de_freeze_comparison_", region, "_", new_cell_type, ".svg"
  )

  save_plot_svg(
    page,
    out_file = fileStr,
    out_dir = out_dir,
    width = 12,
    height = 6
  )

  fit <- lm(new_logFC ~ old_logFC, data = merged)

  data.table::data.table(
    old_cell_type = old_cell_type,
    new_cell_type = new_cell_type,
    n_genes = nrow(merged),
    cor_logFC = cor(merged$old_logFC, merged$new_logFC, use = "complete.obs"),
    cor_neglog10_fdr = cor(merged$old_neglog10_fdr, merged$new_neglog10_fdr, use = "complete.obs"),
    slope_logFC = unname(coef(fit)[["old_logFC"]]),
    mean_abs_diff_logFC = mean(abs(merged$new_logFC - merged$old_logFC), na.rm = TRUE),
    max_abs_diff_logFC = max(abs(merged$new_logFC - merged$old_logFC), na.rm = TRUE)
  )
}

#' One CAP 3 vs CAP 3.1 scatter panel with a y = x reference line and lm fit
#'
#' Mirrors bican.mccarroll.differentialexpression::plot_effect_comparison
#' (age_vs_pmi_rqs.R), adapted for a fixed x/y column pair and a
#' caller-supplied reference line value (0 for effect size, -log10(0.05) for
#' FDR, since -log10(FDR) has no informative zero line).
#'
#' @noRd
.plot_freeze_comparison_panel <- function(dt, x_col, y_col, x_lab, y_lab, title, ref_line) {
  plot_data <- data.frame(
    x_value = dt[[x_col]],
    y_value = dt[[y_col]]
  )

  fit <- lm(y_value ~ x_value, data = plot_data)
  slope <- unname(coef(fit)[["x_value"]])
  correlation <- cor(plot_data$x_value, plot_data$y_value, method = "pearson")

  annotation <- sprintf(
    "Pearson r = %.3f\nSlope = %.3f\nGenes = %s",
    correlation,
    slope,
    format(nrow(plot_data), big.mark = ",")
  )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = x_value, y = y_value)
  ) +
    ggplot2::geom_hline(
      yintercept = ref_line,
      linewidth = 0.3,
      color = "grey75"
    ) +
    ggplot2::geom_vline(
      xintercept = ref_line,
      linewidth = 0.3,
      color = "grey75"
    ) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      color = "grey50",
      linewidth = 0.5
    ) +
    ggrastr::geom_point_rast(
      alpha = 0.5,
      size = 1,
      raster.dpi = 600
    ) +
    ggplot2::geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      linewidth = 0.7
    ) +
    ggplot2::annotate(
      "text",
      x = Inf,
      y = -Inf,
      label = annotation,
      hjust = 1.05,
      vjust = -0.4,
      size = 3.5
    ) +
    ggplot2::labs(
      title = title,
      x = x_lab,
      y = y_lab
    ) +
    cowplot::theme_cowplot() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.margin = ggplot2::margin(8, 12, 12, 8)
    )
}

#' Format a cell type key into a manuscript display label
#'
#' Display-label formatting only; the underlying cell_type value used to
#' subset/join data is never mutated.
#'
#' @noRd
.format_cell_type_label <- function(x) {
  x <- gsub("_", " ", x)
  paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
}

.resolve_de_freeze_comparison_paths <- function(outDir = NULL) {
  out <- .resolve_out_dir(outDir)
  .ensure_dir(out)

  dataset_out_dir <- file.path(out, "pmi_comparison")
  .ensure_dir(dataset_out_dir)

  list(
    outDir = out,
    dataset_out_dir = dataset_out_dir
  )
}
