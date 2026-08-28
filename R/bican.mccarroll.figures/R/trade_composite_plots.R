## R/trade_composite_plots.R
##
## Multi-panel TRADE figures that place this work's barplot next to barplots
## from published datasets, with rows aligned across panels by a layout file
## (row order/subset/color are shared; row *labels* are panel-specific since
## panels may come from different cell type ontologies).
##
## Each exported function below hard-codes the data sources for one figure
## (cache names, DE directories, layout file, titles) and delegates alignment
## and plotting to the private engine .plot_trade_composite(). This mirrors
## the pattern in R/trade_plots.R, whose private helpers (resolve_trade_paths,
## .run_trade_autosomes, save_plot_svg) are reused here without modification.


#' Generate the aligned TRADE composite figure: Burger / Ling / Fröhlich
#'
#' Three barplot panels (this work, Ling et al., Fröhlich et al.) sharing a
#' common row order, subset, and color per row, taken from the layout file
#' \code{ling_frohlich_trade_alignment.txt}. Each panel is labeled with its
#' own cell type names (the layout's per-source columns), since the three
#' datasets use different cell type ontologies. TRADE result tables are
#' loaded from cache when available; \code{force_recompute = TRUE} recomputes
#' and overwrites the caches used by \code{plot_trade_analysis_bican()},
#' \code{plot_trade_analysis_Ling()}, and
#' \code{plot_trade_analysis_PMID_39227716()}.
#'
#' @param force_recompute Logical scalar. If \code{TRUE}, ignore existing
#'   cache files, recompute TRADE results, and overwrite the caches. Defaults
#'   to \code{FALSE}.
#' @param data_cache_dir Directory used to store/read cached TRADE results. If
#'   \code{NULL}, the directory is resolved from
#'   \code{options("bican.mccarroll.figures.cache_dir")}.
#' @param outDir Output directory used to write the SVG file. If \code{NULL},
#'   the directory is resolved from
#'   \code{options("bican.mccarroll.figures.out_dir")}.
#'
#' @return Invisibly returns \code{NULL}. This function is called for its
#'   side effects.
#'
#' @seealso
#'   \code{\link{plot_trade_analysis_bican}},
#'   \code{\link{plot_trade_analysis_Ling}},
#'   \code{\link{plot_trade_analysis_PMID_39227716}}
#' @export
plot_trade_composite_ling_frohlich <- function(force_recompute = FALSE,
                                               data_cache_dir = NULL,
                                               outDir = NULL) {
  panels <- list(
    list(
      layout_col = "this_work",
      title = "Burger et al.",
      cache_name = "trade_BICAN_age_autosomes.tsv",
      de_dir = NULL,
      ct_file = NULL,
      gene_to_chr_path = NULL,
      contrast = "age",
      x_breaks = c(0, 0.001, 0.002),
      x_labels = c("0.000", "0.001", "0.002")
    ),
    list(
      layout_col = "ling",
      title = "Ling et al.",
      cache_name = "trade_Ling_age_autosomes.tsv",
      de_dir = "/broad/mccarroll/dropulation/analysis/SNAP200_controls/differential_expression/results/sex_age/cell_type",
      ct_file = NULL,
      gene_to_chr_path = NULL,
      contrast = "age",
      x_breaks = c(0, 0.004, 0.006),
      x_labels = c("0.000", "0.003", "0.006")
    ),
    list(
      layout_col = "frohlich",
      title = "Fröhlich et al.",
      cache_name = "trade_PMID_39227716_age_autosomes.tsv",
      de_dir = "differential_expression/external_comparison_PMID_39227716/voom-like",
      ct_file = NULL,
      gene_to_chr_path = NULL,
      contrast = "age",
      x_breaks = c(0, 0.001, 0.002),
      x_labels = c("0.00", "0.001", "0.002")
    )
  )

  .plot_trade_composite(
    layout_file = "differential_expression/metadata/ling_frohlich_trade_alignment.txt",
    panels = panels,
    out_file = "trade_composite_ling_frohlich.svg",
    width = 12,
    height = 6,
    data_cache_dir = data_cache_dir,
    outDir = outDir,
    force_recompute = force_recompute
  )
}


#' Generate the aligned TRADE composite figure: Burger / Gabitto
#'
#' Two barplot panels (this work, Gabitto et al., \code{ad_cps} contrast)
#' sharing a common row order, subset, and color per row, taken from the
#' layout file \code{gabitto_trade_alignment.txt}. Both panels use SEA-AD
#' supertype cell type names (the layout's \code{this_work} and
#' \code{gabitto} columns are identical ontologies), since this work's
#' \code{LEVEL_6_sea_ad_mtg_mmc} mapping was generated to match Gabitto et
#' al.'s cell types directly. TRADE result tables are loaded from cache when
#' available; \code{force_recompute = TRUE} recomputes and overwrites the
#' caches used by \code{plot_trade_analysis_bican_sea_ad_mtg_mmc()} and
#' \code{plot_trade_analysis_PMID_39402379()}.
#'
#' @param force_recompute Logical scalar. If \code{TRUE}, ignore existing
#'   cache files, recompute TRADE results, and overwrite the caches. Defaults
#'   to \code{FALSE}.
#' @param data_cache_dir Directory used to store/read cached TRADE results. If
#'   \code{NULL}, the directory is resolved from
#'   \code{options("bican.mccarroll.figures.cache_dir")}.
#' @param outDir Output directory used to write the SVG file. If \code{NULL},
#'   the directory is resolved from
#'   \code{options("bican.mccarroll.figures.out_dir")}.
#'
#' @return Invisibly returns \code{NULL}. This function is called for its
#'   side effects.
#'
#' @seealso
#'   \code{\link{plot_trade_analysis_bican_sea_ad_mtg_mmc}},
#'   \code{\link{plot_trade_analysis_PMID_39402379}}
#' @export
plot_trade_composite_gabitto <- function(force_recompute = FALSE,
                                         data_cache_dir = NULL,
                                         outDir = NULL) {
  ct_file <- "differential_expression/metadata/cell_types_for_sea_ad_mtg_mmc_plots.txt"

  panels <- list(
    list(
      layout_col = "this_work",
      title = "Burger et al.",
      cache_name = "trade_BICAN_sea_ad_mtg_mmc_age_autosomes.tsv",
      de_dir = "differential_expression/results/LEVEL_6_sea_ad_mtg_mmc/sex_age/cell_type",
      ct_file = ct_file,
      gene_to_chr_path = NULL,
      contrast = "__age_DE_results\\.txt$",
      x_breaks = NULL,
      x_labels = NULL
    ),
    list(
      layout_col = "gabitto",
      title = "Gabitto et al.",
      cache_name = "trade_PMID_39402379_ad_cps_autosomes.tsv",
      de_dir = "differential_expression/external_comparison_PMID_39402379/voom-like",
      ct_file = ct_file,
      gene_to_chr_path = NULL,
      contrast = "__MTG__ad_cps_DE_results\\.txt$",
      x_breaks = NULL,
      x_labels = NULL
    )
  )

  .plot_trade_composite(
    layout_file = "differential_expression/metadata/gabitto_trade_alignment.txt",
    panels = panels,
    out_file = "trade_composite_gabitto_ad_cps.svg",
    width = 8,
    height = 9,
    data_cache_dir = data_cache_dir,
    outDir = outDir,
    force_recompute = force_recompute
  )
}


#' Generate the aligned TRADE composite figure: Burger / Kana
#'
#' Two barplot panels (this work, Kana et al.) sharing a common row order,
#' subset, and color per row, taken from the layout file
#' \code{kana_trade_alignment.txt}: the six striatal cell types for which a
#' Kana-to-BICAN cell type name mapping was defined. The \code{this_work}
#' panel uses this work's age contrast; the \code{kana} panel uses the Kana
#' \code{ad_cps} contrast (there is no BICAN vs. Kana age comparison; see
#' \code{\link{plot_trade_analysis_kana_2026}}). TRADE result tables are
#' loaded from cache when available; \code{force_recompute = TRUE} recomputes
#' and overwrites the caches used by \code{plot_trade_analysis_bican()} and
#' \code{plot_trade_analysis_kana_2026()}.
#'
#' @param force_recompute Logical scalar. If \code{TRUE}, ignore existing
#'   cache files, recompute TRADE results, and overwrite the caches. Defaults
#'   to \code{FALSE}.
#' @param data_cache_dir Directory used to store/read cached TRADE results. If
#'   \code{NULL}, the directory is resolved from
#'   \code{options("bican.mccarroll.figures.cache_dir")}.
#' @param outDir Output directory used to write the SVG file. If \code{NULL},
#'   the directory is resolved from
#'   \code{options("bican.mccarroll.figures.out_dir")}.
#'
#' @return Invisibly returns \code{NULL}. This function is called for its
#'   side effects.
#'
#' @seealso
#'   \code{\link{plot_trade_analysis_bican}},
#'   \code{\link{plot_trade_analysis_kana_2026}}
#' @export
plot_trade_composite_kana <- function(force_recompute = FALSE,
                                      data_cache_dir = NULL,
                                      outDir = NULL) {
  panels <- list(
    list(
      layout_col = "this_work",
      title = "Burger et al.",
      cache_name = "trade_BICAN_age_autosomes.tsv",
      de_dir = NULL,
      ct_file = NULL,
      gene_to_chr_path = NULL,
      contrast = "age",
      x_breaks = c(0, 0.001, 0.002),
      x_labels = c("0.000", "0.001", "0.002")
    ),
    list(
      layout_col = "kana",
      title = "Kana et al.",
      cache_name = "trade_KANA_2026_ad_cps_autosomes.tsv",
      de_dir = "differential_expression/external_comparison_kana_2026_biorxiv/voom-like",
      ct_file = "differential_expression/metadata/cell_types_for_kana_2026_plots.txt",
      gene_to_chr_path = NULL,
      contrast = "__CaH__ad_cps_DE_results\\.txt$",
      x_breaks = NULL,
      x_labels = NULL
    )
  )

  .plot_trade_composite(
    layout_file = "differential_expression/metadata/kana_trade_alignment.txt",
    panels = panels,
    out_file = "trade_composite_kana.svg",
    width = 8,
    height = 4,
    data_cache_dir = data_cache_dir,
    outDir = outDir,
    force_recompute = force_recompute
  )
}


## -----------------------------------------------------------------------
## Private engine
## -----------------------------------------------------------------------

## Align one or more TRADE data sources to a shared row layout (order,
## subset, and per-row color) and render them as a single multi-panel SVG.
##
## layout_file: path (relative to the data root, or absolute) to a
##   tab-delimited file with one row per aligned cell type, a "color"
##   column, and one column per panel (see `panels`). "NA" cells mean that
##   panel has no counterpart for that row.
## panels: list of panel specs, each a list with:
##   layout_col       column name in the layout file for this panel
##   title             panel title
##   cache_name        TRADE cache file name (read if present, else computed)
##   de_dir, ct_file, gene_to_chr_path  passed to resolve_trade_paths()/
##                     load_trade_data() when the cache must be (re)computed
##   contrast          contrast regex passed to load_trade_data()
##   x_breaks, x_labels  optional; NULL lets ggplot2 pick breaks
.plot_trade_composite <- function(layout_file,
                                  panels,
                                  out_file,
                                  width,
                                  height,
                                  x_label = "Transcriptome-wide impact (TRADE)",
                                  rel_widths = NULL,
                                  data_cache_dir = NULL,
                                  outDir = NULL,
                                  force_recompute = FALSE) {
  root <- .resolve_data_root_dir(NULL)
  layout_path <- .resolve_under_root(root, layout_file)

  layout <- data.table::fread(
    layout_path,
    colClasses = "character",
    na.strings = c("NA", "")
  )

  if (!("color" %in% names(layout))) {
    stop(
      ".plot_trade_composite(): layout file must contain a 'color' column: ",
      layout_path
    )
  }

  layout_cols <- vapply(panels, function(spec) spec$layout_col, character(1))
  missing_cols <- setdiff(layout_cols, names(layout))
  if (length(missing_cols) > 0L) {
    stop(
      ".plot_trade_composite(): layout_col(s) not found in layout file '",
      layout_path, "': ", paste(missing_cols, collapse = ", ")
    )
  }

  row_id <- sprintf("r%03d", seq_len(nrow(layout)))
  row_levels <- rev(row_id)
  color_vec <- stats::setNames(layout[["color"]], row_id)

  panel_plots <- lapply(panels, function(spec) {
    .trade_composite_panel(
      layout = layout,
      row_id = row_id,
      row_levels = row_levels,
      color_vec = color_vec,
      spec = spec,
      x_label = x_label,
      data_cache_dir = data_cache_dir,
      outDir = outDir,
      force_recompute = force_recompute
    )
  })

  plot_grid_args <- list(
    plotlist = panel_plots,
    nrow = 1,
    align = "h",
    axis = "tb"
  )
  if (!is.null(rel_widths)) {
    plot_grid_args$rel_widths <- rel_widths
  }
  composite <- do.call(cowplot::plot_grid, plot_grid_args)

  top_paths <- resolve_trade_paths(data_cache_dir = data_cache_dir, outDir = outDir)

  save_plot_svg(
    composite,
    out_file = out_file,
    out_dir = top_paths$outDir,
    width = width,
    height = height
  )

  logger::log_info("DONE plotting composite TRADE figure: {out_file}")
  invisible(NULL)
}


## Build one panel: load/compute the panel's TRADE table, align it to the
## shared row layout, and return a ggplot barplot.
.trade_composite_panel <- function(layout,
                                   row_id,
                                   row_levels,
                                   color_vec,
                                   spec,
                                   x_label,
                                   data_cache_dir,
                                   outDir,
                                   force_recompute) {
  paths <- resolve_trade_paths(
    de_dir = spec$de_dir,
    gene_to_chr_path = spec$gene_to_chr_path,
    ct_file = spec$ct_file,
    data_cache_dir = data_cache_dir,
    outDir = outDir
  )

  cache_file <- file.path(paths$data_cache_dir, spec$cache_name)

  trade_auto <- .load_trade_cached(
    cache_file = cache_file,
    de_dir = paths$de_dir,
    gene_to_chr_path = paths$gene_to_chr_path,
    ct_file = paths$ct_file,
    contrast = spec$contrast,
    force_recompute = force_recompute
  )

  dt <- data.table::as.data.table(trade_auto)

  dup_ct <- dt[, .N, by = "cell_type"][N > 1L]
  if (nrow(dup_ct) > 0L) {
    stop(
      ".trade_composite_panel(): source '", spec$title,
      "' has more than one row for cell type(s): ",
      paste(dup_ct$cell_type, collapse = ", ")
    )
  }

  value_by_ct <- stats::setNames(dt$trade_twi, dt$cell_type)
  layout_names <- layout[[spec$layout_col]]

  unmatched <- layout_names[!is.na(layout_names) & !(layout_names %in% names(value_by_ct))]
  if (length(unmatched) > 0L) {
    logger::log_warn(
      "{spec$title}: layout name(s) not found in source data: {paste(unmatched, collapse = ', ')}"
    )
  }

  unused <- setdiff(names(value_by_ct), layout_names[!is.na(layout_names)])
  if (length(unused) > 0L) {
    logger::log_warn(
      "{spec$title}: source cell type(s) not used by the layout: {paste(unused, collapse = ', ')}"
    )
  }

  plot_dt <- data.table::data.table(
    row_id = row_id,
    label = ifelse(is.na(layout_names), "", layout_names),
    value = as.numeric(value_by_ct[layout_names])
  )
  plot_dt$row_id <- factor(plot_dt$row_id, levels = row_levels)

  dt_plot <- plot_dt[!is.na(value)]

  p <- ggplot2::ggplot(
    dt_plot,
    ggplot2::aes(x = value, y = row_id, fill = row_id)
  ) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = color_vec, guide = "none") +
    ggplot2::scale_y_discrete(
      limits = row_levels,
      breaks = row_id,
      labels = plot_dt$label,
      drop = FALSE
    ) +
    ggplot2::labs(title = spec$title, x = x_label, y = NULL) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0),
      axis.title.x = ggplot2::element_text(size = ggplot2::rel(0.8))
    )

  if (!is.null(spec$x_breaks)) {
    xlim_range <- range(c(spec$x_breaks, dt_plot$value), na.rm = TRUE)
    xlim_range[2] <- xlim_range[2] + diff(xlim_range) * 0.08

    scale_args <- list(breaks = spec$x_breaks)
    if (!is.null(spec$x_labels)) {
      scale_args$labels <- spec$x_labels
    }

    p <- p +
      do.call(ggplot2::scale_x_continuous, scale_args) +
      ggplot2::coord_cartesian(xlim = xlim_range)
  }

  p
}


## Load (or compute and cache) an autosomal TRADE result table. Local
## counterpart of the cache-or-compute block in .plot_trade_analysis()
## (R/trade_plots.R); duplicated rather than shared so this file makes no
## changes to trade_plots.R.
.load_trade_cached <- function(cache_file,
                               de_dir,
                               gene_to_chr_path,
                               ct_file,
                               contrast,
                               force_recompute = FALSE) {
  if (!force_recompute && file.exists(cache_file)) {
    return(data.table::fread(cache_file))
  }

  de_dt <- bican.mccarroll.differentialexpression::load_trade_data(
    data_path = de_dir,
    contrast = contrast,
    gene_to_chr_path = gene_to_chr_path,
    cellTypeListFile = ct_file,
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

  trade_auto
}
