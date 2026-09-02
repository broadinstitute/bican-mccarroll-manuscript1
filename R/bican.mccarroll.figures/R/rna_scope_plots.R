# source("R/paths.R")
#
# options(
#     bican.mccarroll.figures.data_root_dir =
#         "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis",
#
#     bican.mccarroll.figures.out_dir =
#         "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/figure_repository"
# )
#
# plot_rna_scope_validation()

#' Plot snRNA-seq vs RNAscope cell type proportion validation
#'
#' Cell type proportions (CTP) are normally estimated from single-nucleus
#' RNA-seq: run the single cell data, transfer cell type labels, count
#' per-donor cells of each type, and convert to a percentage of the donor's
#' total. That approach is scalable but indirect. RNAscope counts cells
#' directly and is far less scalable, so it is used here to validate the
#' snRNA-seq CTP estimates.
#'
#' Reads \code{RNAScope/rna_scope_results.txt} under the configured data
#' root directory (see \code{\link{get_data_root_dir}}), builds one scatter
#' plot per region/cell-type pair present in the file, and assembles them
#' left to right in a single row, in the order CaH OPCs, Pu OPCs, CaH
#' Astrocytes, with shared, outer x/y axis titles rather than per-panel
#' titles. Each panel uses matching x- and y-axis scaling and
#' is annotated with the (non-parametric) Spearman correlation and its
#' p-value, following the same convention as
#' \code{bican.mccarroll.celltypeproportions::plot_metric_correlation()}.
#' Writes the composite as \code{RNAScope_results.svg} to the configured
#' figure output directory (see \code{\link{get_out_dir}}).
#'
#' The input table is a small, hand-curated 12-donor spreadsheet, so it is
#' read directly on every call; no caching is performed.
#'
#' @return Invisibly returns a list with elements \code{data} (the tidy
#'   long-format data frame used for plotting), \code{plots} (the three
#'   individual ggplot objects, named by panel), and \code{composite} (the
#'   assembled cowplot grid).
#'
#' @export
plot_rna_scope_validation <- function() {
  paths <- .resolve_rna_scope_paths()
  panel_spec <- .rna_scope_panel_spec()

  plot_data <- .read_rna_scope_results(paths$rna_scope_file, panel_spec)

  panels <- lapply(seq_len(nrow(panel_spec)), function(i) {
    .build_rna_scope_scatter(
      plot_data[plot_data$panel == panel_spec$panel[i], ],
      title = panel_spec$title[i]
    )
  })
  names(panels) <- panel_spec$panel

  grid <- cowplot::plot_grid(
    plotlist = panels,
    nrow = 1,
    align = "hv"
  )

  x_title <- cowplot::ggdraw() + cowplot::draw_label("snRNA-seq (% of cells)", size = 14)
  y_title <- cowplot::ggdraw() + cowplot::draw_label("RNAscope (% of cells)", angle = 90, size = 14)

  grid_with_x <- cowplot::plot_grid(grid, x_title, ncol = 1, rel_heights = c(1, 0.06))
  composite <- cowplot::plot_grid(y_title, grid_with_x, nrow = 1, rel_widths = c(0.06, 1))

  save_plot_svg(
    composite,
    out_file = "RNAScope_results.svg",
    out_dir = paths$outDir,
    width = 12,
    height = 4
  )

  logger::log_info("DONE plotting RNAscope validation results")

  invisible(list(data = plot_data, plots = panels, composite = composite))
}

## -------------------------------------------------------------------------
## Internal helpers
## -------------------------------------------------------------------------

#' Panel specification for the RNAscope validation figure
#'
#' One row per region/cell-type pair: a short panel key, a display title,
#' and the snRNA-seq / RNAscope column names as they appear in the source
#' spreadsheet.
#'
#' @noRd
.rna_scope_panel_spec <- function() {
  data.frame(
    panel = c("CaH_OPC", "Pu_OPC", "CaH_Astrocyte"),
    title = c("CaH OPCs", "Pu OPCs", "CaH Astrocytes"),
    snrnaseq_col = c(
      "CaH OPCs (Village)",
      "Pu OPCs (Village)",
      "CaH Astrocytes (Village)"
    ),
    rnascope_col = c(
      "CaH OPCs (RNAscope)",
      "Pu OPCs (RNAscope)",
      "CaH Astrocytes (RNAscope)"
    ),
    stringsAsFactors = FALSE
  )
}

#' Read and tidy the RNAscope validation results table
#'
#' The source file interleaves each measurement pair with a blank spacer
#' column, and at least one RNAscope column contains the literal string
#' "pending" for donors not yet scored, so columns are selected by name
#' (never position) and coerced to numeric, treating unparseable values
#' (including "pending") as missing.
#'
#' @param file Path to \code{rna_scope_results.txt}.
#' @param panel_spec Data frame as returned by \code{.rna_scope_panel_spec()}.
#' @return A tidy data frame with columns \code{donor}, \code{panel},
#'   \code{title}, \code{snrnaseq}, \code{rnascope}; rows missing either
#'   measurement are dropped.
#' @noRd
.read_rna_scope_results <- function(file, panel_spec = .rna_scope_panel_spec()) {
  raw <- as.data.frame(
    data.table::fread(file, sep = "\t", header = TRUE, check.names = FALSE, na.strings = "NA")
  )

  missing_cols <- setdiff(c(panel_spec$snrnaseq_col, panel_spec$rnascope_col), names(raw))
  if (length(missing_cols) > 0) {
    stop(
      "RNAscope results file is missing expected column(s): ",
      paste(missing_cols, collapse = ", "),
      ". File: ", file
    )
  }

  rows <- lapply(seq_len(nrow(panel_spec)), function(i) {
    data.frame(
      donor = raw[["Donor"]],
      panel = panel_spec$panel[i],
      title = panel_spec$title[i],
      snrnaseq = suppressWarnings(as.numeric(raw[[panel_spec$snrnaseq_col[i]]])),
      rnascope = suppressWarnings(as.numeric(raw[[panel_spec$rnascope_col[i]]])),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  out[stats::complete.cases(out[, c("snrnaseq", "rnascope")]), ]
}

#' Build one snRNA-seq vs RNAscope scatter panel
#'
#' Uses matching x- and y-axis limits (the combined range of both
#' measurements) so the y=x reference line reads consistently across
#' panels, and annotates the non-parametric Spearman correlation and its
#' p-value only (no Pearson/parametric statistic), matching the convention
#' in \code{bican.mccarroll.celltypeproportions::plot_metric_correlation()}.
#'
#' @param panel_df Data frame subset for a single panel, with columns
#'   \code{snrnaseq} and \code{rnascope}.
#' @param title Panel title.
#' @return A ggplot object.
#' @noRd
.build_rna_scope_scatter <- function(panel_df, title) {
  axis_limits <- range(c(panel_df$snrnaseq, panel_df$rnascope))

  cor_test <- stats::cor.test(
    panel_df$snrnaseq, panel_df$rnascope,
    method = "spearman", exact = FALSE
  )
  annotation <- paste(
    sprintf("ρ = %+.2f", unname(cor_test$estimate)),
    paste0("p = ", format.pval(cor_test$p.value, digits = 2)),
    sep = "\n"
  )

  ggplot2::ggplot(panel_df, ggplot2::aes(x = snrnaseq, y = rnascope)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    ggplot2::geom_point(size = 2) +
    ggplot2::xlim(axis_limits[1], axis_limits[2]) +
    ggplot2::ylim(axis_limits[1], axis_limits[2]) +
    ggplot2::annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = annotation,
      hjust = -0.15,
      vjust = 1.3,
      size = 3.5
    ) +
    ggplot2::labs(title = title) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.title = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(8, 12, 12, 8)
    )
}

#' Resolve input/output paths for the RNAscope validation figure
#' @noRd
.resolve_rna_scope_paths <- function() {
  root <- .resolve_data_root_dir(NULL)
  out <- .resolve_out_dir(NULL)
  .ensure_dir(out)

  list(
    data_root_dir = root,
    rna_scope_file = file.path(root, "RNAScope/rna_scope_results.txt"),
    outDir = out
  )
}
