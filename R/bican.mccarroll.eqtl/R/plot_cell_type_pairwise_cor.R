# r_squared_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3/cell_type_pairwise_r_squared.tsv"
# output_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3/cell_type_cor_plot.png"
# bican.mccarroll.eqtl::plot_cell_type_pairwise_cor(r_squared_path, output_path)


#' Plot heatmap of pairwise R-squared of eQTL effect sizes across cell types
#'
#' Reads the pairwise R-squared matrix produced by
#' \code{\link{get_cell_type_pairwise_cor_matrix}} and draws a clustered
#' heatmap using \pkg{ComplexHeatmap}.
#'
#' @param r_squared_path Character scalar.  Path to the R-squared matrix TSV
#'   (output of \code{\link{get_cell_type_pairwise_cor_matrix}}).  The first
#'   column should be named \code{cell_type} and contain row names.
#' @param output_path Character scalar or \code{NULL}.  If non-NULL, the
#'   heatmap is saved to this path as a PNG file.
#' @param width Numeric.  PNG width in pixels.  Default 3400.
#' @param height Numeric.  PNG height in pixels.  Default 3000.
#' @param res Numeric.  PNG resolution in DPI.  Default 300.
#' @param title Character scalar.  Plot title.
#'
#' @return The \code{ComplexHeatmap::Heatmap} object (invisibly).
#'
#' @export
#' @importFrom data.table fread
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar unit grid.text
#' @importFrom grDevices png dev.off
#' @importFrom logger log_info
plot_cell_type_pairwise_cor <- function(r_squared_path,
                                        output_path = NULL,
                                        width = 3400,
                                        height = 3000,
                                        res = 300,
                                        title = "Pairwise R\u00B2 of eQTL effect sizes across cell types") {

    dt <- data.table::fread(r_squared_path)
    rn <- dt[["cell_type"]]
    dt[, cell_type := NULL]
    cor_matrix_r_squared <- as.matrix(dt)
    rownames(cor_matrix_r_squared) <- rn

    logger::log_info("Building heatmap for {nrow(cor_matrix_r_squared)} cell types")

    col_fun <- circlize::colorRamp2(
        c(0, 0.2, 0.4, 0.6, 0.8, 1),
        c("white", "#fff7bc", "#fee391", "#fdae61", "#f46d43", "#d73027")
    )

    ht <- ComplexHeatmap::Heatmap(
        cor_matrix_r_squared,
        name = "R_squared",
        col = col_fun,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        row_names_side = "left",
        row_names_gp = grid::gpar(fontsize = 12, fontface = "bold"),
        column_names_gp = grid::gpar(fontsize = 12, fontface = "bold"),
        column_names_rot = 45,
        heatmap_legend_param = list(
            at = c(0, 0.25, 0.5, 0.75, 1),
            title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
            labels_gp = grid::gpar(fontsize = 12),
            grid_width  = grid::unit(6, "mm"),
            grid_height = grid::unit(6, "mm")
        )
    )

    if (!is.null(output_path)) {
        grDevices::png(output_path, width = width, height = height, res = res)
        ComplexHeatmap::draw(
            ht,
            heatmap_legend_side = "right",
            padding = grid::unit(c(20, 4, 4, 4), "mm")
        )
        grid::grid.text(
            title,
            x = 0.5,
            y = grid::unit(1, "npc") - grid::unit(0.8, "mm"),
            gp = grid::gpar(fontsize = 16, fontface = "bold")
        )
        grDevices::dev.off()
        logger::log_info("Saved to: {output_path}")
    }

    invisible(ht)
}
