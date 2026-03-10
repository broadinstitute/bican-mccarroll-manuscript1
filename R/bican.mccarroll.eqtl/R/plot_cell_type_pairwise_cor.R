# r_squared_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/eqtl_analysis_pipeline_run_jim/cell_type_pairwise_r_squared_qval_0.01.tsv"
# output_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/eqtl_analysis_pipeline_run_jim/cell_type_cor_plot_qval_0.01.svg"
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
#'   heatmap is saved to this path as an SVG file.
#' @param width Numeric.  SVG width in pixels (divided by res for inches).  Default 3150.
#' @param height Numeric.  SVG height in pixels (divided by res for inches).  Default 3000.
#' @param res Numeric.  Scaling factor for width/height.  Default 300.
#' @param title Character scalar.  Plot title.
#'
#' @return The \code{ComplexHeatmap::Heatmap} object (invisibly).
#'
#' @export
#' @importFrom data.table fread
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar unit grid.text
#' @importFrom grDevices svg dev.off
#' @importFrom logger log_info
plot_cell_type_pairwise_cor <- function(r_squared_path,
                                        output_path = NULL,
                                        width = 10.25,
                                        height = 9,
                                        #title = "Pairwise R\u00B2 of eQTL effect sizes across cell types",
                                        title=NULL) {

    celltype_label_map <- c(
        "MSN_D2_matrix__CaH"      = "MSN D2 matrix (CaH)",
        "MSN_D1_matrix__CaH"      = "MSN D1 matrix (CaH)",
        "MSN_D1_striosome__CaH"   = "MSN D1 striosome (CaH)",
        "MSN_D2_striosome__CaH"   = "MSN D2 striosome (CaH)",
        "GABA_MGE_CAP__CaH"       = "MGE-derived GABAergic (CaH)",
        "GABA_MGE_DFC__DFC"        = "MGE-derived GABAergic (DFC)",
        "GABA_CGE_DFC__DFC"        = "CGE-derived GABAergic (DFC)",
        "glutamatergic_L23IT__DFC" = "Glutamatergic L2/3 IT (DFC)",
        "glutamatergic_L5IT__DFC"  = "Glutamatergic L5 IT (DFC)",
        "astrocyte__DFC"           = "Astrocyte (DFC)",
        "astrocyte__CaH"           = "Astrocyte (CaH)",
        "oligodendrocyte__DFC"     = "Oligodendrocyte (DFC)",
        "oligodendrocyte__CaH"     = "Oligodendrocyte (CaH)",
        "OPC__DFC"                 = "OPC (DFC)",
        "OPC__CaH"                 = "OPC (CaH)",
        "microglia__DFC"           = "Microglia (DFC)",
        "microglia__CaH"           = "Microglia (CaH)"
    )

    dt <- data.table::fread(r_squared_path)
    rn <- dt[["cell_type"]]
    dt[, cell_type := NULL]
    cor_matrix_r_squared <- as.matrix(dt)

    # Map row and column names to human-readable labels
    labels <- ifelse(rn %in% names(celltype_label_map),
                     celltype_label_map[rn], rn)
    rownames(cor_matrix_r_squared) <- labels
    colnames(cor_matrix_r_squared) <- labels

    logger::log_info("Building heatmap for {nrow(cor_matrix_r_squared)} cell types")

    col_fun <- circlize::colorRamp2(
        c(0, 0.2, 0.4, 0.6, 0.8, 1),
        c("white", "#fff7bc", "#fee391", "#fdae61", "#f46d43", "#d73027")
    )

    ht <- ComplexHeatmap::Heatmap(
        cor_matrix_r_squared,
        name = "R_sq",
        col = col_fun,
        width  = grid::unit(ncol(cor_matrix_r_squared) * 0.85, "cm"),
        height = grid::unit(nrow(cor_matrix_r_squared) * 0.85, "cm"),
        column_title = title,
        column_title_gp = grid::gpar(fontsize = 16),
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        row_names_side = "left",
        row_names_gp = grid::gpar(fontsize = 18),
        column_names_gp = grid::gpar(fontsize = 18),
        row_names_max_width = grid::unit(12, "cm"),
        column_names_max_height = grid::unit(12, "cm"),
        column_names_rot = 45,
        heatmap_legend_param = list(
            title = expression(R^2),
            at = c(0, 0.25, 0.5, 0.75, 1),
            title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
            labels_gp = grid::gpar(fontsize = 12),
            grid_width  = grid::unit(6, "mm"),
            grid_height = grid::unit(6, "mm")
        )
    )

    if (!is.null(output_path)) {
        grDevices::svg(output_path, width = width, height = height)
        ComplexHeatmap::draw(
            ht,
            heatmap_legend_side = "right",
            padding = grid::unit(c(10, 4, 4, 4), "mm")
        )
        grDevices::dev.off()
        logger::log_info("Saved to: {output_path}")
    }

    invisible(ht)
}
