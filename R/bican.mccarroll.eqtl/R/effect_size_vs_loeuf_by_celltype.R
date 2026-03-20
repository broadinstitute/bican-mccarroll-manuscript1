#' Effect size vs LOEUF stratified by cell type
#'
#' For each cell type, reads the tensorQTL \code{cis_qtl} output, filters to
#' eGenes (qval < 0.01), restricts to protein-coding genes with gnomAD LOEUF
#' scores, and computes the Spearman correlation between |slope| and LOEUF.
#' Produces a faceted scatterplot and a TSV of per-cell-type results.
#'
#' @param tensorqtl_dir Character scalar.  Path to the tensorQTL results
#'   directory (e.g. \code{.../LEVEL_3}).  Each cell type subdirectory should
#'   contain \code{<ct>.cis_qtl.txt.gz}.
#' @param gnomad_path Character scalar.  Path to gnomAD constraint metrics TSV.
#' @param gene_biotype_path Character scalar.  Path to gene biotype TSV (output
#'   of \code{\link{get_gene_biotype_from_gtf}}) with columns
#'   \code{phenotype_id} and \code{gene_type}.
#' @param plot_output_path Character scalar.  Path for the output SVG plot.
#' @param cell_types Character vector.  The cell type identifiers to analyse.
#'   Defaults to the 17 cell types used in the manuscript.
#' @param qval_threshold Numeric.  FDR threshold for defining eGenes.
#'   Default 0.01.
#' @param plot_width Numeric.  Plot width in inches.  Default 13.4.
#' @param plot_height Numeric or \code{NULL}.  Plot height in inches.  If
#'   \code{NULL}, computed automatically from the number of cell types.
#' @param ncol Integer.  Number of facet columns.  Default 5.
#'
#' @return A \code{data.table} with per-cell-type Spearman results.  Called
#'   for side effects (writes TSV and SVG plot).
#'
#' @export
#' @importFrom data.table fread fwrite data.table rbindlist setnames
#' @importFrom stats cor.test p.adjust
#' @importFrom logger log_info
effect_size_vs_loeuf_by_celltype <- function(
    tensorqtl_dir,
    gnomad_path,
    gene_biotype_path,
    plot_output_path,
    cell_types = c(
        "MSN_D2_matrix__CaH",
        "MSN_D1_matrix__CaH",
        "MSN_D1_striosome__CaH",
        "MSN_D2_striosome__CaH",
        "GABA_MGE_CAP__CaH",
        "GABA_MGE_DFC__DFC",
        "GABA_CGE_DFC__DFC",
        "glutamatergic_L23IT__DFC",
        "glutamatergic_L5IT__DFC",
        "astrocyte__DFC",
        "astrocyte__CaH",
        "oligodendrocyte__DFC",
        "oligodendrocyte__CaH",
        "OPC__DFC",
        "OPC__CaH",
        "microglia__DFC",
        "microglia__CaH"
    ),
    qval_threshold = 0.01,
    plot_width = 13.4,
    plot_height = NULL,
    ncol = 5
) {

    # ---- 1. Filter to protein-coding genes ----
    logger::log_info("Reading gene biotype (filtering to protein-coding)...")
    map_dt <- data.table::fread(gene_biotype_path)
    pc_genes <- map_dt[gene_type == "protein_coding"]$phenotype_id
    logger::log_info("{length(pc_genes)} protein-coding genes")

    # ---- 2. Read gnomAD LOEUF ----
    logger::log_info("Reading gnomAD constraint metrics...")
    gnomad_dt <- data.table::fread(gnomad_path)
    loeuf_dt <- gnomad_dt[canonical == TRUE & mane_select == TRUE &
                           !is.na(gene) & !is.na(lof.oe_ci.upper),
                          .(loeuf = unique(lof.oe_ci.upper)), by = gene]
    data.table::setnames(loeuf_dt, "gene", "phenotype_id")
    logger::log_info("{nrow(loeuf_dt)} genes with LOEUF (MANE Select canonical)")

    # ---- 3. Per-cell-type analysis ----
    logger::log_info("Running per-cell-type Spearman tests...")

    results <- list()
    plot_data <- list()

    for (ct in cell_types) {
        fpath <- file.path(tensorqtl_dir, ct, paste0(ct, ".cis_qtl.txt.gz"))
        logger::log_info("  {ct} ...")

        cis_dt <- data.table::fread(fpath, select = c("phenotype_id", "slope", "qval"))

        # Filter to eGenes
        egenes <- cis_dt[qval < qval_threshold]
        n_egenes <- nrow(egenes)
        egenes[, abs_slope := abs(slope)]

        # Filter to protein-coding
        egenes_pc <- egenes[phenotype_id %in% pc_genes]
        n_pc <- nrow(egenes_pc)

        # Filter to those with LOEUF
        egenes <- merge(egenes_pc, loeuf_dt, by = "phenotype_id")
        n_sig <- nrow(egenes)

        logger::log_info("    eGenes: {n_egenes}, protein-coding: {n_pc}/{n_egenes}, with LOEUF: {n_sig}/{n_pc}")

        if (n_sig >= 10) {
            sp <- stats::cor.test(egenes$loeuf, egenes$abs_slope,
                                  method = "spearman", exact = FALSE)
            spearman_rho <- round(sp$estimate, 4)
            spearman_p   <- sp$p.value
        } else {
            spearman_rho <- NA
            spearman_p   <- NA
        }

        results[[length(results) + 1]] <- data.table::data.table(
            cell_type    = ct,
            n_egenes     = n_egenes,
            n_pc_loeuf   = n_sig,
            spearman_rho = spearman_rho,
            spearman_p   = spearman_p
        )

        if (n_sig >= 10) {
            egenes[, cell_type := ct]
            plot_data[[length(plot_data) + 1]] <- egenes[, .(phenotype_id, abs_slope, loeuf, cell_type)]
        }
    }

    res_dt <- data.table::rbindlist(results)
    res_dt[, spearman_p_BH := stats::p.adjust(spearman_p, method = "BH")]

    logger::log_info("Per-cell-type results:\n{paste(capture.output(print(res_dt)), collapse = '\n')}")

    # ---- 4. Scatterplot faceted by cell type ----
    logger::log_info("Generating scatterplot...")

    plot_dt <- data.table::rbindlist(plot_data)

    celltype_label_map <- c(
        "GABA_CGE_DFC__DFC"      = "CGE-derived GABAergic (DFC)",
        "GABA_MGE_CAP__CaH"      = "MGE-derived GABAergic (CaH)",
        "GABA_MGE_DFC__DFC"      = "MGE-derived GABAergic (DFC)",
        "MSN_D1_matrix__CaH"     = "MSN D1 matrix (CaH)",
        "MSN_D1_striosome__CaH"  = "MSN D1 striosome (CaH)",
        "MSN_D2_matrix__CaH"     = "MSN D2 matrix (CaH)",
        "MSN_D2_striosome__CaH"  = "MSN D2 striosome (CaH)",
        "OPC__CaH"               = "OPC (CaH)",
        "OPC__DFC"               = "OPC (DFC)",
        "astrocyte__CaH"         = "Astrocyte (CaH)",
        "astrocyte__DFC"         = "Astrocyte (DFC)",
        "glutamatergic_L23IT__DFC" = "Glutamatergic L2/3 IT (DFC)",
        "glutamatergic_L5IT__DFC"  = "Glutamatergic L5 IT (DFC)",
        "microglia__CaH"         = "Microglia (CaH)",
        "microglia__DFC"         = "Microglia (DFC)",
        "oligodendrocyte__CaH"   = "Oligodendrocyte (CaH)",
        "oligodendrocyte__DFC"   = "Oligodendrocyte (DFC)"
    )
    celltype_label_levels <- celltype_label_map[cell_types]

    plot_dt[, cell_type_level := factor(celltype_label_map[as.character(cell_type)],
                                        levels = celltype_label_levels)]

    anno_dt <- res_dt[!is.na(spearman_rho), .(cell_type, spearman_rho, spearman_p_BH)]
    anno_dt[, label := sprintf("rho == %.3f * ',' ~ italic(p) == '%.2g'",
                               spearman_rho, spearman_p_BH)]
    anno_dt[, cell_type_level := factor(celltype_label_map[as.character(cell_type)],
                                        levels = celltype_label_levels)]

    n_ct <- length(unique(plot_dt$cell_type))
    nrow_plot <- ceiling(n_ct / ncol)
    if (is.null(plot_height)) {
        plot_height <- 2.2 * nrow_plot + 1
    }

    p <- ggplot2::ggplot(plot_dt, ggplot2::aes(x = loeuf, y = abs_slope)) +
        ggplot2::geom_point(size = 0.5, alpha = 0.3, color = "black") +
        ggplot2::geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.6) +
        ggplot2::geom_text(data = anno_dt, ggplot2::aes(label = label),
                  x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3,
                  size = 4, color = "red", parse = TRUE) +
        ggplot2::facet_wrap(~ cell_type_level, ncol = ncol, scales = "fixed") +
        ggplot2::labs(x = "LOEUF", y = "|eQTL effect size|") +
        ggplot2::theme_classic(base_size = 11, base_family = "Arial") +
        ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", size = 8))

    ggplot2::ggsave(plot_output_path, plot = p, width = plot_width,
                    height = plot_height, device = "svg")
    logger::log_info("Plot saved to: {plot_output_path}")

    invisible(res_dt)
}
