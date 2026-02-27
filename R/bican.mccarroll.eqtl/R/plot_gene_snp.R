# Example usage:
#   devtools::load_all()
#   plot_gene_snp(
#       gene = "CEP112", chr = "chr17", pos = 66192315,
#       vcf_path = "/broad/bican_um1_mccarroll/vcfs/2025-05-05/gvs_concat_outputs_2025-05-05T14-10-02.donors_renamed_filtered_norm.vcf.gz",
#       expression_path = "/broad/.../combined_tpm_expression_across_cell_types.tsv",
#       output_path = "/broad/.../manuscript_data/CEP112_chr17_66192315.png"
#   )


#' Plot gene expression by SNP genotype across cell types
#'
#' Reads a VCF at a single genomic position, extracts genotypes, joins with
#' per-donor gene expression, computes per-cell-type linear-model p-values,
#' and draws a violin + beeswarm plot faceted by cell type (all panels in a
#' single row).
#'
#' By default 10 cell types are plotted.  When \code{gene == "NPAS3"},
#' \code{microglia__DFC} is appended for 11 panels.
#'
#' @param gene Character scalar.  Gene symbol (e.g. \code{"CEP112"}).
#' @param chr Character scalar.  Chromosome (e.g. \code{"chr17"}).
#' @param pos Integer scalar.  Genomic position (1-based).
#' @param vcf_path Character scalar.  Path to the indexed VCF/BCF file.
#' @param expression_path Character scalar.  Path to the combined gene
#'   expression TPM TSV (output of
#'   \code{\link{combine_expression_across_cell_types}}).  First column is
#'   \code{pid} (gene ID); remaining columns are
#'   \code{<donor>_<celltype>__<region>}.
#' @param output_path Character scalar or \code{NULL}.  If non-NULL, the
#'   plot is saved to this path as a PNG file.
#' @param width Numeric.  PNG width in pixels.  Default 2600.
#' @param height Numeric.  PNG height in pixels.  Default 1000.
#' @param res Numeric.  PNG resolution in DPI.  Default 150.
#'
#' @return The \code{ggplot} object (invisibly).
#'
#' @export
#' @importFrom VariantAnnotation readVcf geno
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom data.table fread melt tstrsplit :=
#' @importFrom ggplot2 ggplot aes geom_violin geom_text facet_wrap labs
#'   theme_bw theme element_text element_blank unit
#' @importFrom ggbeeswarm geom_beeswarm
#' @importFrom grDevices png dev.off
#' @importFrom stats lm sd p.adjust
#' @importFrom logger log_info
plot_gene_snp <- function(gene,
                          chr,
                          pos,
                          vcf_path,
                          expression_path,
                          output_path = NULL,
                          width = 2600,
                          height = 1000,
                          res = 150) {

    celltype_order <- c(
        "MSN_D1_matrix__CaH",
        "MSN_D2_matrix__CaH",
        "GABA_MGE_DFC__DFC",
        "GABA_CGE_DFC__DFC",
        "GABA_MGE_CAP__CaH",
        "glutamatergic_L23IT__DFC",
        "astrocyte__CaH",
        "oligodendrocyte__CaH",
        "OPC__CaH",
        "microglia__CaH"
    )

    if (gene == "NPAS3") {
        celltype_order <- c(celltype_order, "microglia__DFC")
    }

    celltype_label_map <- c(
        "MSN_D1_matrix__CaH"        = "MSN D1 matrix (CaH)",
        "MSN_D2_matrix__CaH"        = "MSN D2 matrix (CaH)",
        "GABA_MGE_DFC__DFC"          = "MGE-derived GABAergic (DFC)",
        "GABA_CGE_DFC__DFC"          = "CGE-derived GABAergic (DFC)",
        "GABA_MGE_CAP__CaH"          = "MGE-derived GABAergic (CaH)",
        "glutamatergic_L23IT__DFC"   = "Glutamatergic L2/3 IT (DFC)",
        "astrocyte__CaH"             = "Astrocyte (CaH)",
        "oligodendrocyte__CaH"       = "Oligodendrocyte (CaH)",
        "OPC__CaH"                   = "OPC (CaH)",
        "microglia__CaH"             = "Microglia (CaH)",
        "microglia__DFC"             = "Microglia (DFC)"
    )

    # --- Read VCF at position ---
    gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(pos, pos))
    vcf_snp <- VariantAnnotation::readVcf(vcf_path, "hg38", param = gr)

    if (length(vcf_snp) > 1L) {
        warning("More than one variant found at ", chr, ":", pos)
    }

    genotype_matrix <- VariantAnnotation::geno(vcf_snp)$GT
    genotype_matrix <- gsub("\\|", "/", genotype_matrix)
    genotype_matrix <- gsub("^1/0$", "0/1", genotype_matrix)

    # Extract ref / alt alleles from variant ID (e.g. "chr17:66192315_A/G")
    variant_id <- rownames(genotype_matrix)[1]
    allele_str <- strsplit(variant_id, "_")[[1]][2]
    allele_parts <- strsplit(allele_str, "/")[[1]]
    ref_allele <- allele_parts[1]
    alt_allele <- allele_parts[2]

    genotype_v <- genotype_matrix[1, ]
    names(genotype_v) <- colnames(genotype_matrix)

    # --- Read expression & reshape to long format ---
    logger::log_info("Reading expression from: {expression_path}")
    expr_dt <- data.table::fread(expression_path)
    gene_dt <- expr_dt[expr_dt$pid == gene]
    gene_dt[, pid := NULL]

    long_dt <- data.table::melt(
        gene_dt, measure.vars = names(gene_dt),
        variable.name = "donor_cell_type", value.name = "expression"
    )

    # Split "donorID_celltype__region" into donor and cell_type
    long_dt[, c("donor", "cell_type") := data.table::tstrsplit(
        donor_cell_type, "_", fixed = TRUE, keep = c(1L, NA)
    )]
    # tstrsplit with keep=c(1, NA) only gets the first piece; reconstruct
    # the cell_type by removing the donor prefix
    long_dt[, cell_type := sub("^[^_]+_", "", donor_cell_type)]
    long_dt[, donor_cell_type := NULL]

    # --- Join genotypes ---
    long_dt[, genotype := genotype_v[donor]]

    # Remove unascertained genotypes, compute log expression and dosage
    long_dt <- long_dt[genotype != "./."]
    long_dt[, log_expression := log2(expression + 1)]
    long_dt[, dosage := fifelse(
        genotype == "0/0", 0,
        fifelse(genotype %in% c("0/1", "1/0"), 1,
                fifelse(genotype == "1/1", 2, NA_real_))
    )]

    # --- Compute per-cell-type p-values (linear model) ---
    pval_dt <- long_dt[
        !is.na(log_expression) & !is.na(dosage),
        {
            n_unique <- length(unique(dosage))
            if (n_unique < 2L || .N < 3L || sd(dosage) == 0) {
                list(pval = NA_real_, beta = NA_real_)
            } else {
                model <- lm(log_expression ~ dosage)
                coefs <- summary(model)$coefficients
                if (nrow(coefs) >= 2L) {
                    list(pval = coefs[2, 4], beta = coefs[2, 1])
                } else {
                    list(pval = NA_real_, beta = NA_real_)
                }
            }
        },
        by = cell_type
    ]

    pval_dt[, pval_adj := stats::p.adjust(pval, method = "BH")]
    pval_dt[, p_label := paste0("p = ", signif(pval_adj, digits = 3))]

    # --- Merge p-values back and prepare plot data ---
    plot_dt <- merge(long_dt, pval_dt[, .(cell_type, pval_adj, p_label)],
                     by = "cell_type", all.x = TRUE)

    # Map cell type labels
    plot_dt[, cell_type_label := celltype_label_map[cell_type]]

    # Filter to desired cell types, non-zero expression
    plot_dt <- plot_dt[cell_type %in% celltype_order & log_expression != 0]
    plot_dt[, cell_type := factor(cell_type, levels = celltype_order)]
    plot_dt[, cell_type_label := factor(
        cell_type_label, levels = celltype_label_map[celltype_order]
    )]

    # Allele-based genotype labels
    plot_dt[, allele_genotype := fifelse(
        genotype == "0/0", paste0(ref_allele, "/", ref_allele),
        fifelse(genotype == "0/1", paste0(ref_allele, "/", alt_allele),
                paste0(alt_allele, "/", alt_allele))
    )]

    # --- P-value label positions ---
    global_y_pos <- max(plot_dt$log_expression, na.rm = TRUE) + 0.3

    label_dt <- plot_dt[, .(pval_adj = pval_adj[1], p_label = p_label[1]),
                        by = cell_type_label]
    label_dt[, y_pos := global_y_pos]
    label_dt <- label_dt[!is.na(pval_adj) & pval_adj < 0.05]

    # --- Build plot ---
    n_panels <- length(celltype_order)

    p <- ggplot2::ggplot(plot_dt,
                         ggplot2::aes(x = allele_genotype,
                                      y = log_expression)) +
        ggbeeswarm::geom_beeswarm(size = 0.3, alpha = 0.8) +
        ggplot2::geom_violin(
            trim = FALSE, scale = "width",
            fill = "#9bd3dd", color = NA, alpha = 0.6
        ) +
        ggplot2::facet_wrap(
            ~ cell_type_label, ncol = n_panels, drop = FALSE,
            labeller = ggplot2::label_wrap_gen(width = 18)
        ) +
        ggplot2::geom_text(
            data = label_dt,
            ggplot2::aes(y = y_pos, label = p_label),
            color = "#d62728",
            x = Inf, hjust = 1, vjust = 1,
            nudge_x = -0.15, nudge_y = -0.05,
            size = 6, fontface = "bold", inherit.aes = FALSE
        ) +
        ggplot2::labs(
            title = paste0("Expression of ", gene, " by Genotype at ",
                           chr, ":", pos),
            x = "Genotype",
            y = paste0("Log2(", gene, " expression)")
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(
                hjust = 0.5, size = 24, face = "bold"),
            strip.text = ggplot2::element_text(size = 14, face = "bold"),
            axis.title.x = ggplot2::element_text(size = 20),
            axis.title.y = ggplot2::element_text(size = 20),
            panel.spacing.x = ggplot2::unit(0.18, "lines"),
            panel.spacing.y = ggplot2::unit(0.25, "lines"),
            strip.background = ggplot2::element_blank()
        )

    if (!is.null(output_path)) {
        grDevices::png(output_path, width = width, height = height, res = res)
        print(p)
        grDevices::dev.off()
        logger::log_info("Saved to: {output_path}")
    }

    invisible(p)
}
