# median_expression_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/manuscript_data/test/heatmap_index_snp_median_expression_qval_0.01.tsv"
# cluster_assignments_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/manuscript_data/test/cluster_assignments_qval_0.01_k13.tsv"
# output_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/manuscript_data/test/cluster_expression_fraction_in_expected_celltype.tsv"
# bican.mccarroll.eqtl::get_cluster_expression_fraction(median_expression_path, cluster_assignments_path, output_path)


#' Compute fraction of genes per cluster whose highest expression is in the expected cell type
#'
#' For each cluster with a known expected cell-type family, computes the
#' fraction of genes whose argmax expression (highest median expression
#' across all cell types) falls within the expected cell-type family.
#'
#' @param median_expression_path Character scalar.  Path to the median
#'   expression matrix TSV (output of
#'   \code{\link{get_median_expression_matrix}}).  First column should be
#'   \code{Gene}; remaining columns are cell type / region groups.
#' @param cluster_assignments_path Character scalar.  Path to the cluster
#'   assignments TSV with columns \code{gene} and \code{cluster}.
#' @param output_path Character scalar or \code{NULL}.  If non-NULL, the
#'   result is written to this path as a tab-delimited file.
#' @param cluster_to_expected A named list mapping cluster IDs (as
#'   character) to character vectors of expected cell type column names.
#'   If \code{NULL}, uses the default mapping for K=13 cell-type-specific
#'   clusters.
#'
#' @return A \code{data.table} with one row per cluster containing the
#'   fraction of genes whose highest expression is in the expected cell type.
#'
#' @export
#' @importFrom data.table fread fwrite data.table rbindlist setnames
#' @importFrom logger log_info
get_cluster_expression_fraction <- function(median_expression_path,
                                             cluster_assignments_path,
                                             output_path = NULL,
                                             cluster_to_expected = NULL) {

    # Default mapping: cluster ID -> expected cell types
    if (is.null(cluster_to_expected)) {
        cluster_to_expected <- list(
            "3"  = c("astrocyte__DFC", "astrocyte__CaH"),
            "7"  = c("oligodendrocyte__DFC", "oligodendrocyte__CaH"),
            "10"  = c("OPC__DFC", "OPC__CaH"),
            "8"  = c("microglia__DFC", "microglia__CaH"),
            "12"  = c("MSN_D1_matrix__CaH", "MSN_D2_matrix__CaH",
                      "MSN_D1_striosome__CaH", "MSN_D2_striosome__CaH"),
            "9"  = c("MSN_D1_matrix__CaH", "MSN_D2_matrix__CaH",
                      "MSN_D1_striosome__CaH", "MSN_D2_striosome__CaH"),
            "6" = c("glutamatergic_L23IT__DFC", "glutamatergic_L5IT__DFC"),
            "1" = c("GABA_MGE_CAP__CaH", "GABA_MGE_DFC__DFC", "GABA_CGE_DFC__DFC")
        )
    }

    # 1. Read data
    cluster_dt <- data.table::fread(cluster_assignments_path)
    if ("gene" %in% names(cluster_dt)) data.table::setnames(cluster_dt, "gene", "Gene")
    cluster_dt[, cluster := as.character(cluster)]
    logger::log_info("Cluster assignments: {nrow(cluster_dt)} genes")

    expr_dt <- data.table::fread(median_expression_path)
    expr_genes <- expr_dt[[1]]
    expr_m <- as.matrix(expr_dt[, -1, with = FALSE])
    rownames(expr_m) <- expr_genes
    logger::log_info("Median expression: {nrow(expr_m)} genes x {ncol(expr_m)} cell types")

    # 2. Compute fraction per cluster
    results <- list()

    for (cl in names(cluster_to_expected)) {
        expected_cts <- cluster_to_expected[[cl]]
        genes <- cluster_dt[cluster == cl, Gene]
        genes <- genes[genes %in% expr_genes]

        expected_cts_present <- expected_cts[expected_cts %in% colnames(expr_m)]

        if (length(genes) == 0 || length(expected_cts_present) == 0) {
            results[[length(results) + 1]] <- data.table::data.table(
                cluster = as.integer(cl),
                expected_celltypes = paste(expected_cts, collapse = ","),
                n_genes = length(genes),
                n_nonzero = 0L,
                n_max_in_expected = 0L,
                fraction_max_in_expected = NA_real_
            )
            next
        }

        sub_m <- expr_m[genes, , drop = FALSE]

        # For each gene, find the cell type with highest expression
        max_ct <- colnames(sub_m)[apply(sub_m, 1, which.max)]

        # Exclude genes with all-zero expression (which.max returns 1 arbitrarily)
        gene_max <- apply(sub_m, 1, max, na.rm = TRUE)
        valid <- gene_max > 0
        n_valid <- sum(valid)

        # Count genes whose argmax cell type is in the expected subset
        n_hits <- sum(max_ct[valid] %in% expected_cts_present)
        frac <- n_hits / length(genes)

        results[[length(results) + 1]] <- data.table::data.table(
            cluster = as.integer(cl),
            expected_celltypes = paste(expected_cts_present, collapse = ","),
            n_genes = length(genes),
            n_nonzero = n_valid,
            n_max_in_expected = n_hits,
            fraction_max_in_expected = round(frac, 4)
        )
    }

    result_dt <- data.table::rbindlist(results)
    result_dt <- result_dt[order(cluster)]
    logger::log_info("Results:\n{paste(capture.output(print(result_dt)), collapse = '\n')}")

    if (!is.null(output_path)) {
        data.table::fwrite(result_dt, output_path, sep = "\t")
        logger::log_info("Written to: {output_path}")
    }

    invisible(result_dt)
}
