# slope_matrix_path="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3/slope_matrix_qval_0.01.tsv"
# pval_nominal_matrix_path="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3/pval_nominal_matrix_qval_0.01.tsv"
# pval_nominal_threshold_matrix_path="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3/pval_nominal_threshold_matrix_qval_0.01.tsv"
# egene_union_pairs_path="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3/egene_union_pairs_qval_0.01.tsv"
# region_cell_type_path="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/region_cell_type.tsv"
# output_path="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3/cell_type_pairwise_r_squared.tsv"
# bican.mccarroll.eqtl::get_cell_type_pairwise_cor_matrix(slope_matrix_path, pval_nominal_matrix_path, pval_nominal_threshold_matrix_path, egene_union_pairs_path, region_cell_type_path, output_path)


#' Compute pairwise Spearman correlation matrix of eQTL effect sizes across cell types
#'
#' For each pair of cell type / region groups, identifies eGene-variant pairs
#' that are nominally significant in at least one of the two groups (using
#' \code{pval_nominal < pval_nominal_threshold}), then computes the Spearman
#' correlation of slopes for those pairs.  Returns the R-squared matrix.
#'
#' @param slope_matrix_path Character scalar.  Path to the slope matrix TSV
#'   (output of \code{\link{get_slope_matrix}}).
#' @param pval_nominal_matrix_path Character scalar.  Path to the pval_nominal matrix TSV
#'   (output of \code{\link{get_pval_nominal_matrix}}).
#' @param pval_nominal_threshold_matrix_path Character scalar.  Path to the pval_nominal_threshold
#'   matrix TSV (output of \code{\link{get_pval_nominal_threshold_matrix}}).
#' @param egene_union_pairs_path Character scalar.  Path to the eGene union pairs TSV
#'   (output of \code{\link{get_egene_union_pairs}}).
#' @param region_cell_type_path Character scalar.  Path to a tab-delimited file
#'   with columns \code{cell_type} and \code{region}.
#' @param output_path Character scalar or \code{NULL}.  If non-NULL, the
#'   R-squared matrix is written to this path as a tab-delimited file.
#'
#' @return A matrix of pairwise Spearman R-squared values (cell types x cell types).
#'
#' @export
#' @importFrom data.table fread fwrite
#' @importFrom logger log_info
get_cell_type_pairwise_cor_matrix <- function(slope_matrix_path,
                                              pval_nominal_matrix_path,
                                              pval_nominal_threshold_matrix_path,
                                              egene_union_pairs_path,
                                              region_cell_type_path,
                                              output_path = NULL) {

    slope_dt <- data.table::fread(slope_matrix_path)
    pval_dt <- data.table::fread(pval_nominal_matrix_path)
    pval_threshold_dt <- data.table::fread(pval_nominal_threshold_matrix_path)
    egene_dt <- data.table::fread(egene_union_pairs_path, select = c("phenotype_id", "variant_id"))

    # Filter to eGene union pairs
    slope_dt <- merge(slope_dt, egene_dt, by = c("phenotype_id", "variant_id"))
    pval_dt <- merge(pval_dt, egene_dt, by = c("phenotype_id", "variant_id"))
    pval_threshold_dt <- merge(pval_threshold_dt, egene_dt, by = c("phenotype_id", "variant_id"))

    # Filter to cell type/region columns
    region_cell_type_dt <- data.table::fread(region_cell_type_path)
    ct_cols <- paste0(region_cell_type_dt$cell_type, "__", region_cell_type_dt$region)
    ct_cols <- intersect(ct_cols, names(slope_dt))

    slope_m <- as.matrix(slope_dt[, ct_cols, with = FALSE])
    pval_m <- as.matrix(pval_dt[, ct_cols, with = FALSE])
    pval_thresh_m <- as.matrix(pval_threshold_dt[, ct_cols, with = FALSE])

    n <- length(ct_cols)
    logger::log_info("Computing pairwise Spearman correlations for {n} cell type/regions")

    cor_matrix <- matrix(NA_real_, nrow = n, ncol = n,
                         dimnames = list(ct_cols, ct_cols))

    for (i in seq_len(n)) {
        for (j in seq_len(n)) {
            slope1 <- slope_m[, i]
            slope2 <- slope_m[, j]
            pval1 <- pval_m[, i]
            pval2 <- pval_m[, j]
            thresh1 <- pval_thresh_m[, i]
            thresh2 <- pval_thresh_m[, j]

            sig_idx <- which(
                (!is.na(thresh1) & !is.na(pval1) & pval1 < thresh1) |
                (!is.na(thresh2) & !is.na(pval2) & pval2 < thresh2)
            )

            if (length(sig_idx) == 0) {
                cor_matrix[i, j] <- 0
                next
            }

            cor_matrix[i, j] <- stats::cor(
                slope1[sig_idx], slope2[sig_idx],
                use = "pairwise.complete.obs", method = "spearman"
            )
        }
    }

    r_squared <- cor_matrix^2

    logger::log_info("Correlation matrix complete")

    if (!is.null(output_path)) {
        out_dt <- data.table::as.data.table(r_squared, keep.rownames = "cell_type")
        data.table::fwrite(out_dt, output_path, sep = "\t")
        logger::log_info("Written to: {output_path}")
    }

    return(r_squared)
}
