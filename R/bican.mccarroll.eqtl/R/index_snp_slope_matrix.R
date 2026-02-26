# slope_matrix_path="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3/slope_matrix_qval_0.01.tsv"
# output_path="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3/index_snp_slope_matrix_with_median_impute_qval_0.01.tsv"
# bican.mccarroll.eqtl::get_index_snp_slope_matrix_with_median_impute(slope_matrix_path, output_path=output_path)


#' Select index SNPs per gene and build a sign-adjusted, median-imputed slope matrix
#'
#' Starting from the full slope matrix (eGene-variant pairs x cell type/region),
#' this function:
#' \enumerate{
#'   \item Selects one index variant per gene — the variant with the largest
#'         absolute slope across all cell types.
#'   \item Adjusts signs so the cell type with the largest absolute slope is
#'         always positive (ensures consistent orientation across genes).
#'   \item Filters to genes with at least \code{min_non_na} non-NA slope values.
#'   \item Imputes remaining NA values with the row (gene) median.
#' }
#'
#' @param slope_matrix_path Character scalar.  Path to the slope matrix file
#'   (output of \code{\link{get_slope_matrix}}), with columns
#'   \code{phenotype_id}, \code{variant_id}, and one slope column per
#'   cell type / region.
#' @param min_non_na Integer scalar.  Minimum number of non-NA slope values
#'   required to retain a gene.  Default \code{3}.
#' @param output_path Character scalar or \code{NULL}.  If non-NULL, the result
#'   table is written to this path as a tab-delimited file.
#'
#' @return A \code{data.table} with columns \code{phenotype_id},
#'   \code{variant_id}, and one sign-adjusted, imputed slope column per
#'   cell type / region.
#'
#' @export
#' @importFrom data.table fread fwrite as.data.table
#' @importFrom logger log_info
get_index_snp_slope_matrix_with_median_impute <- function(slope_matrix_path,
                                                          min_non_na = 3,
                                                          output_path = NULL) {

    slope_dt <- data.table::fread(slope_matrix_path)
    slope_cols <- setdiff(names(slope_dt), c("phenotype_id", "variant_id"))

    logger::log_info("Input: {nrow(slope_dt)} eGene-variant pairs x {length(slope_cols)} cell type/regions")

    # --- 1. Select index SNP per gene (variant with max abs slope across all cell types) ---
    slope_m <- as.matrix(slope_dt[, slope_cols, with = FALSE])
    max_abs_slope <- apply(abs(slope_m), 1, max, na.rm = TRUE)
    slope_dt[, max_abs_slope := max_abs_slope]

    # For each gene, keep the variant with the largest max_abs_slope
    slope_dt[, rank := frank(-max_abs_slope, ties.method = "first"), by = "phenotype_id"]
    index_dt <- slope_dt[rank == 1]
    index_dt[, c("max_abs_slope", "rank") := NULL]

    logger::log_info("Selected {nrow(index_dt)} index SNPs (one per gene)")

    # --- 2. Sign-adjust so largest absolute value is always positive ---
    index_m <- as.matrix(index_dt[, slope_cols, with = FALSE])

    for (i in seq_len(nrow(index_m))) {
        row <- index_m[i, ]
        max_idx <- which.max(abs(row))
        index_m[i, ] <- row * sign(row[max_idx])
    }

    # --- 3. Filter to genes with enough non-NA values ---
    non_na_counts <- rowSums(!is.na(index_m))
    keep <- non_na_counts >= min_non_na
    index_m <- index_m[keep, , drop = FALSE]
    index_ids <- index_dt[keep, .(phenotype_id, variant_id)]

    logger::log_info("Retained {nrow(index_m)} genes with >= {min_non_na} non-NA values (dropped {sum(!keep)})")

    # --- 4. Median impute remaining NAs ---
    for (i in seq_len(nrow(index_m))) {
        row <- index_m[i, ]
        na_mask <- is.na(row)
        if (any(na_mask)) {
            index_m[i, na_mask] <- median(row, na.rm = TRUE)
        }
    }

    result_dt <- cbind(index_ids, data.table::as.data.table(index_m))

    logger::log_info("Output: {nrow(result_dt)} genes x {length(slope_cols)} cell type/regions (fully imputed)")

    if (!is.null(output_path)) {
        data.table::fwrite(result_dt, output_path, sep = "\t")
        logger::log_info("Written to: {output_path}")
    }

    return(result_dt)
}
