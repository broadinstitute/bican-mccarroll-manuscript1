# slope_matrix_path="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3/slope_matrix_qval_0.01.tsv"
# output_path="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3/index_snp_slope_matrix_with_zero_impute_qval_0.01.tsv"
# bican.mccarroll.eqtl::get_index_snp_slope_matrix_with_impute(slope_matrix_path, output_path=output_path)


#' Select index SNPs per gene and build a sign-adjusted, zero-imputed slope matrix
#'
#' Starting from the full slope matrix (eGene-variant pairs x cell type/region),
#' this function:
#' \enumerate{
#'   \item Selects one index variant per gene - the variant with the largest
#'         absolute slope across all cell types.
#'   \item Adjusts signs so the cell type with the largest absolute slope is
#'         always positive (ensures consistent orientation across genes).
#'   \item Imputes missing values with zero (NA -> 0), reflecting the
#'         assumption that untested cell types have no detectable effect.
#' }
#'
#' @param slope_matrix_path Character scalar. Path to the slope matrix file.
#' @param output_path Character scalar or \code{NULL}. If non-NULL, the result
#'   table is written to this path as a tab-delimited file.
#'
#' @return A \code{data.table} with identifier columns and one sign-adjusted,
#'   zero-imputed slope column per cell type / region.
#'
#' @export
#' @importFrom data.table fread fwrite as.data.table frank
#' @importFrom logger log_info
get_index_snp_slope_matrix_with_impute <- function(slope_matrix_path,
                                                   output_path = NULL) {

    slope_dt <- data.table::fread(slope_matrix_path)
    slope_cols <- setdiff(names(slope_dt),
                          c("phenotype_id", "variant_id"))

    logger::log_info(
        "Input: {nrow(slope_dt)} eGene-variant pairs x {length(slope_cols)} cell type/regions"
    )

    # --- 1. Select index SNP per gene ---
    slope_m <- as.matrix(slope_dt[, slope_cols, with = FALSE])
    max_abs_slope <- apply(abs(slope_m), 1, max, na.rm = TRUE)
    slope_dt[, max_abs_slope := max_abs_slope]

    slope_dt[, rank := data.table::frank(-max_abs_slope,
                                         ties.method = "first"),
             by = "phenotype_id"]

    index_dt <- slope_dt[rank == 1]
    index_dt[, c("max_abs_slope", "rank") := NULL]

    logger::log_info("Selected {nrow(index_dt)} index SNPs (one per gene)")

    # --- 2. Sign-adjust ---
    index_m <- as.matrix(index_dt[, slope_cols, with = FALSE])

    for (i in seq_len(nrow(index_m))) {
        row <- index_m[i, ]
        max_idx <- which.max(abs(row))
        index_m[i, ] <- row * sign(row[max_idx])
    }

    index_dt[, (slope_cols) := data.table::as.data.table(index_m)]

    # --- 3. Zero-impute missing values ---
    m <- as.matrix(index_dt[, slope_cols, with = FALSE])
    na_mask <- is.na(m)
    if (any(na_mask)) {
        m[na_mask] <- 0
    }
    index_dt[, (slope_cols) := data.table::as.data.table(m)]

    logger::log_info(
        "Output: {nrow(index_dt)} genes x {length(slope_cols)} cell type/regions (zero-imputed)"
    )

    # Sort by gene for reproducible row order (k-means is order-sensitive)
    #Make R CMD CHECK Happy
    phenotype_id <- NULL
    data.table::setorder(index_dt, phenotype_id)

    if (!is.null(output_path)) {
        data.table::fwrite(index_dt, output_path, sep = "\t")
        logger::log_info("Written to: {output_path}")
    }

    index_dt
}
