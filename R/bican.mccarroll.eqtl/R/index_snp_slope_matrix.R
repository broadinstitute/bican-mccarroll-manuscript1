# slope_matrix_path="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3/slope_matrix_qval_0.01.tsv"
# output_path="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3/index_snp_slope_matrix_with_median_impute_qval_0.01.tsv"
# bican.mccarroll.eqtl::get_index_snp_slope_matrix_with_median_impute(slope_matrix_path, output_path=output_path)


#' Select index SNPs per gene and build a sign-adjusted, imputed slope matrix
#'
#' Starting from the full slope matrix (eGene-variant pairs x cell type/region),
#' this function:
#' \enumerate{
#'   \item Selects one index variant per gene - the variant with the largest
#'         absolute slope across all cell types.
#'   \item Adjusts signs so the cell type with the largest absolute slope is
#'         always positive (ensures consistent orientation across genes).
#'   \item Imputes missing values using \code{imputation_method}. When
#'         \code{"median"} is used, genes with fewer than \code{min_non_na}
#'         non-NA values are dropped.
#' }
#'
#' @param slope_matrix_path Character scalar. Path to the slope matrix file.
#' @param min_non_na Integer scalar. Minimum number of non-NA slope values
#'   required to retain a gene when \code{imputation_method = "median"}.
#'   Default \code{3}.
#' @param imputation_method Character scalar. One of \code{"median"} or
#'   \code{"zero"}. Default \code{"median"}.
#' @param output_path Character scalar or \code{NULL}. If non-NULL, the result
#'   table is written to this path as a tab-delimited file.
#'
#' @return A \code{data.table} with identifier columns and one sign-adjusted,
#'   imputed slope column per cell type / region.
#'
#' @export
#' @importFrom data.table fread fwrite as.data.table frank
#' @importFrom logger log_info
get_index_snp_slope_matrix_with_impute <- function(slope_matrix_path,
                                                   min_non_na = 3,
                                                   imputation_method = "median",
                                                   output_path = NULL) {

    if (!imputation_method %in% c("median", "zero")) {
        stop("imputation_method must be one of: 'median', 'zero'.")
    }

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

    # --- 3. Impute (fully encapsulated) ---
    index_dt <- .impute_slope_dt(
        dt = index_dt,
        slope_cols = slope_cols,
        method = imputation_method,
        min_non_na = min_non_na
    )

    logger::log_info(
        "Output: {nrow(index_dt)} genes x {length(slope_cols)} cell type/regions (imputed via '{imputation_method}')"
    )

    # Sort by gene for reproducible row order (k-means is order-sensitive)
    data.table::setorder(index_dt, phenotype_id)

    if (!is.null(output_path)) {
        data.table::fwrite(index_dt, output_path, sep = "\t")
        logger::log_info("Written to: {output_path}")
    }

    index_dt
}

#' Impute slope columns in a data.table
#'
#' Applies row-wise imputation to slope columns in a \code{data.table}.
#'
#' \itemize{
#'   \item \code{"median"}: rows with fewer than \code{min_non_na}
#'         non-NA slope values are dropped. Remaining \code{NA}s are
#'         replaced with the row median.
#'   \item \code{"zero"}: no rows are dropped. All \code{NA}s are
#'         replaced with 0.
#' }
#'
#' The number of dropped rows (if any) is logged.
#'
#' @param dt A \code{data.table} containing identifier columns and slope columns.
#' @param slope_cols Character vector of column names corresponding to slopes.
#' @param method Character scalar. One of \code{"median"} or \code{"zero"}.
#' @param min_non_na Integer scalar. Minimum number of non-NA slope values
#'   required to retain a row when \code{method = "median"}.
#'
#' @return A \code{data.table} with imputed slope columns.
#' @noRd
#' @keywords internal
#' @importFrom logger log_info
.impute_slope_dt <- function(dt,
                             slope_cols,
                             method,
                             min_non_na = 3) {

    if (!method %in% c("median", "zero")) {
        stop("method must be one of: 'median', 'zero'.")
    }

    m <- as.matrix(dt[, slope_cols, with = FALSE])

    if (method == "median") {

        non_na_counts <- rowSums(!is.na(m))
        keep <- non_na_counts >= min_non_na
        n_dropped <- sum(!keep)

        if (n_dropped > 0) {
            logger::log_info(
                "Dropping {n_dropped} rows with < {min_non_na} non-NA values prior to median imputation"
            )
        }

        dt <- dt[keep]
        m <- m[keep, , drop = FALSE]

        for (i in seq_len(nrow(m))) {
            row <- m[i, ]
            na_mask <- is.na(row)
            if (any(na_mask)) {
                m[i, na_mask] <- median(row, na.rm = TRUE)
            }
        }

    } else {

        na_mask <- is.na(m)
        if (any(na_mask)) {
            m[na_mask] <- 0
        }
    }

    dt[, (slope_cols) := data.table::as.data.table(m)]
    dt
}
