# index_snp_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3/index_snp_slope_matrix_with_median_impute_qval_0.01.tsv"
# region_cell_type_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/region_cell_type.tsv"
# expression_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3/combined_gene_expression_tpm.tsv"
# output_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3/heatmap_index_snp_median_expression_qval_0.01.tsv"
# bican.mccarroll.eqtl::get_heatmap_index_snp_median_expression(index_snp_path, region_cell_type_path, expression_path, output_path)


#' Get median gene expression per cell type for heatmap index-SNP genes
#'
#' For each cell type / region group, computes the median TPM expression
#' across samples, then filters to the genes present in the index-SNP
#' slope matrix (i.e. heatmap genes).  NA values within each gene row
#' are replaced with the row median.
#'
#' @param index_snp_path Character scalar.  Path to the index-SNP slope
#'   matrix TSV (output of
#'   \code{\link{get_index_snp_slope_matrix_with_median_impute}}).
#'   Only the \code{phenotype_id} column is used.
#' @param region_cell_type_path Character scalar.  Path to a tab-delimited
#'   file with columns \code{cell_type} and \code{region}.
#' @param expression_path Character scalar.  Path to the combined gene
#'   expression TPM matrix TSV.  First column should be gene IDs
#'   (\code{pid}); remaining columns are sample-level expression values
#'   with names formatted as \code{<sample>_<cell_type>__<region>}.
#' @param output_path Character scalar or \code{NULL}.  If non-NULL, the
#'   result is written to this path as a tab-delimited file.
#'
#' @return A \code{data.table} with columns \code{Gene} and one column
#'   per cell type / region group containing median expression values.
#'
#' @export
#' @importFrom data.table fread fwrite as.data.table
#' @importFrom stats median
#' @importFrom logger log_info
get_heatmap_index_snp_median_expression <- function(index_snp_path,
                                                     region_cell_type_path,
                                                     expression_path,
                                                     output_path = NULL) {

    # Get heatmap gene list from index-SNP matrix
    index_snp_dt <- data.table::fread(index_snp_path, select = "phenotype_id")
    heatmap_genes <- unique(index_snp_dt[["phenotype_id"]])
    logger::log_info("Loaded {length(heatmap_genes)} heatmap genes")

    # Get cell type keys
    region_cell_type_dt <- data.table::fread(region_cell_type_path)
    ct_keys <- paste0(region_cell_type_dt$cell_type, "__", region_cell_type_dt$region)
    logger::log_info("Cell type/region groups: {length(ct_keys)}")

    # Read expression matrix
    expression_dt <- data.table::fread(expression_path)
    gene_ids <- expression_dt[[1]]
    expression_m <- as.matrix(expression_dt[, -1, with = FALSE])
    logger::log_info("Expression matrix: {nrow(expression_m)} genes x {ncol(expression_m)} columns")

    # Extract cell type from column names (strip sample prefix before first _)
    col_cell_types <- sub("^[^_]*_", "", colnames(expression_m))

    # Compute median expression per cell type
    median_m <- sapply(ct_keys, function(ct) {
        cols <- which(col_cell_types == ct)
        if (length(cols) == 0) return(rep(NA_real_, nrow(expression_m)))
        apply(expression_m[, cols, drop = FALSE], 1, stats::median, na.rm = TRUE)
    })

    # Median-impute NAs within each row
    median_m <- t(apply(median_m, 1, function(x) {
        x[is.na(x)] <- stats::median(x, na.rm = TRUE)
        x
    }))

    # Build output data.table
    result_dt <- data.table::as.data.table(median_m)
    result_dt[, Gene := gene_ids]
    data.table::setcolorder(result_dt, c("Gene", ct_keys))

    # Filter to heatmap genes
    result_dt <- result_dt[Gene %in% heatmap_genes]
    logger::log_info("Filtered to {nrow(result_dt)} heatmap genes")

    if (!is.null(output_path)) {
        data.table::fwrite(result_dt, output_path, sep = "\t")
        logger::log_info("Written to: {output_path}")
    }

    return(result_dt)
}
