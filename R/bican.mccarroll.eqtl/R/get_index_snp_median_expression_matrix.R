#' Compute median TPM expression per cell type for index-SNP genes
#'
#' Takes the combined gene expression TPM matrix (output of
#' \code{\link{combine_expression_across_cell_types}}), computes the median
#' expression across donors for each gene in each cell type / region group,
#' then filters to genes present in the index-SNP slope matrix.  NA values
#' are zero-imputed, consistent with the zero-imputed slope matrix.
#'
#' @param expression_path Character scalar.  Path to the combined gene
#'   expression TPM matrix TSV (output of
#'   \code{\link{combine_expression_across_cell_types}}).  First column
#'   should be gene IDs (\code{pid}); remaining columns are sample-level
#'   expression values named \code{<sample>_<cell_type>__<region>}.
#' @param region_cell_type_path Character scalar.  Path to a tab-delimited
#'   file with columns \code{cell_type} and \code{region}.
#' @param index_snp_path Character scalar.  Path to the index-SNP slope
#'   matrix TSV (output of
#'   \code{\link{get_index_snp_slope_matrix_with_impute}}).
#'   Only the \code{phenotype_id} column is used to filter genes.
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
get_index_snp_median_expression_matrix <- function(expression_path,
                                            region_cell_type_path,
                                            index_snp_path,
                                            output_path = NULL) {

    # 1. Get heatmap gene list from index-SNP matrix
    index_snp_dt <- data.table::fread(index_snp_path, select = "phenotype_id")
    heatmap_genes <- unique(index_snp_dt[["phenotype_id"]])
    logger::log_info("Loaded {length(heatmap_genes)} heatmap genes")

    # 2. Get cell type keys
    region_cell_type_dt <- data.table::fread(region_cell_type_path)
    ct_keys <- paste0(region_cell_type_dt$cell_type, "__", region_cell_type_dt$region)
    logger::log_info("Cell type/region groups: {length(ct_keys)}")

    # 3. Read combined expression matrix
    expression_dt <- data.table::fread(expression_path)
    gene_ids <- expression_dt[[1]]
    expression_m <- as.matrix(expression_dt[, -1, with = FALSE])
    logger::log_info("Expression matrix: {nrow(expression_m)} genes x {ncol(expression_m)} columns")

    # 4. Extract cell type from column names (strip sample prefix before first _)
    col_cell_types <- sub("^[^_]*_", "", colnames(expression_m))

    # 5. Compute median expression per cell type
    median_m <- sapply(ct_keys, function(ct) {
        cols <- which(col_cell_types == ct)
        if (length(cols) == 0) return(rep(NA_real_, nrow(expression_m)))
        apply(expression_m[, cols, drop = FALSE], 1, stats::median, na.rm = TRUE)
    })

    # 6. Zero-impute NAs
    median_m[is.na(median_m)] <- 0

    # 7. Build output data.table
    result_dt <- data.table::as.data.table(median_m)
    result_dt[, Gene := gene_ids]
    data.table::setcolorder(result_dt, c("Gene", ct_keys))

    # 8. Filter to heatmap genes
    result_dt <- result_dt[Gene %in% heatmap_genes]
    logger::log_info("Filtered to {nrow(result_dt)} heatmap genes")

    # 9. Write output
    if (!is.null(output_path)) {
        data.table::fwrite(result_dt, output_path, sep = "\t")
        logger::log_info("Written to: {output_path}")
    }

    invisible(result_dt)
}
