# eqtl_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/LEVEL_3"
# region_cell_type_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/manuscript_data/region_cell_type.tsv"
# output_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/manuscript_data/combined_tpm_expression_across_cell_types.tsv"
# bican.mccarroll.eqtl::combine_expression_across_cell_types(eqtl_dir, region_cell_type_path, output_path)


#' Combine gene expression TPM across cell types into a single matrix
#'
#' For each cell type / region group, reads the per-sample gene expression
#' TPM BED file, appends the cell type name to sample column headers, and
#' joins all cell types into one wide matrix keyed by gene ID (\code{pid}).
#'
#' @param eqtl_dir Character scalar.  Path to the eQTL results directory
#'   (e.g. \code{.../LEVEL_3}).  Each cell type subdirectory should contain
#'   a file named \code{<ct>__<reg>.gene_expression_tpm.bed.gz}.
#' @param region_cell_type_path Character scalar.  Path to a tab-delimited
#'   file with columns \code{cell_type} and \code{region}.
#' @param output_path Character scalar or \code{NULL}.  If non-NULL, the
#'   combined matrix is written to this path as a tab-delimited file.
#'
#' @return A \code{data.table} with column \code{pid} (gene ID) and one
#'   column per sample per cell type, named \code{<sample>_<ct>__<reg>}.
#'
#' @export
#' @importFrom data.table fread fwrite
#' @importFrom logger log_info
combine_expression_across_cell_types <- function(eqtl_dir,
                                                  region_cell_type_path,
                                                  output_path = NULL) {

    region_cell_type_dt <- data.table::fread(region_cell_type_path)
    ct_keys <- paste0(region_cell_type_dt$cell_type, "__", region_cell_type_dt$region)
    logger::log_info("Combining expression for {length(ct_keys)} cell type/region groups")

    read_and_format <- function(ct_key) {
        bed_path <- file.path(eqtl_dir, ct_key,
                              paste0(ct_key, ".gene_expression_tpm.bed.gz"))
        logger::log_info("Reading: {ct_key}")
        dt <- data.table::fread(bed_path)

        # Remove chr, start, end columns (first 3)
        dt[, c(1, 2, 3) := NULL]

        # Rename sample columns: append _<ct_key> (skip first col which is pid)
        old_names <- names(dt)[-1]
        new_names <- paste0(old_names, "_", ct_key)
        data.table::setnames(dt, old_names, new_names)

        return(dt)
    }

    expression_list <- lapply(ct_keys, read_and_format)

    # Merge all cell types by pid (full outer join)
    combined_dt <- expression_list[[1]]
    for (i in seq(2, length(expression_list))) {
        combined_dt <- merge(combined_dt, expression_list[[i]],
                             by = "pid", all = TRUE)
    }

    logger::log_info("Combined: {nrow(combined_dt)} genes x {ncol(combined_dt) - 1} sample columns")

    if (!is.null(output_path)) {
        data.table::fwrite(combined_dt, output_path, sep = "\t")
        logger::log_info("Written to: {output_path}")
    }

    return(combined_dt)
}
