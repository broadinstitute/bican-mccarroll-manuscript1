# eqtl_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/LEVEL_3"
# region_cell_type_path="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/region_cell_type.tsv"
# egene_union_pairs_path="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3/egene_union_pairs_qval_0.01.tsv"
# output_path="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3/pval_nominal_threshold_matrix_qval_0.01.tsv"
# bican.mccarroll.eqtl::get_pval_nominal_threshold_matrix(eqtl_dir, region_cell_type_path, egene_union_pairs_path, output_path)


#' Build a p-value nominal threshold matrix across cell types and regions
#'
#' For each cell type / region in the region-cell-type table, reads the
#' tensorQTL cis-eQTL index file and extracts the \code{pval_nominal_threshold}
#' for each eGene-variant pair.  The result is a matrix with one row per
#' eGene-variant pair and one column per cell type / region.
#'
#' Cell type / region combinations whose index file lacks a
#' \code{pval_nominal_threshold} column are skipped with a warning.
#'
#' @param eqtl_dir Character scalar.  Base directory containing per-cell-type
#'   eQTL result subdirectories (\code{<cell_type>__<region>/}).
#' @param region_cell_type_path Character scalar.  Path to a tab-delimited file
#'   with columns \code{cell_type} and \code{region}.
#' @param egene_union_pairs_path Character scalar.  Path to the eGene union
#'   pairs file (output of \code{\link{get_egene_union_pairs}}), with columns
#'   \code{phenotype_id} and \code{variant_id}.
#' @param output_path Character scalar or \code{NULL}.  If non-NULL, the result
#'   table is written to this path as a tab-delimited file.
#'
#' @return A \code{data.table} with columns \code{phenotype_id},
#'   \code{variant_id}, and one \code{pval_nominal_threshold} column per
#'   cell type / region (named \code{<cell_type>__<region>}).
#'
#' @export
#' @importFrom data.table fread fwrite
#' @importFrom logger log_info log_warn
get_pval_nominal_threshold_matrix <- function(eqtl_dir,
                                              region_cell_type_path,
                                              egene_union_pairs_path,
                                              output_path = NULL) {

    region_cell_type_dt <- data.table::fread(region_cell_type_path)

    result_dt <- data.table::fread(egene_union_pairs_path,
                                   select = c("phenotype_id", "variant_id"))

    for (i in seq_len(nrow(region_cell_type_dt))) {
        cell_type <- region_cell_type_dt$cell_type[i]
        region    <- region_cell_type_dt$region[i]
        subdir    <- paste0(cell_type, "__", region)
        eqtl_file <- file.path(eqtl_dir, subdir, paste0(subdir, ".cis_qtl.txt.gz"))

        logger::log_info("Processing: {subdir}")

        index_dt <- data.table::fread(eqtl_file)

        if (!("pval_nominal_threshold" %in% names(index_dt))) {
            logger::log_warn("Skipping (no pval_nominal_threshold): {subdir}")
            next
        }

        index_dt <- index_dt[, .(phenotype_id, variant_id, pval_nominal_threshold)]

        col_name <- subdir
        result_dt[index_dt, (col_name) := i.pval_nominal_threshold,
                  on = .(phenotype_id, variant_id)]

        logger::log_info("  Matched {sum(!is.na(result_dt[[col_name]]))} of {nrow(result_dt)} pairs")
    }

    if (!is.null(output_path)) {
        data.table::fwrite(result_dt, output_path, sep = "\t")
        logger::log_info("Written to: {output_path}")
    }

    return(result_dt)
}
