# eqtl_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/LEVEL_3"
# region_cell_type_path="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/region_cell_type.tsv"
# egene_union_pairs_path="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3/egene_union_pairs_qval_0.01.tsv"
# output_path="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3/pval_nominal_matrix_qval_0.01.tsv"
# bican.mccarroll.eqtl::get_pval_nominal_matrix(eqtl_dir, region_cell_type_path, egene_union_pairs_path, output_path)


#' Build a pval_nominal matrix for eGene-variant pairs across cell types and regions
#'
#' For each cell type / region in the region-cell-type table, reads the
#' all-pairs tensorQTL output and extracts \code{pval_nominal} values for the
#' eGene-variant pairs provided in \code{egene_union_pairs_path}.  The result
#' is a matrix with one row per eGene-variant pair and one column per
#' cell type / region.
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
#'   \code{phenotype_id}, \code{variant_id}, and \code{qval}.
#' @param output_path Character scalar or \code{NULL}.  If non-NULL, the result
#'   table is written to this path as a tab-delimited file.
#'
#' @return A \code{data.table} with columns \code{phenotype_id},
#'   \code{variant_id}, and one pval_nominal column per cell type / region
#'   (named \code{<cell_type>__<region>}).
#'
#' @export
#' @importFrom data.table fread fwrite
#' @importFrom logger log_info log_warn
get_pval_nominal_matrix <- function(eqtl_dir,
                                    region_cell_type_path,
                                    egene_union_pairs_path,
                                    output_path = NULL) {

    region_cell_type_dt <- data.table::fread(region_cell_type_path)

    result_dt <- data.table::fread(egene_union_pairs_path)
    result_dt[, pair_key := paste0(phenotype_id, "_", variant_id)]
    result_dt[, qval := NULL]

    for (i in seq_len(nrow(region_cell_type_dt))) {
        cell_type <- region_cell_type_dt$cell_type[i]
        region    <- region_cell_type_dt$region[i]
        subdir    <- paste0(cell_type, "__", region)

        eqtl_file     <- file.path(eqtl_dir, subdir, paste0(subdir, ".cis_qtl.txt.gz"))
        allpairs_file <- file.path(eqtl_dir, subdir, paste0(subdir, ".cis_qtl_pairs.txt.gz"))

        logger::log_info("Processing: {subdir}")

        # Check that index file has pval_nominal_threshold
        index_dt <- data.table::fread(eqtl_file, nrows = 0)
        if (!("pval_nominal_threshold" %in% names(index_dt))) {
            logger::log_warn("Skipping (no pval_nominal_threshold): {subdir}")
            next
        }

        # Read all-pairs pval_nominal and join to eGene union
        allpairs_dt <- data.table::fread(
            allpairs_file,
            select = c("phenotype_id", "variant_id", "pval_nominal"),
            showProgress = TRUE
        )
        allpairs_dt[, pair_key := paste0(phenotype_id, "_", variant_id)]

        col_name <- subdir
        result_dt[allpairs_dt, (col_name) := i.pval_nominal, on = "pair_key"]

        logger::log_info("  Matched {sum(!is.na(result_dt[[col_name]]))} of {nrow(result_dt)} pairs")
    }

    result_dt[, pair_key := NULL]

    if (!is.null(output_path)) {
        data.table::fwrite(result_dt, output_path, sep = "\t")
        logger::log_info("Written to: {output_path}")
    }

    return(result_dt)
}
