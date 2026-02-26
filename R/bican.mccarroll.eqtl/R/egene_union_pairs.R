# eqtl_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/LEVEL_3"
# region_cell_type_path="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/region_cell_type.tsv"
# qval_threshold=0.05
# output_path="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3/egene_union_pairs_qval_0.05.tsv"
# bican.mccarroll.eqtl::get_egene_union_pairs(eqtl_dir, region_cell_type_path, qval_threshold, output_path)


#' Get union of eGene-variant pairs across cell types and regions
#'
#' For each cell type / region combination in the region-cell-type table,
#' reads the tensorQTL cis-eQTL index file, filters to eGenes at
#' \code{qval < qval_threshold}, and returns the union of unique
#' (phenotype_id, variant_id) pairs across all groups.
#'
#' Each subdirectory under \code{eqtl_dir} is expected to follow the naming
#' convention \code{<cell_type>__<region>/}, containing a file named
#' \code{<cell_type>__<region>.cis_qtl.txt.gz}.
#'
#' @param eqtl_dir Character scalar.  Base directory containing per-cell-type
#'   eQTL result subdirectories.
#' @param region_cell_type_path Character scalar.  Path to a tab-delimited file
#'   with columns \code{cell_type} and \code{region}.
#' @param qval_threshold Numeric scalar.  q-value threshold for calling an
#'   eGene significant.  Default \code{0.05}.
#' @param output_path Character scalar or \code{NULL}.  If non-NULL, the result
#'   table is written to this path as a tab-delimited file.
#'
#' @return A \code{data.table} with columns \code{phenotype_id}, \code{variant_id},
#'   and \code{qval}, containing one row per unique eGene-variant pair.
#'
#' @export
#' @importFrom data.table fread fwrite rbindlist
#' @importFrom logger log_info
get_egene_union_pairs <- function(eqtl_dir,
                                  region_cell_type_path,
                                  qval_threshold = 0.05,
                                  output_path = NULL) {

    region_cell_type_dt <- data.table::fread(region_cell_type_path)

    results <- vector("list", nrow(region_cell_type_dt))

    for (i in seq_len(nrow(region_cell_type_dt))) {
        cell_type <- region_cell_type_dt$cell_type[i]
        region    <- region_cell_type_dt$region[i]
        subdir    <- paste0(cell_type, "__", region)
        filename  <- paste0(subdir, ".cis_qtl.txt.gz")
        eqtl_file <- file.path(eqtl_dir, subdir, filename)

        logger::log_info("Reading eQTL index: {eqtl_file}")

        dt <- data.table::fread(eqtl_file, select = c("phenotype_id", "variant_id", "qval"))
        dt <- dt[dt$qval < qval_threshold, ]

        logger::log_info("  {cell_type} / {region}: {nrow(dt)} eGenes at qval < {qval_threshold}")
        results[[i]] <- dt
    }

    combined <- data.table::rbindlist(results)

    # Deduplicate to unique (phenotype_id, variant_id) pairs
    combined[, pair_key := paste0(phenotype_id, "_", variant_id)]
    result <- combined[!duplicated(pair_key), .(phenotype_id, variant_id, qval)]

    logger::log_info("Total unique eGene-variant pairs: {nrow(result)}")

    if (!is.null(output_path)) {
        data.table::fwrite(result, output_path, sep = "\t")
        logger::log_info("Written to: {output_path}")
    }

    return(result)
}
