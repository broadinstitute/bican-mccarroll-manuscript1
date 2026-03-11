# Example usage:
#   devtools::load_all()
#   get_index_snp_start_distance(
#       eqtl_dir = "/broad/.../results/LEVEL_3",
#       region_cell_type_path = "/broad/.../manuscript_data/region_cell_type.tsv",
#       index_snp_matrix_path = "/broad/.../manuscript_data/index_snp_slope_matrix_with_zero_impute_qval_0.01.tsv",
#       output_path = "/broad/.../manuscript_data/index_snp_start_distance_qval_0.01.tsv"
#   )

#' Get start_distance for each index SNP gene-variant pair
#'
#' Reads \code{start_distance} (variant position minus TSS, in bp) from the
#' tensorQTL index files for each cell type / region, deduplicates, and joins
#' to the index SNP list.
#'
#' @param eqtl_dir Character scalar.  Base directory containing per-cell-type
#'   eQTL result subdirectories (\code{<cell_type>__<region>/}).
#' @param region_cell_type_path Character scalar.  Path to a tab-delimited file
#'   with columns \code{cell_type} and \code{region}.
#' @param index_snp_matrix_path Character scalar.  Path to the index SNP slope
#'   matrix TSV (output of \code{\link{get_index_snp_slope_matrix_with_impute}}),
#'   must have columns \code{phenotype_id} and \code{variant_id}.
#' @param output_path Character scalar or \code{NULL}.  If non-NULL, the result
#'   table is written to this path as a tab-delimited file.
#'
#' @return A \code{data.table} with columns \code{phenotype_id},
#'   \code{variant_id}, and \code{start_distance}.
#'
#' @export
#' @importFrom data.table fread fwrite rbindlist
#' @importFrom logger log_info
get_index_snp_start_distance <- function(eqtl_dir,
                                          region_cell_type_path,
                                          index_snp_matrix_path,
                                          output_path = NULL) {

    index_dt <- data.table::fread(index_snp_matrix_path,
                                  select = c("phenotype_id", "variant_id"))
    logger::log_info("Index SNPs: {nrow(index_dt)} gene-variant pairs")

    region_cell_type_dt <- data.table::fread(region_cell_type_path)

    dist_list <- vector("list", nrow(region_cell_type_dt))
    for (i in seq_len(nrow(region_cell_type_dt))) {
        ct  <- region_cell_type_dt$cell_type[i]
        reg <- region_cell_type_dt$region[i]
        subdir <- paste0(ct, "__", reg)
        eqtl_file <- file.path(eqtl_dir, subdir,
                               paste0(subdir, ".cis_qtl.txt.gz"))

        logger::log_info("Reading: {subdir}")
        dt <- data.table::fread(eqtl_file,
                                select = c("phenotype_id", "variant_id",
                                           "start_distance"))
        dist_list[[i]] <- dt
    }

    all_dist <- data.table::rbindlist(dist_list)
    all_dist <- unique(all_dist, by = c("phenotype_id", "variant_id"))
    logger::log_info("Unique gene-variant pairs with start_distance: {nrow(all_dist)}")

    result <- merge(index_dt, all_dist,
                    by = c("phenotype_id", "variant_id"), all.x = TRUE)

    n_matched <- sum(!is.na(result$start_distance))
    n_missing <- sum(is.na(result$start_distance))
    logger::log_info("Matched: {n_matched} / {nrow(result)} (missing: {n_missing})")

    if (!is.null(output_path)) {
        data.table::fwrite(result, output_path, sep = "\t")
        logger::log_info("Written to: {output_path}")
    }

    invisible(result)
}
