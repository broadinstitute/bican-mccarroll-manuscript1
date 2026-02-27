# Example usage:
#   devtools::load_all()
#   get_sig_coloc(
#       coloc_dir = "/broad/.../coloc/results/LEVEL_3_EUR",
#       gwas_trait = "SCZ_eur",
#       pp_h4_threshold = 0.9,
#       region_cell_type_path = "/broad/.../manuscript_data/region_cell_type.tsv",
#       output_path = "/broad/.../manuscript_data/SCZ_eur_coloc_genes_pp_h4_0.9.tsv"
#   )

#' Get significant colocalization genes
#'
#' For each cell type / region combination, reads the coloc summary file
#' and filters to genes with posterior probability of shared causal variant
#' (PP.H4) at or above a threshold.
#'
#' @param coloc_dir Path to the base coloc results directory
#'   (e.g., \code{.../coloc/results/LEVEL_3_EUR}).  The GWAS trait
#'   subdirectory is appended automatically.
#' @param gwas_trait Name of the GWAS trait subdirectory
#'   (e.g., \code{"SCZ_eur"}, \code{"AD_2022"}).
#' @param pp_h4_threshold Minimum PP.H4 value to retain a gene (numeric).
#' @param region_cell_type_path Path to the region / cell-type TSV
#'   (columns: \code{cell_type}, \code{region}).
#' @param output_path If not \code{NULL}, writes the result to this path.
#'
#' @return A \code{data.table} with columns \code{region}, \code{cell_type},
#'   \code{phenotype_id}, and \code{pp_h4}.
#'
#' @importFrom data.table fread fwrite rbindlist
#' @importFrom logger log_info
#' @export
get_sig_coloc <- function(coloc_dir,
                          gwas_trait,
                          pp_h4_threshold,
                          region_cell_type_path,
                          output_path = NULL) {

    trait_dir <- file.path(coloc_dir, gwas_trait)
    region_cell_type_dt <- fread(region_cell_type_path)
    pp_h4_threshold <- as.numeric(pp_h4_threshold)

    log_info("Collecting coloc genes for {gwas_trait} (PP.H4 >= {pp_h4_threshold})")

    results <- list()

    for (i in seq_len(nrow(region_cell_type_dt))) {
        ct  <- region_cell_type_dt[["cell_type"]][i]
        reg <- region_cell_type_dt[["region"]][i]
        fname <- paste0(ct, "__", reg, ".summary_pp_h4.txt")
        fpath <- file.path(trait_dir, fname)

        dt <- fread(fpath)
        dt <- dt[dt[["pp_h4"]] >= pp_h4_threshold, ]

        if (nrow(dt) > 0L) {
            results[[length(results) + 1L]] <- data.table::data.table(
                region       = reg,
                cell_type    = ct,
                phenotype_id = dt[["phenotype_id"]],
                pp_h4        = dt[["pp_h4"]]
            )
        }
    }

    coloc_dt <- rbindlist(results)
    log_info("Found {nrow(coloc_dt)} coloc entries across {length(results)} cell type/regions")

    if (!is.null(output_path)) {
        fwrite(coloc_dt, output_path, sep = "\t")
        log_info("Written to: {output_path}")
    }

    return(coloc_dt)
}
