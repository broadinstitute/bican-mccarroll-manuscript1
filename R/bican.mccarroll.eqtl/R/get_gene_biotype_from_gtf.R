#' Extract gene biotype from GTF
#'
#' Parses the GTF used for eQTL discovery and matches phenotype IDs to gene
#' records (by gene_name, then gene_id with version, then gene_id without
#' version) to extract \code{gene_type} (e.g., protein_coding, lncRNA).
#'
#' @param gtf_path Character scalar. Path to the GTF file used for eQTL
#'   discovery (e.g., \code{GRCh38_ensembl_v43.gtf}).
#' @param index_snp_matrix_path Character scalar. Path to the index SNP slope
#'   matrix TSV (must have a \code{phenotype_id} column).
#' @param output_path Character scalar or \code{NULL}. If non-NULL, the result
#'   table is written to this path as a tab-delimited file.
#'
#' @return A \code{data.table} with columns \code{phenotype_id},
#'   \code{gene_id}, \code{gene_id_noversion}, and \code{gene_type}.
#'
#' @export
#' @importFrom data.table fread fwrite data.table
#' @importFrom logger log_info
get_gene_biotype_from_gtf <- function(gtf_path,
                                       index_snp_matrix_path,
                                       output_path = NULL) {

    # ---- 1. Parse GTF ----
    log_info("Parsing GTF: {gtf_path}")
    gtf <- rtracklayer::import(gtf_path)
    genes_gr <- gtf[gtf$type == "gene"]

    gtf_dt <- data.table::data.table(
        gene_id   = genes_gr$gene_id,
        gene_name = genes_gr$gene_name,
        gene_type = genes_gr$gene_type
    )
    gtf_dt[, gene_id_noversion := sub("\\..*", "", gene_id)]

    log_info("GTF: {nrow(gtf_dt)} genes")

    # ---- 2. Load phenotype IDs ----
    index_dt <- data.table::fread(index_snp_matrix_path, select = "phenotype_id")
    pheno_ids <- unique(index_dt$phenotype_id)
    log_info("Unique phenotype_ids: {length(pheno_ids)}")

    # ---- 3. Match phenotype_ids to GTF ----
    result_dt <- data.table::data.table(phenotype_id = pheno_ids)

    match_name <- gtf_dt[, .(gene_id, gene_name, gene_type, gene_id_noversion)]
    result_dt <- merge(result_dt, match_name,
                       by.x = "phenotype_id", by.y = "gene_name", all.x = TRUE)
    n_by_name <- sum(!is.na(result_dt$gene_id))
    log_info("Matched by gene_name: {n_by_name} / {length(pheno_ids)}")

    n_matched <- sum(!is.na(result_dt$gene_type))
    log_info("Total GTF match: {n_matched} / {length(pheno_ids)} ({round(100 * n_matched / length(pheno_ids), 1)}%)")

    # ---- 4. Output ----
    out <- result_dt[, .(phenotype_id, gene_id, gene_id_noversion, gene_type)]

    if (!is.null(output_path)) {
        data.table::fwrite(out, output_path, sep = "\t")
        log_info("Written to: {output_path}")
    }

    invisible(out)
}
