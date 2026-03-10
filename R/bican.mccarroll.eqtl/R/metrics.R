# Example usage:
#   devtools::load_all()
#
#   # --- Cluster proportions ---
#   get_all_cluster_proportions(
#       cluster_label_path = "/broad/.../cluster_assignments_qval_0.01_k13.tsv",
#       cluster_group_shared = c(11, 0, 5, 4, 2),
#       cluster_group_specificity = c(5, 4, 2, 12, 9, 1, 6, 3, 7, 10, 8)
#   )
#
#   # --- eGene unions ---
#   get_all_egene_unions(
#       eqtl_dir = "/broad/.../results/LEVEL_3"
#   )


# ============================================================
# Cluster proportions
# ============================================================

#' Compute the proportion of genes in a cluster group
#'
#' Reads a cluster assignment file and computes the number and fraction of
#' genes belonging to the specified cluster IDs.
#'
#' @param cluster_label_path Character scalar.  Path to a tab-delimited cluster
#'   assignments file with (at minimum) a \code{cluster} column.
#' @param cluster_group Integer vector.  Cluster IDs to count.
#' @param cluster_name Character scalar.  Label for this group
#'   (e.g. \code{"shared"} or \code{"specificity"}).
#'
#' @return A \code{data.frame} with one row and columns \code{name},
#'   \code{n}, \code{total}, and \code{frac}.
#'
#' @export
#' @importFrom data.table fread
get_cluster_proportions <- function(cluster_label_path, cluster_group,
                                    cluster_name) {
    df <- data.table::fread(cluster_label_path)
    n     <- sum(df$cluster %in% cluster_group)
    total <- nrow(df)
    frac  <- n / total
    return(data.frame(name = cluster_name, n = n, total = total, frac = frac))
}


#' Compute shared and specificity cluster proportions
#'
#' Reads a cluster assignment file and computes the proportion of genes in
#' two cluster groups: shared (non-cell-type-specific) and specificity
#' (clusters exhibiting some degree of cell-type specificity).
#'
#' @param cluster_label_path Character scalar.  Path to a tab-delimited cluster
#'   assignments file with (at minimum) a \code{cluster} column.
#' @param cluster_group_shared Integer vector.  Cluster IDs for the shared /
#'   non-cell-type-specific group.
#' @param cluster_group_specificity Integer vector.  Cluster IDs for the
#'   specificity group.
#'
#' @return A \code{data.frame} with two rows (\code{"shared"} and
#'   \code{"specificity"}) and columns \code{name}, \code{n}, \code{total},
#'   and \code{frac}.
#'
#' @export
get_all_cluster_proportions <- function(cluster_label_path,
                                        cluster_group_shared = c(11, 0, 5, 4, 2),
                                        cluster_group_specificity = c(5, 4, 2, 12, 9, 1, 6, 3, 7, 10, 8)) {
    s <- get_cluster_proportions(cluster_label_path, cluster_group_shared, "shared")
    p <- get_cluster_proportions(cluster_label_path, cluster_group_specificity, "specificity")
    r <- rbind(s, p)
    return(r)
}


# ============================================================
# eGene unions
# ============================================================

#' Compute the union of eGenes across a set of cell types
#'
#' Reads tensorQTL \code{cis_qtl.txt.gz} output files for each specified cell
#' type, filters to eGenes with \code{qval < qval_threshold}, and returns the
#' number of unique genes in the union across all cell types in the group.
#'
#' @param eqtl_dir Character scalar.  Base directory containing per-cell-type
#'   tensorQTL result subdirectories (e.g. \code{MSN_D1_matrix__CaH/}).
#' @param cell_types Character vector.  Names of cell type subdirectories to
#'   include in the union.
#' @param group_name Character scalar.  Label for this group.
#' @param qval_threshold Numeric scalar.  Maximum q-value for a gene to be
#'   considered an eGene (default 0.01).
#'
#' @return A \code{data.frame} with one row and columns \code{name},
#'   \code{n_egenes}, and per-cell-type eGene counts.
#'
#' @export
#' @importFrom data.table fread
get_egene_union <- function(eqtl_dir, cell_types, group_name,
                            qval_threshold = 0.01) {
    all_egenes <- character(0)
    per_ct <- list()

    for (ct in cell_types) {
        f <- file.path(eqtl_dir, ct, paste0(ct, ".cis_qtl.txt.gz"))
        dt <- data.table::fread(f)
        egenes <- dt$phenotype_id[dt$qval < qval_threshold]
        all_egenes <- union(all_egenes, egenes)
        per_ct[[ct]] <- length(egenes)
    }

    return(data.frame(name = group_name,
                      n_egenes = length(all_egenes)))
}


#' Compute eGene unions for MSN subtypes and microglia regions
#'
#' Reads tensorQTL output and computes the union of eGenes (qval < 0.01) for
#' two predefined cell type groups: (1) four MSN subtypes in caudate and
#' (2) microglia in caudate and DFC.
#'
#' @param eqtl_dir Character scalar.  Base directory containing per-cell-type
#'   tensorQTL result subdirectories.
#' @param msn_cell_types Character vector.  Subdirectory names for MSN subtypes.
#' @param microglia_cell_types Character vector.  Subdirectory names for
#'   microglia regions.
#' @param qval_threshold Numeric scalar.  Maximum q-value (default 0.01).
#'
#' @return A \code{data.frame} with two rows (\code{"MSN"} and
#'   \code{"microglia"}) and columns \code{name}, \code{n_egenes}, and
#'   per-cell-type eGene counts.
#'
#' @export
get_all_egene_unions <- function(eqtl_dir,
                                 msn_cell_types = c("MSN_D1_matrix__CaH",
                                                    "MSN_D1_striosome__CaH",
                                                    "MSN_D2_matrix__CaH",
                                                    "MSN_D2_striosome__CaH"),
                                 microglia_cell_types = c("microglia__CaH",
                                                          "microglia__DFC"),
                                 qval_threshold = 0.01) {
    m <- get_egene_union(eqtl_dir, msn_cell_types, "MSN", qval_threshold)
    g <- get_egene_union(eqtl_dir, microglia_cell_types, "microglia", qval_threshold)
    r <- rbind(m, g)
    return(r)
}
