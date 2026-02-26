# source("R/paths.R")
#
# options(
#     bican.mccarroll.figures.data_root_dir =
#         "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis",
#
#     bican.mccarroll.figures.out_dir =
#         "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/figure_repository",
#
#     bican.mccarroll.figures.cache_dir =
#         "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/figure_repository/data_cache"
# )

#ct_file <- gene_to_chr_file <- cell_metadata_file <- metacells_file <- de_dir <- de_region_interaction_dir <- data_cache_dir <- outDir <- NULL

#' Generate K-means QC heatmaps for age differential expression matrices
#'
#' Produces two manuscript-ready K-means clustering heatmaps based on age
#' differential expression results: (1) a region-combined heatmap used for
#' gene clustering and (2) a region-specific heatmap using the same gene
#' ordering. Intermediate objects required for plotting are cached as RDS
#' files (one per plot) to avoid recomputation.
#'
#' Column names in the matrices are cleaned prior to plotting by replacing
#' double underscores ("__") and single underscores ("_") with spaces for
#' presentation.
#'
#' All file path arguments may be \code{NULL}, in which case they are resolved
#' using \code{resolve_kmeans_qc_paths()} and the configured data root,
#' cache, and output directory options.
#'
#' @param ct_file Path to the cell type list file. If \code{NULL}, resolved
#'   under the data root.
#' @param gene_to_chr_file Path to the gene-to-chromosome mapping file. If
#'   \code{NULL}, resolved under the data root.
#' @param cell_metadata_file Path to the annotated cell metadata file. If
#'   \code{NULL}, resolved under the data root.
#' @param metacells_file Path to the metacell counts file. If \code{NULL},
#'   resolved under the data root.
#' @param de_dir Directory containing region-combined differential expression
#'   results. If \code{NULL}, resolved under the data root.
#' @param de_region_interaction_dir Directory containing region-interaction
#'   differential expression results. If \code{NULL}, resolved under the data
#'   root.
#' @param outDir Output directory for generated SVG plots. If \code{NULL},
#'   resolved via configured output directory options.
#' @param data_cache_dir Directory used to store cached RDS objects for each
#'   heatmap. If \code{NULL}, resolved via configured cache directory options.
#'
#' @export
plot_kmeans_age <- function(
        ct_file = NULL,
        gene_to_chr_file = NULL,
        cell_metadata_file = NULL,
        metacells_file = NULL,
        de_dir = NULL,
        de_region_interaction_dir = NULL,
        outDir = NULL,
        data_cache_dir = NULL) {

    paths <- resolve_kmeans_paths(
        ct_file = ct_file,
        gene_to_chr_file = gene_to_chr_file,
        cell_metadata_file = cell_metadata_file,
        metacells_file = metacells_file,
        de_dir = de_dir,
        de_region_interaction_dir = de_region_interaction_dir,
        outDir = outDir,
        data_cache_dir = data_cache_dir)

    test <- "age"

    scaling_factor <- 5
    k_use <- 19
    cluster_level_order <- c(2, 10, 6, 3, 5, 14, 9, 13, 4, 1, 15, 19, 18, 7, 17, 8, 11, 12)

    cache_combined <- file.path(paths$data_cache_dir, "kmeans_qc_age_heatmap_region_combined.rds")
    cache_region <- file.path(paths$data_cache_dir, "kmeans_qc_age_heatmap_region_specific.rds")

    # K-means clustering plot with region-combined DE results

    combined_obj <- get_or_build_kmeans_combined_cache_age(
        cache_file = cache_combined,
        paths = paths
    )

    #clean up the names
    combined_obj$lfc_mat_z <- clean_matrix_colnames(combined_obj$lfc_mat_z)
    combined_obj$lfc_mat   <- clean_matrix_colnames(combined_obj$lfc_mat)

    out_file <- file.path(paths$outDir, "kmeans_qc_age_heatmap_region_combined.svg")
    grDevices::svg(out_file, width = 14, height = 7)

    gene_clusters <- bican.mccarroll.de.analysis::plot_kmeans_heatmap_with_cluster_labels(
        combined_obj$lfc_mat_z,
        combined_obj$lfc_mat,
        scaling_factor = scaling_factor,
        k = k_use,
        cluster_level_order = cluster_level_order
    )

    grDevices::dev.off()

    # K-means clustering plot with region-specific DE results (same gene order as region-combined)
    region_obj <- get_or_build_kmeans_region_cache_age(
        cache_file = cache_region,
        combined_cache_file = cache_combined,
        paths = paths
    )

    #clean up the names
    region_obj$lfc_mat_z <- clean_matrix_colnames(region_obj$lfc_mat_z)
    region_obj$lfc_mat   <- clean_matrix_colnames(region_obj$lfc_mat)

    out_file <- file.path(paths$outDir, "kmeans_qc_age_heatmap_region_specific.svg")
    grDevices::svg(out_file, width = 14, height = 7)

    z=bican.mccarroll.de.analysis::plot_kmeans_heatmap_with_cluster_labels(
        region_obj$lfc_mat_z,
        region_obj$lfc_mat,
        scaling_factor = scaling_factor
    )

    grDevices::dev.off()

    invisible(gene_clusters)
}

#clean up the names a bit.
clean_matrix_colnames <- function(mat) {

    if (!is.matrix(mat)) {
        stop("clean_matrix_colnames expects a matrix.")
    }

    cn <- colnames(mat)

    if (!is.null(cn)) {

        ## First replace double underscores
        cn <- gsub("__", " ", cn, fixed = TRUE)

        ## Then replace single underscores
        cn <- gsub("_", " ", cn, fixed = TRUE)

        colnames(mat) <- cn
    }

    mat
}

get_or_build_kmeans_combined_cache_age <- function(cache_file, paths) {

    if (file.exists(cache_file)) {
        return(readRDS(cache_file))
    }

    test <- "age"

    regions_use_metacells <- c("CaH", "Pu", "NAC", "ic", "DFC")

    fdr_cutoff <- 0.01
    abs_lfc_cutoff <- log2(1.05)
    min_tpm <- 10
    regions_use_de_mats <- c("CaH", "DFC")

    cell_types_use <- bican.mccarroll.de.analysis::read_cell_types(paths$ct_file)

    cell_metadata <- bican.mccarroll.de.analysis::read_cell_metadata(paths$cell_metadata_file)
    donor_ages <- bican.mccarroll.de.analysis::extract_donor_ages(cell_metadata)

    gene_to_chr <- bican.mccarroll.de.analysis::read_gene_to_chr(paths$gene_to_chr_file)

    de_age <- bican.mccarroll.de.analysis::read_de_results(
        paths$de_dir,
        test,
        paths$ct_file,
        gene_to_chr
    )

    tmp <- bican.mccarroll.de.analysis::read_metacells(
        paths$metacells_file,
        cell_types_use = cell_types_use,
        regions_use = regions_use_metacells
    )

    metacell_summary <- bican.mccarroll.de.analysis::summarize_metacells(
        tmp$metacells,
        tmp$col_metadata,
        donor_ages
    )

    de_age_mat_list <- bican.mccarroll.de.analysis::prep_de_matrices(
        de_age,
        metacell_summary,
        cell_types_use,
        fdr_cutoff = fdr_cutoff,
        abs_lfc_cutoff = abs_lfc_cutoff,
        min_tpm = min_tpm,
        regions_use = regions_use_de_mats
    )

    obj <- list(
        lfc_mat_z = de_age_mat_list$lfc_mat_z,
        lfc_mat = de_age_mat_list$lfc_mat
    )

    saveRDS(obj, cache_file)

    obj
}


get_or_build_kmeans_region_cache_age <- function(cache_file, combined_cache_file, paths) {

    if (file.exists(cache_file)) {
        return(readRDS(cache_file))
    }

    test <- "age"
    regions_use_region_lfc <- c("CaH", "DFC")

    combined_obj <- readRDS(combined_cache_file)

    cell_types_use <- bican.mccarroll.de.analysis::read_cell_types(paths$ct_file)

    gene_to_chr <- bican.mccarroll.de.analysis::read_gene_to_chr(paths$gene_to_chr_file)

    de_ri_age <- bican.mccarroll.de.analysis::read_de_results(
        paths$de_region_interaction_dir,
        test,
        paths$ct_file,
        gene_to_chr
    )

    de_ri_age_lfc_mat <- bican.mccarroll.de.analysis::prep_region_lfc_matrix(
        de_dt = de_ri_age,
        genes_use = rownames(combined_obj$lfc_mat),
        cell_types_use = cell_types_use,
        regions_use = regions_use_region_lfc
    )

    obj <- list(
        lfc_mat_z = combined_obj$lfc_mat_z,
        lfc_mat = de_ri_age_lfc_mat
    )

    saveRDS(obj, cache_file)

    obj
}


resolve_kmeans_paths <- function(
        ct_file = NULL,
        gene_to_chr_file = NULL,
        cell_metadata_file = NULL,
        metacells_file = NULL,
        de_dir = NULL,
        de_region_interaction_dir = NULL,
        outDir = NULL,
        data_cache_dir = NULL) {

    root <- .resolve_data_root_dir(NULL)

    rel <- list(
        ct_file =
            "sburger_tmp/cell_types_use.txt",

        gene_to_chr_file =
            "sburger_tmp/gene_to_chromosome.txt",

        cell_metadata_file =
            "metadata/CAP_cell_metadata.annotated.txt.gz",

        metacells_file =
            "metacells/LEVEL_3/donor_rxn_DGEList_counts.tsv.gz",

        de_dir =
            "differential_expression/results/LEVEL_3/sex_age/cell_type",

        de_region_interaction_dir =
            "differential_expression/results/LEVEL_3/sex_age/cell_type_region_interaction_absolute_effects"
    )

    pick_in <- function(x, key) {
        if (is.null(x)) {
            return(file.path(root, rel[[key]]))
        }
        .resolve_under_root(root, x)
    }

    out <- .resolve_out_dir(outDir)
    cache <- .resolve_cache_dir(data_cache_dir)

    .ensure_dir(out)
    .ensure_dir(cache)

    list(
        data_root_dir             = root,
        ct_file                   = pick_in(ct_file, "ct_file"),
        gene_to_chr_file          = pick_in(gene_to_chr_file, "gene_to_chr_file"),
        cell_metadata_file        = pick_in(cell_metadata_file, "cell_metadata_file"),
        metacells_file            = pick_in(metacells_file, "metacells_file"),
        de_dir                    = pick_in(de_dir, "de_dir"),
        de_region_interaction_dir = pick_in(de_region_interaction_dir, "de_region_interaction_dir"),
        outDir                    = out,
        data_cache_dir            = cache
    )
}
