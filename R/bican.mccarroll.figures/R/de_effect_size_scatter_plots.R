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

#ct_file <- gene_to_chr_file <- de_results_dir<- de_region_interaction_dir<- data_cache_dir<- outDir<- NULL

#' Generate sex and age differential expression scatter plots for manuscript figures
#'
#' @param ct_file Path to the cell type list file. If \code{NULL}, resolved
#'   under the data root.
#' @param gene_to_chr_file Path to the gene-to-chromosome mapping file. If
#'   \code{NULL}, resolved under the data root.
#' @param de_results_dir Base directory for differential expression results.
#'   If \code{NULL}, resolved under the data root.
#' @param de_region_interaction_dir Directory containing region-interaction
#'   differential expression results. If \code{NULL}, resolved under
#'   \code{de_results_dir}.
#' @param data_cache_dir Directory used to store cached DE tables as TSV
#'   files. If \code{NULL}, resolved via configured cache directory options.
#' @param outDir Output directory for generated SVG plots. If \code{NULL},
#'   resolved via configured output directory options.
#'
#' @export
de_sex_age_scatter_plots <- function(
        ct_file = NULL,
        gene_to_chr_file = NULL,
        de_results_dir = NULL,
        de_region_interaction_dir = NULL,
        data_cache_dir = NULL,
        outDir = NULL) {


    paths <- .resolve_sex_age_scatter_plot_paths(
        ct_file = ct_file,
        gene_to_chr_file = gene_to_chr_file,
        de_results_dir = de_results_dir,
        de_region_interaction_dir = de_region_interaction_dir,
        data_cache_dir = data_cache_dir,
        outDir = outDir)

    cache_file_age <- file.path(paths$data_cache_dir, "de_age.tsv")
    cache_file_sex <- file.path(paths$data_cache_dir, "de_sex.tsv")

    age_df <- get_or_build_de_cache(
        test = "age",
        cache_file = cache_file_age,
        de_region_interaction_dir = paths$de_region_interaction_dir,
        ct_file = paths$ct_file,
        gene_to_chr_file = paths$gene_to_chr_file)

    sex_df <- get_or_build_de_cache(
        test = "female_vs_male",
        cache_file = cache_file_sex,
        de_region_interaction_dir = paths$de_region_interaction_dir,
        ct_file = paths$ct_file,
        gene_to_chr_file = paths$gene_to_chr_file)

    #########################
    #Main figure plot group
    #########################

    plot_de_scatter_svg(
        df = sex_df,
        test = "sex",
        cell_type1 = "OPC",
        cell_type2 = "astrocyte",
        region1 = "CaH",
        region2 = "CaH",
        xlab_prefix="Sex DE, ",
        outDir = paths$outDir,
        fdr_cutoff = 0.05,
        add_fit = TRUE)

    plot_de_scatter_svg(
        df = age_df,
        test = "age",
        cell_type1 = "MSN_D1_matrix",
        cell_type2 = "MSN_D2_matrix",
        region1 = "CaH",
        region2 = "CaH",
        xlab_prefix="Age DE, ",
        outDir = paths$outDir,
        fdr_cutoff = 0.05,
        add_fit = TRUE)

    plot_de_scatter_svg(
        df = age_df,
        test = "age",
        cell_type1 = "astrocyte",
        cell_type2 = "astrocyte",
        region1 = "CaH",
        region2 = "DFC",
        xlab_prefix="Age DE, ",
        outDir = paths$outDir,
        fdr_cutoff = 0.05,
        add_fit = TRUE)

    plot_de_scatter_svg(
        df = age_df,
        test = "age",
        cell_type1 = "OPC",
        cell_type2 = "astrocyte",
        region1 = "CaH",
        region2 = "CaH",
        xlab_prefix="Age DE, ",
        outDir = paths$outDir,
        fdr_cutoff = 0.05,
        add_fit = TRUE)

    ####################
    #Supplemental 1:
    ####################

    plot_de_scatter_svg(
        df = age_df,
        test = "age",
        cell_type1 = "MSN_D1_matrix",
        cell_type2 = "MSN_D1_matrix",
        region1 = "CaH",
        region2 = "Pu",
        outDir = paths$outDir,
        fdr_cutoff = 0.05,
        add_fit = TRUE)

    plot_de_scatter_svg(
        df = age_df,
        test = "age",
        cell_type1 = "MSN_D1_matrix",
        cell_type2 = "MSN_D1_matrix",
        region1 = "CaH",
        region2 = "NAC",
        outDir = paths$outDir,
        fdr_cutoff = 0.05,
        add_fit = TRUE)

    plot_de_scatter_svg(
        df = age_df,
        test = "age",
        cell_type1 = "MSN_D1_matrix",
        cell_type2 = "MSN_D1_matrix",
        region1 = "Pu",
        region2 = "NAC",
        outDir = paths$outDir,
        fdr_cutoff = 0.05,
        add_fit = TRUE)

    #Supplemental 2:
    plot_de_scatter_svg(
        df = age_df,
        test = "age",
        cell_type1 = "astrocyte",
        cell_type2 = "astrocyte",
        region1 = "CaH",
        region2 = "Pu",
        outDir = paths$outDir,
        fdr_cutoff = 0.05,
        add_fit = TRUE)

    plot_de_scatter_svg(
        df = age_df,
        test = "age",
        cell_type1 = "astrocyte",
        cell_type2 = "astrocyte",
        region1 = "CaH",
        region2 = "NAC",
        outDir = paths$outDir,
        fdr_cutoff = 0.05,
        add_fit = TRUE)

    plot_de_scatter_svg(
        df = age_df,
        test = "age",
        cell_type1 = "astrocyte",
        cell_type2 = "astrocyte",
        region1 = "Pu",
        region2 = "NAC",
        outDir = paths$outDir,
        fdr_cutoff = 0.05,
        add_fit = TRUE)


}


plot_de_scatter_svg <- function(
        df,
        test,
        cell_type1,
        cell_type2,
        region1,
        region2,
        xlab_prefix=NULL,
        outDir,
        fdr_cutoff = 0.05,
        add_fit = TRUE,
        width = 7,
        height = 7) {

    fileStr <- paste(
        "de_scatter_plot_",
        test,
        "_", cell_type1, "_", cell_type2,
        "_", region1, "_vs_", region2,
        ".svg",
        sep = "")

    if (!is.null(outDir)) {
        out_file <- file.path(outDir, fileStr)
        grDevices::svg(out_file, width = width, height = height)
        on.exit(grDevices::dev.off(), add = TRUE)
    }

    bican.mccarroll.de.analysis::plot_de_scatter(
        df,
        cell_type1,
        cell_type2,
        region1,
        region2,
        fdr_cutoff = fdr_cutoff,
        add_fit = add_fit,
        xlab_prefix=xlab_prefix)

    invisible()
}
get_or_build_de_cache <- function(
        test,
        cache_file,
        de_region_interaction_dir,
        ct_file,
        gene_to_chr_file) {

    if (file.exists(cache_file)) {
        logger::log_info("Using cached data from {cache_file}")
        de_ri <- utils::read.table(
            cache_file, header = TRUE, sep = "\t",
            stringsAsFactors = FALSE, check.names = FALSE)
        data.table::setDT(de_ri)

    } else {
        logger::log_info(
            "No cached data from {cache_file} regenerating data from sources.  This can take a few minutes"
        )
        gene_to_chr <- bican.mccarroll.de.analysis::read_gene_to_chr(
            gene_to_chr_file)

        de_ri <- bican.mccarroll.de.analysis::read_de_results(
            de_region_interaction_dir,
            test,
            ct_file,
            gene_to_chr)

        utils::write.table(
            de_ri,
            file = cache_file,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
    }

    de_ri
}

.resolve_sex_age_scatter_plot_paths <- function(
        de_results_dir = NULL,
        de_region_interaction_dir = NULL,
        gene_to_chr_file = NULL,
        ct_file = NULL,
        outDir = NULL,
        data_cache_dir = NULL) {

    root <- .resolve_data_root_dir(NULL)

    rel <- list(
        de_results_dir =
            "differential_expression/results",

        de_region_interaction_dir =
            "differential_expression/results/LEVEL_3/sex_age/cell_type_region_interaction_absolute_effects",

        gene_to_chr_file =
            "sburger_tmp/gene_to_chromosome.txt",

        ct_file =
            "sburger_tmp/cell_types_use.txt"
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
        de_results_dir            = pick_in(de_results_dir, "de_results_dir"),
        de_region_interaction_dir = pick_in(de_region_interaction_dir, "de_region_interaction_dir"),
        gene_to_chr_file          = pick_in(gene_to_chr_file, "gene_to_chr_file"),
        ct_file                   = pick_in(ct_file, "ct_file"),
        outDir                    = out,
        data_cache_dir            = cache
    )
}
