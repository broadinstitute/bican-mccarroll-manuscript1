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

# de_dir<-gene_to_chr_file<-ct_file<-outDir<-data_cache_dir<- NULL

#' Generate volcano plots for sex and age differential expression
#'
#' @param de_dir Directory containing region-combined differential expression
#'   results. If \code{NULL}, resolved under the data root.
#' @param gene_to_chr_file Path to the gene-to-chromosome mapping file. If
#'   \code{NULL}, resolved under the data root.
#' @param ct_file Path to the cell type list file. If \code{NULL}, resolved
#'   under the data root.
#' @param outDir Output directory for generated SVG plots. If \code{NULL},
#'   resolved via configured output directory options.
#' @param data_cache_dir Directory used to store cached DE tables as TSV
#'   files. If \code{NULL}, resolved via configured cache directory options.
#'
#' @export
plot_de_volcano <- function(
        de_dir = NULL,
        gene_to_chr_file = NULL,
        ct_file = NULL,
        outDir = NULL,
        data_cache_dir = NULL) {

    paths <- resolve_de_volcano_paths(
        de_dir = de_dir,
        gene_to_chr_file = gene_to_chr_file,
        ct_file = ct_file,
        outDir = outDir,
        data_cache_dir = data_cache_dir)

    region_use <- NA
    ct <- "microglia"

    fdr_cutoff <- 0.05
    abs_log_fc_cutoff <- log2(1.05)

    chr_color_map=c("X"="cornflowerblue", "Y"="tomato", "autosome"="black")

    ###############################################
    #Gather the sex DE data averaged across regions
    ###############################################
    test <- "female_vs_male"
    cache_file <- file.path(paths$data_cache_dir, paste0("volcano_", test, "_", ct, ".tsv"))

    sex_df <- get_or_build_de_volcano_cache(
        test = test,
        cache_file = cache_file,
        de_dir = paths$de_dir,
        ct_file = paths$ct_file,
        gene_to_chr_file = paths$gene_to_chr_file)

    test="female_vs_male"
    sex_df <- modify_chromosome_labels(sex_df)

    p1<-bican.mccarroll.de.analysis::plot_de_volcano_gg(
        sex_df,
        cell_type_use = ct,
        region_use = region_use,
        fdr_cutoff = fdr_cutoff,
        abs_log_fc_cutoff = abs_log_fc_cutoff,
        show_title = FALSE,
        chr_color_map=chr_color_map)

    p1<-add_style_volcano(p1)
    p1 <- p1 + ggplot2::labs(x = "Sex DE log2 fold change")

    fileStr <- paste("de_volcano_", test, "_", ct, ".svg", sep = "")
    out_file <- file.path(paths$outDir, fileStr)
    ggplot2::ggsave(
        filename = out_file,
        plot = p1,
        width = 7,
        height = 5,
        units = "in"
    )

    #Gather the age DE data averaged across regions
    test <- "age"
    cache_file <- file.path(paths$data_cache_dir, paste0("volcano_", test, "_", ct, ".tsv"))
    age_df <- get_or_build_de_volcano_cache(
        test = test,
        cache_file = cache_file,
        de_dir = paths$de_dir,
        ct_file = paths$ct_file,
        gene_to_chr_file = paths$gene_to_chr_file)

    age_df=modify_chromosome_labels(age_df)

    p2<-bican.mccarroll.de.analysis::plot_de_volcano_gg(
        age_df,
        cell_type_use = ct,
        region_use = region_use,
        fdr_cutoff = fdr_cutoff,
        abs_log_fc_cutoff = abs_log_fc_cutoff,
        show_title = FALSE,
        chr_color_map=chr_color_map)

    p2<-add_style_volcano(p2)
    p2 <- p2 + ggplot2::labs(x = "Age DE log2 fold change")
    fileStr <- paste("de_volcano_", test, "_", ct, ".svg", sep = "")
    out_file <- file.path(paths$outDir, fileStr)

    ggplot2::ggsave(
        filename = out_file,
        plot = p2,
        width = 7,
        height = 5,
        units = "in"
    )

    invisible(out_file)
}

add_style_volcano<-function (p) {
    p <- p +
        ggplot2::theme_classic() +
        ggplot2::theme(
            axis.title.x = ggplot2::element_text(size = ggplot2::rel(2)),
            axis.title.y = ggplot2::element_text(size = ggplot2::rel(2)),
            legend.text = ggplot2::element_text(size = ggplot2::rel(1.5))
        ) +
        ggplot2::guides(
            color = ggplot2::guide_legend(override.aes = list(size = 2.5))
        )
    return (p)
}

#Drop mitochondria, then map to "X", "Y", or "autosome"
modify_chromosome_labels <- function(df) {
    #explicitly remove mitochondria
    df=df[df$chr!="M",]
    #map remaining genes to "X", "Y", or "autosome"
    df$chr <- ifelse(df$chr %in% c("X", "Y"), df$chr, "autosome")
    return (df)
}

get_or_build_de_volcano_cache <- function(
        test,
        cache_file,
        de_dir,
        ct_file,
        gene_to_chr_file) {

    if (file.exists(cache_file)) {
        logger::log_info("Using cached data from {cache_file}")
        df=data.table::fread(cache_file, header = TRUE, sep = "\t")
        return (df)
    }

    logger::log_info(
        "No cached data from {cache_file} regenerating data from sources.  This can take a few minutes"
    )
    gene_to_chr <- bican.mccarroll.de.analysis::read_gene_to_chr(gene_to_chr_file)

    de_dt <- bican.mccarroll.de.analysis::read_de_results(
        de_dir,
        test,
        ct_file,
        gene_to_chr)

    utils::write.table(
        de_dt,
        file = cache_file,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE)

    de_dt
}

resolve_de_volcano_paths <- function(
        de_dir = NULL,
        gene_to_chr_file = NULL,
        ct_file = NULL,
        outDir = NULL,
        data_cache_dir = NULL) {

    root <- .resolve_data_root_dir(NULL)

    rel <- list(
        de_dir =
            "differential_expression/results/LEVEL_3/sex_age/cell_type",

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
        data_root_dir    = root,
        de_dir           = pick_in(de_dir, "de_dir"),
        gene_to_chr_file = pick_in(gene_to_chr_file, "gene_to_chr_file"),
        ct_file          = pick_in(ct_file, "ct_file"),
        outDir           = out,
        data_cache_dir   = cache
    )
}
