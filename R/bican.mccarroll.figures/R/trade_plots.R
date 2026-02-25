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

#' Generate the TRADE manuscript figure
#'
#' Plots TRADE figures for analysis of autosome age effects across cell types and regions.
#' Results are written to SVG files in \code{outDir}. Intermediate TRADE results are cached
#' as tab-delimited text files in \code{data_cache_dir}.
#'
#' @param de_dir Path to differential expression results for the regions-combined dataset.
#'   If \code{NULL}, a default relative path under
#'   \code{options("bican.mccarroll.figures.data_root_dir")} is used via
#'   \code{resolve_trade_paths()}.
#' @param de_region_subset_dir Path to differential expression results for the region-subset dataset.
#'   If \code{NULL}, a default relative path under
#'   \code{options("bican.mccarroll.figures.data_root_dir")} is used via
#'   \code{resolve_trade_paths()}.
#' @param de_region_interaction_dir Path to region-interaction differential expression results.
#'   This function accepts and resolves the path for completeness even if the current plot does not use it.
#' @param gene_to_chr_path Path to a two-column mapping from gene identifier to chromosome.
#'   If \code{NULL}, a default relative path under
#'   \code{options("bican.mccarroll.figures.data_root_dir")} is used via
#'   \code{resolve_trade_paths()}.
#' @param ct_file Path to a text file containing one cell type per line. The file is read with
#'   \code{scan()} to define the displayed cell type subset and order. If \code{NULL}, a default
#'   relative path under \code{options("bican.mccarroll.figures.data_root_dir")} is used via
#'   \code{resolve_trade_paths()}.
#' @param data_cache_dir Directory used to store cached TRADE results as tab-delimited text.
#'   If \code{NULL}, the directory is resolved from
#'   \code{options("bican.mccarroll.figures.cache_dir")} via \code{resolve_trade_paths()}.
#' @param outDir Output directory used to write SVG files. If \code{NULL}, the directory is resolved
#'   from \code{options("bican.mccarroll.figures.out_dir")} via \code{resolve_trade_paths()}.
#'
#' @return Invisibly returns \code{NULL}. This function is called for its side effects.
#'
#' @seealso
#'   \code{bican.mccarroll.differentialexpression::load_trade_data},
#'   \code{bican.mccarroll.differentialexpression::run_trade},
#'   \code{bican.mccarroll.differentialexpression::trade_barplot},
#'   \code{bican.mccarroll.differentialexpression::trade_heatmap}
#' @export
plot_trade_analysis <- function(
        de_dir = NULL,
        de_region_subset_dir = NULL,
        de_region_interaction_dir = NULL,
        gene_to_chr_path = NULL,
        ct_file = NULL,
        data_cache_dir = NULL,
        outDir = NULL) {

    paths <- resolve_trade_paths(
        de_dir = de_dir,
        de_region_subset_dir = de_region_subset_dir,
        de_region_interaction_dir = de_region_interaction_dir,
        gene_to_chr_path = gene_to_chr_path,
        ct_file = ct_file,
        data_cache_dir = data_cache_dir,
        outDir = outDir)

    region_order <- c("CaH", "Pu", "NAC", "ic", "DFC")
    cell_types_use <- scan(paths$ct_file, what = character(), quiet = TRUE)

    run_trade_autosomes <- function(de_dt) {
        chr <- NULL  # R CMD CHECK
        de_auto <- de_dt[chr %in% 1:22]
        bican.mccarroll.differentialexpression::run_trade(de_auto)
    }

    filter_ic_to_non_neurons <- function(de_dt) {
        region <- cell_type <- NULL  # R CMD CHECK
        non_neuron_types <- c("astrocyte", "OPC", "oligodendrocyte", "microglia")
        de_dt[!(region == "ic" & !(cell_type %in% non_neuron_types))]
    }

    cell_types_use <- c("MSN_D1_matrix","MSN_D1_striosome","MSN_D2_matrix","MSN_D2_striosome",
                        "glutamatergic_L23IT","glutamatergic_L4IT","glutamatergic_L5IT","glutamatergic_L6IT",
                        "GABA_TAC3-PLPP4","GABA_PTHLH-PVALB","GABA_PVALB","GABA_SST","GABA_VIP","GABA_LAMP5",
                        "astrocyte","OPC","oligodendrocyte","microglia")
    # --------------------------------------------------------------------------
    # Dataset 1: regions combined, age (AUTOSOMES ONLY)
    # --------------------------------------------------------------------------

    cache_file <- file.path(paths$data_cache_dir,
                            "trade_dataset1_age_autosomes.tsv")

    if (file.exists(cache_file)) {
        #the trade functions assume data.table not data frame.
        trade_auto <- data.table::fread(cache_file)
    } else {
        de_dt <- bican.mccarroll.differentialexpression::load_trade_data(
            data_path = paths$de_dir, contrast = "age",
            gene_to_chr_path = paths$gene_to_chr_path,
            cellTypeListFile = paths$ct_file, regions_use = NULL)

        trade_auto <- run_trade_autosomes(de_dt)

        utils::write.table(trade_auto, file = cache_file, sep = "\t",
                           row.names = FALSE, col.names = TRUE, quote = FALSE)
    }

    p_bar_age <- bican.mccarroll.differentialexpression::trade_barplot(
        trade_auto, cell_types_use = cell_types_use,
        value_var = "trade_twi")
    p_bar_age <- p_bar_age + ggplot2::geom_col(fill = "black")

    save_plot_svg(p_bar_age,
                  out_file = "trade_dataset1_age_autosomes_barplot.svg",
                  out_dir = paths$outDir, width=5, height=9)

    # --------------------------------------------------------------------------
    # Dataset 2: region subset, age (AUTOSOMES ONLY)
    # --------------------------------------------------------------------------

    cache_file <- file.path(paths$data_cache_dir,
                            "trade_dataset2_age_subset_region_autosomes.tsv")

    if (file.exists(cache_file)) {
        trade_auto <- data.table::fread(cache_file)
    } else {

        de_dt <- bican.mccarroll.differentialexpression::load_trade_data(
            data_path = paths$de_region_subset_dir, contrast = "age",
            gene_to_chr_path = paths$gene_to_chr_path,
            cellTypeListFile = paths$ct_file, regions_use = NULL)

        de_dt <- filter_ic_to_non_neurons(de_dt)
        trade_auto <- run_trade_autosomes(de_dt)

        utils::write.table(trade_auto, file = cache_file, sep = "\t",
                           row.names = FALSE, col.names = TRUE, quote = FALSE)
    }

    cell_types <- c("MSN_D1_matrix","MSN_D1_striosome","MSN_D2_matrix","MSN_D2_striosome",
                    "glutamatergic_L23IT","glutamatergic_L4IT","glutamatergic_L5IT","glutamatergic_L6IT",
                    "GABA_TAC3-PLPP4","GABA_PTHLH-PVALB","GABA_PVALB","GABA_SST","GABA_VIP","GABA_LAMP5",
                    "astrocyte","OPC","oligodendrocyte","microglia")

    p_heat_age <- bican.mccarroll.differentialexpression::trade_heatmap(
        trade_auto, cell_types_use = cell_types,
        region_order = region_order, value_var = "trade_twi")

    save_plot_svg(p_heat_age,
                  out_file = "trade_dataset2_age_subset_region_autosomes_heatmap.svg",
                  out_dir = paths$outDir, width=5, height=8)

    # --------------------------------------------------------------------------
    # Dataset 3: region interaction, age (AUTOSOMES ONLY)
    # --------------------------------------------------------------------------

    cache_file <- file.path(paths$data_cache_dir,
                            "trade_dataset3_age_interaction_region_autosomes.tsv")

    if (file.exists(cache_file)) {
        trade_auto <- data.table::fread(cache_file)
    } else {

        de_dt <- bican.mccarroll.differentialexpression::load_trade_data(
            data_path = paths$de_region_interaction_dir, contrast = "age",
            gene_to_chr_path = paths$gene_to_chr_path,
            cellTypeListFile = paths$ct_file, regions_use = NULL)

        de_dt <- filter_ic_to_non_neurons(de_dt)
        trade_auto <- run_trade_autosomes(de_dt)

        utils::write.table(trade_auto, file = cache_file, sep = "\t",
                           row.names = FALSE, col.names = TRUE, quote = FALSE)
    }

    p_heat_age <- bican.mccarroll.differentialexpression::trade_heatmap(
        trade_auto, cell_types_use = NULL,
        region_order = region_order, value_var = "trade_twi")

    save_plot_svg(p_heat_age,
                  out_file = "trade_dataset3_age_interaction_region_autosomes_heatmap.svg",
                  out_dir = paths$outDir, width=5, height=8)

    invisible(NULL)
}



# plot_trade_analysis_barplots <- function(
#         de_dir = NULL,
#         de_region_subset_dir = NULL,
#         de_region_interaction_dir = NULL,
#         gene_to_chr_path = NULL,
#         ct_file = NULL,
#         data_cache_dir = NULL,
#         outDir = NULL) {
#
#     paths <- resolve_trade_paths(
#         de_dir = de_dir,
#         de_region_subset_dir = de_region_subset_dir,
#         de_region_interaction_dir = de_region_interaction_dir,
#         gene_to_chr_path = gene_to_chr_path,
#         ct_file = ct_file,
#         data_cache_dir = data_cache_dir,
#         outDir = outDir)
#
#     region_order <- c("CaH", "Pu", "NAC", "ic", "DFC")
#     cell_types_use <- scan(paths$ct_file, what = character(), quiet = TRUE)
#
#     run_trade_autosomes <- function(de_dt) {
#         chr <- NULL  # R CMD CHECK
#         de_auto <- de_dt[chr %in% 1:22]
#         bican.mccarroll.differentialexpression::run_trade(de_auto)
#     }
#
#     filter_ic_to_non_neurons <- function(de_dt) {
#         region <- cell_type <- NULL  # R CMD CHECK
#         non_neuron_types <- c("astrocyte", "OPC", "oligodendrocyte", "microglia")
#         de_dt[!(region == "ic" & !(cell_type %in% non_neuron_types))]
#     }
#
#
#     # --------------------------------------------------------------------------
#     # Dataset 3: regions interaction, age (AUTOSOMES ONLY)
#     # --------------------------------------------------------------------------
#
#     cache_file <- file.path(paths$data_cache_dir,
#                             "trade_dataset3_age_interaction_region_autosomes.tsv")
#
#     if (file.exists(cache_file)) {
#         trade_auto <- data.table::fread(cache_file)
#     } else {
#
#         de_dt <- bican.mccarroll.differentialexpression::load_trade_data(
#             data_path = paths$de_region_interaction_dir, contrast = "age",
#             gene_to_chr_path = paths$gene_to_chr_path,
#             cellTypeListFile = paths$ct_file, regions_use = NULL)
#
#         de_dt <- filter_ic_to_non_neurons(de_dt)
#         trade_auto <- run_trade_autosomes(de_dt)
#
#         utils::write.table(trade_auto, file = cache_file, sep = "\t",
#                            row.names = FALSE, col.names = TRUE, quote = FALSE)
#     }
#
#
#     hard_coded_cell_type_order<- c(
#         "MSN_D1_matrix",
#         "MSN_D1_striosome",
#         "MSN_D2_matrix",
#         "MSN_D2_striosome",
#         "glutamatergic_L23IT",
#         "glutamatergic_L4IT",
#         "glutamatergic_L5IT",
#         "glutamatergic_L6IT",
#         "GABA_PTHLH-PVALB",
#         "GABA_TAC3-PLPP4",
#         "GABA_PVALB",
#         "GABA_SST",
#         "GABA_VIP",
#         "GABA_LAMP5",
#         "astrocyte",
#         "OPC",
#         "oligodendrocyte",
#         "microglia"
#     )
#
#     regions <- unique(trade_auto$region)
#     result <- list()
#
#     for (r in regions) {
#         d <- data.table::copy(trade_auto)
#         d <- d[region == r, ]
#
#         p <- bican.mccarroll.differentialexpression::trade_barplot(
#             d, cell_types_use = hard_coded_cell_type_order, value_var = "trade_twi")
#
#         # Add title = region name
#         p <- p + ggplot2::ggtitle(r) +
#             ggplot2::theme(
#                 plot.title = ggplot2::element_text(hjust = 0.5)
#             )
#
#         result[[as.character(r)]] <- p
#     }
#
#     # Arrange into 3 rows x 2 cols (fills by row by default)
#     big_plot <- cowplot::plot_grid(
#         plotlist = unname(result),
#         nrow = 3,
#         ncol = 2,
#         align = "hv"
#     )
#
#     # One cell type at a time, for cell types that are in at least 2 regions.
#     cell_types <- unique(trade_auto$cell_type)
#
#     result <- list()
#
#     for (ct in cell_types) {
#         d <- data.table::copy(trade_auto)
#         d <- d[cell_type == ct, ]
#         if (dim (d)[1] <2) {
#             next
#         }
#         p <- trade_barplot_regions(d, regions_use = region_order, value_var = "trade_twi")
#
#         # Add title = cell type
#         p <- p +
#             ggplot2::ggtitle(ct) +
#             ggplot2::theme(
#                 plot.title = ggplot2::element_text(hjust = 0.5)
#             )
#
#         result[[as.character(ct)]] <- p
#     }
#
#     big_plot <- cowplot::plot_grid(
#         plotlist = unname(result),
#         nrow = 4,
#         ncol = 3,
#         align = "hv"
#     )
#
#     big_plot
#     invisible (big_plot)
#
# }

resolve_trade_paths <- function(
        de_dir = NULL,
        de_region_subset_dir = NULL,
        de_region_interaction_dir = NULL,
        gene_to_chr_path = NULL,
        ct_file = NULL,
        outDir = NULL,
        data_cache_dir = NULL) {

    root <- .resolve_data_root_dir(NULL)

    rel <- list(
        de_dir =
            "differential_expression/results/LEVEL_3/sex_age/cell_type",

        de_region_subset_dir =
            "differential_expression/results/LEVEL_3/sex_age/cell_type_subset_region",

        de_region_interaction_dir =
            "differential_expression/results/LEVEL_3/sex_age/cell_type_region_interaction_absolute_effects",

        gene_to_chr_path =
            "metadata/gene_to_chromosome.txt",

        ct_file =
            "differential_expression/metadata/cell_types_for_de_filtering_plot.txt"
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
        de_dir                    = pick_in(de_dir, "de_dir"),
        de_region_subset_dir      = pick_in(de_region_subset_dir, "de_region_subset_dir"),
        de_region_interaction_dir = pick_in(de_region_interaction_dir, "de_region_interaction_dir"),
        gene_to_chr_path          = pick_in(gene_to_chr_path, "gene_to_chr_path"),
        ct_file                   = pick_in(ct_file, "ct_file"),
        outDir                    = out,
        data_cache_dir            = cache
    )
}




save_plot_svg <- function(plot, out_file, out_dir = ".", width = 14, height = 7) {
    out_svg <- file.path(out_dir, out_file)

    svglite::svglite(file = out_svg, width = width, height = height)
    on.exit(grDevices::dev.off(), add = TRUE)

    print(plot)

    invisible(out_svg)
}


#' Plot TRADE results as a barplot across regions (expects input already filtered)
#'
#' @param trade_results A data.table produced by run_trade(), already filtered to a single cell type.
#' @param regions_use Optional character vector specifying the order and subset of regions to display.
#' @param region_var Character scalar naming the region column.
#' @param value_var Character scalar naming the TRADE statistic to plot.
#'
#' @return A ggplot object.
#' @export
trade_barplot_regions <- function(trade_results,
                                  regions_use = NULL,
                                  region_var = "region",
                                  value_var = "trade_twi") {

    # Make R CMD CHECK happy
    region <- value <- N <- .N <- NULL

    dt <- data.table::as.data.table(trade_results)

    if (!(region_var %in% names(dt))) {
        stop("trade_barplot_regions(): region_var not found in trade_results.")
    }
    if (!(value_var %in% names(dt))) {
        stop("trade_barplot_regions(): value_var not found in trade_results.")
    }

    # Sanity: require exactly one cell type in input
    if ("cell_type" %in% names(dt)) {
        n_ct <- data.table::uniqueN(dt[["cell_type"]])
        if (n_ct != 1L) {
            stop("trade_barplot_regions(): requires trade_results to contain exactly one cell_type. Filter before calling.")
        }
    }

    # Require one row per region
    dup_rg <- dt[, .N, by = region_var][N > 1L]
    if (nrow(dup_rg) > 0L) {
        stop("trade_barplot_regions(): requires one row per region. Filter/aggregate first.")
    }

    dt <- dt[, list(region = get(region_var), value = get(value_var))]

    if (!is.null(regions_use)) {
        dt <- dt[region %in% regions_use]
        dt[, region := factor(region, levels = rev(regions_use))]
    } else {
        dt <- dt[order(value)]
        dt[, region := factor(region, levels = region)]
    }

    ggplot2::ggplot(dt, ggplot2::aes(x = value, y = region)) +
        ggplot2::geom_col() +
        ggplot2::labs(x = paste0("TRADE ", value_var), y = NULL) +
        ggplot2::theme_classic()
}
