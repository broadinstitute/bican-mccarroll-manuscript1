source("R/paths.R")

options(
    bican.mccarroll.figures.data_root_dir =
        "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis",

    bican.mccarroll.figures.out_dir =
        "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/figure_repository",

    bican.mccarroll.figures.cache_dir =
        "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/figure_repository/data_cache"
)



#' Generate donor-level GEX plots vs age
#'
#' This function computes all required intermediate objects (with RDS caching)
#' and then generates the manuscript plots (heatmap, scatters, and the 5x6 panel grid).
#' It does not save files; it prints plots to the active device.
#'
#' @param ct_file Path to cell types list. If NULL, resolved by resolver.
#' @param cell_metadata_file Path to annotated cell metadata. If NULL, resolved by resolver.
#' @param metacells_file Path to metacells DGEList counts. If NULL, resolved by resolver.
#' @param outDir Output directory. If NULL, resolved by resolver.
#' @param data_cache_dir Cache directory. If NULL, resolved by resolver.
#'
#' @return Invisibly NULL.
plot_donor_gex_vs_age <- function(ct_file = NULL,
                                  cell_metadata_file = NULL,
                                  metacells_file = NULL,
                                  outDir = NULL,
                                  data_cache_dir = NULL) {

    paths <- .resolve_donor_gex_age_paths(
        ct_file = ct_file,
        cell_metadata_file = cell_metadata_file,
        metacells_file = metacells_file,
        outDir = outDir,
        data_cache_dir = data_cache_dir
    )

    cache <- load_or_build_donor_gex_age_cache(paths)
    donor_ages <- cache$donor_ages
    metacell_cr_list <- cache$metacell_cr_list

    ## -----------------------
    ## Plots
    ## -----------------------

    gene_list <- c("GPNMB","CD163","CD163L1","MS4A4A","MS4A6A","CYP27A1","CLEC5A","CLEC2B","FPR3","GPR141","DOCK5","IQGAP2","IL15","CYTL1","MAFB","FOXP1","ESR1","THRB")
    heatmap_key <- "microglia__CaH"
    gs_gaps <- c(6,12,14)

    #This uses pheatmap, so needs to captured by the current plotting device.
    # save to SVG donor_age_microglia_heatmap.svg
    p_heatmap <- bican.mccarroll.de.analysis::plot_donor_gex_age_heatmap(
        metacell_cr_list[[heatmap_key]],
        gs = gene_list,
        donor_ages = donor_ages,
        gs_gaps = gs_gaps,
        cluster_gs = FALSE,
        transpose = FALSE
    )

    gene_list <- c("PCLO","CAPS2","STX3","MCTP2","JPH1","SLIT1","SEMA5A","TSHZ3","IQGAP1","SIPA1L2","CTNNAL1","DNMBP","KIF13B","EPB41L4A","LAMB1","FBLN5","TLL1","ABI3BP","PTPN13","PTPRK","DENND3","LRRK1","MKNK1","SGSM1","ABCG1","LEPR","C2","PPP1R1B","GPR88","RGS9","RGS12","ADCY3","PDE1B","PKIA","KCNQ3","SCN4B","SYT6","SYT12","LRRTM3","SYNDIG1","ASTN2","PCDH9","EPHA5","SLITRK2","ADGRB2","CSMD2","HMGCS1","MSMO1","MICOS10","PEX5","ABCD2","SAT1","DCX","FOXG1")
    heatmap_key <- "MSN_D1_matrix__CaH"
    gs_gaps <- c(5,8,14,18,24,27,34,36,46,52)

    # save to SVG donor_age_D1_matrix_heatmap.svg
    bican.mccarroll.de.analysis::plot_donor_gex_age_heatmap(
        metacell_cr_list[[heatmap_key]],
        gs = gene_list,
        donor_ages = donor_ages,
        gs_gaps = gs_gaps,
        cluster_gs = FALSE,
        transpose = FALSE
    )

    #############################################
    # Big scatter plot/table of genes/cell types
    #############################################

    gene_list <- c("FKBP5", "RGS9", "RYR3", "GRIA1", "CLEC2B")

    keys <- c("MSN_D1_matrix__CaH", "glutamatergic_L23IT__DFC", "astrocyte__CaH",
        "OPC__CaH", "oligodendrocyte__CaH", "microglia__CaH")

    nice_names <- c("D1 matrix MSN", "L23IT glutamatergic\nneuron", "Astrocyte",
        "OPC", "Oligodendrocyte", "Microglia")

    plot_list <- list()
    idx <- 0L
    n_row <- length(gene_list)
    n_col <- length(keys)

    for (r in seq_len(n_row)) {
        for (c in seq_len(n_col)) {
            idx <- idx + 1L
            plot_list[[idx]] <- make_panel(
                metacell_cr_list = metacell_cr_list,
                key = keys[c],
                gene = gene_list[r],
                donor_ages = donor_ages,
                show_x = (r == n_row),
                show_spearman = TRUE,
                spearman_text_size = 4
            )
        }
    }

    grid <- cowplot::plot_grid(
        plotlist = plot_list,
        nrow = n_row,
        align = "hv"
    )

    col_label_grobs <- lapply(
        nice_names,
        function(x) cowplot::ggdraw() + cowplot::draw_label(x, size = 12)
    )
    col_labels <- cowplot::plot_grid(plotlist = col_label_grobs, nrow = 1)

    grid_with_cols <- cowplot::plot_grid(
        grid,
        col_labels,
        ncol = 1,
        rel_heights = c(1, 0.08)
    )

    row_label_grobs <- lapply(
        gene_list,
        function(g) cowplot::ggdraw() + cowplot::draw_label(g, angle = 90, size = 14)
    )
    row_labels <- cowplot::plot_grid(plotlist = row_label_grobs, ncol = 1)

    #save to SVG donor_age_scatterplot_grid.svg
    all_plot <- cowplot::plot_grid(
        row_labels,
        grid_with_cols,
        nrow = 1,
        rel_widths = c(0.06, 1)
    )


    invisible(NULL)
}

#' Generate donor-level GEX plots vs age
#'
#' This function computes all required intermediate objects (with RDS caching)
#' and then generates the manuscript plots and saves them as SVG files.
#'
#' @param ct_file Path to cell types list. If NULL, resolved by resolver.
#' @param cell_metadata_file Path to annotated cell metadata. If NULL, resolved by resolver.
#' @param metacells_file Path to metacells DGEList counts. If NULL, resolved by resolver.
#' @param outDir Output directory. If NULL, resolved by resolver.
#' @param data_cache_dir Cache directory. If NULL, resolved by resolver.
#' @param width SVG width (inches). Default 8.
#' @param height SVG height (inches). Default 8.
#'
#' @return Invisibly NULL.
plot_donor_gex_vs_age <- function(ct_file = NULL,
                                  cell_metadata_file = NULL,
                                  metacells_file = NULL,
                                  outDir = NULL,
                                  data_cache_dir = NULL,
                                  width = 8,
                                  height = 8) {

    paths <- .resolve_donor_gex_age_paths(
        ct_file = ct_file,
        cell_metadata_file = cell_metadata_file,
        metacells_file = metacells_file,
        outDir = outDir,
        data_cache_dir = data_cache_dir
    )

    cache <- load_or_build_donor_gex_age_cache(paths)
    donor_ages <- cache$donor_ages
    metacell_cr_list <- cache$metacell_cr_list

    ## ============================================================
    ## Heatmap 1: Microglia
    ## ============================================================

    gene_list <- c("GPNMB","CD163","CD163L1","MS4A4A","MS4A6A","CYP27A1","CLEC5A","CLEC2B","FPR3","GPR141","DOCK5","IQGAP2","IL15","CYTL1","MAFB","FOXP1","ESR1","THRB")
    heatmap_key <- "microglia__CaH"
    gs_gaps <- c(6,12,14)

    out_file <- file.path(paths$outDir, "donor_age_microglia_heatmap.svg")

    svglite::svglite(out_file, width = 9, height = 6)
    bican.mccarroll.de.analysis::plot_donor_gex_age_heatmap(
        metacell_cr_list[[heatmap_key]],
        gs = gene_list,
        donor_ages = donor_ages,
        gs_gaps = gs_gaps,
        cluster_gs = FALSE,
        transpose = FALSE
    )
    grDevices::dev.off()

    ## ============================================================
    ## Heatmap 2: D1 matrix MSN
    ## ============================================================

    gene_list <- c("PCLO","CAPS2","STX3","MCTP2","JPH1","SLIT1","SEMA5A","TSHZ3","IQGAP1","SIPA1L2","CTNNAL1","DNMBP","KIF13B","EPB41L4A","LAMB1","FBLN5","TLL1","ABI3BP","PTPN13","PTPRK","DENND3","LRRK1","MKNK1","SGSM1","ABCG1","LEPR","C2","PPP1R1B","GPR88","RGS9","RGS12","ADCY3","PDE1B","PKIA","KCNQ3","SCN4B","SYT6","SYT12","LRRTM3","SYNDIG1","ASTN2","PCDH9","EPHA5","SLITRK2","ADGRB2","CSMD2","HMGCS1","MSMO1","MICOS10","PEX5","ABCD2","SAT1","DCX","FOXG1")
    heatmap_key <- "MSN_D1_matrix__CaH"
    gs_gaps <- c(5,8,14,18,24,27,34,36,46,52)

    out_file <- file.path(paths$outDir, "donor_age_D1_matrix_heatmap.svg")

    svglite::svglite(out_file, width = 7, height = 8)
    bican.mccarroll.de.analysis::plot_donor_gex_age_heatmap(
        metacell_cr_list[[heatmap_key]],
        gs = gene_list,
        donor_ages = donor_ages,
        gs_gaps = gs_gaps,
        cluster_gs = FALSE,
        transpose = FALSE
    )
    grDevices::dev.off()

    ## ============================================================
    ## Scatterplot grid
    ## ============================================================

    gene_list <- c("FKBP5", "RGS9", "RYR3", "GRIA1", "CLEC2B")

    keys <- c("MSN_D1_matrix__CaH", "glutamatergic_L23IT__DFC", "astrocyte__CaH",
              "OPC__CaH", "oligodendrocyte__CaH", "microglia__CaH")

    nice_names <- c("D1 matrix MSN", "L23IT glutamatergic\nneuron", "Astrocyte",
                    "OPC", "Oligodendrocyte", "Microglia")

    plot_list <- list()
    idx <- 0L
    n_row <- length(gene_list)
    n_col <- length(keys)

    for (r in seq_len(n_row)) {
        for (c in seq_len(n_col)) {
            idx <- idx + 1L
            plot_list[[idx]] <- make_panel(
                metacell_cr_list = metacell_cr_list,
                key = keys[c],
                gene = gene_list[r],
                donor_ages = donor_ages,
                show_x = (r == n_row),
                show_spearman = TRUE,
                spearman_text_size = 4
            )
        }
    }

    grid <- cowplot::plot_grid(
        plotlist = plot_list,
        nrow = n_row,
        align = "hv"
    )

    col_label_grobs <- lapply(
        nice_names,
        function(x) cowplot::ggdraw() + cowplot::draw_label(x, size = 10)
    )
    col_labels <- cowplot::plot_grid(plotlist = col_label_grobs, nrow = 1)

    col_strip_height <- 0.08

    grid_with_cols <- cowplot::plot_grid(
        grid,
        col_labels,
        ncol = 1,
        rel_heights = c(1, col_strip_height)
    )

    row_label_grobs <- lapply(
        gene_list,
        function(g) cowplot::ggdraw() + cowplot::draw_label(g, angle = 90, size = 14)
    )
    row_labels <- cowplot::plot_grid(plotlist = row_label_grobs, ncol = 1)

    # Add a bottom spacer so the left label column has the same total height as grid_with_cols.
    row_labels_padded <- cowplot::plot_grid(
        row_labels,
        cowplot::ggdraw(),          # blank spacer
        ncol = 1,
        rel_heights = c(1, col_strip_height)
    )

    all_plot <- cowplot::plot_grid(
        row_labels_padded,
        grid_with_cols,
        nrow = 1,
        rel_widths = c(0.06, 1)
    )

    out_file <- file.path(paths$outDir, "donor_age_scatterplot_grid.svg")

    svglite::svglite(out_file, width = 8, height = 8)
    print(all_plot)
    grDevices::dev.off()

    invisible(NULL)
}

#for each individual plot, apply some modifications.
make_panel <- function(metacell_cr_list,
                       key,
                       gene,
                       donor_ages,
                       show_x = FALSE,
                       show_spearman = TRUE,
                       spearman_text_size = 4,
                       rho_threshold = 0.2,
                       y_axis_floor = 10) {

    exp_vec <- metacell_cr_list[[key]][gene, ]

    p <- bican.mccarroll.de.analysis::plot_donor_gex_age_scatterplot(
        exp_vec,
        donor_ages,
        main = "",
        show_spearman = show_spearman,
        size = spearman_text_size,
        rho_threshold = rho_threshold,
        y_axis_floor = y_axis_floor
    )

    p +
        ggplot2::theme_classic() +
        ggplot2::theme(
            axis.title = ggplot2::element_blank(),

            ## X axis: only show on bottom row
            axis.text.x  = if (show_x) ggplot2::element_text(size = 8) else ggplot2::element_blank(),
            axis.ticks.x = if (show_x) ggplot2::element_line() else ggplot2::element_blank(),

            ## Y axis: ALWAYS show for every panel
            axis.text.y  = ggplot2::element_text(size = 8),
            axis.ticks.y = ggplot2::element_line(),

            plot.margin = ggplot2::margin(2, 2, 2, 2)
        )
}

#' Load or build cached data for donor-level GEX vs age plots
#'
#' Caching policy:
#' - If the cache RDS exists, it is loaded and trusted without validation.
#' - Otherwise, inputs are read and the bundled list is saved to cache.
#'
#' @param paths A list produced by resolve_donor_gex_age_paths().
#'
#' @return A named list with donor_ages and metacell_cr_list.
load_or_build_donor_gex_age_cache <- function(paths) {

    if (!is.list(paths)) {
        stop("paths must be a list produced by resolve_donor_gex_age_paths().", call. = FALSE)
    }
    required <- c("ct_file", "cell_metadata_file", "metacells_file", "cache_rds")
    missing <- setdiff(required, names(paths))
    if (length(missing) > 0L) {
        stop("paths is missing required fields: ", paste(missing, collapse = ", "), call. = FALSE)
    }

    if (file.exists(paths$cache_rds)) {
        cache <- readRDS(paths$cache_rds)
        return(cache)
    }

    regions_use_metacells <- c("CaH", "Pu", "NAC", "ic", "DFC")

    cell_types_use <- bican.mccarroll.de.analysis::read_cell_types(paths$ct_file)

    cell_metadata <- bican.mccarroll.de.analysis::read_cell_metadata(paths$cell_metadata_file)
    donor_ages <- bican.mccarroll.de.analysis::extract_donor_ages(cell_metadata)

    tmp <- bican.mccarroll.de.analysis::read_metacells(
        paths$metacells_file,
        cell_types_use = cell_types_use,
        regions_use = regions_use_metacells
    )

    metacell_cr_list <- bican.mccarroll.de.analysis::split_metacells_by_cell_type_region(
        tmp$metacells,
        tmp$col_metadata,
        donor_ages
    )

    cache <- list(
        donor_ages = donor_ages,
        metacell_cr_list = metacell_cr_list
    )

    saveRDS(cache, paths$cache_rds)

    cache
}


#' Resolve input/output paths for donor-level GEX vs age plots
#'
#' This resolver maps NULL path parameters onto defaults derived from:
#'   - getOption("bican.mccarroll.figures.data_root_dir")
#'   - getOption("bican.mccarroll.figures.out_dir")
#'   - getOption("bican.mccarroll.figures.cache_dir")
#'
#' It also creates outDir and data_cache_dir if they do not exist.
#'
#' @param ct_file Path to cell types list. If NULL, resolved relative to data_root_dir.
#' @param cell_metadata_file Path to annotated cell metadata. If NULL, resolved relative
#'   to dirname(data_root_dir).
#' @param metacells_file Path to metacells DGEList counts. If NULL, resolved relative to data_root_dir.
#' @param outDir Output directory. If NULL, taken from option bican.mccarroll.figures.out_dir.
#' @param data_cache_dir Cache directory. If NULL, taken from option bican.mccarroll.figures.cache_dir.
#'
#' @return A named list of resolved paths.
.resolve_donor_gex_age_paths <- function(ct_file = NULL,
                                         cell_metadata_file = NULL,
                                         metacells_file = NULL,
                                         outDir = NULL,
                                         data_cache_dir = NULL) {

    get_opt_or_stop <- function(opt_name) {
        x <- getOption(opt_name)
        if (is.null(x) || !nzchar(x)) {
            stop("Required option is not set: ", opt_name, call. = FALSE)
        }
        x
    }

    data_root_dir <- get_opt_or_stop("bican.mccarroll.figures.data_root_dir")

    if (is.null(outDir)) {
        outDir <- get_opt_or_stop("bican.mccarroll.figures.out_dir")
    }
    if (is.null(data_cache_dir)) {
        data_cache_dir <- get_opt_or_stop("bican.mccarroll.figures.cache_dir")
    }

    if (!dir.exists(outDir)) {
        dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
    }
    if (!dir.exists(data_cache_dir)) {
        dir.create(data_cache_dir, recursive = TRUE, showWarnings = FALSE)
    }

    if (is.null(ct_file)) {
        ct_file <- file.path(data_root_dir, "sburger_tmp", "cell_types_use.txt")
    }

    if (is.null(cell_metadata_file)) {
        analysis_dir <- dirname(data_root_dir)
        cell_metadata_file <- file.path(
            analysis_dir,
            "cellarium_upload",
            "CAP_freeze_3",
            "CAP_cell_metadata.annotated.txt.gz"
        )
    }

    if (is.null(metacells_file)) {
        metacells_file <- file.path(
            data_root_dir,
            "metacells",
            "LEVEL_3",
            "donor_rxn_DGEList_counts.tsv.gz"
        )
    }

    list(
        data_root_dir = data_root_dir,
        outDir = outDir,
        data_cache_dir = data_cache_dir,
        ct_file = ct_file,
        cell_metadata_file = cell_metadata_file,
        metacells_file = metacells_file,
        cache_rds = file.path(data_cache_dir, "donor_gex_vs_age_cache.rds")
    )
}
