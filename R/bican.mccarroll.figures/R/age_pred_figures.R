#' Manuscript figure: residual correlation heatmap + Jaccard overlap heatmap
#'
#' Create a single SVG with two heatmaps on the same row:
#' (1) Jaccard overlap of aging programs (ordered to match the residual correlation heatmap)
#' (2) Residual correlation heatmap (clustered; provides the ordering)
#'
#' Cell types are restricted to those listed in `cellTypeListFile` (one cell type per line,
#' no header). Data are fetched via `get_age_prediction_results()` and cached as usual.
#'
#' @param region Character scalar region name (e.g., "CaH").
#' @param cellTypeListFile Path to a text file containing one cell type per line (no header).
#' @param data_dir See `get_age_prediction_results()`.
#' @param data_name See `get_age_prediction_results()`.
#' @param age_de_results_dir See `get_age_prediction_results()`.
#' @param contig_yaml_file See `get_age_prediction_results()`.
#' @param reduced_gtf_file See `get_age_prediction_results()`.
#' @param n_cores See `get_age_prediction_results()`.
#' @param data_cache_dir Cache root directory used by `get_age_prediction_results()`.
#' @param outDir Output directory for the SVG.
#'
#' @return Invisibly returns a list with components `corr_out`, `jaccard_out`, and `out_svg`.
#' @export
age_prediction_residual_corr_and_jaccard_heatmaps_region <- function(
        region = "CaH",
        cellTypeListFile = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/metadata/cell_types_for_de_filtering_plot.txt",
        data_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/metacells/LEVEL_3",
        data_name = "donor_rxn_DGEList",
        age_de_results_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results/old/LEVEL_3/sex_age/cell_type_region_interaction_absolute_effects",
        contig_yaml_file = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/metadata/GRCh38_ensembl_v43.contig_groups.yaml",
        reduced_gtf_file = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/metadata/GRCh38_ensembl_v43.reduced.gtf.gz",
        n_cores = 14,
        data_cache_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/figure_repository/data_cache",
        outDir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/figure_repository") {

    cell_types <- read.table(cellTypeListFile, header=F)$V1

    results <- get_age_prediction_results(
        data_dir = data_dir,
        data_name = data_name,
        age_de_results_dir = age_de_results_dir,
        contig_yaml_file = contig_yaml_file,
        reduced_gtf_file = reduced_gtf_file,
        n_cores = n_cores,
        data_cache_dir = data_cache_dir
    )

    model_predictions <- results$donor_predictions
    all_models <- results$model_coefficients

    model_predictions <- model_predictions[
        model_predictions$region == region &
            model_predictions$cell_type %in% cell_types,
    ]

    all_models <- all_models[
        all_models$region == region &
            all_models$cell_type %in% cell_types,
    ]

    corr_title <- sprintf(
        "Age Prediction residuals (predicted - actual)\nregion [%s]",
        region
    )

    jac_title <- sprintf(
        "Cell-type gene overlap in aging programs\nregion [%s]",
        region
    )

    corr_out <- bican.mccarroll.differentialexpression::plot_residual_corr_heatmap(
        model_predictions = model_predictions,
        mode = "within_region",
        region = region,
        value_var = "resid_mean_corrected",
        title = corr_title,
        annotate_cells = TRUE,
        row_fontsize = 10,
        col_fontsize = 10,
        cell_fontsize = 9
    )

    jaccard_out <- bican.mccarroll.differentialexpression::plot_jaccard_overlap_heatmap(
        all_models = all_models,
        mode = "within_region",
        region = region,
        title = jac_title,
        coef_thresh = 0,
        annotate_cells = TRUE,
        row_fontsize = 10,
        col_fontsize = 10,
        cell_fontsize = 9,
        row_order_names = corr_out$row_order_names,
        column_order_names = corr_out$column_order_names
    )

    #Yucky way to plot both together.
    # heatmaps are objects:
    #   corr_out$heatmap
    #   jaccard_out$heatmap

    g_corr <- grid::grid.grabExpr(
        ComplexHeatmap::draw(
            corr_out$heatmap,
            newpage = FALSE,
            heatmap_legend_side = "right"
        )
    )

    g_jac <- grid::grid.grabExpr(
        ComplexHeatmap::draw(
            jaccard_out$heatmap,
            newpage = FALSE,
            heatmap_legend_side = "right"
        )
    )

    p_corr <- cowplot::ggdraw(g_corr)
    p_jac  <- cowplot::ggdraw(g_jac)

    final <- cowplot::plot_grid(
        p_jac, p_corr,
        nrow = 1,
        rel_widths = c(1, 1)
    )

    # Then ggsave(final) or draw on svg device

    out_svg <- file.path(
        outDir,
        sprintf("age_prediction_residual_corr_and_jaccard_region_%s.svg", region)
    )

    grDevices::svg(filename = out_svg, width = 14, height = 7)
    print(final)
    grDevices::dev.off()
}




#' Manuscript figure: pairwise corrected residual scatterplots within a region
#'
#' Generate a single SVG containing all pairwise scatterplots of corrected donor
#' residuals across a specified list of cell types within a fixed region.
#' Points correspond to donors present in both groups and are colored by donor age.
#'
#' @param cell_type_list Character vector of cell types to include.
#' @param region Character scalar region name (e.g., "CaH").
#' @param data_dir See `get_age_prediction_results()`.
#' @param data_name See `get_age_prediction_results()`.
#' @param age_de_results_dir See `get_age_prediction_results()`.
#' @param contig_yaml_file See `get_age_prediction_results()`.
#' @param reduced_gtf_file See `get_age_prediction_results()`.
#' @param n_cores See `get_age_prediction_results()`.
#' @param data_cache_dir Cache root directory used by `get_age_prediction_results()`.
#' @param outDir Output directory for the SVG.
#' @param ncol Number of columns in the panel grid.
#'
#' @return Invisibly returns the assembled ggplot object.
#' @export
age_prediction_corrected_residual_pairwise_scatter_region <- function(
        cell_type_list = c("astrocyte", "OPC", "microglia", "MSN_D1"),
        region = "CaH",
        data_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/metacells/LEVEL_3",
        data_name = "donor_rxn_DGEList",
        age_de_results_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results/old/LEVEL_3/sex_age/cell_type_region_interaction_absolute_effects",
        contig_yaml_file = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/metadata/GRCh38_ensembl_v43.contig_groups.yaml",
        reduced_gtf_file = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/metadata/GRCh38_ensembl_v43.reduced.gtf.gz",
        n_cores = 14,
        data_cache_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/figure_repository/data_cache",
        outDir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/figure_repository",
        ncol = 3) {

    results <- get_age_prediction_results(
        data_dir = data_dir,
        data_name = data_name,
        age_de_results_dir = age_de_results_dir,
        contig_yaml_file = contig_yaml_file,
        reduced_gtf_file = reduced_gtf_file,
        n_cores = n_cores,
        data_cache_dir = data_cache_dir
    )

    donor_pred <- results$donor_predictions

    donor_pred <- donor_pred[
        donor_pred$region == region &
            donor_pred$cell_type %in% cell_type_list,
    ]

    prs <- utils::combn(cell_type_list, 2, simplify = FALSE)

    plots <- list()
    legend_plot <- NULL

    for (pair in prs) {

        x_group <- pair[1]
        y_group <- pair[2]

        p <- bican.mccarroll.differentialexpression::plot_residual_pair_scatter_one(
            model_predictions = donor_pred,
            mode = "within_region",
            region = region,
            x_group = x_group,
            y_group = y_group,
            value_var = "resid_mean_corrected",
            color_var = "age",
            color_title = "Donor age"
        )

        p <- p +
            ggplot2::labs(
                title = paste0(x_group, " vs ", y_group),
                x = NULL,
                y = NULL
            ) +
            ggplot2::theme(
                legend.position = "none",
                plot.title = ggplot2::element_text(hjust = 0, size = 10),
                axis.title.x = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_blank()
            )

        if (is.null(legend_plot)) {
            legend_plot <- cowplot::get_legend(
                p + ggplot2::theme(legend.position = "right")
            )
        }

        plots[[length(plots) + 1]] <- p
    }

    grid <- cowplot::plot_grid(plotlist = plots, ncol = ncol, align = "hv")
    grid <- cowplot::plot_grid(grid, legend_plot, ncol = 2, rel_widths = c(1, 0.18))

    header_text <- sprintf(
        "Corrected residual pairwise scatterplots | Region: %s",
        region
    )

    core <- cowplot::plot_grid(
        cowplot::ggdraw() +
            cowplot::draw_label(header_text, x = 0, hjust = 0, fontface = "bold", size = 12),
        grid,
        cowplot::ggdraw() +
            cowplot::draw_label("Corrected residual (predicted - actual)", size = 11),
        ncol = 1,
        rel_heights = c(0.12, 1, 0.12)
    )

    # Fix A: reserve space on left for shared Y label
    left_pad <- 0.03

    final <- cowplot::ggdraw() +
        cowplot::draw_plot(core, x = left_pad, y = 0, width = 1 - left_pad, height = 1) +
        cowplot::draw_label(
            "Corrected residual (predicted - actual)",
            angle = 90,
            x = left_pad * 0.35,
            y = 0.5,
            vjust = 0.5,
            size = 11
        )

    left_pad <- 0.03

    final_padded <- cowplot::ggdraw() +
        cowplot::draw_plot(
            final,
            x = left_pad,
            y = 0,
            width = 1 - left_pad,
            height = 1
        )

    out_svg <- file.path(
        outDir,
        sprintf("age_prediction_corrected_residual_pairwise_scatter_region_%s.svg", region)
    )

    ggplot2::ggsave(
        filename = out_svg,
        plot = final_padded,
        device = "svg",
        width = 14,
        height = 8
    )

    invisible(final)
}

#' Manuscript figure: pairwise uncorrected residual scatterplots within a region
#'
#' Generate a single SVG containing all pairwise scatterplots of corrected donor
#' residuals across a specified list of cell types within a fixed region.
#' Points correspond to donors present in both groups and are colored by donor age.
#'
#' @param cell_type_list Character vector of cell types to include.
#' @param region Character scalar region name (e.g., "CaH").
#' @param data_dir See `get_age_prediction_results()`.
#' @param data_name See `get_age_prediction_results()`.
#' @param age_de_results_dir See `get_age_prediction_results()`.
#' @param contig_yaml_file See `get_age_prediction_results()`.
#' @param reduced_gtf_file See `get_age_prediction_results()`.
#' @param n_cores See `get_age_prediction_results()`.
#' @param data_cache_dir Cache root directory used by `get_age_prediction_results()`.
#' @param outDir Output directory for the SVG.
#' @param ncol Number of columns in the panel grid.
#'
#' @return Invisibly returns the assembled ggplot object.
#' @export
age_prediction_uncorrected_residual_pairwise_scatter_region <- function(
        cell_type_list = c("astrocyte", "OPC", "microglia", "MSN_D1"),
        region = "CaH",
        data_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/metacells/LEVEL_3",
        data_name = "donor_rxn_DGEList",
        age_de_results_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results/old/LEVEL_3/sex_age/cell_type_region_interaction_absolute_effects",
        contig_yaml_file = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/metadata/GRCh38_ensembl_v43.contig_groups.yaml",
        reduced_gtf_file = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/metadata/GRCh38_ensembl_v43.reduced.gtf.gz",
        n_cores = 14,
        data_cache_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/figure_repository/data_cache",
        outDir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/figure_repository",
        ncol = 3) {

    results <- get_age_prediction_results(
        data_dir = data_dir,
        data_name = data_name,
        age_de_results_dir = age_de_results_dir,
        contig_yaml_file = contig_yaml_file,
        reduced_gtf_file = reduced_gtf_file,
        n_cores = n_cores,
        data_cache_dir = data_cache_dir
    )

    donor_pred <- results$donor_predictions

    donor_pred <- donor_pred[
        donor_pred$region == region &
            donor_pred$cell_type %in% cell_type_list,
    ]

    prs <- utils::combn(cell_type_list, 2, simplify = FALSE)

    plots <- list()
    legend_plot <- NULL

    for (pair in prs) {

        x_group <- pair[1]
        y_group <- pair[2]

        p <- bican.mccarroll.differentialexpression::plot_residual_pair_scatter_one(
            model_predictions = donor_pred,
            mode = "within_region",
            region = region,
            x_group = x_group,
            y_group = y_group,
            value_var = "resid_mean",
            color_var = "age",
            color_title = "Donor age"
        )

        p <- p +
            ggplot2::labs(
                title = paste0(x_group, " vs ", y_group),
                x = NULL,
                y = NULL
            ) +
            ggplot2::theme(
                legend.position = "none",
                plot.title = ggplot2::element_text(hjust = 0, size = 10),
                axis.title.x = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_blank()
            )

        if (is.null(legend_plot)) {
            legend_plot <- cowplot::get_legend(
                p + ggplot2::theme(legend.position = "right")
            )
        }

        plots[[length(plots) + 1]] <- p
    }

    grid <- cowplot::plot_grid(plotlist = plots, ncol = ncol, align = "hv")
    grid <- cowplot::plot_grid(grid, legend_plot, ncol = 2, rel_widths = c(1, 0.18))

    header_text <- sprintf(
        "Uncorrected residual pairwise scatterplots | Region: %s",
        region
    )

    core <- cowplot::plot_grid(
        cowplot::ggdraw() +
            cowplot::draw_label(header_text, x = 0, hjust = 0, fontface = "bold", size = 12),
        grid,
        cowplot::ggdraw() +
            cowplot::draw_label("Uncorrected residual (predicted - actual)", size = 11),
        ncol = 1,
        rel_heights = c(0.12, 1, 0.12)
    )

    # Fix A: reserve space on left for shared Y label
    left_pad <- 0.03

    final <- cowplot::ggdraw() +
        cowplot::draw_plot(core, x = left_pad, y = 0, width = 1 - left_pad, height = 1) +
        cowplot::draw_label(
            "Uncorrected residual (predicted - actual)",
            angle = 90,
            x = left_pad * 0.35,
            y = 0.5,
            vjust = 0.5,
            size = 11
        )

    left_pad <- 0.03

    final_padded <- cowplot::ggdraw() +
        cowplot::draw_plot(
            final,
            x = left_pad,
            y = 0,
            width = 1 - left_pad,
            height = 1
        )

    out_svg <- file.path(
        outDir,
        sprintf("age_prediction_uncorrected_residual_pairwise_scatter_region_%s.svg", region)
    )

    ggplot2::ggsave(
        filename = out_svg,
        plot = final_padded,
        device = "svg",
        width = 14,
        height = 8
    )

    invisible(final)
}





age_prediction_examples <- function(
        cell_type_list = c("astrocyte", "OPC", "microglia", "MSN_D1"),
        region = "CaH",
        data_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/metacells/LEVEL_3",
        data_name = "donor_rxn_DGEList",
        age_de_results_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results/old/LEVEL_3/sex_age/cell_type_region_interaction_absolute_effects",
        contig_yaml_file = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/metadata/GRCh38_ensembl_v43.contig_groups.yaml",
        reduced_gtf_file = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/metadata/GRCh38_ensembl_v43.reduced.gtf.gz",
        n_cores = 14,
        data_cache_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/figure_repository/data_cache",
        outDir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/figure_repository"
) {

    results <- get_age_prediction_results(
        data_dir = data_dir,
        data_name = data_name,
        age_de_results_dir = age_de_results_dir,
        contig_yaml_file = contig_yaml_file,
        reduced_gtf_file = reduced_gtf_file,
        n_cores = n_cores,
        data_cache_dir = data_cache_dir
    )

    donor_pred <- results$donor_predictions
    gam_fit <- results$gam_fit

    make_one_plot <- function(cell_type) {

        dp <- donor_pred[donor_pred$cell_type == cell_type & donor_pred$region == region, ]
        gf <- gam_fit[gam_fit$cell_type == cell_type & gam_fit$region == region, ]

        bican.mccarroll.differentialexpression::plot_mc_donor_predictions(
            donor_predictions = dp,
            gam_fit_df = gf,
            y_var = "pred_mean",
            color_var = "resid_mean",
            legend_title = "Residuals",
            alpha_points = 0.8,
            errorbar_width = 0.1
        ) +
            ggplot2::labs(title = cell_type, subtitle = NULL, x = NULL, y = NULL) +
            ggplot2::theme(
                plot.title = ggplot2::element_text(hjust = 0, size = 11),
                axis.title.x = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_blank()
            )
    }

    plot_list <- lapply(cell_type_list, make_one_plot)

    grid <- cowplot::plot_grid(plotlist = plot_list, ncol = 4, align = "hv")

    header_text <- sprintf("Monte Carlo cross fold donor age predictions | Region: %s", region)

    core <- cowplot::plot_grid(
        cowplot::ggdraw() +
            cowplot::draw_label(header_text, x = 0, hjust = 0, fontface = "bold", size = 12),
        grid,
        cowplot::ggdraw() +
            cowplot::draw_label("Chronological age", size = 11),
        ncol = 1,
        rel_heights = c(0.12, 1, 0.12)
    )

    # Fix A: reserve space on the left for the shared y label
    left_pad <- 0.025

    final <- cowplot::ggdraw() +
        cowplot::draw_plot(core, x = left_pad, y = 0, width = 1 - left_pad, height = 1) +
        cowplot::draw_label(
            "Predicted age (MC mean)",
            angle = 90,
            x = left_pad * 0.35,
            y = 0.5,
            vjust = 0.5,
            size = 11
        )

    if (!is.null(outDir)) {
        output_svg <- file.path(outDir, "age_prediction_cell_type_examples.svg")
        ggplot2::ggsave(filename = output_svg, plot = final, device = "svg", width = 16, height = 4)
    }

    invisible(final)
}


# Generation of raw data for all plots for age prediction.
# These results will be cached.
# This returns a list of dataframes that can be used for plotting.
get_age_prediction_results<-function (data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/metacells/LEVEL_3",
                              data_name="donor_rxn_DGEList",
                              age_de_results_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results/old/LEVEL_3/sex_age/cell_type_region_interaction_absolute_effects",
                              contig_yaml_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/metadata/GRCh38_ensembl_v43.contig_groups.yaml",
                              reduced_gtf_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/metadata/GRCh38_ensembl_v43.reduced.gtf.gz",
                              n_cores=14,
                              data_cache_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/figure_repository/data_cache"
                              ) {

    cache_dir <- file.path(data_cache_dir, "age_prediction_cache")
    if (!dir.exists(cache_dir)) {
        dir.create(cache_dir, recursive = TRUE)
    }

    #get files.  If they are null, then regenerate them and save.
    results <- read_age_prediction_results(cache_dir)

    if (!is.null(results)) {
        logger::log_info("Using cached data from {cache_dir}")
        return (results)
        #read in all the files and return
    } else { # process the data as usual, write to the cache
        #clustering_min_genes and num_clusters don't do anything when outDir is NULL
        logger::log_info("No cached data from {cache_dir} regenerating data from sources.  This can take a while")

        #this will write files to the cache dir.
        bican.mccarroll.differentialexpression::predict_age_by_celltype_region(
            data_dir = data_dir,
            data_name = data_name,
            age_de_results_dir = age_de_results_dir,
            outPDFFile=NULL,
            result_dir = cache_dir,
            contig_yaml_file = contig_yaml_file,
            reduced_gtf_file = reduced_gtf_file,
            optimize_alpha = FALSE,
            alpha_fixed = 0,
            n_cores = n_cores)

        #read in all the files and return
        return (read_age_prediction_results(cache_dir))
    }




}


read_age_prediction_results <- function(cache_dir) {

    if (!dir.exists(cache_dir)) {
        stop(sprintf(
            "Cache directory does not exist: %s. Run run_age_prediction() to generate data.",
            cache_dir
        ), call. = FALSE)
    }

    files <- list(
        donor_predictions = "age_prediction_results_donor_predictions.txt",
        gam_fit           = "age_prediction_results_gam_fit.txt",
        model_metrics     = "age_prediction_results_model_metrics.txt",
        model_coefficients     = "age_prediction_results_model_coefficients.txt"
    )

    paths <- file.path(cache_dir, unlist(files))

    missing <- !file.exists(paths)

    if (any(missing)) {
        logger::log_warn(sprintf(
            "Missing required files in %s:\n%s",
            cache_dir,
            paste(paths[missing], collapse = "\n")
        ))
        return(NULL)
    }

    # If we reach here, all files exist.
    # Any read error should propagate as a real error.
    result <- list(
        donor_predictions = data.table::fread(paths[1], sep = "\t", header = TRUE, data.table = FALSE),
        gam_fit           = data.table::fread(paths[2], sep = "\t", header = TRUE, data.table = FALSE),
        model_metrics     = data.table::fread(paths[3], sep = "\t", header = TRUE, data.table = FALSE),
        model_coefficients   = data.table::fread(paths[4], sep = "\t", header = TRUE, data.table = FALSE)
    )

    result
}
