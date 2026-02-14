


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


age_prediction_examples<-function (
            cell_type_list=c("astrocyte", "OPC", "microglia", "MSN_D1"),
            region="CaH",
            data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/metacells/LEVEL_3",
            data_name="donor_rxn_DGEList",
            age_de_results_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results/old/LEVEL_3/sex_age/cell_type_region_interaction_absolute_effects",
            contig_yaml_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/metadata/GRCh38_ensembl_v43.contig_groups.yaml",
            reduced_gtf_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/metadata/GRCh38_ensembl_v43.reduced.gtf.gz",
            n_cores=14,
            data_cache_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/figure_repository/data_cache",
            outDir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/figure_repository") {


    results <- get_age_prediction_results(data_dir = data_dir,
                              data_name = data_name,
                              age_de_results_dir = age_de_results_dir,
                              contig_yaml_file = contig_yaml_file,
                              reduced_gtf_file = reduced_gtf_file,
                              n_cores = n_cores,
                              data_cache_dir = data_cache_dir)

    donor_pred<-results$donor_predictions
    gam_fit<-results$gam_fit

    plot_list=list()

    for (cell_type in cell_type_list) {
        dp=donor_pred[donor_pred$cell_type == cell_type & donor_pred$region == region,]
        gf=gam_fit[gam_fit$cell_type == cell_type & gam_fit$region == region,]

        p<-bican.mccarroll.differentialexpression::plot_mc_donor_predictions(
            donor_predictions=dp,
            gam_fit_df = gf,
            y_var = "pred_mean",
            color_var = "resid_mean",
            legend_title = "Residuals",
            alpha_points = 0.8,
            errorbar_width = 0.1)

        p <- p +
            ggplot2::labs(title = cell_type, subtitle = NULL, x = NULL) +
            ggplot2::theme(
                plot.title = ggplot2::element_text(hjust = 0, size = 11),
                axis.title.x = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_blank()
            )

        plot_list[[cell_type]]<-p

    }

    p<-cowplot::plot_grid(plotlist = plot_list, ncol=4)

    grid <- cowplot::plot_grid(plotlist = plot_list, ncol = 4, align = "hv")

    # Shared header + shared x label (and optional shared y label)
    header_text <- sprintf("Monte Carlo cross fold donor age predictions | Region: %s", region)

    final <- cowplot::plot_grid(
        cowplot::ggdraw() +
            cowplot::draw_label(header_text, x = 0, hjust = 0, fontface = "bold", size = 12) +
            ggplot2::theme(plot.margin = ggplot2::margin(0, 5.5, 0, 5.5)),
        grid,
        cowplot::ggdraw() +
            cowplot::draw_label("Chronological age", size = 11) +
            ggplot2::theme(plot.margin = ggplot2::margin(0, 5.5, 0, 5.5)),
        ncol = 1,
        rel_heights = c(0.12, 1, 0.12)
    )

    final <- cowplot::ggdraw() +
        cowplot::draw_plot(final, x = 0.04, y = 0, width = 0.96, height = 1) +
        cowplot::draw_label("Predicted age (MC mean)", angle = 90, x = 0, y = 0.5, vjust = 0.5)

    if (!is.null(outDir)) {
        output_svg <- file.path(outDir, "age_prediction_cell_type_examples.svg")
        ggplot2::ggsave(filename = output_svg, plot = final, device = "svg", width = 16, height = 4)
    }


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
        model_metrics     = "age_prediction_results_model_metrics.txt"
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
        model_metrics     = data.table::fread(paths[3], sep = "\t", header = TRUE, data.table = FALSE)
    )

    result
}
