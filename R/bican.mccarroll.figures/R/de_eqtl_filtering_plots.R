library (ggplot2)


de_filtering_plot<-function (
        data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results",
        data_name="donor_rxn_DGEList",
        data_cache_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/figure_repository/data_cache",
        cellTypeListFile="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/metadata/cell_types_for_de_filtering_plot.txt",
        outDir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/figure_repository",
        clustering_min_genes=150, num_clusters=4) {

    #if the expected cache file is present, use it instead of rebuilding the data.
    cache_file <- file.path(data_cache_dir, "de_filtering_plot_cache.txt")
    if (file.exists(cache_file)) {
        logger::log_info("Using cached data from {cache_file}")
        df=read.table(cache_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    } else { # process the data as usual, write to the cache
        #clustering_min_genes and num_clusters don't do anything when outDir is NULL
        logger::log_info("No cached data from {cache_file} regenerating data from sources")
        df=bican.mccarroll.differentialexpression::compare_all_age_de_runs (data_dir, outDir=NULL, filter_levels=c(0,1,2,3,4), fdr_cutoff=0.05, clustering_min_genes=clustering_min_genes, num_clusters=num_clusters)
        write.table(df, file=cache_file, sep="\t", row.names=FALSE, quote=FALSE)
    }

    #filter by the cell type list before clustering
    cell_types_to_use <- read.table(cellTypeListFile, header=FALSE, stringsAsFactors=FALSE)[,1]
    df <- df[df$cell_type %in% cell_types_to_use, ]
    missing=setdiff(cell_types_to_use, df$cell_type)
    if (length(missing)>0) {
        missing_list=paste(missing, collapse=', ')
        logger::log_warn("The following cell types were requested but not found in the data: {missing_list}")
    }

    df_filtered=df[df$num_genes_significant_old>=clustering_min_genes, ]
    #df_filtered=df

    # Steve has asked to include level 0 - this is a temporary hack until that's approved as the
    # best way to visualize data.
    df2=.add_baseline_comparison_level (df_filtered)
    res <- bican.mccarroll.differentialexpression::cluster_filtering_trajectories(df2, K = num_clusters)
    p<- make_filtering_cluster_panels(df2, res$clusters, legend_scale=0.8)

    if (!is.null(outDir)) {
        output_svg <- file.path(outDir, "de_filtering_plot.svg")
        ggplot2::ggsave(filename = output_svg, plot = p, device = "svg", width = 8, height = 8)
    }

    #I need a plot of the initial number of discoveries.
    z=df_filtered[df_filtered$base_level==0 & df_filtered$comparison_level==1, c("cell_type", "num_genes_significant_old")]
    p_num_de_genes<-plot_num_de_genes(z)


    p1<-res$plot_trajectories
    p1 <- p1 + ggplot2::labs(title = "Change in DE results with filtering") +
        ggplot2::xlab("Filtering level")

    p2<-res$plot_mapping

    #additional idea: instead of plotting all of the curves at once, use the K-means
    #to generate facets and plot each set of cell types within that cluster as a group, and
    #label the cell types within the facet.

    heatmap_plot<-bican.mccarroll.differentialexpression::plot_filtering_trajectories_heatmap(df2)



}

plot_num_de_genes<-function (df) {
    cell_type <- num_genes_significant_old <- NULL

    df$cell_type <- factor(
        z$cell_type,
        levels = z$cell_type[order(z$num_genes_significant_old)]
    )

    ggplot(df, aes(x = cell_type, y = num_genes_significant_old)) +
        geom_col() +
        geom_text(
            aes(label = num_genes_significant_old),
            hjust = -0.1,
            size = 3
        ) +
        scale_y_log10() +
        coord_flip() +
        labs(
            x = "cell type",
            y = "Number of age DE genes"
        ) +
        theme_bw()

}

make_filtering_cluster_panels <- function(
        df,
        clusters_df,
        ncol = 2,
        legend_position = c(0.02, 0.02),   # default: lower left
        legend_justification = c(0, 0),
        panel_titles = NULL,
        legend_scale = 1                   # NEW: scale factor for inset legend
) {
    cell_type <- comparison_level <- frac_genes_discovered <- base_level <- cluster <- NULL

    df2 <- merge(df, clusters_df, by = "cell_type", all.x = TRUE, sort = FALSE)
    df2 <- df2[df2$base_level == 0 & !is.na(df2$cluster), , drop = FALSE]
    df2$cluster <- as.factor(df2$cluster)

    y_lim <- range(df2$frac_genes_discovered, na.rm = TRUE)
    ks <- sort(unique(df2$cluster))

    # Titles
    if (is.null(panel_titles)) {
        title_map <- setNames(paste0("Cluster ", ks), ks)
    } else {
        if (!is.null(names(panel_titles)) && all(as.character(ks) %in% names(panel_titles))) {
            title_map <- panel_titles[as.character(ks)]
        } else {
            if (length(panel_titles) != length(ks)) {
                stop("panel_titles must be NULL, named for all clusters, or same length as clusters.")
            }
            title_map <- setNames(panel_titles, ks)
        }
    }

    make_one <- function(k) {
        sub <- df2[df2$cluster == k, , drop = FALSE]

        max_x <- max(sub$comparison_level, na.rm = TRUE)
        end_df <- sub[sub$comparison_level == max_x,
                      c("cell_type", "frac_genes_discovered")]
        end_df <- end_df[!is.na(end_df$frac_genes_discovered), , drop = FALSE]
        end_df <- end_df[order(-end_df$frac_genes_discovered, end_df$cell_type), , drop = FALSE]

        sub$cell_type <- factor(sub$cell_type, levels = end_df$cell_type)

        ggplot2::ggplot(
            sub,
            ggplot2::aes(
                x = comparison_level,
                y = frac_genes_discovered,
                group = cell_type,
                color = cell_type
            )
        ) +
            ggplot2::geom_hline(yintercept = 1, linewidth = 0.3) +
            ggplot2::geom_line(alpha = 0.7, linewidth = 0.6) +
            ggplot2::geom_point(alpha = 0.7, size = 1.2) +
            ggplot2::scale_x_continuous(breaks = sort(unique(df2$comparison_level))) +
            ggplot2::coord_cartesian(ylim = y_lim) +
            ggplot2::labs(
                title = title_map[[as.character(k)]],
                x = "Filtering level",
                y = "Fraction of genes discovered"
            ) +
            ggplot2::theme_bw() +
            ggplot2::theme(
                plot.title = ggplot2::element_text(hjust = 0.5),

                legend.title = ggplot2::element_blank(),
                legend.position = legend_position,
                legend.justification = legend_justification,

                legend.text = ggplot2::element_text(size = 9 * legend_scale),
                legend.key.height = ggplot2::unit(0.35 * legend_scale, "cm"),
                legend.key.width  = ggplot2::unit(0.35 * legend_scale, "cm"),
                legend.spacing.y  = ggplot2::unit(0.2 * legend_scale, "cm"),

                legend.background = ggplot2::element_rect(
                    fill = "white",
                    color = "grey70"
                )
            )
    }

    plots <- lapply(ks, make_one)
    cowplot::plot_grid(plotlist = plots, ncol = ncol, align = "hv")
}


.add_baseline_comparison_level <- function(df) {
    ## Expect exactly one non-zero comparison level per (cell_type, base_level)
    ## Baseline rows are copied from comparison_level == 1

    df_lvl1 <- df[df$comparison_level == 1, , drop = FALSE]

    baseline <- df_lvl1

    baseline$comparison_level <- 0

    baseline$logFC_correlation <- 1
    baseline$logFC_sign_agreement <- 1
    baseline$FDR_correlation <- 1
    baseline$frac_genes_discovered <- 1

    baseline$n_dropped_significant <- 0
    baseline$n_dropped_non_significant <- 0

    baseline$num_genes_significant_new <- baseline$num_genes_significant_old

    rbind(baseline, df)
}

