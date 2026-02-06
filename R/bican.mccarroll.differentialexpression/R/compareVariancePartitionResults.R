# library(ggplot2)
# library(cowplot)

# data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metacells"
# data_name="donor_rxn_DGEList"
# cell_type="microglia"
# runOneRDS="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/variance_partition_results/microglia_variance_partition.rds"
# runTwoRDS="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/variance_partition_results_no_village/microglia_variance_partition.rds"
# runOneName="With Village"
# runTwoName="No Village"
# outPDF="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/variance_partition_results/microglia_variance_partition_with_village_comparison.pdf"
#
# data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metacells"
# data_name="donor_rxn_DGEList"
# cell_type="microglia"
# runOneRDS="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/variance_partition_results/microglia_variance_partition.rds"
# runTwoRDS="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/variance_partition_results_no_hbcac_biobank/microglia_variance_partition.rds"
# runOneName="With HBCAC/BioBank"
# runTwoName="Without"
# outPDF="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/variance_partition_results/microglia_variance_partition_with_hbcac_biobank_comparison.pdf"

# data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metacells"
# data_name="donor_rxn_DGEList"
# cell_type="microglia"
# runOneRDS="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/variance_partition_results/microglia_variance_partition.rds"
# runTwoRDS="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/variance_partition_results_no_hbcac/microglia_variance_partition.rds"
# runOneName="With HBCAC"
# runTwoName="Without"
# outPDF="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/variance_partition_results/microglia_variance_partition_no_hbcac_comparison.pdf"

# data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metacells"
# data_name="donor_rxn_DGEList"
# cell_type="microglia"
# runOneRDS="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/variance_partition_results/microglia_variance_partition.rds"
# runTwoRDS="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/variance_partition_results_no_hbcac_no_pmi/microglia_variance_partition.rds"
# runOneName="Default"
# runTwoName="Without HBCAC and PMI"
# outPDF="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/variance_partition_results/microglia_variance_partition_no_hbcac_pmi_comparison.pdf"

#' Compare two variance partition runs where different sets of covariates were used.
#'
#' @param runOneRDS Path to first variance partition RDS file.
#' @param runTwoRDS Path to second variance partition RDS file.
#' @param runOneName Name for the first run, used in plots.
#' @param runTwoName Name for the second run, used in plots.
#' @param outPDF Optional path to output PDF file for plots. If NULL, plots are printed to the console.
#' @return None. Generates and optionally saves plots comparing the two variance partition runs.
#' @noRd
compareTwoVariancePartitionRuns<-function (data_dir, data_name, cell_type, runOneRDS, runTwoRDS, runOneName, runTwoName, outPDF=NULL) {
    a=readRDS(runOneRDS)
    b=readRDS(runTwoRDS)
    both=intersect(rownames(a), rownames(b))
    a=a[both,]
    b=b[both,]

    corPlot=plotCorrelationBarPlot(a, b, runOneName, runTwoName)
    scatterPlots=plot_scatter_correlations(a, b, runOneName, runTwoName, plots_per_page = 4)

    if (!is.null(outPDF)) {
        pdf(outPDF)

        print(corPlot)
        for (page in scatterPlots) {
            print(page)
        }
        dev.off()
    } else {
        print(corPlot)
        for (page in scatterPlots) {
            print(page)
        }
    }

    #adhoc
    genes=rownames(a)[which(a$imputed_sex>0.05 & b$imputed_sex<0.05)]
    excludeColumns=c("PC1", "PC2", "PC3", "PC4", "PC5", "pct_intronic", "pmi_hr", "region")

    all_vars <- union(colnames(a), colnames(b))
    all_vars <- setdiff(all_vars, excludeColumns)

    # Generate a color palette (e.g., Glasbey or RColorBrewer or manually)
    palette_colors = c(variancePartition::ggColorHue(length(all_vars) - 1), "grey85")
    names(palette_colors) <- c(all_vars)

    # Variables used in the plots
    vars_a <- setdiff(colnames(a), excludeColumns)
    vars_b <- setdiff(colnames(b), excludeColumns)

    # Subset color map to match columns in each dataset
    colors_a <- palette_colors[vars_a]
    colors_b <- palette_colors[vars_b]

    #plot with consistent colors
    variancePartition::plotPercentBars(a[genes, vars_a], col = colors_a) +
        ggplot2::ggtitle("Include BioBank and HBCAC")

    variancePartition::plotPercentBars(b[genes, vars_b], col = colors_b) +
        ggplot2::ggtitle("Exclude BioBank and HBCAC")


}

# plotVariableCorrelation<-function (data_dir, data_name, cell_type, a, b) {
#     if (is.null(data_dir) | is.null(data_name) | is.null (cell_type)) {
#         logger::log_warn("data_dir, data_name, or cell_type is NULL. Can't generate feature correlation heatmap.")
#         return (NULL)
#     }
#
#     dge=bican.mccarroll.differentialexpression::loadDGEList(data_dir, prefix = data_name)
#     dge=dge[,dge$samples$cell_type == cell_type]
#     required_vars=union(colnames(a), colnames(b))
#     required_vars=intersect(required_vars, colnames(dge$samples))
#
#     fString=paste(required_vars, collapse = " + ")
#     f<-stats::formula(paste("~", fString, sep=" "))
#     C<- variancePartition::canCorPairs(formula=f, data=dge$samples)
#
#     base_plot_function <- function() {
#         variancePartition::plotCorrMatrix(C)
#     }
#
#     # Capture the Base R plot using ggdraw
#     p_base <- cowplot::ggdraw(base_plot_function)
#
#     return (p_base)
#
# }

plotCorrelationBarPlot<-function (a, b, runOneName, runTwoName) {
    # Get the list of common columns
    common_cols <- intersect(colnames(a), colnames(b))

    # Compute correlations between corresponding columns
    correlations <- sapply(common_cols, function(col) {
        stats::cor(a[[col]], b[[col]], use = "complete.obs", method = "pearson")
    })

    # Convert to a dataframe for plotting
    cor_df <- data.frame(Column = common_cols, Correlation = correlations)

    # Sort columns by correlation for visual clarity
    cor_df$Column <- factor(cor_df$Column, levels = cor_df$Column[order(cor_df$Correlation)])

    strTitle=paste("Variance Partition correlation [", runOneName, "] vs [", runTwoName, "]", sep="")

    # Make R CMD CHECK Happy
    Column<-Correlation<-NULL

    p<- ggplot(cor_df, aes(x = Column, y = Correlation)) +
        geom_col(fill = "grey70") +
        geom_text(aes(label = round(Correlation, 3)), vjust = -0.5, size = 3) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(size = 14, face = "bold")) +
        labs(title = strTitle)
    return (p)
}


plot_scatter_correlations <- function(a, b, runOneName, runTwoName, plots_per_page = 4) {
    common_cols <- intersect(colnames(a), colnames(b))

    scatter_plots <- list()

    # Generate ggplot scatter plots with title as column name + correlation
    for (col in common_cols) {
        corr_val <- stats::cor(a[[col]], b[[col]], use = "complete.obs")
        p <- ggplot(data.frame(a = a[[col]], b = b[[col]]), aes(x = a, y = b)) +
            geom_point(alpha = 0.5) +
            geom_abline(intercept = 0, slope = 1, color = "red") +
            labs(title = paste0(col, " (r = ", round(corr_val, 3), ")"),
                 x = runOneName, y = runTwoName) +
            theme_minimal()
        scatter_plots[[col]] <- p
    }

    # Paginate into cowplot grids
    paginated_pages <- list()
    total_plots <- length(scatter_plots)
    num_pages <- ceiling(total_plots / plots_per_page)

    plot_names <- names(scatter_plots)

    for (i in seq_len(num_pages)) {
        start_idx <- (i - 1) * plots_per_page + 1
        end_idx <- min(i * plots_per_page, total_plots)
        page_plots <- scatter_plots[start_idx:end_idx]

        # Pad with blank plots if needed
        while (length(page_plots) < plots_per_page) {
            page_plots <- c(page_plots, list(ggplot() + theme_void()))
        }

        # Create the aggregated page layout
        page_plot <- cowplot::plot_grid(plotlist = page_plots, ncol = 2)
        paginated_pages[[i]] <- page_plot
    }

    return(paginated_pages)
}

# Usage:
# cor_summary <- plot_scatter_correlations(a, b)



