# Compare old vs new CPM filtering of differential expression results

# library(ggplot2)
# library(cowplot)
# library(ggrepel)
# library(logger)
# library (ComplexHeatmap)

#
# old_data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression_old_gene_filtering/sex_age/cell_type"
# new_data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age/cell_type"
# outPDF="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age/cell_type/compare_feature_selection.pdf"
# fdr_cutoff=0.05
# cell_type="microglia"

#compare_age_de_runs(old_data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression_old_gene_filtering/sex_age/cell_type", new_data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age/cell_type", outPDF="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age/cell_type/compare_feature_selection.pdf", fdr_cutoff=0.05)

#compare all levels
#compare_all_age_de_runs(data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results", outDir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/compare_results_by_level")


#' Compare age differential expression results across filtering trajectories
#'
#' Compares age-related differential expression (DE) results across a sequence
#' of filtering stringency levels, treating the first entry of
#' \code{filter_levels} as the baseline and comparing it to each subsequent
#' level. For each adjacent comparison \eqn{(LEVEL_{base}, LEVEL_{k})}, the
#' function calls \code{compare_age_de_runs()} to compute concordance metrics
#' across cell types (e.g., logFC and FDR correlations and the fraction of
#' discoveries retained), then aggregates results across comparisons.
#'
#' After aggregating, the function:
#' \itemize{
#'   \item Plots per-cell-type discovery trajectories across comparison levels.
#'   \item Filters to cell types with at least \code{clustering_min_genes}
#'         significant genes in the baseline run.
#'   \item Clusters discovery trajectories using \code{cluster_filtering_trajectories()}
#'         (currently with \code{K = 4}).
#'   \item Writes a PDF summary containing the trajectory plots and additional
#'         diagnostics relating discovery reduction to baseline signal.
#' }
#'
#' The function assumes a standardized directory layout for each filtering level:
#' \preformatted{
#'   data_dir/LEVEL_<level>/sex_age/cell_type/
#' }
#'
#' Output is written to \code{outDir} as \code{"compare_de_age_all_levels_summary.pdf"}.
#'
#' @param data_dir Character scalar. Base directory containing DE results
#'   organized by filtering level (e.g., \code{LEVEL_0}, \code{LEVEL_1}, ...).
#' @param outDir Character scalar. Output directory where the summary PDF will
#'   be written. If \code{NULL}, no output file is written and the aggregated
#'   data frame is returned.
#' @param filter_levels Integer vector. Filtering levels to compare. The first
#'   value is treated as the baseline level, and each subsequent value is
#'   compared against it. Default is \code{c(0, 1, 2, 3, 4)}.
#' @param fdr_cutoff Numeric scalar. False discovery rate (FDR) threshold used
#'   within \code{compare_age_de_runs()} when computing per-cell-type summaries.
#'   Default is \code{0.05}.
#' @param clustering_min_genes Integer scalar. Minimum number of significant
#'   genes in the baseline run (\code{num_genes_significant_old}) required for a
#'   cell type to be included in trajectory clustering. Default is \code{100}.
#'   Only used if \code{outDir} is not \code{NULL}.
#' @param num_clusters Integer scalar. Number of clusters to use in
#'  \code{cluster_filtering_trajectories()}. Default is \code{4}.  Only used
#'  if \code{outDir} is not \code{NULL}.
#'
#' @export
compare_all_age_de_runs<-function (data_dir, outDir=NULL, filter_levels=c(0,1,2,3,4), fdr_cutoff=0.05, clustering_min_genes=100, num_clusters=4) {
    base_level=filter_levels[1]
    results=list()
    for (i in 1:(length(filter_levels)-1)) {
        comparison_level=filter_levels[i+1]

        baseline_name=paste0("LEVEL_", base_level)
        comparison_name=paste0("LEVEL_", comparison_level)

        #only emit the per-comparison PDF if outDir is specified, otherwise just return the data frame at the end.
        outPDF=NULL
        if (!is.null(outDir))
            outPDF=paste(outDir, "/compare_age_DE_", baseline_name, "_vs_", comparison_name, ".pdf", sep="")

        old_data_dir=paste(data_dir,"/", baseline_name, "/sex_age/cell_type", sep="")
        new_data_dir=paste(data_dir,"/", comparison_name, "/sex_age/cell_type", sep="")
        logger::log_info(paste0("Comparing age DE results between LEVEL ", base_level, " and LEVEL ", comparison_level, "\n"))

        z=compare_age_de_runs(old_data_dir=old_data_dir, new_data_dir=new_data_dir,
                              baseline_name=baseline_name, comparison_name=comparison_name,
                              outPDF=outPDF, outFile=NULL, fdr_cutoff=fdr_cutoff)
        df=z$df
        df$base_level=base_level
        df$comparison_level=comparison_level
        results[[i]]=df
    }

    df=do.call(rbind, results)

    if (is.null(outDir)) {
        return (df)
    }

    #outDir isn't null, write the output
    dataSummaryFile=paste(outDir,"/compare_de_age_all_levels_summary.txt", sep="")
    write.table (df, file=dataSummaryFile, sep="\t", quote=FALSE, row.names=FALSE)

    #unknown how many clusters to pick.
    df_filtered=df[df$num_genes_significant_old>=clustering_min_genes, ]
    res <- cluster_filtering_trajectories(df_filtered, K = num_clusters)

    p1<-res$plot_trajectories
    p2<-res$plot_mapping

    combined <- cowplot::plot_grid(
        p1,
        p2,
        nrow = 2
    )
    heatmap_plot<-plot_filtering_trajectories_heatmap(df_filtered)

    p2<-plot_frac_lines(df_filtered)

    # is the dropoff a result of the number of initial observations?
    # IE: we just barely had power before, and now we don't?
    # only compare the maximum level.
    level=max (df$comparison_level)
    p3<-plot_reduction_vs_initial(df, cluster_df=res$clusters, comparison_level=level)
    p4<- plot_reduction_by_parent_type(df, comparison_level=level)


    pdfSummaryFile=paste(outDir,"/compare_de_age_all_levels_summary.pdf", sep="")
    grDevices::pdf(pdfSummaryFile, width = 11, height = 11)
    on.exit(grDevices::dev.off(), add = TRUE)

    #print all the plots to the PDF, which auto-closes.
    print (combined)
    print (heatmap_plot)
    print (p3)
    print (p4)
    print (p2)

}


#' Compare age differential expression results across runs
#'
#' Compare age-related differential expression results between an old and a new
#' analysis run across multiple cell types. Generates summary plots and an
#' aggregated results table, with optional PDF and tab-delimited output.
#'
#' @param old_data_dir Directory containing age DE result files from the old run.
#' @param new_data_dir Directory containing age DE result files from the new run.
#' @param baseline_name Name of the baseline group used in the DE analysis.
#' @param comparison_name Name of the comparison group used in the DE analysis.
#' @param outPDF Path to a PDF file for saving plots. If NULL, no PDF is written.
#' @param outFile Path to a tab-delimited file for saving the summary table.
#'   If NULL, no file is written.
#' @param fdr_cutoff FDR threshold used to define significant genes.
#'
#' @return A list containing summary statistics and plots produced by
#'   \code{plot_summary()}.
#'
#' @export
compare_age_de_runs<-function (old_data_dir, new_data_dir, baseline_name, comparison_name, outPDF, outFile=NULL, fdr_cutoff=0.05) {
    suffix="_age_DE_results.txt"
    f=list.files(old_data_dir, pattern=suffix, full.names = F)
    cell_types_list=sub(suffix, "", f)

    plot_list=list()
    df_list=list()
    merged_list=list()

    for (cell_type in cell_types_list){
        logger::log_info(paste0("Comparing age DE results for cell type: ", cell_type, "\n"))
        p=compare_age_de_run(cell_type=cell_type,
                            old_data_dir=old_data_dir,
                            new_data_dir=new_data_dir,
                            baseline_name=baseline_name,
                            comparison_name=comparison_name,
                            fdr_cutoff=fdr_cutoff)

        plot_page=compose_de_comparison_plot(
            scatter_effect = p$scatter_effect,
            scatter_fdr = p$scatter_fdr,
            venn_plot_genes_tested = p$venn_plot_genes_tested,
            venn_plot_de_genes = p$venn_plot_de_genes,
            cell_type = cell_type
        )

        plot_list[[cell_type]]=plot_page
        df_list[[cell_type]]=p$df
        merged_list[[cell_type]]=p$merged
    }

    df=do.call(rbind, df_list)

    #this generates plots and the summary dataframe
    z=plot_summary(df)

    #plot the fraction discovery vs initial

    d3=df[, c("cell_type", "num_genes_significant_old", "num_genes_significant_new")]
    d3$log10_num_genes_old <- log10(d3$num_genes_significant_old + 1)
    d3$frac_genes_discovered <- d3$num_genes_significant_new /
        (d3$num_genes_significant_old + 1)
    d3$parent_type <- factor(sub("(_.*)$", "", d3$cell_type))

    #Make R CMD CHECK happy
    log10_num_genes_old<-frac_genes_discovered<-parent_type<-NULL

    p<-ggplot(d3,
           aes(x = log10_num_genes_old,
               y = frac_genes_discovered,
               color = parent_type)) +
        geom_point(size = 2) +
        geom_smooth(method = "loess",
                    se = FALSE,
                    color = "black",
                    linetype = "dashed",
                    linewidth = 0.8) +
        labs(
            x = expression(log[10]("Number of significant genes (before filtering)")),
            y = "Fraction of genes discovered (level 3)",
            title = "Reduction vs initial number of discoveries"
        ) +
        theme_bw()


    #plots
    if (!is.null(outPDF)) {
        pdf(outPDF, width=11, height=11)
        print (z$plot)
        print (p)
        for (pp in plot_list){
            print(pp)
        }
        dev.off()
    }

    if (!is.null(outFile)) {
        write.table (z$df, file=outFile, sep="\t", quote=FALSE, row.names=FALSE)
    }

    # test_de_inflation(merged_list[["microglia"]], main_prefix = "DE Analysis Age Microglia", pval_col = "P.Value_old", t_col = "t_old")

    return (z)
}


compose_de_comparison_plot<-function (scatter_effect, scatter_fdr, venn_plot_genes_tested, venn_plot_de_genes, cell_type) {
    p1=cowplot::plot_grid(
        venn_plot_genes_tested,
        venn_plot_de_genes,
        scatter_effect,
        scatter_fdr,
        rel_heights = c(0.4, 0.6),
        ncol=2
    )

    final_plot <- add_supertitle_cowplot(p = p1,
                                       title = paste("DE Gene Comparison:", cell_type),
                                       title_size = 16,
                                       rel_height_title = 0.08
    )

    return (final_plot)
}


plot_summary <- function(df) {

    ## add a column for the fraction of genes discovered in the new filtering
    # to avoid division by zero add one to denominator
    df$frac_genes_discovered <- df$num_genes_significant_new /
        (df$num_genes_significant_old + 1)

    metrics <- c(
        "logFC_correlation",
        "FDR_correlation",
        "frac_genes_discovered"
    )

    ## reshape to long format using base R
    df_long <- reshape(
        df,
        direction = "long",
        varying  = list(metrics),
        v.names  = "value",
        timevar  = "metric",
        times    = metrics
    )

    row.names(df_long) <- NULL
    df_long$metric <- factor(df_long$metric, levels = metrics)

    ## split into main metrics and frac_genes_discovered
    df_main <- df_long[df_long$metric != "frac_genes_discovered", ]
    df_frac <- df_long[df_long$metric == "frac_genes_discovered", ]

    ## z-scores for main metrics
    df_main$z <- ave(
        df_main$value,
        df_main$metric,
        FUN = function(x) {
            m  <- mean(x, na.rm = TRUE)
            sd <- stats::sd(x, na.rm = TRUE)
            if (is.na(sd) || sd == 0) rep(0, length(x)) else (x - m) / sd
        }
    )
    df_main$label <- ifelse(abs(df_main$z) > 2, df_main$cell_type, "")

    ## z-scores and labels for frac_genes_discovered
    if (nrow(df_frac) > 0) {
        m_frac  <- mean(df_frac$value, na.rm = TRUE)
        sd_frac <- stats::sd(df_frac$value, na.rm = TRUE)

        if (is.na(sd_frac) || sd_frac == 0) {
            df_frac$z <- 0
        } else {
            df_frac$z <- (df_frac$value - m_frac) / sd_frac
        }

        ## start with no labels
        df_frac$label <- ""

        ## optional: z-score based labels
        df_frac$label[df_frac$z > 2 | df_frac$z < -2] <-
            df_frac$cell_type[df_frac$z > 2 | df_frac$z < -2]

        ## always label top 3
        top_idx <- order(df_frac$value, decreasing = TRUE)[seq_len(min(3, nrow(df_frac)))]
        df_frac$label[top_idx] <- df_frac$cell_type[top_idx]

        ## always label bottom 3
        bottom_idx <- order(df_frac$value, decreasing = FALSE)[seq_len(min(3, nrow(df_frac)))]
        df_frac$label[bottom_idx] <- df_frac$cell_type[bottom_idx]
    }

    ## main plot
    #Make R CMD CHECK happy
    metric<-value<-label<-NULL

    p_main <- ggplot(df_main, aes(x = metric, y = value)) +
        geom_violin(trim = FALSE) +
        geom_jitter(width = 0.08, height = 0, size = 1, alpha = 0.6) +
        ggrepel::geom_text_repel(
            data = df_main[df_main$label != "", ],
            aes(label = label),
            max.overlaps = Inf,
            box.padding = 0.3,
            point.padding = 0.1,
            min.segment.length = 0,
            size = 2.5
        ) +
        labs(
            x = NULL,
            y = "Value",
            title = "Concordance metrics across cell types"
        ) +
        scale_x_discrete(labels = function(x) gsub("_", " ", x))+

        theme_bw()

    ## frac_genes_discovered plot
    p_frac <- ggplot(df_frac, aes(x = metric, y = value)) +
        geom_violin(trim = FALSE) +
        geom_jitter(width = 0.08, height = 0, size = 1, alpha = 0.6) +
        ggrepel::geom_text_repel(
            data = df_frac[df_frac$label != "", ],
            aes(label = label),
            max.overlaps = Inf,
            box.padding = 0.3,
            point.padding = 0.1,
            min.segment.length = 0,
            size = 2.5
        ) +
        labs(
            x = NULL,
            y = "Value",
            title = ""
        ) +
        scale_x_discrete(labels = function(x) gsub("_", " ", x)) +
        theme_bw()

    combined <- cowplot::plot_grid(
        p_main,
        p_frac,
        ncol       = 2,
        rel_widths = c(3, 1)
    )

    result <- list(df = df, plot = combined)
    return(result)
}


test_de_inflation <- function(
        de_dataframe,
        main_prefix = "DE Analysis",
        pval_col = "P.Value",
        t_col = "t"
) {
    stopifnot(is.data.frame(de_dataframe))
    stopifnot(all(c(pval_col, t_col) %in% names(de_dataframe)))

    pvals <- de_dataframe[[pval_col]]
    tvals <- de_dataframe[[t_col]]
    n <- length(pvals)

    ## 1) Genomic inflation factor (lambda)
    chi2 <- tvals^2
    lambda_gc <- median(chi2, na.rm = TRUE) / qchisq(0.5, df = 1)

    ## 2) KS test vs Uniform(0,1)
    ks_res <- ks.test(pvals, "punif")

    ## 3) Counts of small p-values
    thresholds <- c(1e-2, 1e-3, 1e-4, 1e-5)
    small_counts <- sapply(thresholds, function(thr) sum(pvals <= thr, na.rm = TRUE))

    ## -------- Log output --------
    log_info("Genomic inflation factor lambda_gc: {sprintf('%.3f', lambda_gc)}")
    log_info("KS test D={sprintf('%.3f', ks_res$statistic)}, p={formatC(ks_res$p.value, format='e', digits=2)}")

    for (i in seq_along(thresholds)) {
        thr <- thresholds[i]
        log_info("Number of genes with p <= {thr}: {small_counts[i]}")
    }

    ## -------- QQ plot data --------
    expected <- -log10((seq_len(n) - 0.5) / n)
    observed <- -log10(sort(pvals))
    dfqq <- data.frame(expected = expected, observed = observed)

    ## Annotation coordinates (upper-left)
    x_range <- range(dfqq$expected, finite = TRUE)
    y_range <- range(dfqq$observed, finite = TRUE)

    x_annot <- x_range[1] + 0.05 * diff(x_range)
    y_annot <- y_range[2] - 0.05 * diff(y_range)

    ## Annotation text
    label_lines <- c(
        sprintf("lambda = %.3f", lambda_gc),
        sprintf("KS p = %s", formatC(ks_res$p.value, format = "e", digits = 2)),
        sprintf("p<=1e-3: %d", small_counts[thresholds == 1e-3]),
        sprintf("p<=1e-4: %d", small_counts[thresholds == 1e-4]),
        sprintf("p<=1e-5: %d", small_counts[thresholds == 1e-5])
    )
    annot_label <- paste(label_lines, collapse = "\n")

    ## -------- QQ plot --------
    qqplot <- ggplot(dfqq, aes(x = expected, y = observed)) +
        geom_point(alpha = 0.6, size = 1) +
        geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 0.7) +
        annotate(
            "text",
            x = x_annot, y = y_annot,
            label = annot_label,
            hjust = 0, vjust = 1,
            size = 3
        ) +
        labs(
            title = paste(main_prefix, "QQ-plot of p-values"),
            x = "Expected -log10(p)",
            y = "Observed -log10(p)"
        ) +
        theme_bw()

    ## -------- Histogram --------
    dfhist <- data.frame(pvals = pvals)

    histplot <- ggplot(dfhist, aes(x = pvals)) +
        geom_histogram(bins = 40, color = "black", fill = "grey80") +
        labs(
            title = paste(main_prefix, "P-value histogram"),
            x = "P-value",
            y = "Count"
        ) +
        theme_bw()

    ## -------- Combined panel --------
    combined <- cowplot::plot_grid(
        qqplot,
        histplot,
        ncol = 1,
        rel_heights = c(2, 1),
        align = "v"
    )

    print(combined)

    invisible(list(
        lambda_gc = lambda_gc,
        ks = ks_res,
        small_counts = small_counts
    ))
}

#' Compare age DE results between old and strict filtering (ad hoc)
#'
#' This internal helper compares differential expression results for age
#' between an older, more permissive run and a stricter filtering run for
#' a single cell type. It:
#' \itemize{
#'   \item Reads limma/voom/DREAM result tables from the two directories.
#'   \item Merges overlapping genes and compares log-fold changes and FDR.
#'   \item Computes correlation and sign agreement for log-fold changes.
#'   \item Computes correlation for \eqn{-\log_{10}(\mathrm{FDR})}.
#'   \item Summarizes genes dropped by the strict filtering, split by
#'         significance in the old run.
#'   \item Produces three ggplot2 plots and a combined cowplot layout.
#' }
#'
#' This function is intended for interactive, ad-hoc inspection and is not
#' part of the public API.
#'
#' @param old_data_dir Character scalar. Directory containing the older,
#'   more permissive age DE result files.
#' @param new_data_dir Character scalar. Directory containing the stricter
#'   filtering age DE result files.
#' @param cell_type Character scalar. Cell type prefix used to identify
#'   the DE results file in each directory. The function looks for a file
#'   starting with \code{cell_type} and ending in
#'   \code{"_age_DE_results.txt"}.
#' @param baseline_name Character scalar. Name of the baseline filtering level
#'  (e.g., \code{"LEVEL_0"}).
#' @param comparison_name Character scalar. Name of the comparison filtering level
#'  (e.g., \code{"LEVEL_3"}).
#' @param fdr_cutoff Numeric scalar. FDR threshold used to classify genes
#'   as significant or non-significant in the old run when summarizing
#'   dropped genes. Default is \code{0.05}.
#'
#' @return A combined plot.
#' @export
compare_age_de_run <- function (cell_type,
                                old_data_dir,
                                new_data_dir,
                                baseline_name=baseline_name,
                                comparison_name=comparison_name,
                                fdr_cutoff = 0.05) {

    find_de_file <- function(data_dir, cell_type) {
        pattern <- paste0("^", cell_type, "_age_DE_results\\.txt$")
        files <- list.files(data_dir, pattern = pattern, full.names = TRUE)
        if (length(files) == 0L) {
            stop("No DE file found in ", data_dir, " for cell_type = ", cell_type)
        }
        if (length(files) > 1L) {
            warning("Multiple files found for cell_type; using first: ", basename(files[1L]))
        }
        files[1L]
    }

    read_de_file <- function(filepath) {
        x <- read.table(filepath,
                        header = TRUE,
                        sep = "\t",
                        row.names = 1,
                        stringsAsFactors = FALSE,
                        check.names = FALSE)
        x$gene <- rownames(x)
        rownames(x) <- NULL
        x
    }

    old_file <- find_de_file(old_data_dir, cell_type)
    new_file <- find_de_file(new_data_dir, cell_type)

    old_de <- read_de_file(old_file)
    new_de <- read_de_file(new_file)

    ## Merge overlapping genes for effect/FDR comparison
    merged <- merge(old_de,
                    new_de,
                    by = "gene",
                    suffixes = c("_old", "_new"),
                    all = FALSE)

    #add the cell type
    merged<-cbind(cell_type=cell_type, merged)

    ## Dropped genes and their old significance status
    genes_old <- old_de$gene
    genes_new <- new_de$gene
    genes_dropped <- setdiff(genes_old, genes_new)

    old_sig <- old_de$adj.P.Val < fdr_cutoff

    sig_old_genes     <- old_de$gene[old_sig]
    non_sig_old_genes <- old_de$gene[!old_sig]

    sig_old_dropped     <- intersect(sig_old_genes, genes_dropped)
    non_sig_old_dropped <- intersect(non_sig_old_genes, genes_dropped)

    n_sig_dropped     <- length(sig_old_dropped)
    n_non_sig_dropped <- length(non_sig_old_dropped)

    ## ---------- stats for logFC plot ----------
    logfc_cor <- cor(merged$logFC_old, merged$logFC_new, use = "complete.obs")
    sign_agree <- mean(sign(merged$logFC_old) == sign(merged$logFC_new))

    ann_text_effect <- paste0(
        "r = ", sprintf("%.3f", logfc_cor), "\n",
        "sign agree = ", sprintf("%.1f", 100 * sign_agree), "%"
    )

    x_range_logfc <- range(merged$logFC_old, na.rm = TRUE)
    y_range_logfc <- range(merged$logFC_new, na.rm = TRUE)

    x_pos_logfc <- x_range_logfc[1L] + 0.02 * diff(x_range_logfc)
    y_pos_logfc <- y_range_logfc[2L] - 0.02 * diff(y_range_logfc)

    ## Plot 1: effect size comparison
    lim_logfc <- range(c(merged$logFC_old, merged$logFC_new), na.rm = TRUE)

    #Make R CMD CHECK happy
    logFC_old<-logFC_new<-NULL

    p_effect <- ggplot(merged,
                       aes(x = logFC_old, y = logFC_new)) +
        geom_point(alpha = 0.4, size = 0.7) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        annotate("text",
                 x = x_pos_logfc,
                 y = y_pos_logfc,
                 hjust = 0,
                 vjust = 1,
                 label = ann_text_effect,
                 size = 4) +
        labs(
            title = paste0("Age logFC: ", cell_type),
            x = paste0("logFC ", baseline_name, ""),
            y = paste0("logFC ", comparison_name, "")
        ) +
        coord_fixed(xlim = lim_logfc, ylim = lim_logfc) +
        theme_bw()

    ## ---------- stats and plot for FDR comparison ----------
    merged$neg_log10_fdr_old <- -log10(merged$adj.P.Val_old)
    merged$neg_log10_fdr_new <- -log10(merged$adj.P.Val_new)

    num_genes_significant_old=sum(merged$adj.P.Val_old < fdr_cutoff, na.rm=TRUE)
    num_genes_significant_new=sum(merged$adj.P.Val_new < fdr_cutoff, na.rm=TRUE)
    strTitleSuffix=paste(num_genes_significant_old, " significant ", baseline_name, "; ", num_genes_significant_new, " ", comparison_name, " ", sep="")

    fdr_cor <- cor(merged$neg_log10_fdr_old,
                   merged$neg_log10_fdr_new,
                   use = "complete.obs")

    ann_text_fdr <- paste0(
        "r = ", sprintf("%.3f", fdr_cor)
    )

    x_range_fdr <- range(merged$neg_log10_fdr_old, na.rm = TRUE)
    y_range_fdr <- range(merged$neg_log10_fdr_new, na.rm = TRUE)

    lim_fdr <- range(c(merged$neg_log10_fdr_old, merged$neg_log10_fdr_new),
                     na.rm = TRUE)

    x_pos_fdr <- lim_fdr[1L] + 0.02 * diff(lim_fdr)
    y_pos_fdr <- lim_fdr[2L] - 0.02 * diff(lim_fdr)

    #Make R CMD CHECK happy
    neg_log10_fdr_old<-neg_log10_fdr_new<-NULL

    p_fdr <- ggplot(merged,
                    aes(x = neg_log10_fdr_old, y = neg_log10_fdr_new)) +
        geom_point(alpha = 0.4, size = 0.7) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        annotate("text",
                 x = x_pos_fdr,
                 y = y_pos_fdr,
                 hjust = 0,
                 vjust = 1,
                 label = ann_text_fdr,
                 size = 4) +
        labs(
            title = paste0("Age FDR: ", cell_type, "\n", strTitleSuffix),
            x = paste0("log10(BH-adjusted p-value) ", baseline_name, ""),
            y = paste0("log10(BH-adjusted p-value) ", comparison_name, "")
        ) +
        coord_fixed(xlim = lim_fdr, ylim = lim_fdr) +
        theme_bw()

    ## Plot 3: barplot of dropped genes by old significance status
    drop_counts <- data.frame(
        status = c("non-significant_old", "significant_old"),
        count  = c(n_non_sig_dropped, n_sig_dropped),
        stringsAsFactors = FALSE
    )

    drop_counts$percent <- with(drop_counts, count / c(
        sum(old_de$adj.P.Val >= fdr_cutoff),
        sum(old_de$adj.P.Val <  fdr_cutoff)
    ))

    max_y <- max(drop_counts$count)

    #Make R CMD CHECK happy
    status<-count<-percent<-NULL

    p_drop_bar <- ggplot(drop_counts,
                         aes(x = status, y = count)) +
        geom_col(fill = "lightblue") +
        geom_text(
            aes(
                y = max_y * 0.5,
                label = paste0(count, " (", round(100 * percent, 1), "%)")
            ),
            size = 6
        ) +
        labs(
            title = "Genes removed by strict filtering",
            x = "Status in old run",
            y = "Number of genes \n not tested in strict filtering"
        ) +
        theme_bw(base_size = 14) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

    # add in the venn diagrams to replicate the eQTL analysis.
    z=plot_gene_venn(old_de, new_de, text_size=6, title_size=16,
                     baseline_name=baseline_name, comparison_name=comparison_name,
                     title="Number of genes tested")

    venn_plot_genes_tested=z$plot

    #how many eQTLs are found in both levels (FDR<0.05)
    old_de_fdr=old_de[old_de$adj.P.Val < fdr_cutoff, ]
    new_de_fdr=new_de[new_de$adj.P.Val < fdr_cutoff, ]

    z=plot_gene_venn(old_de_fdr, new_de_fdr, text_size=6, title_size=16,
                     baseline_name=baseline_name, comparison_name=comparison_name,
                     title="Number of DE genes discovered")

    venn_plot_de_genes=z$plot


    df=data.frame(cell_type=cell_type,
               logFC_correlation=logfc_cor,
               logFC_sign_agreement=sign_agree,
               FDR_correlation=fdr_cor,
               n_dropped_significant=n_sig_dropped,
               n_dropped_non_significant=n_non_sig_dropped,
               num_genes_significant_old=num_genes_significant_old,
               num_genes_significant_new=num_genes_significant_new)


    result=list(df=df, scatter_effect=p_effect, scatter_fdr=p_fdr, venn_plot_genes_tested=venn_plot_genes_tested, venn_plot_de_genes=venn_plot_de_genes)
    return (result)
}


#' Venn diagram comparing overlap between two eQTL result sets
#'
#' Creates a two-set Venn diagram showing overlap of items defined by
#' \code{gene_name} (or any pre-filtered subset of rows, e.g. FDR < 0.05).
#' Also computes union and intersection counts for tracking across runs.
#'
#' @param old_de data.frame for the baseline run.
#'   Must include a \code{gene_name} column.
#' @param new_de data.frame for the comparison run.
#'   Must include a \code{gene_name} column.
#' @param baseline_name character scalar. Label for the baseline dataset.
#' @param comparison_name character scalar. Label for the comparison dataset.
#' @param title character scalar. Plot title.
#' @param text_size numeric scalar. Text size for counts and percentages in the Venn.
#' @param title_size numeric scalar. Title font size.
#'
#' @return A list with:
#' \describe{
#' \item{\code{plot}}{A ggplot object (the Venn diagram).}
#' \item{\code{stats}}{A one-row data.frame with union and intersection counts.}
#' }
#'
#' @importFrom ggvenn ggvenn
#' @importFrom ggplot2 ggtitle theme element_text
plot_gene_venn <- function(old_de,
                           new_de,
                           baseline_name = "baseline",
                           comparison_name = "comparison",
                           title = "Number of genes tested",
                           text_size = 6,
                           title_size = 16) {
    gene_name <- NULL

    genes_baseline <- unique(old_de[["gene"]])
    genes_comp <- unique(new_de[["gene"]])

    n_intersect <- length(intersect(genes_baseline, genes_comp))
    n_union <- length(union(genes_baseline, genes_comp))

    stats <- data.frame(
        baseline_name = baseline_name,
        comparison_name = comparison_name,
        n_union = n_union,
        n_intersect = n_intersect,
        stringsAsFactors = FALSE
    )

    sets <- list(
        baseline = genes_baseline,
        comparison = genes_comp
    )
    names(sets) <- c(baseline_name, comparison_name)

    p <- ggvenn::ggvenn(
        sets,
        fill_color = c("#0072B2", "#009E73"),
        text_size = text_size
    ) +
        ggplot2::ggtitle(title) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(
                size = title_size,
                hjust = 0.5
            )
        )

    list(plot = p, stats = stats)
}

#Can we cluster cell types by how they lose data as we increase filtering stringency?
plot_frac_lines <- function(df) {

    ## Sort by comparison_level for correct line drawing
    df_sorted <- df[order(df$cell_type, df$comparison_level), ]

    #Make R CMD CHECK Happy
    comparison_level<-frac_genes_discovered<-cell_type<-NULL

    ggplot(df_sorted,
           aes(x = comparison_level,
               y = frac_genes_discovered,
               group = cell_type)) +
        geom_line(alpha = 0.6) +
        geom_point(size = 1) +
        scale_x_continuous(breaks = sort(unique(df$comparison_level))) +
        labs(
            x = "Comparison level",
            y = "Fraction of genes discovered",
            title = "Filtering response per cell type"
        ) +
        theme_bw()
}

#' Cluster filtering-response trajectories
#'
#' Given per-cell-type discovery fractions across filtering comparison levels,
#' cluster the discovery trajectories using k-means and return convenience objects
#' for plotting and downstream summaries.
#'
#' @param df Data frame with columns `cell_type`, `comparison_level`, and
#'   `frac_genes_discovered`. Each row represents one cell type at one comparison level.
#' @param K Integer scalar. Number of k-means clusters.
#'
#' @return A list with elements:
#'   \itemize{
#'     \item `plot_trajectories`: ggplot of per-cell-type trajectories with cluster mean
#'       trajectories overlaid.
#'     \item `plot_mapping`: ggplot mapping each cell type to its cluster.
#'     \item `clusters`: data frame mapping `cell_type` to `cluster`.
#'     \item `mean_trajectories`: data frame of mean trajectories per cluster.
#'     \item `wide_matrix`: numeric matrix used for clustering (rows = cell types,
#'       columns = comparison levels).
#'   }
#'
#' @export
cluster_filtering_trajectories <- function(df, K = 4) {
    ## df must have: cell_type, comparison_level, frac_genes_discovered

    ## Wide matrix: rows = cell types, columns = comparison levels
    wide <- reshape(
        df[, c("cell_type", "comparison_level", "frac_genes_discovered")],
        idvar   = "cell_type",
        timevar = "comparison_level",
        direction = "wide"
    )

    ## Order columns by numeric comparison_level
    col_idx <- order(as.numeric(sub("frac_genes_discovered\\.", "", names(wide)[-1])))
    wide <- wide[, c(1, 1 + col_idx)]

    mat_numeric <- as.matrix(wide[, -1])
    rownames(mat_numeric) <- wide$cell_type

    ## K-means clustering on trajectories
    set.seed(1)
    km <- kmeans(mat_numeric, centers = K)

    ## Mean trajectory per cluster (rows = clusters, cols = comparison levels)
    cl_means <- rowsum(mat_numeric, group = km$cluster) /
        as.vector(table(km$cluster))

    ## Order clusters by average mean fraction across levels (descending)
    cl_order <- order(rowMeans(cl_means), decreasing = TRUE)

    ## Remap old cluster labels to ordered 1..K
    old_to_new <- setNames(seq_along(cl_order), cl_order)
    new_cluster <- old_to_new[as.character(km$cluster)]

    ## Cluster df
    cluster_df <- data.frame(
        cell_type = wide$cell_type,
        cluster   = new_cluster,
        stringsAsFactors = FALSE
    )
    cluster_df$cluster <- factor(cluster_df$cluster, levels = seq_len(K))

    ## Reorder cl_means to match new cluster labels
    cl_means_ord <- cl_means[cl_order, , drop = FALSE]

    ## Long df of cluster mean trajectories for plotting
    comp_levels <- as.numeric(sub("frac_genes_discovered\\.", "", colnames(cl_means_ord)))
    mean_long <- data.frame(
        cluster             = rep(seq_len(K), each = ncol(cl_means_ord)),
        comparison_level    = rep(comp_levels, times = K),
        frac_genes_mean     = as.vector(t(cl_means_ord))
    )
    mean_long$cluster <- factor(mean_long$cluster, levels = seq_len(K))

    ## Merge cluster assignment back into original df
    df2 <- merge(df, cluster_df, by = "cell_type")

    ## Order for nice lines
    df2 <- df2[order(df2$cluster, df2$cell_type, df2$comparison_level), ]

    ## Trajectory plot with cluster means (thick dashed)

    #Make R CMD CHECK Happy
    comparison_level <- frac_genes_discovered <- cluster <- cell_type <- frac_genes_mean <- NULL

    p_trajectories <- ggplot(
        df2,
        aes(
            x = comparison_level,
            y = frac_genes_discovered,
            group = cell_type,
            color = cluster
        )
    ) +
        geom_line(alpha = 0.7) +
        geom_point(size = 1) +
        ## cluster mean trajectories
        geom_line(
            data = mean_long,
            aes(
                x = comparison_level,
                y = frac_genes_mean,
                group = cluster,
                color = cluster
            ),
            inherit.aes = FALSE,
            linetype = "dashed",
            linewidth = 1.1
        ) +
        scale_x_continuous(breaks = sort(unique(df$comparison_level))) +
        labs(
            x = "Comparison level",
            y = "Fraction of genes discovered",
            color = "Cluster",
            title = paste0("Filtering response clusters (K = ", K, ")")
        ) +
        theme_bw()

    ## Plot mapping clusters to cell types
    ## Order cell types by (cluster, then name)
    ord_idx <- order(cluster_df$cluster, cluster_df$cell_type)
    cluster_df$cell_type <- factor(cluster_df$cell_type,
                                   levels = cluster_df$cell_type[ord_idx])

    p_map <- ggplot(cluster_df,
                    aes(x = cluster, y = cell_type, color = cluster)) +
        geom_point(size = 2) +
        scale_y_discrete(name = "Cell type") +
        scale_x_discrete(name = "Cluster") +
        guides(color = "none") +
        theme_bw()

    list(
        plot_trajectories = p_trajectories,
        plot_mapping      = p_map,
        clusters          = cluster_df,
        mean_trajectories = mean_long,
        wide_matrix       = mat_numeric
    )
}

#' Heatmap of filtering-response trajectories
#'
#' Convert per-cell-type discovery fractions across comparison levels into a wide matrix
#' and draw a heatmap with hierarchical clustering of cell types.
#'
#' @param df Data frame with columns `cell_type`, `comparison_level`, and
#'   `frac_genes_discovered`. Each row represents one cell type at one comparison level.
#'
#' @return A `ComplexHeatmap::Heatmap` object.
#'
#' @export
plot_filtering_trajectories_heatmap <- function(df) {
    ## df must have: cell_type, comparison_level, frac_genes_discovered

    wide <- reshape(
        df[, c("cell_type", "comparison_level", "frac_genes_discovered")],
        idvar     = "cell_type",
        timevar   = "comparison_level",
        direction = "wide"
    )

    ## order columns by numeric comparison_level (1, 2, 3, 4, ...)
    comp_levels <- as.numeric(sub("frac_genes_discovered\\.", "", names(wide)[-1]))
    col_idx <- order(comp_levels)
    wide <- wide[, c(1, 1 + col_idx)]

    mat <- as.matrix(wide[, -1])
    rownames(mat) <- wide$cell_type

    ## hierarchical clustering on rows only
    row_hc <- stats::hclust(stats::dist(mat))

    ## diverging color scale centered at 1
    col_fun <- circlize::colorRamp2(
        c(min(mat, na.rm = TRUE), 1, max(mat, na.rm = TRUE)),
        c("blue", "white", "red")
    )

    ComplexHeatmap::Heatmap(
        mat,
        name = "frac_genes_discovered",
        col = col_fun,
        cluster_rows = row_hc,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE
    )
}

plot_filtering_trajectories_heatmap2 <- function(df) {
    ## Required columns:
    ##   cell_type
    ##   comparison_level
    ##   num_genes_significant_old   (defines level 0)
    ##   num_genes_significant_new   (levels 1..)

    df$comparison_level <- as.numeric(df$comparison_level)

    ## Wide matrix for "new" counts
    wide_new <- reshape(
        df[, c("cell_type", "comparison_level", "num_genes_significant_new")],
        idvar     = "cell_type",
        timevar   = "comparison_level",
        direction = "wide"
    )

    ## Order new-level columns numerically
    levs <- as.numeric(sub("num_genes_significant_new\\.", "", names(wide_new)[-1]))
    col_idx <- order(levs)
    wide_new <- wide_new[, c(1, 1 + col_idx)]

    mat_new <- as.matrix(wide_new[, -1, drop = FALSE])
    rownames(mat_new) <- wide_new$cell_type
    colnames(mat_new) <- as.character(sort(levs))

    ## Level 0 column from old counts
    old_by_ct <- tapply(df$num_genes_significant_old, df$cell_type, function(x) x[1])
    old_by_ct <- old_by_ct[rownames(mat_new)]
    mat0 <- matrix(old_by_ct, ncol = 1)
    colnames(mat0) <- "0"
    rownames(mat0) <- rownames(mat_new)

    ## Combine and normalize per row to [0, 1]
    mat_counts <- cbind(mat0, mat_new)
    row_max <- apply(mat_counts, 1, max, na.rm = TRUE)
    row_max[!is.finite(row_max) | row_max <= 0] <- NA_real_
    mat_norm <- mat_counts / row_max

    ## Cluster rows only
    row_hc <- stats::hclust(stats::dist(mat_norm))

    ## Single-hue sequential color scale (low = white, high = blue)
    col_fun <- circlize::colorRamp2(
        c(0, 1),
        c("white", "blue")
    )

    ComplexHeatmap::Heatmap(
        mat_norm,
        name = "frac_of_max",
        col = col_fun,
        cluster_rows = row_hc,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        column_title = "Comparison level (0 = old)"
    )
}



plot_reduction_vs_initial <- function(df, cluster_df, comparison_level=4) {

    ## attach cluster to main df
    df_cl <- merge(df, cluster_df, by = "cell_type")

    ## keep only comparison_level maximum
    d3 <- df_cl[df_cl$comparison_level == comparison_level, ]

    ## log10 of initial number of significant genes
    d3$log10_num_genes_old <- log10(d3$num_genes_significant_old)

    #Make R CMD CHECK happy
    log10_num_genes_old<-frac_genes_discovered<-parent_type<-NULL

    #Make R CMD CHECK happy
    cluster<-NULL

    ggplot(d3,
           aes(x = log10_num_genes_old,
               y = frac_genes_discovered,
               color = cluster)) +
        geom_point(size = 2) +
        geom_smooth(method = "loess",
                    se = FALSE,
                    color = "black",
                    linetype = "dashed",
                    linewidth = 0.8) +
        labs(
            x = expression(log[10]("Number of significant genes (before filtering)")),
            y = paste0("Fraction of genes discovered (level ",comparison_level, " )"),
            color = "Cluster",
            title = "Reduction vs initial number of discoveries"
        ) +
        theme_bw()
}


plot_reduction_by_parent_type <- function(df, comparison_level=4) {

    ## Keep only comparison_level = 3
    d3 <- df[df$comparison_level == comparison_level, ]

    ## log10 of original # significant genes
    d3$log10_num_genes_old <- log10(d3$num_genes_significant_old)

    ## Create parent cell type label (truncate at first underscore)
    ## If no underscore exists, keep original
    d3$parent_type <- sub("(_.*)$", "", d3$cell_type)

    ## Convert to factor for stable legend ordering
    d3$parent_type <- factor(d3$parent_type)

    #Make R CMD CHECK happy
    log10_num_genes_old<-frac_genes_discovered<-parent_type<-NULL

    ggplot(d3,
           aes(x = log10_num_genes_old,
               y = frac_genes_discovered,
               color = parent_type)) +
        geom_point(size = 2) +
        geom_smooth(method = "loess",
                    se = FALSE,
                    color = "black",
                    linetype = "dashed",
                    linewidth = 0.8) +
        labs(
            x = expression(log[10]("Number of significant genes (old filtering)")),
            y = paste0("Fraction of genes discovered (level ",comparison_level, " )"),
            color = "Parent cell type",
            title = "Reduction vs initial discoveries by parent cell group"
        ) +
        theme_bw()
}

add_supertitle_cowplot <- function(p,
                                   title,
                                   title_size = 16,
                                   rel_height_title = 0.08) {
    title_grob <- cowplot::ggdraw() +
        cowplot::draw_label(
            title,
            fontface = "bold",
            size = title_size,
            x = 0.5,
            hjust = 0.5
        )

    cowplot::plot_grid(
        title_grob,
        p,
        ncol = 1,
        rel_heights = c(rel_height_title, 1)
    )
}


