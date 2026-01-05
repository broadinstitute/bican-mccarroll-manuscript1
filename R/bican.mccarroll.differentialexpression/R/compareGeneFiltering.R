# Compare old vs new CPM filtering of differential expression results

#library(ggplot2)
#library(cowplot)
#library(ggrepel)
#library(logger)

#
# old_data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression_old_gene_filtering/sex_age/cell_type"
# new_data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age/cell_type"
# # outPDF="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age/cell_type/compare_feature_selection.pdf"
# fdr_cutoff=0.05
# cell_type="microglia"

#compare_age_de_runs(old_data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression_old_gene_filtering/sex_age/cell_type", new_data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age/cell_type", outPDF="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age/cell_type/compare_feature_selection.pdf", fdr_cutoff=0.05)

#compare all levels
#compare_all_age_de_runs(data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results", outDir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/compare_results_by_level")


#compare all levels to each other.
compare_all_age_de_runs<-function (data_dir, outDir, filter_levels=c(0,1,2,3), fdr_cutoff=0.05) {
    for (base_level in filter_levels) {
        for (comparison_level in filter_levels) {
            if (base_level < comparison_level) {
                old_data_dir=paste(data_dir,"/LEVEL_", base_level, "/sex_age/cell_type", sep="")
                new_data_dir=paste(data_dir,"/LEVEL_", comparison_level, "/sex_age/cell_type", sep="")
                outPDF=paste(outDir, "/compare_age_DE_LEVEL_", base_level, "_vs_LEVEL_", comparison_level, ".pdf", sep="")
                outFile=paste(outDir, "/compare_age_DE_LEVEL_", base_level, "_vs_LEVEL_", comparison_level, ".txt", sep="")
                logger::log_info(paste0("Comparing age DE results between LEVEL ", base_level, " and LEVEL ", comparison_level, "\n"))
                z=compare_age_de_runs(old_data_dir=old_data_dir,
                                    new_data_dir=new_data_dir,
                                    outPDF=outPDF,
                                    outFile=outFile,
                                    fdr_cutoff=fdr_cutoff)
            }
        }
    }
}

#compute all the level to level +1 comparisons, merge into a single result dataframe
compare_all_age_de_runs_trajectories<-function (data_dir, outDir, filter_levels=c(0,1,2,3), fdr_cutoff=0.05) {
    base_level=filter_levels[1]
    results=list()
    for (i in 1:(length(filter_levels)-1)) {
        comparison_level=filter_levels[i+1]
        old_data_dir=paste(data_dir,"/LEVEL_", base_level, "/sex_age/cell_type", sep="")
        new_data_dir=paste(data_dir,"/LEVEL_", comparison_level, "/sex_age/cell_type", sep="")
        logger::log_info(paste0("Comparing age DE results between LEVEL ", base_level, " and LEVEL ", comparison_level, "\n"))
        z=compare_age_de_runs(old_data_dir=old_data_dir, new_data_dir=new_data_dir, outPDF=NULL, outFile=NULL, fdr_cutoff=fdr_cutoff)
        df=z$df
        df$base_level=base_level
        df$comparison_level=comparison_level
        results[[i]]=df
    }
    df=do.call(rbind, results)

    plot_frac_lines(df)

    #unknown how many clusters to pick.
    res <- cluster_filtering_trajectories(df, K = 4)
    p1<-res$plot_trajectories
    p2<-res$plot_mapping

    combined <- cowplot::plot_grid(
        p1,
        p2,
        nrow       = 2
    )

    # is the dropoff a result of the number of initial observations?
    # IE: we just barely had power before, and now we don't?
    plot_reduction_vs_initial(df, cluster_df=res$clusters)
    plot_reduction_by_parent_type(df)

}



compare_age_de_runs<-function (old_data_dir, new_data_dir, outPDF, outFile=NULL, fdr_cutoff=0.05) {
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
                            fdr_cutoff=fdr_cutoff)
        plot_list[[cell_type]]=p$plot
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
        pdf(outPDF)
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
#' @param fdr_cutoff Numeric scalar. FDR threshold used to classify genes
#'   as significant or non-significant in the old run when summarizing
#'   dropped genes. Default is \code{0.05}.
#'
#' @return A combined plot.
#' @keywords internal
#' @noRd
compare_age_de_run <- function (cell_type,
                                old_data_dir,
                                new_data_dir,
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
            x = "logFC (old filtering)",
            y = "logFC (strict filtering)"
        ) +
        coord_fixed(xlim = lim_logfc, ylim = lim_logfc) +
        theme_bw()

    ## ---------- stats and plot for FDR comparison ----------
    merged$neg_log10_fdr_old <- -log10(merged$adj.P.Val_old)
    merged$neg_log10_fdr_new <- -log10(merged$adj.P.Val_new)

    num_genes_significant_old=sum(merged$adj.P.Val_old < fdr_cutoff, na.rm=TRUE)
    num_genes_significant_new=sum(merged$adj.P.Val_new < fdr_cutoff, na.rm=TRUE)
    strTitleSuffix=paste(num_genes_significant_old, " significant old; ", num_genes_significant_new, " new", sep="")

    fdr_cor <- cor(merged$neg_log10_fdr_old,
                   merged$neg_log10_fdr_new,
                   use = "complete.obs")

    ann_text_fdr <- paste0(
        "r = ", sprintf("%.3f", fdr_cor)
    )

    x_range_fdr <- range(merged$neg_log10_fdr_old, na.rm = TRUE)
    y_range_fdr <- range(merged$neg_log10_fdr_new, na.rm = TRUE)

    x_pos_fdr <- x_range_fdr[1L] + 0.02 * diff(x_range_fdr)
    y_pos_fdr <- y_range_fdr[2L] - 0.02 * diff(y_range_fdr)

    lim_fdr <- range(c(merged$neg_log10_fdr_old, merged$neg_log10_fdr_new),
                     na.rm = TRUE)

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
            x = expression(-log[10]("FDR (old) filtering")),
            y = expression(-log[10]("FDR (strict filtering)"))
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

    ## Combined cowplot layout:
    top_row <- cowplot::plot_grid(
        p_effect,
        p_fdr,
        ncol = 2
    )

    bottom_row <- cowplot::plot_grid(
        p_drop_bar,
        ncol = 1
    )

    combined <- cowplot::plot_grid(
        top_row,
        bottom_row,
        ncol = 1,
        rel_heights = c(1.5, 1)
    )

    df=data.frame(cell_type=cell_type,
               logFC_correlation=logfc_cor,
               logFC_sign_agreement=sign_agree,
               FDR_correlation=fdr_cor,
               n_dropped_significant=n_sig_dropped,
               n_dropped_non_significant=n_non_sig_dropped,
               num_genes_significant_old=num_genes_significant_old,
               num_genes_significant_new=num_genes_significant_new)

    result=list(df=df, plot=combined, merged=merged)
    return (result)
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

plot_reduction_vs_initial <- function(df, cluster_df) {

    ## attach cluster to main df
    df_cl <- merge(df, cluster_df, by = "cell_type")

    ## keep only comparison_level == 3
    d3 <- df_cl[df_cl$comparison_level == 3, ]

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
            y = "Fraction of genes discovered (level 3)",
            color = "Cluster",
            title = "Reduction vs initial number of discoveries"
        ) +
        theme_bw()
}


plot_reduction_by_parent_type <- function(df) {

    ## Keep only comparison_level = 3
    d3 <- df[df$comparison_level == 3, ]

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
            y = "Fraction of genes discovered (level 3)",
            color = "Parent cell type",
            title = "Reduction vs initial discoveries by parent cell group"
        ) +
        theme_bw()
}


