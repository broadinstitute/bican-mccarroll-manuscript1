# Compare old vs new CPM filtering of differential expression results

# library(ggplot2)
# library(cowplot)
# library(ggrepel)
#
# old_data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression_old_gene_filtering/sex_age/cell_type"
# new_data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age/cell_type"
# # outPDF="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age/cell_type/compare_feature_selection.pdf"
# fdr_cutoff=0.05
# cell_type="microglia"

#compare_age_de_runs(old_data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression_old_gene_filtering/sex_age/cell_type", new_data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age/cell_type", outPDF="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age/cell_type/compare_feature_selection.pdf", fdr_cutoff=0.05)

compare_age_de_runs<-function (old_data_dir, new_data_dir, outPDF, fdr_cutoff=0.05) {
    suffix="_age_DE_results.txt"
    f=list.files(old_data_dir, pattern=suffix, full.names = F)
    cell_types_list=sub(suffix, "", f)

    plot_list=list()
    df_list=list()

    for (cell_type in cell_types_list){
        logger::log_info(paste0("Comparing age DE results for cell type: ", cell_type, "\n"))
        p=compare_age_de_run(cell_type=cell_type,
                            old_data_dir=old_data_dir,
                            new_data_dir=new_data_dir,
                            fdr_cutoff=fdr_cutoff)
        plot_list[[cell_type]]=p$plot
        df_list[[cell_type]]=p$df
    }

    df=do.call(rbind, df_list)

    pdf(outPDF)
    print (plot_summary(df))
    for (p in plot_list){
        print(p)
    }
    dev.off()

}


plot_summary<-function (df) {
    ## df has columns: cell_type, logFC_correlation, logFC_sign_agreement, FDR_correlation

    ## reshape to long format using base R
    df_long <- reshape(
        df,
        direction = "long",
        varying = list(c("logFC_correlation",
                         "logFC_sign_agreement",
                         "FDR_correlation")),
        v.names = "value",
        timevar = "metric",
        times = c("logFC_correlation",
                  "logFC_sign_agreement",
                  "FDR_correlation")
    )

    ## reshape() creates an id column and rownames; drop what you don't need
    row.names(df_long) <- NULL

    ## ensure metric is a factor with a nice order
    df_long$metric <- factor(
        df_long$metric,
        levels = c("logFC_correlation",
                   "logFC_sign_agreement",
                   "FDR_correlation")
    )

    #for labeling outliers.
    df_long$z <- ave(
        df_long$value,
        df_long$metric,
        FUN = function(x) {
            m  <- mean(x, na.rm = TRUE)
            sd <- stats::sd(x, na.rm = TRUE)
            if (sd == 0) {
                rep(0, length(x))
            } else {
                (x - m) / sd
            }
        }
    )

    ## 3) Label outliers with |z| > 2
    df_long$label <- ifelse(abs(df_long$z) > 2, df_long$cell_type, "")

    ## violin plot

    p=ggplot(df_long, aes(x = metric, y = value)) +
        geom_violin(trim = FALSE) +
        geom_jitter(width = 0.08, height = 0, size = 1, alpha = 0.6) +
        geom_text_repel(
            data = df_long[df_long$label != "", ],
            aes(label = label),
            max.overlaps = Inf,
            box.padding = 0.3,
            point.padding = 0.1,
            min.segment.length = 0
        ) +
        labs(
            x = NULL,
            y = "Value",
            title = "Concordance metrics across cell types"
        ) +
        theme_bw()
    p

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
               n_dropped_non_significant=n_non_sig_dropped)

    result=list(df=df, plot=combined)
    return (result)
}
