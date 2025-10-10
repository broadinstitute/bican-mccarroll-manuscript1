# library(ggplot2)
# library(reshape2)
# library(ggforce)


# in_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age/cell_type/"
# file_pattern="age"
# cellTypeListFile="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metadata/mash_cell_type_list_simple.txt"


#' @title plotPairwiseEffectSizes
#' @description
#' Generate summary and pairwise plots of differential expression (DE) effect sizes
#' and p-values across two experiments.
#'
#' @details
#' Reads and parses DE input files, constructs mash input matrices, visualizes
#' p-value and effect size distributions, and produces pairwise scatterplots
#' of effect sizes across selected reference conditions.
#'
#' @param in_dir Directory containing DE result files.
#' @param file_pattern Pattern to match DE result files in `in_dir`.
#' @param cellTypeListFile Path to a file listing cell types to include in the
#' @return Invisibly returns a list of data frames from the pairwise effect plots.
#'
#' @importFrom stats na.omit
#' @importFrom graphics hist
#' @importFrom grDevices dev.off
#' @noRd
plotPairwiseEffectSizes<-function (in_dir, file_pattern, cellTypeListFile) {

    d=parse_de_inputs(in_dir, file_pattern, cellTypeListFile)

    #make mash inputs
    mash_inputs_union<-make_mash_inputs(d, coef_col = "logFC", t_col = "t", fdr_col = "adj.P.Val", gene_mode="union")

    effect_size=mash_inputs_union$Bhat
    effect_size[mash_inputs_union$missing_mask==T]<-NA

    fdr=mash_inputs_union$FDR
    pval=mash_inputs_union$P_val

    #plot pval histograms by cell type
    plot_pval_histograms(pval)
    plot_effect_histograms(effect_size, bins=100, free_y=FALSE, log_y = FALSE)
    plot_effect_histograms(effect_size, bins=100, free_y=FALSE, log_y = TRUE)

    #pairwise effect size plots
    #plot_effect_pairs(effect_size, fdr, condition_ref = "microglia", min_fdr_threshold = 1e-100, ncol=4, nrow=3)
    df1=plot_effect_pairs(effect_size, fdr, condition_ref = "microglia", min_fdr_threshold = 0.5, ncol=4, nrow=3)

    df2=plot_effect_pairs(effect_size, fdr, condition_ref = "MSN_D1", min_fdr_threshold = 0.5, ncol=4, nrow=3)
    df3=plot_effect_pairs(effect_size, fdr, condition_ref = "MSN_D2", min_fdr_threshold = 0.5, ncol=4, nrow=3)



}

plot_effect_pairs <- function(effect_size, fdr, condition_ref,
                              min_fdr_threshold = 0.2,
                              ncol = 2, nrow = 2) {
    stopifnot(is.matrix(effect_size), is.matrix(fdr))
    stopifnot(identical(colnames(effect_size), colnames(fdr)))
    stopifnot(condition_ref %in% colnames(effect_size))

    ref <- condition_ref
    others <- setdiff(colnames(effect_size), ref)

    # build long df
    df_list <- lapply(others, function(ct) {
        d <- data.frame(
            fc1  = effect_size[, ref],
            fc2  = effect_size[, ct],
            fdr1 = fdr[, ref],
            fdr2 = fdr[, ct]
        )

        d$min_fdr <- pmin(d$fdr1, d$fdr2)
        d <- d[!is.na(d$min_fdr) & d$min_fdr > min_fdr_threshold, , drop = FALSE]

        r <- suppressWarnings(cor(d$fc1, d$fc2, use = "complete.obs"))
        d$panel <- sprintf("%s vs %s\nr = %.3f", ref, ct, r)
        d
    })
    df_all <- do.call(rbind, df_list)
    naRows=apply(df_all, 1, function (x) any(is.na(x)))
    df_all=df_all[!naRows, ]

    if (!nrow(df_all)) {
        message("No rows pass the FDR threshold.")
        return(invisible(0L))
    }

    # R CMD CHECK Happy
    fc1 <- fc2 <- panel <- NULL

    base <- ggplot(df_all, aes(x = fc1, y = fc2)) +
        geom_point(alpha = 0.4, size = 1.2, na.rm = TRUE) +
        geom_abline(slope = 1, intercept = 0,
                    color = "gray50", linetype = "dotted", linewidth = 0.6) +
        # prefilter finite rows for the smoother to avoid warnings
        geom_smooth(
            data = function(d) d[is.finite(d$fc1) & is.finite(d$fc2), ],
            method = "lm", se = FALSE,
            color = "red", linetype = "dashed", linewidth = 0.7
        ) +
        ggforce::facet_wrap_paginate(~ panel, ncol = ncol, nrow = nrow, scales = "free") +
        labs(x = paste(ref, "effect size"),
             y = "other condition effect size",
             title = sprintf("Effect-size comparisons vs %s (min FDR > %.2f)", ref, min_fdr_threshold)) +
        theme_minimal(base_size = 14) +
        theme(plot.title = element_text(hjust = 0.5),
              panel.grid.minor = element_blank())

    n_pages <- ceiling(length(unique(df_all$panel)) / (ncol * nrow))
    for (pg in seq_len(n_pages)) {
        print(base + ggforce::facet_wrap_paginate(~ panel, ncol = ncol, nrow = nrow,
                                         scales = "free", page = pg))
    }
    return (df_all)
}


plot_pval_histograms <- function(pmat, bins = 40, free_y = TRUE) {
    stopifnot(is.matrix(pmat))

    # reshape: genes x conditions -> long
    df <- data.frame(
        pval = as.vector(pmat),
        condition = rep(colnames(pmat), each = nrow(pmat)),
        stringsAsFactors = FALSE
    )
    # keep valid [0,1]
    df$pval <- as.numeric(df$pval)
    df <- df[is.finite(df$pval) & df$pval >= 0 & df$pval <= 1, ]

    # R CMD CHECK Happy
    pval <- condition <- NULL

    p <- ggplot(df, aes(x = pval)) +
        geom_histogram(bins = bins, fill = "steelblue", color = "white") +
        facet_wrap(~ condition, scales = if (free_y) "free_y" else "fixed") +
        labs(x = "p-value", y = "Frequency", title = "P-value distributions per condition") +
        theme_minimal(base_size = 14) +
        theme(strip.text = element_text(face = "bold"),
              plot.title = element_text(hjust = 0.5))

    p
}


plot_effect_histograms <- function(emat, bins = 100, free_y = FALSE, log_y = TRUE) {
    stopifnot(is.matrix(emat))

    df <- data.frame(
        effect = as.vector(emat),
        condition = rep(colnames(emat), each = nrow(emat)),
        stringsAsFactors = FALSE
    )

    # keep only finite values
    df <- df[is.finite(df$effect), ]
    if (!nrow(df))
        stop("No finite effect sizes found.")

    # R CMD CHECK Happy
    effect <- condition <- NULL

    p <- ggplot(df, aes(x = effect)) +
        geom_histogram(
            bins = bins,
            fill = "steelblue",
            color = "white",
            alpha = 0.8,
            na.rm = TRUE
        ) +
        geom_vline(
            xintercept = 0,
            color = "black",
            linetype = "dashed",
            linewidth = 0.3
        ) +
        facet_wrap(~ condition, scales = if (free_y) "free_y" else "fixed") +
        labs(
            x = "Effect size",
            y = if (log_y) "Frequency (log10)" else "Frequency",
            title = "Effect-size distributions per condition"
        ) +
        theme_minimal(base_size = 14) +
        theme(
            strip.text = element_text(face = "bold"),
            plot.title = element_text(hjust = 0.5),
            panel.grid.minor = element_blank()
        )

    # apply log10 scaling safely (ggplot ignores zeros automatically)
    if (log_y) {
        p <- p + scale_y_continuous(
            trans  = "log10",
            limits = c(1, NA),           # drop zero-count bins
            oob    = scales::censor,     # censor values outside limits
            breaks = scales::log_breaks(),
            labels = scales::label_number()
        )
    }

    p
}



