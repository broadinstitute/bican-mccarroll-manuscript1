library(ggplot2)

# This code is slightly more adhoc and uses external data to answer reviewer questions.
# plot_bican_snap_manifest(manifest_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/external_comparison_snap200/metadata/bican_dfc_snap_age_de_overlap_manifest.tsv", out_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/external_comparison_snap200/results")


## ------------------------------------------------------------------
## Top-level functions
## ------------------------------------------------------------------

plot_bican_snap_manifest <- function(
        manifest_file,
        out_dir = NULL,
        alpha = 0.05,
        width = 7,
        height = 7) {

    manifest <- read.table(
        manifest_file,
        header = TRUE,
        stringsAsFactors = FALSE,
        sep = "\t",
        check.names = FALSE
    )

    required_cols <- c("cell_type", "bican", "snap")
    missing_cols <- setdiff(required_cols, colnames(manifest))

    if (length(missing_cols) > 0) {
        stop(sprintf(
            "Manifest is missing required columns: %s",
            paste(missing_cols, collapse = ", ")
        ))
    }

    if (!is.null(out_dir) && !dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }

    results <- vector("list", nrow(manifest))
    names(results) <- manifest$cell_type

    sign_results <- vector("list", nrow(manifest))
    names(sign_results) <- manifest$cell_type

    for (i in seq_len(nrow(manifest))) {
        cell_type <- manifest$cell_type[i]

        res <- plot_bican_snap_cell_type(
            cell_type = cell_type,
            bican_file = manifest$bican[i],
            snap_file = manifest$snap[i],
            alpha = alpha
        )

        results[[cell_type]] <- res
        sign_results[[cell_type]] <- res$sign_test

        if (!is.null(out_dir)) {
            out_file <- file.path(
                out_dir,
                sprintf(
                    "bican_snap_age_effect_scatter_%s.svg",
                    gsub("[^A-Za-z0-9_.-]+", "_", cell_type)
                )
            )

            ggplot2::ggsave(
                filename = out_file,
                plot = res$plot,
                device = "svg",
                width = width,
                height = height
            )
        }
    }

    sign_results <- do.call(rbind, sign_results)

    sign_summary_plot <- plot_sign_test_summary(sign_results)

    if (!is.null(out_dir)) {
        ggplot2::ggsave(
            filename = file.path(out_dir, "bican_snap_sign_test_summary.svg"),
            plot = sign_summary_plot,
            device = "svg",
            width = 8,
            height = max(4, 0.35 * nrow(sign_results) + 2)
        )

        write.table(
            sign_results,
            file = file.path(out_dir, "bican_snap_sign_test_summary.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
        )
    }

    invisible (list(
        per_cell_type = results,
        sign_results = sign_results,
        sign_summary_plot = sign_summary_plot
    ))
}


plot_bican_snap_cell_type <- function(cell_type, bican_file, snap_file, alpha = 0.05) {
    df <- make_bican_snap_df(
        bican_file = bican_file,
        snap_file = snap_file,
        alpha = alpha
    )

    sign_res <- sign_test_fraction(
        x = df$snap_logFC,
        y = df$bican_logFC
    )

    cor_val <- stats::cor(
        df$snap_logFC,
        df$bican_logFC,
        use = "complete.obs"
    )

    lims <- range(c(df$snap_logFC, df$bican_logFC), na.rm = TRUE)

    title <- cell_type

    subtitle <- sprintf(
        "BICAN-significant genes: %.1f%% same sign in SNAP (%d/%d); cor=%.2f",
        100 * sign_res$frac,
        sign_res$k,
        sign_res$n,
        cor_val
    )

    p <- ggplot(df, aes(x = snap_logFC, y = bican_logFC)) +
        geom_point(na.rm = TRUE, size = 1.2, alpha = 0.6) +
        geom_abline(intercept = 0, slope = 1, color = "black") +
        coord_cartesian(xlim = lims, ylim = lims) +
        labs(
            title = title,
            subtitle = subtitle,
            x = "SNAP logFC",
            y = "BICAN logFC"
        ) +
        theme_bw(base_size = 12) +
        theme(
            plot.title = element_text(size = 13, face = "bold"),
            plot.subtitle = element_text(size = 11),
            axis.title = element_text(size = 13),
            axis.text = element_text(size = 11)
        )

    list(
        plot = p,
        sign_test = data.frame(
            cell_type = cell_type,
            k = sign_res$k,
            n = sign_res$n,
            fraction = sign_res$frac,
            percent = 100 * sign_res$frac,
            p_value = sign_res$p_value,
            cor = cor_val,
            stringsAsFactors = FALSE
        ),
        data = df
    )
}


## ------------------------------------------------------------------
## Plotting helpers
## ------------------------------------------------------------------

plot_sign_test_summary <- function(sign_results) {
    sign_results <- sign_results[order(sign_results$percent), ]

    sign_results$cell_type <- factor(
        sign_results$cell_type,
        levels = sign_results$cell_type
    )

    ggplot(sign_results, aes(x = cell_type, y = percent)) +
        geom_col(color = "black") +
        geom_text(
            aes(label = paste0(k, "/", n)),
            hjust = -0.1,
            size = 5
        ) +
        coord_flip(ylim = c(0, 105)) +
        labs(
            x = NULL,
            y = "% same sign",
            title = "BICAN-significant age effects: sign concordance in SNAP"
        ) +
        theme_bw(base_size = 12)
}


## ------------------------------------------------------------------
## Data helpers
## ------------------------------------------------------------------

make_bican_snap_df <- function(bican_file, snap_file, alpha = 0.05) {
    bican <- read_de_result(bican_file)
    snap <- read_de_result(snap_file)

    common_genes <- intersect(rownames(bican), rownames(snap))

    genes <- rownames(bican)[bican$adj.P.Val < alpha]
    genes <- intersect(genes, common_genes)

    data.frame(
        gene = genes,
        snap_logFC = snap[genes, "logFC"],
        snap_adjP = snap[genes, "adj.P.Val"],
        bican_logFC = bican[genes, "logFC"],
        bican_adjP = bican[genes, "adj.P.Val"],
        stringsAsFactors = FALSE
    )
}


read_de_result <- function(file) {
    x <- read.table(
        file,
        header = TRUE,
        stringsAsFactors = FALSE,
        sep = "\t",
        check.names = FALSE
    )

    required_cols <- c("logFC", "adj.P.Val")
    missing_cols <- setdiff(required_cols, colnames(x))

    if (length(missing_cols) > 0) {
        stop(sprintf(
            "Missing required columns in %s: %s",
            file,
            paste(missing_cols, collapse = ", ")
        ))
    }

    x
}


## ------------------------------------------------------------------
## Statistics helpers
## ------------------------------------------------------------------

sign_test_fraction <- function(x, y) {
    ok <- is.finite(x) & is.finite(y) & x != 0 & y != 0

    n <- sum(ok)

    if (n == 0) {
        return(list(
            frac = NA_real_,
            n = 0,
            k = 0,
            p_value = NA_real_
        ))
    }

    k <- sum(sign(x[ok]) == sign(y[ok]))

    bt <- stats::binom.test(
        x = k,
        n = n,
        p = 0.5,
        alternative = "greater"
    )

    list(
        frac = k / n,
        n = n,
        k = k,
        p_value = bt$p.value
    )
}
