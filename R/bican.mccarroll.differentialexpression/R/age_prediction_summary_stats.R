library (ggplot2)
library (cowplot)

age_preds_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/age_prediction/age_prediction_region_alpha_0/age_prediction_results_alpha0_donor_predictions.txt"
cell_type_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/metadata/cell_types_for_prs_test.txt"

# summarize median absolute error, mean absolute error, and median absolute error
# in the youngest 20% of donors for each cell type and region, and save to a summary file

#' Summarize age-prediction model performance and covariate associations
#'
#' This function reads per-donor age prediction results, optionally restricts
#' to a user-specified set of cell types, computes model-level error metrics
#' for each cell type by region, generates summary plots, and fits simple
#' regression models relating prediction error to model size / data size
#' covariates.
#'
#' The expected input file is a tab-delimited table with one row per donor per
#' cell type x region model. At minimum, it should contain the columns needed
#' by \code{\link{compute_age_error_metrics}} (including \code{cell_type},
#' \code{region}, \code{donor}, \code{age}, the prediction column used by
#' \code{compute_age_error_metrics}, and \code{num_nuclei}, \code{num_umis},
#' \code{num_features}).
#'
#' @param age_preds_file Character scalar. Path to a tab-delimited file
#' containing per-donor predictions and metadata.
#' @param cell_type_file Optional character scalar. Path to a text file
#' containing one cell type per line. If provided, \code{age_preds_file} is
#' filtered to those cell types.
#' @param out_summary_file Character scalar. Path for writing a tab-delimited
#' summary table of model-level metrics (one row per cell type x region).
#'
#' @return A named list with:
#' \describe{
#'   \item{age_preds}{data.frame of per-donor predictions (possibly filtered).}
#'   \item{metrics_df}{data.frame of model-level error metrics and covariates.}
#'   \item{fit_all}{List returned by \code{\link{fit_age_mae_regression}} for \code{mae}.}
#'   \item{fit_y20}{List returned by \code{\link{fit_age_mae_regression}} for \code{mae_young20}.}
#'   \item{p_mae}{Heatmap of \code{mae}.}
#'   \item{p_mae_y20}{Heatmap of \code{mae_young20}.}
#'   \item{combined}{Combined (cowplot) panel of MAE vs predictors.}
#'   \item{combined_y20}{Combined (cowplot) panel of MAE in youngest 20\% vs predictors.}
#' }
#'
#' @import ggplot2 cowplot
#' @export
summarize_age_prediction_results<-function (age_preds_file, cell_type_file=NULL, out_summary_file) {

    age_preds=read.table(age_preds_file, header=TRUE, sep="\t", stringsAsFactors = FALSE)

    #filter the age_preds to a restricted set of cell types.
    if (!is.null(cell_type_file)) {
        cell_types <- read.table(cell_type_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)$V1
        age_preds <- age_preds[age_preds$cell_type %in% cell_types, ]
    }

    metrics_df <- compute_age_error_metrics(age_preds)

    # Heatmap for any metric
    p1 <- plot_age_metric_heatmap(metrics_df, metric = "mae", text_size = 4)
    p2 <- plot_age_metric_heatmap(metrics_df, metric = "mae_young20")

    # Regressions
    fit_all <- fit_age_mae_regression(metrics_df, outcome = "mae")
    fit_y20 <- fit_age_mae_regression(metrics_df, outcome = "mae_young20")

    p_feat <- plot_error_vs_predictor(metrics_df, x_col = "num_features", y_col = "mae", x_scale = "identity", annotate_r2_size=5, show_y_label = FALSE)
    p_nuc  <- plot_error_vs_predictor(metrics_df, x_col = "num_nuclei",   y_col = "mae", x_scale = "identity", annotate_r2_size=5, show_y_label = FALSE)
    p_umi  <- plot_error_vs_predictor(metrics_df, x_col = "num_umis",     y_col = "mae", x_scale = "identity", annotate_r2_size=5, show_y_label = FALSE)

    y_axis_label <- ggdraw() +
        draw_label("Mean AE (decades)", angle = 90, size = 13)

    combined <- plot_grid(y_axis_label,
        plot_grid(p_feat, p_nuc, p_umi, ncol = 1),
        ncol = 2, rel_widths = c(0.05, 1)
    )

    p_feat <- plot_error_vs_predictor(metrics_df, x_col = "num_features", y_col = "mae_young20", x_scale = "identity", annotate_r2_size=5, show_y_label = FALSE)
    p_nuc  <- plot_error_vs_predictor(metrics_df, x_col = "num_nuclei",   y_col = "mae_young20", x_scale = "identity", annotate_r2_size=5, show_y_label = FALSE)
    p_umi  <- plot_error_vs_predictor(metrics_df, x_col = "num_umis",     y_col = "mae_young20", x_scale = "identity", annotate_r2_size=5, show_y_label = FALSE)

    y_axis_label <- ggdraw() +
        draw_label("Mean AE in youngest 20% (decades)", angle = 90, size = 13)

    combined_y20 <- plot_grid(y_axis_label,
                          plot_grid(p_feat, p_nuc, p_umi, ncol = 1),
                          ncol = 2, rel_widths = c(0.05, 1)
    )

    result=list(age_preds=age_preds, metrics_df=metrics_df, fit_all=fit_all, fit_y20=fit_y20, p_mae=p1, p_mae_y20=p2, p_feat=p_feat, p_nuc=p_nuc, p_umi=p_umi, combined=combined, combined_y20=combined_y20)
    return (result)

}


#' Summarize age-prediction error metrics per cell_type x region
#'
#' @param age_preds data.frame with per-donor predictions and metadata.
#' Required columns: cell_type, region, donor, age, pred_mean_corrected,
#' num_nuclei, num_umis, num_features.
#' @param group_cols character vector of grouping columns.
#' @param pred_col prediction column (corrected predictions).
#' @param age_col chronological age column.
#' @param young_q quantile defining "youngest" subset within each group.
#'
#' @return data.frame with one row per group and error metrics/covariates.
compute_age_error_metrics <- function(age_preds,
                                      group_cols = c("cell_type", "region"),
                                      pred_col = "pred_mean",
                                      age_col = "age",
                                      young_q = 0.2) {
    stopifnot(is.data.frame(age_preds))
    req <- c(group_cols, "donor", age_col, pred_col,
             "num_nuclei", "num_umis", "num_features")
    miss <- setdiff(req, colnames(age_preds))
    if (length(miss) > 0) {
        stop("Missing required columns: ", paste(miss, collapse = ", "))
    }

    dt <- data.table::as.data.table(age_preds)

    dt[, abs_err := abs(get(pred_col) - get(age_col))]

    dt[, young_cut := stats::quantile(get(age_col), probs = young_q, type = 7, na.rm = TRUE),
       by = group_cols]
    dt[, is_young := get(age_col) <= young_cut]

    summarize_one_group <- function(d) {
        nf_u <- unique(d[!is.na(num_features), num_features])
        nf <- if (length(nf_u) == 1) as.numeric(nf_u) else NA_real_

        if (length(nf_u) > 1) {
            grp_vals <- vapply(group_cols, function(gc) as.character(d[[gc]][1]), character(1))
            warning("num_features not unique within group: ",
                    paste(grp_vals, collapse = " / "))
        }

        # Enforce numeric/double outputs where appropriate to keep data.table happy
        list(
            n_donors = as.integer(data.table::uniqueN(d$donor)),
            mae = as.numeric(mean(d$abs_err, na.rm = TRUE)),
            medae = as.numeric(stats::median(d$abs_err, na.rm = TRUE)),
            mae_young20 = as.numeric(mean(d$abs_err[d$is_young], na.rm = TRUE)),
            medae_young20 = as.numeric(stats::median(d$abs_err[d$is_young], na.rm = TRUE)),
            num_nuclei = as.numeric(sum(d$num_nuclei, na.rm = TRUE)),
            num_umis = as.numeric(sum(d$num_umis, na.rm = TRUE)),
            num_features = as.numeric(nf)
        )
    }

    out <- dt[, summarize_one_group(.SD), by = group_cols]
    data.table::setDF(out)
}


#' Heatmap for a single metric across cell_type x region
#'
#' @param metrics_df output of compute_age_error_metrics()
#' @param metric column name to plot (e.g. "mae" or "mae_young20")
#' @param cell_type_col column name for cell type
#' @param region_col column name for region
#' @param na_color color for missing tiles
#'
#' @return ggplot object
plot_age_metric_heatmap <- function(metrics_df,
                                metric,
                                cell_type_col = "cell_type",
                                region_col = "region",
                                na_color = "grey90",
                                add_text = TRUE,
                                text_digits = 1,
                                text_size = 3,
                                text_threshold = NULL) {
    stopifnot(is.data.frame(metrics_df))
    need <- c(cell_type_col, region_col, metric)
    miss <- setdiff(need, colnames(metrics_df))
    if (length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse = ", "))

    cell_type <- region <- value <- label <- text_col <- NULL

    dt <- data.table::as.data.table(metrics_df)
    dt[, value := get(metric)]
    data.table::setnames(dt, c(cell_type_col, region_col), c("cell_type", "region"))

    # Pretty labels (for axes / legend)
    dt[, cell_type_label := gsub("_", " ", cell_type, fixed = TRUE)]
    dt[, region_label := gsub("_", " ", region, fixed = TRUE)]
    metric_label <- gsub("_", " ", metric, fixed = TRUE)

    # Sort rows by average metric across regions (best = lowest error at top)
    cell_means <- dt[, .(mean_value = mean(value, na.rm = TRUE)), by = cell_type_label]
    cell_means <- cell_means[!is.nan(mean_value)]
    # For errors, lower is better; to put best at top, make it the last factor level
    cell_levels <- cell_means[order(mean_value, decreasing = TRUE), cell_type_label]
    dt[, cell_type_label := factor(cell_type_label, levels = cell_levels)]

    # Choose text color based on fill intensity:
    # With "Blues 3" here, low values are darker -> use white text; high values are lighter -> black text.
    rng <- range(dt$value, na.rm = TRUE)
    if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) {
        # degenerate case
        dt[, text_col := "black"]
    } else {
        if (is.null(text_threshold)) {
            text_threshold <- rng[1] + 0.55 * (rng[2] - rng[1])
        }
        dt[, text_col := ifelse(value <= text_threshold, "white", "black")]
    }

    p <- ggplot2::ggplot(dt, ggplot2::aes(x = region_label, y = cell_type_label, fill = value)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradientn(
            colors = grDevices::hcl.colors(256, "Blues 3"),
            na.value = na_color
        ) +
        ggplot2::labs(x = "Region", y = "Cell type", fill = metric_label) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            panel.grid = ggplot2::element_blank()
        )

    if (isTRUE(add_text)) {
        dt_txt <- dt[!is.na(value)]
        dt_txt[, label := sprintf(paste0("%.", text_digits, "f"), value)]
        p <- p + ggplot2::geom_text(
            data = dt_txt,
            ggplot2::aes(label = label, color = text_col),
            size = text_size,
            show.legend = FALSE
        ) +
            ggplot2::scale_color_identity()
    }

    p
}


fit_age_mae_regression <- function(metrics_df,
                                            outcome,
                                            numeric_covars = c("num_nuclei",
                                                               "num_umis",
                                                               "num_features"),
                                            log_pseudocount = 1) {
    stopifnot(is.data.frame(metrics_df))

    need <- c(outcome, numeric_covars)
    miss <- setdiff(need, colnames(metrics_df))
    if (length(miss) > 0) {
        stop("Missing required columns: ", paste(miss, collapse = ", "))
    }

    df <- metrics_df

    # Transform numeric predictors: scale(log10(x))
    for (nm in numeric_covars) {
        x <- df[[nm]]
        x2 <- log10(pmax(x, log_pseudocount, na.rm = TRUE))
        df[[paste0("zlog10_", nm)]] <- as.numeric(scale(x2))
    }

    rhs <- paste(paste0("zlog10_", numeric_covars), collapse = " + ")
    fml <- stats::as.formula(paste(outcome, "~", rhs))

    fit <- stats::lm(fml, data = df)

    sm <- summary(fit)
    coefs <- sm$coefficients
    ci <- try(stats::confint(fit), silent = TRUE)

    coef_df <- data.frame(
        term = rownames(coefs),
        estimate = coefs[, "Estimate"],
        std.error = coefs[, "Std. Error"],
        statistic = coefs[, "t value"],
        p.value = coefs[, "Pr(>|t|)"],
        stringsAsFactors = FALSE
    )

    if (!inherits(ci, "try-error")) {
        coef_df$conf.low <- ci[coef_df$term, 1]
        coef_df$conf.high <- ci[coef_df$term, 2]
    } else {
        coef_df$conf.low <- NA_real_
        coef_df$conf.high <- NA_real_
    }

    list(
        fit = fit,
        coef_table = coef_df,
        r.squared = sm$r.squared,
        adj.r.squared = sm$adj.r.squared
    )
}


plot_mae_vs_num_features <- function(metrics_df,
                                     x_col = "num_features",
                                     y_col = "mae",
                                     color_by = NULL,
                                     smoother_span = 1,
                                     point_size = 2.6,
                                     alpha = 0.9,
                                     show_legend = TRUE) {

    stopifnot(is.data.frame(metrics_df))

    need <- c(x_col, y_col)
    if (!is.null(color_by)) need <- c(need, color_by)

    miss <- setdiff(need, colnames(metrics_df))
    if (length(miss) > 0) {
        stop("Missing required columns: ", paste(miss, collapse = ", "))
    }

    num_features <- mae <- grp <- NULL

    df <- metrics_df
    df$num_features <- df[[x_col]]
    df$mae <- df[[y_col]]

    if (!is.null(color_by)) {
        raw_grp <- as.factor(df[[color_by]])
        lbls <- gsub("_", " ", levels(raw_grp))
        df$grp <- factor(raw_grp, levels = levels(raw_grp), labels = lbls)
    }

    n_grp <- length(levels(df$grp))
    pal <- .make_palette(n_grp)

    if (is.null(color_by)) {
        p <- ggplot2::ggplot(df, ggplot2::aes(x = num_features, y = mae)) +
            ggplot2::geom_point(
                size = point_size,
                alpha = alpha,
                shape = 21,
                fill = "grey70",
                color = "black",
                stroke = 0.6
            )
    } else {
        n_grp <- length(levels(df$grp))
        if (n_grp > length(pal)) {
            stop("Too many groups in ", color_by, " (", n_grp, "). Increase palette length.")
        }

        p <- ggplot2::ggplot(df, ggplot2::aes(x = num_features, y = mae, fill = grp)) +
            ggplot2::geom_point(
                size = point_size,
                alpha = alpha,
                shape = 21,
                color = "black",
                stroke = 0.5
            ) +
            ggplot2::scale_fill_manual(values = pal[seq_len(n_grp)], name = color_by)
    }

    p <- p +
        ggplot2::geom_smooth(
            method = "loess",
            span = smoother_span,
            se = FALSE,
            color = "grey25",
            linewidth = 1.1,
            inherit.aes = FALSE,
            mapping = ggplot2::aes(x = num_features, y = mae),
            data = df
        ) +
        ggplot2::labs(
            x = "Number of age-informative genes",
            y = "Mean absolute error (years)"
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            axis.title = ggplot2::element_text(size = 13),
            axis.text = ggplot2::element_text(size = 11),
            legend.position = if (show_legend && !is.null(color_by)) "right" else "none",
            legend.title = ggplot2::element_text(size = 12),
            legend.text = ggplot2::element_text(size = 10)
        )

    p
}

.make_palette <- function(n) {
    if (requireNamespace("colorspace", quietly = TRUE)) {
        # "Dark 3" can collide; "Set 3" often includes pale/yellowish tones.
        # "Glasbey" is ideal but not always available. This is a good compromise.
        colorspace::qualitative_hcl(n, palette = "Dark 2")
    } else {
        # Hand-curated, distinct, no yellow, no black/grey.
        c(
            "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
            "#66A61E", "#E6AB02", "#A6761D", "#666666",
            "#1F78B4", "#B2DF8A", "#FB9A99", "#CAB2D6"
        )[seq_len(n)]
    }
}

plot_error_vs_predictor <- function(metrics_df,
                                    x_col,
                                    y_col = "mae",
                                    x_label = NULL,
                                    y_label = NULL,
                                    x_scale = c("identity", "log10"),
                                    smoother_span = 1,
                                    add_smoother = TRUE,
                                    annotate_adj_r2 = TRUE,
                                    point_size = 2.6,
                                    alpha = 0.9,
                                    annotate_r2_size = 3.6,
                                    log_pseudocount = 1,
                                    show_y_label = TRUE) {

    stopifnot(is.data.frame(metrics_df))
    x_scale <- match.arg(x_scale)

    need <- c(x_col, y_col)
    miss <- setdiff(need, colnames(metrics_df))
    if (length(miss) > 0) {
        stop("Missing required columns: ", paste(miss, collapse = ", "))
    }

    num_x <- y_val <- NULL

    df <- metrics_df
    df$num_x <- df[[x_col]]
    df$y_val <- df[[y_col]]

    if (is.null(x_label)) {
        x_label <- gsub("_", " ", x_col)
    } else {
        x_label <- gsub("_", " ", x_label)
    }

    if (is.null(y_label)) {
        if (tolower(y_col) == "mae") {
            y_label <- "Mean AE (years)"
        } else if (tolower(y_col) == "medae") {
            y_label <- "Median AE (years)"
        } else if (tolower(y_col) == "mae_young20") {
            y_label <- "Mean AE in youngest 20% (years)"
        } else if (tolower(y_col) == "medae_young20") {
            y_label <- "Median AE in youngest 20% (years)"
        } else {
            y_label <- gsub("_", " ", y_col)
        }
    } else {
        y_label <- gsub("_", " ", y_label)
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x = num_x, y = y_val)) +
        ggplot2::geom_point(
            size = point_size,
            alpha = alpha,
            shape = 21,
            fill = "grey70",
            color = "black",
            stroke = 0.5
        ) +
        ggplot2::labs(
            x = x_label,
            y = if (isTRUE(show_y_label)) y_label else NULL
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            axis.title = ggplot2::element_text(size = 13),
            axis.text = ggplot2::element_text(size = 11)
        )

    if (x_scale == "log10") {
        p <- p + ggplot2::scale_x_log10()
    }

    if (add_smoother) {
        p <- p + ggplot2::geom_smooth(
            method = "loess",
            span = smoother_span,
            se = FALSE,
            color = "grey25",
            linewidth = 1.1
        )
    }

    if (annotate_adj_r2) {
        x_raw <- df$num_x
        x2 <- log10(pmax(x_raw, log_pseudocount, na.rm = TRUE))
        z_x <- as.numeric(scale(x2))

        fit_df <- data.frame(y = df$y_val, z_x = z_x)
        fit <- stats::lm(y ~ z_x, data = fit_df)
        adj_r2 <- summary(fit)$adj.r.squared

        xr <- range(df$num_x, na.rm = TRUE)
        yr <- range(df$y_val, na.rm = TRUE)
        ann_x <- xr[1] + 0.70 * (xr[2] - xr[1])
        ann_y <- yr[1] + 0.92 * (yr[2] - yr[1])

        p <- p + ggplot2::annotate(
            "text",
            x = ann_x,
            y = ann_y,
            label = sprintf("adj R2 = %.2f", adj_r2),
            hjust = 0,
            size = annotate_r2_size,
            color = "grey20"
        )
    }

    p
}
