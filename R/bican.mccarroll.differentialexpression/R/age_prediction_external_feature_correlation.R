age_preds_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/age_prediction/age_prediction_region_alpha_0/age_prediction_results_alpha0_donor_predictions.txt"
ctp_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/cell_type_proportions/LEVEL_1/donor_region.cell_type_proportions.txt"
cellTypeListFile="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/metadata/mash_cell_type_list_simple.txt"
cell_type="OPC"

# For each cell type and region, does the predicted age correlate with the OPC % in that region?
# Compare the predicted age to the chronological age
# Also compare the z-score (pred age) and z-score (chron age)
# I think the comparisons should remain "within" region, but allow any cell type age prediction to predict OPC.

correlate_ctp<-function (age_preds_file, ctp_file, cell_type="OPC", cellTypeListFile=NULL) {
    age_preds=read.table(age_preds_file, header=TRUE, sep="\t", stringsAsFactors = FALSE)
    ap=unique(age_preds[,c("cell_type", "region", "donor", "age", "ctp_fraction", "pred_mean_corrected")])


    ctp=read.table(ctp_file, header=TRUE, sep="\t", stringsAsFactors = FALSE)

    #filter age models to requested cell types
    if (!is.null(cellTypeListFile)) {
        cell_types = read.table(cellTypeListFile,
                                header = FALSE,
                                stringsAsFactors = FALSE)$V1
        age_preds = age_preds[age_preds$cell_type %in% cell_types, ]
    }

    #make a joint key for each input that is the donor+region
    age_preds$key=paste0(age_preds$donor, "_", age_preds$region)
    ctp$key=paste0(ctp$donor_external_id, "_", ctp$brain_region_abbreviation_simple)

    if (!cell_type %in% ctp$annotation_most_specific)
        stop ("Cell Type Requested [", cell_type, "] not found in CTP data")

    ctp=ctp[ctp$annotation_most_specific==cell_type, ]

    common_keys=intersect(age_preds$key, ctp$key)
    ctp=ctp[match(common_keys, ctp$key), ]

    # add the cell type proportions as a new column
    idx=match(age_preds$key, ctp$key)
    age_preds[,"ctp_fraction"]=ctp$fraction_nuclei[idx]*100

    #add z-score scalings. - these should be by cell type / region
    age_preds=add_zscores(age_preds)

    #get correlations
    cor_df <- calc_group_correlations(age_preds)

    plot_cor_pairs_points(cor_df, title = paste0("Correlation of ", cell_type, " CTP with Age Predictions"))

    pdf ("/downloads/steve_plots.pdf", width=9, height=9)
    for (region in unique(age_preds$region)) {
        p1<-plot_ctp_fraction_vs_age(
            age_preds,
            cell_type_select = "oligodendrocyte",
            region_select = region,
            fraction_label = cell_type
        )

        p2<-plot_ctp_fraction_vs_age(
            age_preds,
            cell_type_select = "OPC",
            region_select = region,
            fraction_label = cell_type
        )
        print (cowplot::plot_grid(p1, p2, nrow=2))
    }

    dev.off()

    plot_ctp_fraction_vs_age(
        age_preds,
        cell_type_select = "oligodendrocyte",
        region_select = "DFC",
        fraction_label = cell_type
    )

    #try to residualize OPC fraction with regard to age, then plot against the predicted age residuals.
    res <- calc_opc_age_resid_cor(age_preds)

    plot_opc_age_resid_cor(res$summary_df, title = "Residualized CTP vs resid_mean_corrected")

    # subset to one group and plot the two vectors that generate the correlation
    plot_one(res, cell_type = "astrocyte", region = "NAC")

    #how are the 90 year old donors being predicted.
    plot_predicted_age_dist_at_age(age_preds, age_value=9)

    return(results)
}

plot_predicted_age_dist_at_age <- function(
        age_preds,
        age_value = 9,
        geom = c("density", "histogram"),
        bins = 30,
        title = NULL,
        alpha = 0.5,
        ref_line = TRUE
) {
    geom <- match.arg(geom)

    required_cols <- c("cell_type", "region", "age", "pred_mean_corrected", "pred_mean")
    missing_cols <- setdiff(required_cols, names(age_preds))
    if (length(missing_cols) > 0) {
        stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
    }

    if (!is.numeric(age_value) || length(age_value) != 1) {
        stop("age_value must be a single numeric value")
    }

    idx <- is.finite(age_preds$age) &
        is.finite(age_preds$pred_mean_corrected) &
        (age_preds$age == age_value)

    df <- age_preds[idx, , drop = FALSE]
    if (nrow(df) == 0) {
        stop("No rows matched age == age_value")
    }

    if (is.null(title)) {
        title <- paste0("Predicted age distribution at age ", age_value)
    }

    # Compute symmetric x-axis limits around age_value
    max_dev <- max(abs(df$pred_mean_corrected - age_value), na.rm = TRUE)
    x_limits <- c(age_value - max_dev, age_value + max_dev)

    cell_type <- region <- NULL

    p <- ggplot2::ggplot(
        df,
        ggplot2::aes(
            x = pred_mean_corrected,
            color = region,
            fill = region
        )
    ) +
        ggplot2::facet_wrap(~ cell_type) +
        ggplot2::labs(
            title = title,
            x = "Predicted age",
            y = if (geom == "density") "Density" else "Count",
            color = "Region",
            fill = "Region"
        ) +
        ggplot2::scale_x_continuous(limits = x_limits) +
        ggplot2::theme_classic()

    if (isTRUE(ref_line)) {
        p <- p +
            ggplot2::geom_vline(
                xintercept = age_value,
                linetype = "dashed",
                linewidth = 0.5
            )
    }

    if (geom == "density") {
        p <- p +
            ggplot2::geom_density(alpha = alpha)
    } else {
        p <- p +
            ggplot2::geom_histogram(
                bins = bins,
                position = "identity",
                alpha = alpha
            )
    }

    p
}



plot_one<-function (res,cell_type,  region) {
    # subset to one group and plot the two vectors that generate the correlation
    one <- res$resid_df[res$resid_df$cell_type == cell_type & res$resid_df$region == region, , drop = FALSE]

    ggplot2::ggplot(
        one,
        ggplot2::aes(x = resid_mean_corrected, y = ctp_age_resid)
    ) +
        ggplot2::geom_point() +
        ggplot2::geom_smooth(method = "lm", se = FALSE) +
        ggplot2::theme_classic() +
        ggplot2::labs(
            x = "resid_mean_corrected",
            y = "residuals of ctp_fraction ~ age",
            title = paste0(cell_type, " / ", region)
        )

}

plot_ctp_fraction_vs_age <- function(
        age_preds,
        cell_type_select,
        region_select,
        fraction_label,
        title = NULL,
        point_size = 2.3,
        add_cor_labels = TRUE,
        cor_method = "pearson"
) {
    required_cols <- c(
        "cell_type", "region",
        "age", "pred_mean_corrected",
        "ctp_fraction"
    )
    missing_cols <- setdiff(required_cols, names(age_preds))
    if (length(missing_cols) > 0) {
        stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
    }

    df <- age_preds[
        age_preds$cell_type == cell_type_select &
            age_preds$region == region_select,
        ,
        drop = FALSE
    ]

    if (nrow(df) == 0) {
        stop("No data for requested cell_type / region combination")
    }

    if (is.null(title)) {
        title <- paste0(
            fraction_label,
            " fraction vs age and predicted age\n",
            cell_type_select,
            " / ",
            region_select
        )
    }

    plot_df <- rbind(
        data.frame(
            x = df$age,
            ctp_fraction = df$ctp_fraction,
            age_type = "chronological age",
            stringsAsFactors = FALSE
        ),
        data.frame(
            x = df$pred_mean_corrected,
            ctp_fraction = df$ctp_fraction,
            age_type = "predicted age",
            stringsAsFactors = FALSE
        )
    )

    cor_age <- NA_real_
    cor_pred <- NA_real_
    n_age <- 0L
    n_pred <- 0L

    ok_age <- is.finite(df$ctp_fraction) & is.finite(df$age)
    n_age <- sum(ok_age)
    if (n_age >= 3) {
        cor_age <- stats::cor(df$ctp_fraction[ok_age], df$age[ok_age], method = cor_method)
    }

    ok_pred <- is.finite(df$ctp_fraction) & is.finite(df$pred_mean_corrected)
    n_pred <- sum(ok_pred)
    if (n_pred >= 3) {
        cor_pred <- stats::cor(df$ctp_fraction[ok_pred], df$pred_mean_corrected[ok_pred], method = cor_method)
    }

    cor_label <- NULL
    if (isTRUE(add_cor_labels)) {
        cor_label <- paste0(
            "cor(CTP, age) = ", sprintf("%.3f", cor_age), " (n=", n_age, ")\n",
            "cor(CTP, predicted) = ", sprintf("%.3f", cor_pred), " (n=", n_pred, ")"
        )
    }

    age_type <- NULL

    p <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(
            x = x,
            y = ctp_fraction,
            color = age_type
        )
    ) +
        ggplot2::geom_point(size = point_size) +
        ggplot2::geom_smooth(
            method = "lm",
            se = FALSE,
            linewidth = 0.6
        ) +
        ggplot2::labs(
            title = title,
            x = "Age",
            y = paste0(fraction_label, " fraction"),
            color = NULL
        ) +
        ggplot2::theme_classic()

    if (isTRUE(add_cor_labels)) {
        p <- p +
            ggplot2::annotate(
                "text",
                x = -Inf,
                y = Inf,
                label = cor_label,
                hjust = -0.05,
                vjust = 1.05,
                size = 3.2
            )
    }

    p
}


calc_opc_age_resid_cor <- function(age_preds) {
    required_cols <- c("cell_type", "region", "ctp_fraction", "age", "resid_mean_corrected")
    missing_cols <- setdiff(required_cols, names(age_preds))
    if (length(missing_cols) > 0) {
        stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
    }

    groups <- interaction(age_preds$cell_type, age_preds$region, drop = TRUE)
    group_levels <- levels(groups)

    # Copy so we can add residual columns aligned to original rows
    df_out <- age_preds
    df_out$ctp_age_resid <- NA_real_

    summary_list <- vector("list", length(group_levels))

    for (i in seq_along(group_levels)) {
        g <- group_levels[i]
        idx <- groups == g

        # Fit on complete cases for the residual model
        ok_fit <- idx & is.finite(df_out$ctp_fraction) & is.finite(df_out$age)
        n_fit <- sum(ok_fit)

        if (n_fit >= 2) {
            fit <- stats::lm(ctp_fraction ~ age, data = df_out[ok_fit, , drop = FALSE])
            df_out$ctp_age_resid[ok_fit] <- stats::residuals(fit)
        }

        # Correlate residuals with resid_mean_corrected on complete pairs
        ok_cor <- idx & is.finite(df_out$ctp_age_resid) & is.finite(df_out$resid_mean_corrected)
        n_cor <- sum(ok_cor)

        ct_name <- df_out$cell_type[which(idx)[1]]
        rg_name <- df_out$region[which(idx)[1]]

        if (n_cor < 3) {
            summary_list[[i]] <- data.frame(
                cell_type = ct_name,
                region = rg_name,
                cor_resid = NA_real_,
                p = NA_real_,
                ci_lo = NA_real_,
                ci_hi = NA_real_,
                n_cor = n_cor,
                n_fit = n_fit,
                stringsAsFactors = FALSE
            )
            next
        }

        ct <- suppressWarnings(stats::cor.test(
            df_out$ctp_age_resid[ok_cor],
            df_out$resid_mean_corrected[ok_cor],
            method = "pearson"
        ))

        summary_list[[i]] <- data.frame(
            cell_type = ct_name,
            region = rg_name,
            cor_resid = unname(ct$estimate),
            p = ct$p.value,
            ci_lo = unname(ct$conf.int[1]),
            ci_hi = unname(ct$conf.int[2]),
            n_cor = n_cor,
            n_fit = n_fit,
            stringsAsFactors = FALSE
        )
    }

    summary_df <- do.call(rbind, summary_list)

    list(
        resid_df = df_out,
        summary_df = summary_df
    )
}


add_zscores<-function (age_preds) {
    z_by_group <- function(x) {
        as.numeric(scale(x))
    }

    age_preds$age_z <-
        ave(
            age_preds$age,
            age_preds$cell_type,
            age_preds$region,
            FUN = z_by_group
        )

    age_preds$pred_mean_corrected_z <-
        ave(
            age_preds$pred_mean_corrected,
            age_preds$cell_type,
            age_preds$region,
            FUN = z_by_group
        )

    return(age_preds)

}

calc_group_correlations <- function(age_preds) {
    required_cols <- c(
        "cell_type", "region", "ctp_fraction",
        "age", "pred_mean_corrected", "resid_mean_corrected"
    )
    missing_cols <- setdiff(required_cols, names(age_preds))
    if (length(missing_cols) > 0) {
        stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
    }

    groups <- interaction(age_preds$cell_type, age_preds$region, drop = TRUE)
    group_levels <- levels(groups)

    cor_test_summary <- function(x, y, conf_level = 0.95) {
        ok <- is.finite(x) & is.finite(y)
        x <- x[ok]
        y <- y[ok]
        n_ok <- length(x)

        if (n_ok < 3) {
            return(list(
                cor = NA_real_,
                p = NA_real_,
                lo = NA_real_,
                hi = NA_real_,
                n = n_ok
            ))
        }

        ct <- suppressWarnings(stats::cor.test(x, y, conf.level = conf_level, method = "pearson"))

        list(
            cor = unname(ct$estimate),
            p = ct$p.value,
            lo = unname(ct$conf.int[1]),
            hi = unname(ct$conf.int[2]),
            n = n_ok
        )
    }

    out_list <- vector("list", length(group_levels))

    for (i in seq_along(group_levels)) {
        g <- group_levels[i]
        idx <- groups == g

        ct_name <- age_preds$cell_type[which(idx)[1]]
        rg_name <- age_preds$region[which(idx)[1]]

        ctp <- age_preds$ctp_fraction[idx]

        s_age <- cor_test_summary(ctp, age_preds$age[idx])
        s_pmc <- cor_test_summary(ctp, age_preds$pred_mean_corrected[idx])
        s_rmc <- cor_test_summary(ctp, age_preds$resid_mean_corrected[idx])

        out_list[[i]] <- data.frame(
            cell_type = ct_name,
            region = rg_name,

            cor_ctp_age = s_age$cor,
            cor_ctp_age_p = s_age$p,
            cor_ctp_age_ci_lo = s_age$lo,
            cor_ctp_age_ci_hi = s_age$hi,
            n_ctp_age = s_age$n,

            cor_ctp_pred_mean_corrected = s_pmc$cor,
            cor_ctp_pred_mean_corrected_p = s_pmc$p,
            cor_ctp_pred_mean_corrected_ci_lo = s_pmc$lo,
            cor_ctp_pred_mean_corrected_ci_hi = s_pmc$hi,
            n_ctp_pred_mean_corrected = s_pmc$n,

            cor_ctp_resid_mean_corrected = s_rmc$cor,
            cor_ctp_resid_mean_corrected_p = s_rmc$p,
            cor_ctp_resid_mean_corrected_ci_lo = s_rmc$lo,
            cor_ctp_resid_mean_corrected_ci_hi = s_rmc$hi,
            n_ctp_resid_mean_corrected = s_rmc$n,

            stringsAsFactors = FALSE
        )
    }

    do.call(rbind, out_list)
}

plot_cor_pairs_points <- function(
        cor_df,
        title = NULL,
        point_size = 2.3,
        add_identity = TRUE,
        add_zero_lines = FALSE) {
    required_cols <- c("cell_type", "region", "cor_ctp_age", "cor_ctp_pred_mean_corrected")
    missing_cols <- setdiff(required_cols, names(cor_df))
    if (length(missing_cols) > 0)
    {
        stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
    }

    cell_type <- region <- NULL

    p <- ggplot2::ggplot(
        cor_df,
        ggplot2::aes(
            x = cor_ctp_age,
            y = cor_ctp_pred_mean_corrected,
            color = cell_type,
            shape = region
        )
    )

    if (isTRUE(add_identity))
    {
        p <- p +
            ggplot2::geom_abline(
                slope = 1,
                intercept = 0,
                linewidth = 0.5
            )
    }

    if (isTRUE(add_zero_lines))
    {
        p <- p +
            ggplot2::geom_hline(
                yintercept = 0,
                linetype = "dashed",
                linewidth = 0.4
            ) +
            ggplot2::geom_vline(
                xintercept = 0,
                linetype = "dashed",
                linewidth = 0.4
            )
    }

    p +
        ggplot2::geom_point(size = point_size) +
        ggplot2::labs(
            title = title,
            x = "cor(ctp_fraction, age)",
            y = "cor(ctp_fraction, pred_mean_corrected)",
            color = "Cell type",
            shape = "Region"
        ) +
        ggplot2::theme_classic()
}

plot_opc_age_resid_cor <- function(
        summary_df,
        title = NULL,
        point_size = 2.3,
        dodge_width = 0.6,
        add_zero_line = TRUE
) {
    required_cols <- c("cell_type", "region", "cor_resid")
    missing_cols <- setdiff(required_cols, names(summary_df))
    if (length(missing_cols) > 0) {
        stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
    }

    cell_type <- region <- NULL

    p <- ggplot2::ggplot(
        summary_df,
        ggplot2::aes(
            x = cell_type,
            y = cor_resid,
            color = region
        )
    )

    if (isTRUE(add_zero_line)) {
        p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4)
    }

    p +
        ggplot2::geom_point(
            position = ggplot2::position_dodge(width = dodge_width),
            size = point_size
        ) +
        ggplot2::theme_classic() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::labs(
            title = title,
            x = "Cell type",
            y = "cor(resid(CTP ~ age), resid_mean_corrected)",
            color = "Region"
        )
}
