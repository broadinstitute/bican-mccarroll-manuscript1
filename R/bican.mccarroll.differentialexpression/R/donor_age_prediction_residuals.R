################################################################################
# Focus on the residuals of donor age prediction acrsso cell types or regions
################################################################################
# library(ggplot2)
# library(ComplexHeatmap)
# library(circlize)
# library(cowplot)
# library(data.table)
# library(logger)
#
# model_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/age_prediction/age_prediction_region_alpha_0/age_prediction_results_alpha0_model_coefficients.txt"
# prediction_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/age_prediction/age_prediction_region_alpha_0/age_prediction_results_alpha0_donor_predictions.txt"
#
# # Restrict the cell correlation analysis to a subset of cell types
# cellTypeListFile="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/metadata/mash_cell_type_list_simple.txt"
# out_pdf_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/age_prediction/age_prediction_region_alpha_0/age_prediction_region_residual_comparison.pdf"
#
# #only look at the extremes of the distribution.
# exclude_donor_age_range=c(1,8.9)

compare_age_residuals_celltype_region <- function (model_file, prediction_file,
                                            cellTypeListFile = NULL, exclude_donor_age_range=c(),
                                            out_pdf_file = NULL) {

    model_predictions = load_model_predictions(prediction_file, cellTypeListFile)
    if (length(exclude_donor_age_range) == 2) {
        age_min = exclude_donor_age_range[1]
        age_max = exclude_donor_age_range[2]
        logger::log_info(paste0("Excluding donors with age in range [", age_min, ", ", age_max, "] from analysis"))
        model_predictions = model_predictions[!(model_predictions$age >= age_min & model_predictions$age <= age_max), ]
    }
    all_models = load_models(model_file, cellTypeListFile)

    if (!is.null(out_pdf_file))
        pdf(out_pdf_file, width = 11, height = 11)

    #for each region, compute the residual matrix and plot
    regions=unique(model_predictions$region)
    for (current_region in regions) {
        model_predictions_region = model_predictions[model_predictions$region == current_region, ]
        res_mat <- donor_celltype_matrix(model_predictions_region, value_var = "resid_mean")
        p2 <- plot_residual_pair_scatter_paged(
            res_mat,
            cellType = current_region,
            per_page = 4,
            facet_font_size = 8
        )

        #heat map of all cell type residual correlations
        p3 <- plot_residual_corr_heatmap(res_mat, cellType = current_region, annotate_cells = TRUE)

        #heat map of jaccard index of aging program overlap
        model<-all_models[all_models$region == current_region, ]
        J <- jaccard_by_celltype(model, coef_thresh = 0)  # or small error, e.g., 1e-8
        strTitle= paste0("Cell-type gene overlap in aging programs", " region [", current_region, "]")
        ht <- plot_jaccard_heatmap(J, title = strTitle, annotate_cells = TRUE)

        ComplexHeatmap::draw(ht, heatmap_legend_side = "right")
        ComplexHeatmap::draw(p3, heatmap_legend_side = "right")

        for (p in p2) {
            print (p)
        }
    }



    if (!is.null(out_pdf_file))
        dev.off()

}


#####################
# PLOTS
#####################

# 2. Scatter plots with annotated correlations, paged if too many pairs
plot_residual_pair_scatter_paged <- function(res_mat,
                                             cellType = NULL,
                                             per_page = 12,
                                             facet_font_size = 10,
                                             ncol = NULL)
{
    stopifnot(is.matrix(res_mat))
    regs <- colnames(res_mat)

    prs <- if (length(regs) >= 2)
        utils::combn(regs, 2, simplify = FALSE)
    else
        list()

    if (!length(prs))
        stop("Need >=2 cell types (columns) in res_mat")

    make_panel <- function(title, x, y, xlab, ylab)
    {
        keep <- is.finite(x) & is.finite(y)
        if (sum(keep) < 2)
            return(NULL)

        x <- x[keep]
        y <- y[keep]
        r <- stats::cor(x, y)

        rng <- range(c(x, y))
        pad <- diff(rng) * 0.1
        lim <- c(rng[1] - pad, rng[2] + pad)

        ggplot2::ggplot(data.frame(x, y), ggplot2::aes(x, y)) +
            ggplot2::geom_abline(
                intercept = 0,
                slope = 1,
                color = "black",
                linetype = "dashed"
            ) +
            ggplot2::geom_point(
                size = 2,
                alpha = 0.7,
                color = "steelblue"
            ) +
            ggplot2::geom_smooth(
                method = "lm",
                formula = y ~ x,
                se = FALSE,
                linewidth = 0.6,
                color = "red"
            ) +
            ggplot2::annotate(
                "text",
                x = -Inf,
                y = Inf,
                label = sprintf("r = %.2f", r),
                hjust = -0.1,
                vjust = 1.2,
                size = 3.2
            ) +
            ggplot2::ggtitle(title) +
            ggplot2::coord_cartesian(xlim = lim, ylim = lim) +
            ggplot2::labs(x = xlab, y = ylab) +
            ggplot2::theme_classic(base_size = 12) +
            ggplot2::theme(plot.title = ggplot2::element_text(
                size = facet_font_size,
                hjust = 0.5
            ))
    }

    plots <- list()
    for (p in prs) {
        x_name <- p[1]
        y_name <- p[2]
        title  <- paste(x_name, "vs", y_name)

        pan <- make_panel(
            title = title,
            x = res_mat[, x_name],
            y = res_mat[, y_name],
            xlab = paste("Residual in", x_name),
            ylab = paste("Residual in", y_name)
        )

        if (!is.null(pan))
            plots[[length(plots) + 1]] <- pan
    }

    if (!length(plots))
        stop("No valid pairs after filtering")

    if (is.null(ncol))
        ncol <- ceiling(sqrt(per_page))
    nrow <- ceiling(per_page / ncol)

    blanks <- function(n)
        replicate(n, ggplot2::ggplot() + ggplot2::theme_void(), simplify = FALSE)

    page_title <- "Age Prediction residuals (predicted - actual)"
    if (!is.null(cellType))
        page_title <- paste(cellType, page_title, sep = "\n")

    pages <- list()
    for (s in seq(1, length(plots), by = per_page)) {
        page_plots <- plots[s:min(s + per_page - 1, length(plots))]
        if (length(page_plots) < per_page)
            page_plots <- c(page_plots, blanks(per_page - length(page_plots)))

        grid <- cowplot::plot_grid(
            plotlist = page_plots,
            ncol = ncol,
            nrow = nrow
        )

        pg <- cowplot::ggdraw() +
            cowplot::draw_label(
                page_title,
                x = 0,
                y = 1,
                hjust = 0,
                vjust = 1,
                size = 14
            ) +
            cowplot::draw_plot(grid, y = 0, height = 0.94)

        pages[[length(pages) + 1]] <- pg
    }

    pages
}


# m: donors x regions residual matrix (from compute_residual_matrix)
# method: "pearson" or "spearman"
# cluster: TRUE = hierarchical clustering rows/cols; FALSE = keep input order
# annotate_cells: TRUE = show correlation values in heatmap cells; FALSE = no annotation
plot_residual_corr_heatmap <- function(res_mat,
                                       cellType = NULL,
                                       annotate_cells = TRUE) {
    stopifnot(is.matrix(res_mat))
    C <- cor(res_mat, use = "pairwise.complete.obs")  # regions x regions

    col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("#3b4cc0", "white", "#b40426"))

    column_title <- "Age Prediction residuals (predicted - actual)"
    if (!is.null(cellType))
        column_title <- paste(cellType, column_title, sep = "\n")

    # optional cell annotation
    cf <- if (isTRUE(annotate_cells)) {
        function(j, i, x, y, width, height, fill) {
            grid::grid.text(sprintf("%.2f", C[i, j]), x, y, gp = grid::gpar(fontsize = 12))
        }
    } else
        NULL

    ComplexHeatmap::Heatmap(
        C,
        name = "r",
        col = col_fun,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_title = NULL,
        column_title = column_title,
        show_row_dend = TRUE,
        show_column_dend = TRUE,
        row_names_gp = grid::gpar(fontsize = 9),
        column_names_gp = grid::gpar(fontsize = 9),
        heatmap_legend_param = list(at = c(-1, -0.5, 0, 0.5, 1), title = "Correlation"),
        cell_fun = cf
    )
}


# Build Jaccard matrix of selected features per cell type
jaccard_by_celltype <- function(all_models, coef_thresh = 0) {
    stopifnot(all(c("cell_type", "feature", "coef") %in% names(all_models)))

    # binary incidence matrix: features x cell_types
    sel <- abs(all_models$coef) > coef_thresh
    df  <- all_models[sel, c("feature", "cell_type")]
    df  <- unique(df)  # one row per (feature, cell_type)

    feats <- sort(unique(df$feature))
    cts   <- sort(unique(df$cell_type))
    A <- matrix(
        FALSE,
        nrow = length(feats),
        ncol = length(cts),
        dimnames = list(feats, cts)
    )
    A[cbind(match(df$feature, feats), match(df$cell_type, cts))] <- TRUE

    # Jaccard: J = (A' A) / (n_i + n_j - A' A)
    XtX <- crossprod(A)                           # |intersection|
    n   <- matrix(colSums(A), ncol = ncol(A), nrow = ncol(A))
    J   <- as.matrix(XtX / (n + t(n) - XtX))
    diag(J) <- 1
    J
}

# Heatmap with optional numeric annotation
plot_jaccard_heatmap <- function(J,
                                 title = "Feature overlap (Jaccard Index)",
                                 annotate_cells = TRUE,
                                 cluster = TRUE) {
    stopifnot(
        requireNamespace("ComplexHeatmap", quietly = TRUE),
        requireNamespace("circlize", quietly = TRUE),
        requireNamespace("grid", quietly = TRUE)
    )

    col_fun <- circlize::colorRamp2(c(0, 0.5, 1), c("#f7fbff", "#6baed6", "#08306b"))

    cf <- if (isTRUE(annotate_cells)) {
        function(j, i, x, y, w, h, fill) {
            grid::grid.text(sprintf("%.2f", J[i, j]), x, y, gp = grid::gpar(fontsize = 10))
        }
    } else
        NULL

    ComplexHeatmap::Heatmap(
        J,
        name = "Jaccard",
        col = col_fun,
        cluster_rows = cluster,
        cluster_columns = cluster,
        show_row_dend = cluster,
        show_column_dend = cluster,
        row_names_gp = grid::gpar(fontsize = 10),
        column_names_gp = grid::gpar(fontsize = 10),
        heatmap_legend_param = list(at = c(0, 0.5, 1)),
        cell_fun = cf,
        column_title = title
    )
}











donor_celltype_matrix <- function(model_predictions, value_var) {
    data.table::setDT(model_predictions)

    wide_dt <- data.table::dcast(
        model_predictions,
        donor ~ cell_type,
        value.var = value_var
    )

    mat <- as.matrix(wide_dt[, -1])
    rownames(mat) <- wide_dt$donor

    mat
}


load_models <- function (model_file, cellTypeListFile=NULL) {

    all_models=data.table::fread(model_file,
          header = TRUE,
          sep = "\t",
          stringsAsFactors = FALSE)

    if (!is.null(cellTypeListFile)) {
        cell_types = read.table(cellTypeListFile,
                                header = FALSE,
                                stringsAsFactors = FALSE)$V1
        all_models = all_models[all_models$cell_type %in% cell_types, ]
    }

    strCell=paste0("[", length(unique(all_models$cell_type)), "] cell types")
    strRegion=paste0("[", length(unique(all_models$region)), "] regions")
    str=paste0("Loaded model coefficients for ", strCell, " ", strRegion)
    logger::log_info(str)
    return(all_models)
}

load_model_predictions <- function (prediction_file, cellTypeListFile=NULL) {


    all_preds=data.table::fread(prediction_file,
                                 header = TRUE,
                                 sep = "\t",
                                 stringsAsFactors = FALSE)

    if (!is.null(cellTypeListFile)) {
        cell_types = read.table(cellTypeListFile,
                                header = FALSE,
                                stringsAsFactors = FALSE)$V1
        all_preds = all_preds[all_preds$cell_type %in% cell_types, ]
    }


    strCell=paste0("[", length(unique(all_preds$cell_type)), "] cell types")
    strRegion=paste0("[", length(unique(all_preds$region)), "] regions")
    strDonor=paste0("[", length(unique(all_preds$donor)), "] donors")
    str=paste0("Loaded donor predictions for ", strCell, " ", strRegion, " ", strDonor)

    logger::log_info(str)
    return(all_preds)
}
