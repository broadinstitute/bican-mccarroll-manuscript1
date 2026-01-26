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
# value_var="resid_mean"
# out_pdf_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/age_prediction/age_prediction_region_alpha_0/age_prediction_region_residual_comparison.pdf"
# exclude_donor_age_range=c()

# #only look at the extremes of the distribution.
# exclude_donor_age_range=c(1,8.9)


#' Compare donor age residuals across cell types and regions
#'
#' Loads saved age-prediction outputs (models and per-donor predictions) and
#' generates diagnostic plots comparing donor-level residuals (or other
#' donor-level summary values) across:
#' \itemize{
#'   \item cell types within each region, and
#'   \item regions within each cell type.
#' }
#'
#' For each region, donor-level values are spread into a donor-by-cell-type
#' matrix and pairwise relationships are plotted (via
#' \code{plot_residuals_for_subset()}). The same analysis is repeated for each
#' cell type across regions.
#'
#' Optionally excludes donors whose ages fall within a specified closed interval
#' \code{[age_min, age_max]} prior to plotting.
#'
#' If \code{out_pdf_file} is provided, plots are written to a PDF device and the
#' device is closed on exit.
#'
#' @param model_file Path to a file containing serialized age models (passed to
#'   \code{load_models()}).
#' @param prediction_file Path to a file containing serialized model predictions
#'   (passed to \code{load_model_predictions()}).
#' @param cellTypeListFile Optional path to a cell type list file used by
#'   \code{load_model_predictions()} and \code{load_models()} to restrict or
#'   order cell types.
#' @param value_var Character scalar naming the column in the predictions table
#'   to compare (default \code{"resid_mean"}). Examples include
#'   \code{"resid_mean"} or \code{"resid_mean_corrected"} if available.
#' @param exclude_donor_age_range Numeric vector of length 2 specifying an age
#'   interval \code{c(age_min, age_max)}. Donors with ages in the closed interval
#'   \code{[age_min, age_max]} are excluded. If not length 2, no exclusion is
#'   performed.
#' @param out_pdf_file Optional output PDF filepath. If provided, opens a PDF
#'   device with \code{width = 11} and \code{height = 11}, writes plots, and
#'   closes the device.
#'
#' @return Invisibly returns \code{NULL}. This function is called for its side
#'   effects (plot generation and optional PDF output).
#'
#' @export
compare_age_residuals_celltype_region <- function(model_file, prediction_file,
                                                  cellTypeListFile = NULL,
                                                  value_var="resid_mean",
                                                  exclude_donor_age_range = c(),
                                                  out_pdf_file = NULL) {

    model_predictions <- load_model_predictions(prediction_file, cellTypeListFile)

    if (length(exclude_donor_age_range) == 2) {
        age_min <- exclude_donor_age_range[1]
        age_max <- exclude_donor_age_range[2]
        logger::log_info(paste0(
            "Excluding donors with age in range [", age_min, ", ", age_max, "] from analysis"
        ))
        model_predictions <- model_predictions[!(
            model_predictions$age >= age_min & model_predictions$age <= age_max
        ), ]
    }

    all_models <- load_models(model_file, cellTypeListFile)

    if (!is.null(out_pdf_file))
        grDevices::pdf(out_pdf_file, width = 11, height = 11)

    # 1) per region: correlate residuals across cell types
    for (current_region in unique(model_predictions$region)) {
        logger::log_info("Processing region: ", current_region)
        preds_sub  <- model_predictions[model_predictions$region == current_region, ]
        models_sub <- all_models[all_models$region == current_region, ]

            plot_residuals_for_subset(
            model_predictions_subset = preds_sub,
            all_models_subset = models_sub,
            value_var=value_var,
            spread_var = "cell_type",
            jaccard_group_var = "cell_type",
            title_label = paste0("region [", current_region, "]"),
            annotate_cells = TRUE,
            per_page = 4,
            facet_font_size = 8,
            coef_thresh = 0
        )
    }

    # 2) per cell type: correlate residuals across regions
    for (current_cell_type in unique(model_predictions$cell_type)) {
        logger::log_info("Processing cell type: ", current_cell_type)
        preds_sub  <- model_predictions[model_predictions$cell_type == current_cell_type, ]
        models_sub <- all_models[all_models$cell_type == current_cell_type, ]

        plot_residuals_for_subset(
            model_predictions_subset = preds_sub,
            all_models_subset = models_sub,
            value_var=value_var,
            spread_var = "region",
            jaccard_group_var = "region",
            title_label = paste0("cell_type [", current_cell_type, "]"),
            annotate_cells = TRUE,
            per_page = 4,
            facet_font_size = 8,
            coef_thresh = 0
        )
    }

    #is there a correlation between the jaccard index of a pair of cell types and their residual correlation?

    #Are residuals more correlated across regions or cell types?  Compute the median for each cell type and median for each region
    #and plot the distributions.

    if (!is.null(out_pdf_file))
        grDevices::dev.off()
}

plot_residuals_for_subset <- function(model_predictions_subset,
                                      all_models_subset,
                                      value_var="resid_mean",
                                      spread_var,
                                      jaccard_group_var,
                                      title_label,
                                      annotate_cells = TRUE,
                                      per_page = 4,
                                      facet_font_size = 8,
                                      coef_thresh = 0) {



    row_fontsize=16
    col_fontsize=16
    cell_fontsize=16

    # Defensive: need at least 2 groups to compare (2+ cell types in a region, or 2+ regions in a cell type)
    n_groups <- length(unique(model_predictions_subset[[spread_var]]))
    if (is.na(n_groups) || n_groups < 2) {
        logger::log_info(paste0(
            "Skipping subset ", title_label, ": only ", n_groups, " unique ", spread_var, " present"
        ))
        return(invisible(NULL))
    }

    res_mat <- donor_wide_matrix(
        model_predictions=model_predictions_subset,
        spread_var = spread_var,
        value_var = value_var
    )

    donor_age <- unique(model_predictions_subset[, c("donor", "age")])

    p2 <- plot_residual_pair_scatter_paged(
        res_mat,
        cellType = title_label,
        per_page = per_page,
        facet_font_size = facet_font_size,
        donor_meta = donor_age,
        donor_id_col = "donor",
        color_var = "age",
        color_title = "Donor age"
    )

    p3 <- plot_residual_corr_heatmap(
        res_mat,
        cellType = title_label,
        annotate_cells = annotate_cells,
        row_fontsize = row_fontsize,
        col_fontsize = col_fontsize,
        cell_fontsize = cell_fontsize
    )

    # Feature overlap across the dimension actually being compared
    J <- jaccard_by_group(all_models_subset, group_var = jaccard_group_var, coef_thresh = coef_thresh)

    overlap_label <- if (identical(jaccard_group_var, "cell_type")) "Cell-type" else jaccard_group_var
    strTitle <- paste0(overlap_label, " gene overlap in aging programs ", title_label)
    ht <- plot_jaccard_heatmap(J, title = strTitle, annotate_cells = annotate_cells,
                               row_fontsize = row_fontsize, col_fontsize = col_fontsize,
                              cell_fontsize = cell_fontsize)

    cc <- coef_corr_on_intersect(all_models_subset, group_var = jaccard_group_var)

    ht_c <- plot_coef_corr_heatmap(
        coef_correlation      = cc$coef_correlation,
        overlap_gene_counts   = cc$overlap_gene_counts,
        title = paste0(jaccard_group_var, " coef correlation (intersect genes) ", title_label),
        annotate_cells = annotate_cells,
        row_fontsize = row_fontsize,
        col_fontsize = col_fontsize,
        cell_fontsize = cell_fontsize
    )

    ComplexHeatmap::draw(ht, heatmap_legend_side = "right")
    ComplexHeatmap::draw(p3, heatmap_legend_side = "right")
    ComplexHeatmap::draw(ht_c, heatmap_legend_side = "right")

    for (p in p2) print(p)

    invisible(list(res_mat = res_mat, corr_heatmap = p3, jaccard = J))
}

# Generalized: donors x <spread_var> matrix for a chosen value column
donor_wide_matrix <- function(model_predictions, spread_var, value_var) {
    data.table::setDT(model_predictions)

    f <- stats::as.formula(paste0("donor ~ ", spread_var))
    wide_dt <- data.table::dcast(
        model_predictions,
        formula = f,
        value.var = value_var
    )

    mat <- as.matrix(wide_dt[, -1])
    rownames(mat) <- wide_dt$donor
    mat
}


# calculate correlation on the model coefficients that overlap in a pair of comparisons.
# feature x <group_var> matrix of coefficients (NA where a feature is absent)
coef_matrix_by_group <- function(all_models, group_var) {
    stopifnot(
        all(c("feature", "coef") %in% names(all_models)),
        group_var %in% names(all_models)
    )

    cols <- c("feature", group_var, "coef")

    if (data.table::is.data.table(all_models)) {
        df <- all_models[, ..cols]
    } else {
        df <- all_models[, cols, drop = FALSE]
    }

    f <- stats::as.formula(paste0("feature ~ ", group_var))
    wide_dt <- data.table::dcast(data.table::as.data.table(df), f, value.var = "coef")

    mat <- as.matrix(wide_dt[, -1, with = FALSE])
    rownames(mat) <- wide_dt[["feature"]]
    mat
}

# Pairwise correlation across groups, using intersected genes for each pair
# Returns list(C=cor_matrix, N=overlap_counts)
coef_corr_on_intersect <- function(all_models, group_var, method = "pearson") {
    stopifnot(
        all(c("feature", "coef") %in% names(all_models)),
        group_var %in% names(all_models)
    )

    # feature x group coefficient matrix (NA if gene absent)
    mat <- coef_matrix_by_group(all_models, group_var = group_var)

    # pairwise correlation on intersecting genes
    coef_correlation <- stats::cor(
        mat,
        use = "pairwise.complete.obs",
        method = method
    )

    # number of intersecting genes per pair
    keep <- is.finite(mat)
    overlap_gene_counts <- as.matrix(t(keep) %*% keep)

    list(
        coef_correlation = coef_correlation,
        overlap_gene_counts = overlap_gene_counts
    )
}



#####################
# PLOTS
#####################

plot_coef_corr_heatmap <- function(coef_correlation,
                                   overlap_gene_counts,
                                   title = "Coefficient correlation (intersect genes)",
                                   annotate_cells = TRUE,
                                   cluster_if_complete = TRUE,
                                   row_fontsize=10,
                                   col_fontsize=10,
                                   cell_fontsize=10) {
    stopifnot(
        requireNamespace("ComplexHeatmap", quietly = TRUE),
        requireNamespace("circlize", quietly = TRUE),
        requireNamespace("grid", quietly = TRUE)
    )

    # Only cluster if the matrix is complete (no NA/NaN/Inf)
    complete_mat <- all(is.finite(coef_correlation))
    do_cluster <- isTRUE(cluster_if_complete) && complete_mat

    col_fun <- circlize::colorRamp2(
        c(-1, 0, 1),
        c("#3b4cc0", "white", "#b40426")
    )

    cell_fun <- if (isTRUE(annotate_cells)) {
        function(j, i, x, y, w, h, fill) {
            r <- coef_correlation[i, j]
            n <- overlap_gene_counts[i, j]
            lab <- if (!is.finite(r)) {
                sprintf("NA\n(n=%d)", n)
            } else {
                sprintf("%.2f\n(n=%d)", r, n)
            }
            grid::grid.text(lab, x, y, gp = grid::gpar(fontsize = cell_fontsize))
        }
    } else {
        NULL
    }

    ComplexHeatmap::Heatmap(
        coef_correlation,
        name = "r",
        col = col_fun,
        na_col = "grey90",
        cluster_rows = do_cluster,
        cluster_columns = do_cluster,
        show_row_dend = do_cluster,
        show_column_dend = do_cluster,
        row_names_gp = grid::gpar(fontsize = row_fontsize),
        column_names_gp = grid::gpar(fontsize = col_fontsize),
        heatmap_legend_param = list(
            at = c(-1, -0.5, 0, 0.5, 1),
            title = "Correlation"
        ),
        cell_fun = cell_fun,
        column_title = title
    )
}


# 2. Scatter plots with annotated correlations, paged if too many pairs
# plot_residual_pair_scatter_paged <- function(res_mat,
#                                              cellType = NULL,
#                                              per_page = 12,
#                                              facet_font_size = 10,
#                                              ncol = NULL)
# {
#     stopifnot(is.matrix(res_mat))
#     regs <- colnames(res_mat)
#
#     prs <- if (length(regs) >= 2)
#         utils::combn(regs, 2, simplify = FALSE)
#     else
#         list()
#
#     if (!length(prs))
#         stop("Need >=2 cell types (columns) in res_mat")
#
#     make_panel <- function(title, x, y, xlab, ylab)
#     {
#         keep <- is.finite(x) & is.finite(y)
#         if (sum(keep) < 2)
#             return(NULL)
#
#         x <- x[keep]
#         y <- y[keep]
#         r <- stats::cor(x, y)
#
#         rng <- range(c(x, y))
#         pad <- diff(rng) * 0.1
#         lim <- c(rng[1] - pad, rng[2] + pad)
#
#         ggplot2::ggplot(data.frame(x, y), ggplot2::aes(x, y)) +
#             ggplot2::geom_abline(
#                 intercept = 0,
#                 slope = 1,
#                 color = "black",
#                 linetype = "dashed"
#             ) +
#             ggplot2::geom_point(
#                 size = 2,
#                 alpha = 0.7,
#                 color = "steelblue"
#             ) +
#             ggplot2::geom_smooth(
#                 method = "lm",
#                 formula = y ~ x,
#                 se = FALSE,
#                 linewidth = 0.6,
#                 color = "red"
#             ) +
#             ggplot2::annotate(
#                 "text",
#                 x = -Inf,
#                 y = Inf,
#                 label = sprintf("r = %.2f", r),
#                 hjust = -0.1,
#                 vjust = 1.2,
#                 size = 3.2
#             ) +
#             ggplot2::ggtitle(title) +
#             ggplot2::coord_cartesian(xlim = lim, ylim = lim) +
#             ggplot2::labs(x = xlab, y = ylab) +
#             ggplot2::theme_classic(base_size = 12) +
#             ggplot2::theme(plot.title = ggplot2::element_text(
#                 size = facet_font_size,
#                 hjust = 0.5
#             ))
#     }
#
#     plots <- list()
#     for (p in prs) {
#         x_name <- p[1]
#         y_name <- p[2]
#         title  <- paste(x_name, "vs", y_name)
#
#         pan <- make_panel(
#             title = title,
#             x = res_mat[, x_name],
#             y = res_mat[, y_name],
#             xlab = paste("Residual in", x_name),
#             ylab = paste("Residual in", y_name)
#         )
#
#         if (!is.null(pan))
#             plots[[length(plots) + 1]] <- pan
#     }
#
#     if (!length(plots))
#         stop("No valid pairs after filtering")
#
#     if (is.null(ncol))
#         ncol <- ceiling(sqrt(per_page))
#     nrow <- ceiling(per_page / ncol)
#
#     blanks <- function(n)
#         replicate(n, ggplot2::ggplot() + ggplot2::theme_void(), simplify = FALSE)
#
#     page_title <- "Age Prediction residuals (predicted - actual)"
#     if (!is.null(cellType))
#         page_title <- paste(cellType, page_title, sep = "\n")
#
#     pages <- list()
#     for (s in seq(1, length(plots), by = per_page)) {
#         page_plots <- plots[s:min(s + per_page - 1, length(plots))]
#         if (length(page_plots) < per_page)
#             page_plots <- c(page_plots, blanks(per_page - length(page_plots)))
#
#         grid <- cowplot::plot_grid(
#             plotlist = page_plots,
#             ncol = ncol,
#             nrow = nrow
#         )
#
#         pg <- cowplot::ggdraw() +
#             cowplot::draw_label(
#                 page_title,
#                 x = 0,
#                 y = 1,
#                 hjust = 0,
#                 vjust = 1,
#                 size = 14
#             ) +
#             cowplot::draw_plot(grid, y = 0, height = 0.94)
#
#         pages[[length(pages) + 1]] <- pg
#     }
#
#     pages
# }

plot_residual_pair_scatter_paged <- function(res_mat,
                                             cellType = NULL,
                                             per_page = 12,
                                             facet_font_size = 10,
                                             ncol = NULL,
                                             donor_meta = NULL,
                                             donor_id_col = "donor",
                                             color_var = "age",
                                             color_title = NULL) {
    stopifnot(is.matrix(res_mat))

    regs <- colnames(res_mat)
    if (length(regs) < 2)
        stop("Need >=2 cell types (columns) in res_mat")

    donors <- rownames(res_mat)
    if (is.null(donors))
        stop("res_mat must have rownames containing donor IDs")

    meta <- NULL
    if (!is.null(donor_meta)) {
        stopifnot(is.data.frame(donor_meta) || data.table::is.data.table(donor_meta))
        stopifnot(donor_id_col %in% names(donor_meta))
        stopifnot(color_var %in% names(donor_meta))

        cols <- c(donor_id_col, color_var)

        if (data.table::is.data.table(donor_meta)) {
            meta <- unique(donor_meta[, ..cols])
        } else {
            meta <- unique(donor_meta[, cols, drop = FALSE])
        }

        data.table::setnames(meta, cols, c("donor", "color_value"))
    }

    prs <- utils::combn(regs, 2, simplify = FALSE)

    make_panel <- function(title, x, y, xlab, ylab)
    {
        keep <- is.finite(x) & is.finite(y)
        if (sum(keep) < 2)
            return(NULL)

        x <- x[keep]
        y <- y[keep]
        d <- donors[keep]

        df <- data.frame(
            donor = d,
            x = x,
            y = y,
            stringsAsFactors = FALSE
        )

        if (!is.null(meta)) {
            df <- merge(df, meta, by = "donor", all.x = TRUE, sort = FALSE)
        } else {
            df$color_value <- NA_real_
        }

        r <- stats::cor(df$x, df$y)

        rng <- range(c(df$x, df$y))
        pad <- diff(rng) * 0.1
        lim <- c(rng[1] - pad, rng[2] + pad)

        use_color <- any(is.finite(df$color_value))
        if (is.null(color_title))
            color_title <- color_var

        p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
            ggplot2::geom_abline(
                intercept = 0,
                slope = 1,
                color = "black",
                linetype = "dashed"
            )

        if (use_color) {
            p <- p +
                ggplot2::geom_point(
                    ggplot2::aes(color = color_value),
                    size = 2,
                    alpha = 0.7
                ) +
                ggplot2::scale_color_viridis_c(name = color_title)
        } else {
            p <- p +
                ggplot2::geom_point(
                    size = 2,
                    alpha = 0.7,
                    color = "steelblue"
                )
        }

        p +
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
                size = 5
            ) +
            ggplot2::coord_cartesian(xlim = lim, ylim = lim) +
            ggplot2::labs(x = xlab, y = ylab) +
            ggplot2::theme_classic(base_size = 12) +
            ggplot2::theme(
                axis.title = ggplot2::element_text(size = 14),
                axis.text  = ggplot2::element_text(size = 12)
            )
    }

    plots <- list()
    for (p in prs) {
        x_name <- p[1]
        y_name <- p[2]

        pan <- make_panel(
            title = NULL,
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
                                       annotate_cells = TRUE,
                                       row_fontsize=9,
                                       col_fontsize=9,
                                       cell_fontsize=9) {
    stopifnot(is.matrix(res_mat))
    C <- cor(res_mat, use = "pairwise.complete.obs")  # regions x regions

    col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("#3b4cc0", "white", "#b40426"))

    column_title <- "Age Prediction residuals (predicted - actual)"
    if (!is.null(cellType))
        column_title <- paste(cellType, column_title, sep = "\n")

    # optional cell annotation
    cf <- if (isTRUE(annotate_cells)) {
        function(j, i, x, y, width, height, fill) {
            grid::grid.text(sprintf("%.2f", C[i, j]), x, y, gp = grid::gpar(fontsize = cell_fontsize))
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
        row_names_gp = grid::gpar(fontsize = row_fontsize),
        column_names_gp = grid::gpar(fontsize = col_fontsize),
        heatmap_legend_param = list(at = c(-1, -0.5, 0, 0.5, 1), title = "Correlation"),
        cell_fun = cf
    )
}

# Build Jaccard matrix of selected features per level of `group_var`
jaccard_by_group <- function(all_models, group_var, coef_thresh = 0) {
    stopifnot(
        all(c("feature", "coef") %in% names(all_models)),
        group_var %in% names(all_models)
    )

    cols <- c("feature", group_var)
    sel  <- abs(all_models[["coef"]]) > coef_thresh

    if (data.table::is.data.table(all_models)) {
        df <- all_models[sel, ..cols]
    } else {
        df <- all_models[sel, cols, drop = FALSE]
    }

    df <- unique(df)  # one row per (feature, group)

    feats <- sort(unique(df[["feature"]]))
    grps  <- sort(unique(df[[group_var]]))

    A <- matrix(
        FALSE,
        nrow = length(feats),
        ncol = length(grps),
        dimnames = list(feats, grps)
    )

    A[cbind(match(df[["feature"]], feats), match(df[[group_var]], grps))] <- TRUE

    XtX <- crossprod(A)
    n   <- matrix(colSums(A), ncol = ncol(A), nrow = ncol(A))
    J   <- as.matrix(XtX / (n + t(n) - XtX))
    diag(J) <- 1
    J
}

# Heatmap with optional numeric annotation
plot_jaccard_heatmap <- function(J,
                                 title = "Feature overlap (Jaccard Index)",
                                 annotate_cells = TRUE,
                                 cluster = TRUE,
                                 row_fontsize=10,
                                 col_fontsize=10,
                                 cell_fontsize=10) {
    stopifnot(
        requireNamespace("ComplexHeatmap", quietly = TRUE),
        requireNamespace("circlize", quietly = TRUE),
        requireNamespace("grid", quietly = TRUE)
    )

    col_fun <- circlize::colorRamp2(c(0, 0.5, 1), c("#f7fbff", "#6baed6", "#08306b"))

    cf <- if (isTRUE(annotate_cells)) {
        function(j, i, x, y, w, h, fill) {
            grid::grid.text(sprintf("%.2f", J[i, j]), x, y, gp = grid::gpar(fontsize = cell_fontsize))
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
        row_names_gp = grid::gpar(fontsize = row_fontsize),
        column_names_gp = grid::gpar(fontsize = col_fontsize),
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

