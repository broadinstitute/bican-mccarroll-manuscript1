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
# out_pdf_file="/downloads/age_prediction_region_residual_comparison.pdf"
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
                                                  value_var = "resid_mean",
                                                  exclude_donor_age_range = c(),
                                                  out_pdf_file = NULL) {

    model_predictions <- load_model_predictions(prediction_file, cellTypeListFile)

    if (length(exclude_donor_age_range) == 2) {
        age_min <- exclude_donor_age_range[1]
        age_max <- exclude_donor_age_range[2]
        logger::log_info(sprintf(
            "Excluding donors with age in range [%s, %s] from analysis",
            age_min, age_max
        ))
        model_predictions <- model_predictions[!(
            model_predictions$age >= age_min & model_predictions$age <= age_max
        ), ]
    }

    all_models <- load_models(model_file, cellTypeListFile)

    if (!is.null(out_pdf_file))
        grDevices::pdf(out_pdf_file, width = 11, height = 11)
    on.exit({
        if (!is.null(out_pdf_file)) grDevices::dev.off()
    }, add = TRUE)

    for (current_region in unique(model_predictions$region)) {
        logger::log_info("Processing region: ", current_region)
        out <- plot_residuals_for_subset(
            model_predictions = model_predictions,
            all_models = all_models,
            mode = "within_region",
            region = current_region,
            value_var = value_var
        )
        if (!is.null(out)) {
            ComplexHeatmap::draw(out$jaccard_heatmap, heatmap_legend_side = "right")
            ComplexHeatmap::draw(out$corr_heatmap, heatmap_legend_side = "right")
            ComplexHeatmap::draw(out$coef_corr_heatmap, heatmap_legend_side = "right")
            for (p in out$scatter_pages) print(p)
        }
    }

    for (current_cell_type in unique(model_predictions$cell_type)) {
        logger::log_info("Processing cell type: ", current_cell_type)
        out <- plot_residuals_for_subset(
            model_predictions = model_predictions,
            all_models = all_models,
            mode = "within_cell_type",
            cell_type = current_cell_type,
            value_var = value_var
        )
        if (!is.null(out)) {
            ComplexHeatmap::draw(out$jaccard_heatmap, heatmap_legend_side = "right")
            ComplexHeatmap::draw(out$corr_heatmap, heatmap_legend_side = "right")
            ComplexHeatmap::draw(out$coef_corr_heatmap, heatmap_legend_side = "right")
            for (p in out$scatter_pages) print(p)
        }
    }

    invisible(NULL)
}

#############################
# PLOTS
#############################


#' Plot residual scatter for one pair of groups within a slice
#'
#' Create a single residual-vs-residual scatter plot comparing two groups
#' (two cell types within a region, or two regions within a cell type).
#' Points correspond to donors present in both groups (intersection of donors).
#' The plot includes a y = x reference line, a linear fit, and an annotation of
#' the Pearson correlation computed on finite donor pairs.
#'
#' @param model_predictions A data.frame or data.table of donor-level predictions.
#'   Must contain columns `donor`, `cell_type`, `region`, `age`, and the column
#'   specified by `value_var`.
#' @param mode Character scalar. Either `"within_region"` (compare cell types
#'   within a fixed region) or `"within_cell_type"` (compare regions within a
#'   fixed cell type).
#' @param cell_type Character scalar giving the cell type slice when
#'   `mode = "within_cell_type"`. Ignored when `mode = "within_region"`.
#' @param region Character scalar giving the region slice when
#'   `mode = "within_region"`. Ignored when `mode = "within_cell_type"`.
#' @param x_group Character scalar giving the first group to plot on the x axis.
#'   This must be a value of the compared dimension (cell type for
#'   `"within_region"`, region for `"within_cell_type"`).
#' @param y_group Character scalar giving the second group to plot on the y axis.
#' @param value_var Character scalar giving the residual column to plot
#'   (e.g., `"resid_mean"` or `"resid_mean_corrected"`).
#' @param donor_id_col Column name for donor IDs in `model_predictions`.
#' @param color_var Column name in `model_predictions` used to color points
#'   (default `"age"`).
#' @param color_title Legend title for the color scale. If `NULL`, defaults to
#'   `color_var`.
#' @param point_size Numeric point size.
#' @param point_alpha Numeric point transparency.
#' @param lm_linewidth Line width for the linear regression fit.
#' @param annotate_r_size Text size for the correlation annotation.
#' @param pad_frac Fraction of the data range used to pad shared x/y limits.
#'
#' @return A `ggplot` object, or `NULL` if fewer than two donors are present in
#'   the intersection after filtering.
#'
#' @export
plot_residual_pair_scatter_one <- function(model_predictions,
                                           mode = c("within_region", "within_cell_type"),
                                           cell_type = NULL,
                                           region = NULL,
                                           x_group,
                                           y_group,
                                           value_var = "resid_mean",
                                           donor_id_col = "donor",
                                           color_var = "age",
                                           color_title = "Donor age",
                                           point_size = 2,
                                           point_alpha = 0.7,
                                           lm_linewidth = 0.6,
                                           annotate_r_size = 5,
                                           pad_frac = 0.1) {

    mode <- match.arg(mode)
    .require_slice_args(mode, cell_type, region)
    dims <- .mode_dims(mode)

    mp <- .slice_predictions(model_predictions, mode, cell_type, region)
    if (.warn_if_fewer_than_2_groups(mp, dims$group_var, .default_title_label(mode, cell_type, region))) {
        return(NULL)
    }

    # Pull the two groups and intersect donors
    data.table::setDT(mp)
    dt_x <- mp[get(dims$group_var) == x_group, c(donor_id_col, value_var, color_var), with = FALSE]
    dt_y <- mp[get(dims$group_var) == y_group, c(donor_id_col, value_var), with = FALSE]

    data.table::setnames(dt_x, c(donor_id_col, value_var, color_var), c("donor", "x", "color_value"))
    data.table::setnames(dt_y, c(donor_id_col, value_var), c("donor", "y"))

    df <- merge(dt_x, dt_y, by = "donor", all = FALSE, sort = FALSE)
    if (nrow(df) < 2) {
        logger::log_warn("Skipping pair plot: <2 donors in intersection for requested groups")
        return(NULL)
    }

    # Build a 2-col matrix and reuse the matrix plotter for consistent styling
    res_mat <- as.matrix(df[, c("x", "y")])
    rownames(res_mat) <- df[["donor"]]
    colnames(res_mat) <- c(x_group, y_group)

    donor_meta <- df[, c("donor", "color_value")]
    data.table::setnames(donor_meta, c("donor", "color_value"), c(donor_id_col, color_var))

    .plot_residual_pair_scatter_one_matrix(
        res_mat = res_mat,
        x_name = x_group,
        y_name = y_group,
        donor_meta = donor_meta,
        donor_id_col = donor_id_col,
        color_var = color_var,
        color_title = color_title,
        point_size = point_size,
        point_alpha = point_alpha,
        lm_linewidth = lm_linewidth,
        annotate_r_size = annotate_r_size,
        pad_frac = pad_frac
    )
}

#' Plot residual correlation heatmap within a slice
#'
#' Compute the pairwise Pearson correlation matrix of donor residuals across
#' groups (cell types within a region, or regions within a cell type) and render
#' the result as a clustered heatmap using `ComplexHeatmap`.
#'
#' Correlations are computed with `stats::cor(..., use = "pairwise.complete.obs")`.
#' The heatmap is clustered on both axes using hierarchical clustering.
#'
#' @param model_predictions A data.frame or data.table of donor-level predictions.
#'   Must contain columns `donor`, `cell_type`, `region`, and the column specified
#'   by `value_var`.
#' @param mode Character scalar. Either `"within_region"` or `"within_cell_type"`.
#' @param cell_type Character scalar giving the cell type slice when
#'   `mode = "within_cell_type"`.
#' @param region Character scalar giving the region slice when
#'   `mode = "within_region"`.
#' @param value_var Character scalar giving the residual column to correlate
#'   (e.g., `"resid_mean"` or `"resid_mean_corrected"`).
#' @param title Optional character scalar to use as the heatmap title. If `NULL`,
#'   a default title derived from `mode` and the slice value is used.
#' @param annotate_cells Logical; if `TRUE`, annotate cells with correlation values.
#' @param row_fontsize Font size for row labels.
#' @param col_fontsize Font size for column labels.
#' @param cell_fontsize Font size for cell annotations.
#'
#' @return A list with components:
#' \describe{
#'   \item{heatmap}{A `ComplexHeatmap::Heatmap` object. Draw with `ComplexHeatmap::draw()`.}
#'   \item{row_order_names}{Character vector giving the clustered row order.}
#'   \item{column_order_names}{Character vector giving the clustered column order.}
#' }
#'
#' @details
#' To render the returned heatmap:
#' \preformatted{
#' out <- plot_residual_corr_heatmap(...)
#' ComplexHeatmap::draw(out$heatmap, heatmap_legend_side = "right")
#' }
#'
#' @export
plot_residual_corr_heatmap <- function(model_predictions,
                                       mode = c("within_region", "within_cell_type"),
                                       cell_type = NULL,
                                       region = NULL,
                                       value_var = "resid_mean",
                                       title = NULL,
                                       annotate_cells = TRUE,
                                       row_fontsize = 9,
                                       col_fontsize = 9,
                                       cell_fontsize = 9) {

    mode <- match.arg(mode)

    # Slice + choose compared dimension
    dims <- .mode_dims(mode)
    .require_slice_args(mode = mode, cell_type = cell_type, region = region)

    mp <- model_predictions
    if (identical(dims$slice_var, "region")) {
        mp <- mp[mp$region == region, , drop = FALSE]
        spread_var <- "cell_type"
        title_label <- if (is.null(title)) sprintf("Age Prediction residuals (predicted - actual)\nregion [%s]", region) else title
    } else {
        mp <- mp[mp$cell_type == cell_type, , drop = FALSE]
        spread_var <- "region"
        title_label <- if (is.null(title)) sprintf("Age Prediction residuals (predicted - actual)\ncell type [%s]", cell_type) else title
    }

    # Need >=2 groups
    n_groups <- length(unique(mp[[spread_var]]))
    if (n_groups < 2) {
        logger::log_warn(sprintf("Only %d unique %s present; returning NULL", n_groups, spread_var))
        return(NULL)
    }

    # donor x group residual matrix
    res_mat <- donor_wide_matrix(
        model_predictions = mp,
        spread_var = spread_var,
        value_var = value_var
    )

    # Heatmap object (clustered)
    ht <- .plot_residual_corr_heatmap_matrix(
        res_mat = res_mat,
        title = title_label,
        annotate_cells = annotate_cells,
        row_fontsize = row_fontsize,
        col_fontsize = col_fontsize,
        cell_fontsize = cell_fontsize
    )

    # Extract clustered order without drawing to the active device
    ord <- .get_heatmap_order_names(ht, colnames(res_mat))

    list(
        heatmap = ht,
        row_order_names = ord$row_order_names,
        column_order_names = ord$column_order_names
    )
}

#' Plot Jaccard overlap heatmap for aging programs within a slice
#'
#' Compute pairwise Jaccard overlap of aging gene programs across groups
#' (cell types within a region, or regions within a cell type) and generate a
#' heatmap visualization.
#'
#' Gene sets are defined per group using features with `abs(coef) > coef_thresh`.
#'
#' @param all_models A data.frame or data.table of model coefficients.
#'   Must contain columns `cell_type`, `region`, `feature`, and `coef`.
#' @param mode Character scalar. Either `"within_region"` or `"within_cell_type"`.
#' @param cell_type Character scalar giving the cell type slice when
#'   `mode = "within_cell_type"`.
#' @param region Character scalar giving the region slice when
#'   `mode = "within_region"`.
#' @param coef_thresh Numeric threshold applied to `abs(coef)` when defining gene sets.
#' @param title Optional character scalar to use as the heatmap title. If `NULL`,
#'   a default title derived from `mode` and the slice value is used.
#' @param annotate_cells Logical; if `TRUE`, annotate heatmap cells with Jaccard values.
#' @param row_fontsize Font size for heatmap row labels.
#' @param col_fontsize Font size for heatmap column labels.
#' @param cell_fontsize Font size for numeric cell annotations.
#' @param row_order_names Optional character vector of row names specifying the desired
#'   row order. Names not present in the Jaccard matrix are ignored. If provided,
#'   row clustering is disabled.
#' @param column_order_names Optional character vector of column names specifying the
#'   desired column order. Names not present in the Jaccard matrix are ignored. If
#'   provided, column clustering is disabled.
#'
#' @return A list with components:
#' \describe{
#'   \item{jaccard}{Numeric matrix of pairwise Jaccard coefficients.}
#'   \item{heatmap}{A `ComplexHeatmap::Heatmap` object. Draw with `ComplexHeatmap::draw()`.}
#' }
#'
#' @details
#' To render the returned heatmap:
#' \preformatted{
#' out <- plot_jaccard_overlap_heatmap(...)
#' ComplexHeatmap::draw(out$heatmap, heatmap_legend_side = "right")
#' }
#'
#' @export
plot_jaccard_overlap_heatmap <- function(all_models,
                                         mode = c("within_region", "within_cell_type"),
                                         cell_type = NULL,
                                         region = NULL,
                                         coef_thresh = 0,
                                         title = NULL,
                                         annotate_cells = TRUE,
                                         row_fontsize = 16,
                                         col_fontsize = 16,
                                         cell_fontsize = 16,
                                         row_order_names = NULL,
                                         column_order_names = NULL) {

    mode <- match.arg(mode)
    .require_slice_args(mode, cell_type, region)
    dims <- .mode_dims(mode)

    am <- .slice_models(all_models, mode, cell_type, region)

    # Default title logic matches plot_residual_corr_heatmap:
    if (is.null(title)) {
        if (identical(dims$slice_var, "region")) {
            title <- sprintf("Cell-type gene overlap in aging programs\nregion [%s]", region)
        } else {
            title <- sprintf("Region gene overlap in aging programs\ncell type [%s]", cell_type)
        }
    }

    J <- jaccard_by_group(
        am,
        group_var = dims$group_var,
        coef_thresh = coef_thresh
    )

    ht <- plot_jaccard_heatmap(
        J,
        title = title,
        annotate_cells = annotate_cells,
        row_fontsize = row_fontsize,
        col_fontsize = col_fontsize,
        cell_fontsize = cell_fontsize,
        row_order_names = row_order_names,
        column_order_names = column_order_names
    )

    list(jaccard = J, heatmap = ht)
}

plot_residuals_for_subset <- function(model_predictions,
                                      all_models,
                                      mode = c("within_region", "within_cell_type"),
                                      cell_type = NULL,
                                      region = NULL,
                                      value_var = "resid_mean",
                                      annotate_cells = TRUE,
                                      per_page = 4,
                                      facet_font_size = 8,
                                      coef_thresh = 0,
                                      corr_title = NULL) {

    mode <- match.arg(mode)
    .require_slice_args(mode, cell_type, region)
    dims <- .mode_dims(mode)

    title_label <- .default_title_label(mode, cell_type, region)

    # Build res_mat once for paged scatter (efficient) and for extracting ordering
    mp_sub <- .slice_predictions(model_predictions, mode, cell_type, region)
    if (.warn_if_fewer_than_2_groups(mp_sub, dims$group_var, title_label)) {
        return(invisible(NULL))
    }

    res_mat <- donor_wide_matrix(
        model_predictions = mp_sub,
        spread_var = dims$group_var,
        value_var = value_var
    )

    donor_age <- unique(mp_sub[, c("donor", "age")])

    p_scatter_pages <- plot_residual_pair_scatter_paged(
        res_mat = res_mat,
        cellType = title_label,
        per_page = per_page,
        facet_font_size = facet_font_size,
        donor_meta = donor_age,
        donor_id_col = "donor",
        color_var = "age",
        color_title = "Donor age"
    )

    ht_corr_out <- plot_residual_corr_heatmap(
        model_predictions = mp_sub,
        mode = mode,
        cell_type = cell_type,
        region = region,
        value_var = value_var,
        title = corr_title,
        annotate_cells = annotate_cells,
        row_fontsize = 16,
        col_fontsize = 16,
        cell_fontsize = 16
    )

    # Pull heatmap + clustered ordering directly from the returned list
    ht_corr <- NULL
    row_order_names <- NULL
    col_order_names <- NULL

    if (!is.null(ht_corr_out)) {
        ht_corr <- ht_corr_out$heatmap
        row_order_names <- ht_corr_out$row_order_names
        col_order_names <- ht_corr_out$column_order_names
    }

    tmp_j <- plot_jaccard_overlap_heatmap(
        all_models = all_models,
        mode = mode,
        cell_type = cell_type,
        region = region,
        coef_thresh = coef_thresh,
        title = sprintf("Cell-type gene overlap in aging programs\nregion [%s]", region),
        annotate_cells = annotate_cells,
        row_fontsize = 16,
        col_fontsize = 16,
        cell_fontsize = 16,
        row_order_names = row_order_names,
        column_order_names = col_order_names
    )


    cc <- coef_corr_on_intersect(
        .slice_models(all_models, mode, cell_type, region),
        group_var = dims$group_var
    )

    ht_coef <- plot_coef_corr_heatmap(
        coef_correlation = cc$coef_correlation,
        overlap_gene_counts = cc$overlap_gene_counts,
        title = paste0(dims$group_var, " coef correlation (intersect genes) ", title_label),
        annotate_cells = annotate_cells,
        row_fontsize = 16,
        col_fontsize = 16,
        cell_fontsize = 16
    )

    # Return objects; caller decides whether to draw/print
    invisible(list(
        res_mat = res_mat,
        scatter_pages = p_scatter_pages,
        corr_heatmap = ht_corr,
        jaccard = tmp_j$jaccard,
        jaccard_heatmap = tmp_j$heatmap,
        coef_corr_heatmap = ht_coef
    ))
}

plot_residual_pair_scatter_paged <- function(res_mat,
                                             cellType = NULL,
                                             per_page = 12,
                                             facet_font_size = 10,
                                             ncol = NULL,
                                             donor_meta = NULL,
                                             donor_id_col = "donor",
                                             color_var = "age",
                                             color_title = NULL,
                                             point_size = 2,
                                             point_alpha = 0.7,
                                             lm_linewidth = 0.6,
                                             annotate_r_size = 5,
                                             pad_frac = 0.1) {
    stopifnot(is.matrix(res_mat))

    regs <- colnames(res_mat)
    if (length(regs) < 2)
        stop("Need >=2 groups (columns) in res_mat")

    prs <- utils::combn(regs, 2, simplify = FALSE)

    plots <- list()
    for (pair in prs) {
        pan <- .plot_residual_pair_scatter_one_matrix(
            res_mat = res_mat,
            x_name = pair[1],
            y_name = pair[2],
            donor_meta = donor_meta,
            donor_id_col = donor_id_col,
            color_var = color_var,
            color_title = color_title,
            point_size = point_size,
            point_alpha = point_alpha,
            lm_linewidth = lm_linewidth,
            annotate_r_size = annotate_r_size,
            pad_frac = pad_frac
        )
        if (!is.null(pan))
            plots[[length(plots) + 1]] <- pan
    }

    if (!length(plots))
        stop("No valid pairs after filtering")

    if (is.null(ncol))
        ncol <- ceiling(sqrt(per_page))
    nrow <- ceiling(per_page / ncol)

    blanks <- function(n) {
        replicate(n, ggplot2::ggplot() + ggplot2::theme_void(), simplify = FALSE)
    }

    page_title <- "Age Prediction residuals (predicted - actual)"
    if (!is.null(cellType))
        page_title <- paste(cellType, page_title, sep = "\n")

    pages <- list()
    for (s in seq(1, length(plots), by = per_page)) {
        page_plots <- plots[s:min(s + per_page - 1, length(plots))]
        if (length(page_plots) < per_page)
            page_plots <- c(page_plots, blanks(per_page - length(page_plots)))

        grid <- cowplot::plot_grid(plotlist = page_plots, ncol = ncol, nrow = nrow)

        pages[[length(pages) + 1]] <- cowplot::ggdraw() +
            cowplot::draw_label(page_title, x = 0, y = 1, hjust = 0, vjust = 1, size = 14) +
            cowplot::draw_plot(grid, y = 0, height = 0.94)
    }

    pages
}

########################
# HELPERS
########################

# Generalized: donors x <spread_var> matrix for a chosen value column
#' Convert donor-level predictions to a wide residual matrix
#'
#' Reshape a long-format donor-level prediction data.table or data.frame into
#' a wide matrix suitable for pairwise residual comparisons. Rows correspond
#' to donors and columns correspond to levels of `spread_var`. Cell values
#' are taken from `value_var`.
#'
#' This function is typically used to prepare input for residual scatter
#' plots or correlation analyses across regions or cell types.
#'
#' @param model_predictions A data.frame or data.table containing at minimum
#'   the columns `donor`, `spread_var`, and `value_var`.
#' @param spread_var A single string giving the column name whose unique values
#'   will become columns in the output matrix (e.g., region or cell type).
#' @param value_var A single string giving the column name containing numeric
#'   values to populate the matrix (e.g., residuals).
#'
#' @return A numeric matrix with donors as rownames and levels of
#'   `spread_var` as column names.
#'
#' @importFrom data.table setDT dcast
#' @importFrom stats as.formula
#' @noRd
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

# Build Jaccard matrix of selected features per level of `group_var`
jaccard_by_group <- function(all_models, group_var, coef_thresh = 0) {
    stopifnot(
        all(c("feature", "coef") %in% names(all_models)),
        group_var %in% names(all_models)
    )

    cols <- c("feature", group_var)
    sel  <- abs(all_models[["coef"]]) > coef_thresh

    if (data.table::is.data.table(all_models)) {
        df <- all_models[sel, cols, with = FALSE]
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

plot_jaccard_heatmap <- function(J,
                                 title = NULL,
                                 annotate_cells = TRUE,
                                 row_fontsize = 16,
                                 col_fontsize = 16,
                                 cell_fontsize = 16,
                                 row_order_names = NULL,
                                 column_order_names = NULL) {
    stopifnot(is.matrix(J))

    # Validate / map optional name-based ordering to indices
    row_order <- NULL
    column_order <- NULL

    if (!is.null(row_order_names)) {
        stopifnot(is.character(row_order_names))
        if (is.null(rownames(J)))
            stop("J must have rownames when row_order_names is provided")

        miss <- setdiff(row_order_names, rownames(J))
        if (length(miss) > 0L) {
            warning(sprintf(
                "row_order_names contains %d names not present in J (showing up to 10): %s",
                length(miss),
                paste(utils::head(miss, 10), collapse = ", ")
            ), call. = FALSE)
        }

        row_order_names <- row_order_names[row_order_names %in% rownames(J)]
        row_order <- match(row_order_names, rownames(J))

        if (length(row_order) == 0L)
            row_order <- NULL
    }

    if (!is.null(column_order_names)) {
        stopifnot(is.character(column_order_names))
        if (is.null(colnames(J)))
            stop("J must have colnames when column_order_names is provided")

        miss <- setdiff(column_order_names, colnames(J))
        if (length(miss) > 0L) {
            warning(sprintf(
                "column_order_names contains %d names not present in J (showing up to 10): %s",
                length(miss),
                paste(utils::head(miss, 10), collapse = ", ")
            ), call. = FALSE)
        }

        column_order_names <- column_order_names[column_order_names %in% colnames(J)]
        column_order <- match(column_order_names, colnames(J))

        if (length(column_order) == 0L)
            column_order <- NULL
    }

    # If an explicit order is provided, do not recluster that axis
    cluster_rows <- is.null(row_order)
    cluster_columns <- is.null(column_order)

    col_fun <- circlize::colorRamp2(c(0, 0.5, 1), c("#f7fbff", "#6baed6", "#08306b"))

    cf <- if (isTRUE(annotate_cells)) {
        function(j, i, x, y, width, height, fill) {
            grid::grid.text(sprintf("%.2f", J[i, j]), x, y, gp = grid::gpar(fontsize = cell_fontsize))
        }
    } else {
        NULL
    }

    ComplexHeatmap::Heatmap(
        J,
        name = "J",
        col = col_fun,
        cluster_rows = cluster_rows,
        cluster_columns = cluster_columns,
        row_order = row_order,
        column_order = column_order,
        row_title = NULL,
        column_title = title,
        show_row_dend = cluster_rows,
        show_column_dend = cluster_columns,
        row_names_gp = grid::gpar(fontsize = row_fontsize),
        column_names_gp = grid::gpar(fontsize = col_fontsize),
        heatmap_legend_param = list(at = c(0, 0.5, 1), title = "Jaccard"),
        cell_fun = cf
    )
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
        df <- all_models[, cols, with = FALSE]
    } else {
        df <- all_models[, cols, drop = FALSE]
    }

    f <- stats::as.formula(paste0("feature ~ ", group_var))
    wide_dt <- data.table::dcast(data.table::as.data.table(df), f, value.var = "coef")

    mat <- as.matrix(wide_dt[, -1, with = FALSE])
    rownames(mat) <- wide_dt[["feature"]]
    mat
}



####################
# FILE I/O
####################


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

#' Plot heatmap of residual correlations
#'
#' Compute the pairwise Pearson correlation matrix from a residual matrix and
#' display it as a clustered heatmap using `ComplexHeatmap`. Rows and columns
#' correspond to the variables (e.g., regions or cell types) in `res_mat`.
#'
#' Optionally annotates each heatmap cell with the numeric correlation value.
#'
#' @param res_mat A numeric matrix of residuals with donors in rows and
#'   variables (e.g., regions or cell types) in columns.
#' @param cellType Optional string used to prefix the heatmap title.
#' @param annotate_cells Logical; if `TRUE`, overlay each heatmap cell with
#'   the corresponding correlation value rounded to two decimals.
#' @param row_fontsize Numeric font size for row labels.
#' @param col_fontsize Numeric font size for column labels.
#' @param cell_fontsize Numeric font size for correlation annotations inside cells.
#'
#' @return A `ComplexHeatmap::Heatmap` object.
#'
#' @details
#' Correlations are computed using `stats::cor()` with
#' `use = "pairwise.complete.obs"`. Rows and columns are hierarchically
#' clustered. The color scale is fixed to the range [-1, 1] with blue for
#' negative correlations, white for zero, and red for positive correlations.
#'
#' @importFrom stats cor
#' @importFrom circlize colorRamp2
#' @importFrom grid grid.text gpar
#' @importFrom ComplexHeatmap Heatmap
#' @noRd
.plot_residual_corr_heatmap_matrix <- function(res_mat,
                                               title = NULL,
                                               annotate_cells = TRUE,
                                               row_fontsize = 9,
                                               col_fontsize = 9,
                                               cell_fontsize = 9) {
    stopifnot(is.matrix(res_mat))

    C <- stats::cor(res_mat, use = "pairwise.complete.obs")

    col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("#3b4cc0", "white", "#b40426"))

    if (is.null(title)) {
        title <- "Age Prediction residuals (predicted - actual)"
    }

    cf <- if (isTRUE(annotate_cells)) {
        function(j, i, x, y, width, height, fill) {
            grid::grid.text(
                sprintf("%.2f", C[i, j]),
                x,
                y,
                gp = grid::gpar(fontsize = cell_fontsize)
            )
        }
    } else {
        NULL
    }

    ComplexHeatmap::Heatmap(
        C,
        name = "r",
        col = col_fun,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_title = NULL,
        column_title = title,
        show_row_dend = TRUE,
        show_column_dend = TRUE,
        row_names_gp = grid::gpar(fontsize = row_fontsize),
        column_names_gp = grid::gpar(fontsize = col_fontsize),
        heatmap_legend_param = list(at = c(-1, -0.5, 0, 0.5, 1), title = "Correlation"),
        cell_fun = cf
    )
}

#' Plot residual scatter for a single pair of columns
#'
#' Create a single residual-vs-residual scatter plot for two columns of a residual
#' matrix, optionally coloring points by a donor-level metadata variable. The plot
#' includes a y=x reference line, a linear regression fit, and an annotation of the
#' Pearson correlation computed on finite x/y pairs. The x and y limits are shared
#' and padded to enforce a square comparison scale.
#'
#' @param res_mat A numeric matrix of residuals with donors in rows and variables
#'   (e.g., regions or cell types) in columns. Must have `rownames(res_mat)` as
#'   donor IDs and `colnames(res_mat)` containing `x_name` and `y_name`.
#' @param x_name A single string giving the column name in `res_mat` to plot on
#'   the x axis.
#' @param y_name A single string giving the column name in `res_mat` to plot on
#'   the y axis.
#' @param donor_meta Optional data.frame or data.table containing donor-level
#'   metadata for coloring points. If provided, must include `donor_id_col` and
#'   `color_var`.
#' @param donor_id_col Column name in `donor_meta` containing donor IDs.
#'   Matched to `rownames(res_mat)`.
#' @param color_var Column name in `donor_meta` used to color points. If no finite
#'   values are present after merging, points are drawn in a single color.
#' @param color_title Optional legend title for the color scale. If `NULL`,
#'   defaults to `color_var`.
#' @param point_size Point size used for the scatter.
#' @param point_alpha Point transparency used for the scatter.
#' @param lm_linewidth Line width for the linear regression fit.
#' @param annotate_r_size Text size for the correlation annotation.
#' @param pad_frac Fraction of the data range used to pad the shared x/y limits.
#'
#' @return A `ggplot` object, or `NULL` if fewer than two finite (x, y) pairs are
#'   available after filtering.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_abline geom_smooth annotate
#' @importFrom ggplot2 coord_cartesian labs theme_classic theme element_text
#' @importFrom stats cor
#' @noRd
.plot_residual_pair_scatter_one_matrix <- function(res_mat,
                                                   x_name,
                                                   y_name,
                                                   donor_meta = NULL,
                                                   donor_id_col = "donor",
                                                   color_var = "age",
                                                   color_title = NULL,
                                                   point_size = 2,
                                                   point_alpha = 0.7,
                                                   lm_linewidth = 0.6,
                                                   annotate_r_size = 5,
                                                   pad_frac = 0.1) {

    stopifnot(is.matrix(res_mat))
    stopifnot(is.character(x_name) && length(x_name) == 1L)
    stopifnot(is.character(y_name) && length(y_name) == 1L)

    regs <- colnames(res_mat)
    if (is.null(regs))
        stop("res_mat must have colnames")

    if (!(x_name %in% regs))
        stop(sprintf("x_name not found in res_mat colnames: %s", x_name))
    if (!(y_name %in% regs))
        stop(sprintf("y_name not found in res_mat colnames: %s", y_name))
    if (x_name == y_name)
        stop("x_name and y_name must be different")

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
            meta <- unique(donor_meta[, cols, with = FALSE])
        } else {
            meta <- unique(donor_meta[, cols, drop = FALSE])
        }

        data.table::setnames(meta, cols, c("donor", "color_value"))
    }

    x <- res_mat[, x_name]
    y <- res_mat[, y_name]

    keep <- is.finite(x) & is.finite(y)
    if (sum(keep) < 2)
        return(NULL)

    df <- data.frame(
        donor = donors[keep],
        x = x[keep],
        y = y[keep],
        stringsAsFactors = FALSE
    )

    if (!is.null(meta)) {
        df <- merge(df, meta, by = "donor", all.x = TRUE, sort = FALSE)
    } else {
        df$color_value <- NA_real_
    }

    use_color <- any(is.finite(df$color_value))
    if (is.null(color_title))
        color_title <- color_var

    r <- stats::cor(df$x, df$y)

    rng <- range(c(df$x, df$y))
    pad <- diff(rng) * pad_frac
    lim <- c(rng[1] - pad, rng[2] + pad)

    # R CMD CHECK
    intercept <- slope <- color_value <- NULL

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
                size = point_size,
                alpha = point_alpha
            ) +
            ggplot2::scale_color_viridis_c(name = color_title)
    } else {
        p <- p +
            ggplot2::geom_point(
                size = point_size,
                alpha = point_alpha,
                color = "steelblue"
            )
    }

    p +
        ggplot2::geom_smooth(
            method = "lm",
            formula = y ~ x,
            se = FALSE,
            linewidth = lm_linewidth,
            color = "red"
        ) +
        ggplot2::annotate(
            "text",
            x = -Inf,
            y = Inf,
            label = sprintf("r = %.2f", r),
            hjust = -0.1,
            vjust = 1.2,
            size = annotate_r_size
        ) +
        ggplot2::coord_cartesian(xlim = lim, ylim = lim) +
        ggplot2::labs(
            x = paste("Residual in", x_name),
            y = paste("Residual in", y_name)
        ) +
        ggplot2::theme_classic(base_size = 12) +
        ggplot2::theme(
            axis.title = ggplot2::element_text(size = 14),
            axis.text  = ggplot2::element_text(size = 12)
        )
}

#################
# OTHER HELPERS
#################

.mode_dims <- function(mode) {
    if (identical(mode, "within_region")) {
        return(list(
            slice_var = "region",
            group_var = "cell_type",
            slice_label = "region"
        ))
    }
    if (identical(mode, "within_cell_type")) {
        return(list(
            slice_var = "cell_type",
            group_var = "region",
            slice_label = "cell_type"
        ))
    }
    stop("mode must be 'within_region' or 'within_cell_type'", call. = FALSE)
}

.require_slice_args <- function(mode, cell_type, region) {
    if (identical(mode, "within_region")) {
        if (is.null(region) || length(region) != 1L) {
            stop("For mode='within_region', region must be a single value", call. = FALSE)
        }
    } else {
        if (is.null(cell_type) || length(cell_type) != 1L) {
            stop("For mode='within_cell_type', cell_type must be a single value", call. = FALSE)
        }
    }
}

.slice_predictions <- function(model_predictions, mode, cell_type, region) {
    dims <- .mode_dims(mode)
    data.table::setDT(model_predictions)

    if (identical(dims$slice_var, "region")) {
        region_value <- region
        model_predictions[region == region_value]
    } else {
        cell_type_value <- cell_type
        model_predictions[cell_type == cell_type_value]
    }
}

.slice_models <- function(all_models, mode, cell_type, region) {
    dims <- .mode_dims(mode)
    data.table::setDT(all_models)

    if (identical(dims$slice_var, "region")) {
        region_value <- region
        all_models[region == region_value]
    } else {
        cell_type_value <- cell_type
        all_models[cell_type == cell_type_value]
    }
}

.warn_if_fewer_than_2_groups <- function(dt, group_var, label) {
    n_groups <- length(unique(dt[[group_var]]))
    if (is.na(n_groups) || n_groups < 2) {
        logger::log_warn(sprintf(
            "Skipping %s: only %d unique %s present",
            label, n_groups, group_var
        ))
        return(TRUE)
    }
    FALSE
}

.default_title_label <- function(mode, cell_type, region) {
    if (identical(mode, "within_region")) {
        paste0("region [", region, "]")
    } else {
        paste0("cell_type [", cell_type, "]")
    }
}

.get_heatmap_order_names <- function(ht, name_vec) {
    tmp_pdf <- tempfile(fileext = ".pdf")
    grDevices::pdf(tmp_pdf, width = 4, height = 4)
    ht_drawn <- ComplexHeatmap::draw(ht, newpage = TRUE)
    grDevices::dev.off()
    unlink(tmp_pdf)

    ro <- ComplexHeatmap::row_order(ht_drawn)
    co <- ComplexHeatmap::column_order(ht_drawn)

    if (is.list(ro)) ro <- ro[[1]]
    if (is.list(co)) co <- co[[1]]

    list(
        row_order_names = name_vec[ro],
        column_order_names = name_vec[co]
    )
}


