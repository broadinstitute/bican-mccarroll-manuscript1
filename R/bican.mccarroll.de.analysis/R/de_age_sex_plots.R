#########################
# Basic DE plots  - volcano plots, scatter plots comparing effect sizes,
# and correlation heatmaps across cell type x region groups.
######################

#' Read cell types from a text file
#'
#' Reads one cell type label per line (or whitespace-separated tokens).
#'
#' @param ct_file Path to a text file of cell types.
#' @return A character vector of cell type names.
#' @export
read_cell_types <- function(ct_file) {
    scan(ct_file, what = character(), quiet = TRUE)
}

#' Read gene-to-chromosome mapping
#'
#' The input file must contain at least columns "gene" and "chr".
#' Chromosomes are filtered to 1:22, X, Y, and M by default.
#'
#' @param gene_to_chr_file Path to a tabular file readable by data.table::fread.
#' @return A data.table with columns including gene and chr.
#' @export
#' @import data.table
read_gene_to_chr <- function(gene_to_chr_file) {
  chr <- gene <- NULL

  dt <- data.table::fread(gene_to_chr_file)

  keep <- c(as.character(1:22), "X", "Y", "M")
  dt <- dt[dt[["chr"]] %in% keep,]

  dt[]
}

#' Read and format differential expression results
#'
#' This reads DE results via bican.mccarroll.differentialexpression::parse_de_inputs
#' and standardizes column names, merges chromosome annotations, and computes
#' log_fc standard errors from the t-statistic.
#'
#' @param de_dir Directory containing DE results.
#' @param test Test name passed to parse_de_inputs (e.g., "age", "female_vs_male").
#' @param ct_file Path to cell types file.
#' @param gene_to_chr A data.table with at least "gene" and "chr".
#' @return A data.table with standardized DE results.
#' @export
read_de_results <- function(de_dir, test, ct_file, gene_to_chr) {
    df <- bican.mccarroll.differentialexpression::parse_de_inputs(de_dir, test, ct_file)
    prep_de(df, gene_to_chr)
}

#' Prepare DE results for plotting and downstream analysis
#' (originally format_de_results)
#' @param df A data.frame-like object returned by parse_de_inputs.
#' @param gene_to_chr A data.table with at least "gene" and "chr".
#' @return A data.table with standardized columns and annotations.
#' @export
prep_de <- function(df, gene_to_chr) {
    gene <- chr <- log_fc <- t <- log_fc_se<- NULL

    dt <- data.table::as.data.table(df)

    data.table::setnames(
        dt,
        c("gene", "cell_type", "region", "test", "log_fc", "ave_expr",
          "t", "p_value", "adj_p_val", "b", "z_std")
    )

    dt <- merge(gene_to_chr, dt, by = "gene", all.y = TRUE)
    dt[, log_fc_se := log_fc / t]

    return (dt)
}

#' Volcano plot for DE results
#'
#' Produces a volcano plot for a single cell type (and optionally region).
#'
#' Overplotting order:
#' 1) Non-significant points (adj_p_val >= fdr_cutoff) are plotted first in light grey.
#' 2) Significant points (adj_p_val < fdr_cutoff) are then overplotted:
#'    - If chr_color_map is NULL: all significant points are plotted together using significant_color.
#'    - If chr_color_map is provided: significant points are plotted in groups in the REVERSE order
#'      of names(chr_color_map), so the LAST name in chr_color_map is drawn last and appears on top.
#'      The special key "default" is treated as "all significant points whose chr does not match any
#'      non-default key in chr_color_map", and it participates in that same reversed drawing order.
#'
#' @param de_dt A prepared DE data.table from prep_de/read_de_results.
#' @param cell_type_use Cell type label to plot.
#' @param region_use Region label to plot, or NA for region-combined results.
#' @param fdr_cutoff Adjusted p-value threshold.
#' @param abs_log_fc_cutoff Absolute log2 fold-change threshold.
#' @param show_title Whether to show the cell type as the plot title.
#' @param significant_color Color used for points passing the FDR threshold when chr_color_map is NULL.
#' @param chr_color_map Optional named character vector or list mapping chromosome -> color.
#'   Must include "default" for unmatched chromosomes (e.g., c(X="orange", Y="green", default="blue")).
#'   If NULL, uses the original scheme (lightgrey background + significant_color for FDR-pass).
#' @return Invisibly returns NULL.
#' @export
plot_de_volcano <- function(de_dt,
                            cell_type_use,
                            region_use,
                            fdr_cutoff = 0.05,
                            abs_log_fc_cutoff = log2(1.05),
                            show_title = TRUE,
                            significant_color = "cornflowerblue",
                            chr_color_map = NULL) {

  cell_type <- region <- log_fc <- adj_p_val <- chr <- NULL

  if (is.na(region_use)) {
    dt <- de_dt[cell_type == cell_type_use & is.na(region)]
  } else {
    dt <- de_dt[cell_type == cell_type_use & region == region_use]
  }

  rng <- dt[, max(abs(log_fc), na.rm = TRUE)]
  rng <- c(-rng, rng)

  # Precompute y for ALL points so ylim includes significant points.
  # Guard against adj_p_val == 0 (would yield Inf).
  adj_p_floor <- dt[, min(adj_p_val[adj_p_val > 0], na.rm = TRUE)]
  if (!is.finite(adj_p_floor)) {
    adj_p_floor <- .Machine$double.xmin
  }
  y_all <- -log10(pmax(dt$adj_p_val, adj_p_floor))
  ylim <- c(0, max(y_all[is.finite(y_all)], na.rm = TRUE) * 1.02)

  graphics::par(pty = "s", xpd = FALSE)

  non_sig <- dt[adj_p_val >= fdr_cutoff]
  sig_dt <- dt[adj_p_val < fdr_cutoff]

  graphics::plot(
    non_sig$log_fc,
    -log10(pmax(non_sig$adj_p_val, adj_p_floor)),
    pch = 20,
    xlim = rng,
    ylim = ylim,
    col = "lightgrey",
    xlab = "Effect size, log2",
    ylab = "Adjusted p-value, -log10",
    main = ""
  )

  if (nrow(sig_dt) > 0) {

    if (is.null(chr_color_map)) {

      graphics::points(
        sig_dt$log_fc,
        -log10(pmax(sig_dt$adj_p_val, adj_p_floor)),
        pch = 20,
        col = significant_color
      )

    } else {

      if (is.list(chr_color_map)) {
        chr_color_map <- unlist(chr_color_map, use.names = TRUE)
      }

      if (!is.character(chr_color_map) || is.null(names(chr_color_map))) {
        stop("chr_color_map must be a named character vector (or list coercible to one).")
      }

      if (!("default" %in% names(chr_color_map))) {
        stop('chr_color_map must include a "default" entry.')
      }

      map_names <- names(chr_color_map)
      non_default_names <- map_names[map_names != "default"]
      plot_order <- rev(map_names)

      sig_y <- -log10(pmax(sig_dt$adj_p_val, adj_p_floor))

      for (chr_name in plot_order) {

        if (chr_name == "default") {
          idx <- !(sig_dt$chr %in% non_default_names)
        } else {
          idx <- sig_dt$chr == chr_name
        }

        if (any(idx)) {
          graphics::points(
            sig_dt$log_fc[idx],
            sig_y[idx],
            pch = 20,
            col = chr_color_map[[chr_name]]
          )
        }
      }
    }
  }

  if (show_title)
    graphics::title(main = cell_type_use, adj = 0)

  graphics::abline(h = -log10(fdr_cutoff), lty = 2)
  graphics::abline(v = c(-abs_log_fc_cutoff, abs_log_fc_cutoff), lty = 2)

  # Make R CMD CHECK Happy
  .N <- NULL

  up <- dt[log_fc > abs_log_fc_cutoff & adj_p_val < fdr_cutoff, .N]
  down <- dt[log_fc < -abs_log_fc_cutoff & adj_p_val < fdr_cutoff, .N]

  p <- stats::binom.test(up, up + down, 0.5)$p.value

  graphics::par(xpd = NA)
  if (p < 0.05) {
    p_txt <- formatC(p, format = "e", digits = 1)
    graphics::legend(
      "topright",
      inset = c(0, -0.16),
      c(
        paste("proportion up =", round(up / (up + down), 2)),
        paste("p-value =", p_txt)
      ),
      bty = "n"
    )
  }

  invisible(NULL)
}

#' Volcano plot for DE results (ggplot2)
#'
#' Produces a volcano plot for a single cell type (and optionally region) and returns a ggplot object.
#'
#' Coloring:
#' If chr_color_map is provided, points are colored by chromosome and chr_color_map must include
#' a color for every chromosome observed in the plotted data (full coverage is enforced).
#' If chr_color_map is NULL, all points are black.
#'
#' Overplotting order:
#' 1) Non-significant points (adj_p_val >= fdr_cutoff) are drawn first.
#' 2) Significant points (adj_p_val < fdr_cutoff) are drawn second.
#' Within each significance layer, chromosomes follow the reverse order of names(chr_color_map),
#' so the last entry in chr_color_map appears on top.
#'
#' @param de_dt A prepared DE data.table from prep_de/read_de_results.
#' @param cell_type_use Cell type label to plot.
#' @param region_use Region label to plot, or NA for region-combined results.
#' @param fdr_cutoff Adjusted p-value threshold.
#' @param abs_log_fc_cutoff Absolute log2 fold-change threshold.
#' @param show_title Whether to show the cell type as the plot title.
#' @param chr_color_map Optional named character vector (or list) mapping chromosome -> color.
#'   If provided, it must cover every chromosome observed in the plotted data.
#'
#' @return A ggplot object.
#' @export
plot_de_volcano_gg <- function(de_dt,
                               cell_type_use,
                               region_use,
                               fdr_cutoff = 0.05,
                               abs_log_fc_cutoff = log2(1.05),
                               show_title = TRUE,
                               chr_color_map = NULL) {

  cell_type <- region <- log_fc <- adj_p_val <- chr <- NULL
  y <- is_sig <- col_group <- chr_draw <- draw_order <- NULL

  if (is.na(region_use)) {
    dt <- de_dt[cell_type == cell_type_use & is.na(region)]
  } else {
    dt <- de_dt[cell_type == cell_type_use & region == region_use]
  }

  rng <- dt[, max(abs(log_fc), na.rm = TRUE)]
  x_limits <- c(-rng, rng)

  adj_p_floor <- dt[, min(adj_p_val[adj_p_val > 0], na.rm = TRUE)]
  if (!is.finite(adj_p_floor)) {
    adj_p_floor <- .Machine$double.xmin
  }

  dt_plot <- data.table::copy(dt)
  dt_plot[, y := -log10(pmax(adj_p_val, adj_p_floor))]
  dt_plot[, is_sig := adj_p_val < fdr_cutoff]

  if (!is.null(chr_color_map)) {

    if (is.list(chr_color_map)) {
      chr_color_map <- unlist(chr_color_map, use.names = TRUE)
    }

    chroms_obs <- sort(unique(dt_plot[!is.na(chr), chr]))
    missing_chroms <- setdiff(chroms_obs, names(chr_color_map))
    if (length(missing_chroms) > 0) {
      stop(
        paste0(
          "chr_color_map does not cover all chromosomes in the data. Missing: ",
          paste(missing_chroms, collapse = ", ")
        )
      )
    }

    dt_plot[, col_group := chr]

    plot_order <- rev(names(chr_color_map))
    draw_map <- stats::setNames(seq_along(plot_order), plot_order)

    dt_plot[, chr_draw := draw_map[col_group]]
    dt_plot[, draw_order := (is_sig * 1000L) + chr_draw]

    color_values <- chr_color_map

  } else {

    dt_plot[, col_group := "all"]
    dt_plot[, draw_order := ifelse(is_sig, 2L, 1L)]

    color_values <- c(all = "black")
  }

  data.table::setorder(dt_plot, draw_order)

  p <- ggplot2::ggplot(
    dt_plot,
    ggplot2::aes(x = log_fc, y = y)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(color = col_group),
      size = 1.2
    ) +
    ggplot2::scale_x_continuous(limits = x_limits) +
    ggplot2::scale_color_manual(values = color_values, name = NULL) +
    ggplot2::geom_hline(yintercept = -log10(fdr_cutoff), linetype = 2) +
    ggplot2::geom_vline(xintercept = c(-abs_log_fc_cutoff, abs_log_fc_cutoff), linetype = 2) +
    ggplot2::labs(
      x = "Effect size, log2",
      y = "Adjusted p-value, -log10",
      title = if (show_title) cell_type_use else NULL
    )

  p
}

#' Scatter plot comparing DE effect sizes
#'
#' Produces an effect-size scatter plot for two DE results defined by cell type and
#' region pairs, and returns a ggplot object.
#'
#' All genes are plotted as black points. The plotting range is determined from the
#' genes that pass the adjusted p-value cutoff in at least one of the two DE results.
#'
#' Spearman rho^2 is computed on the same subset of genes, namely those significant
#' in at least one comparison, and is annotated in the top-left corner of the plot.
#'
#' @param de_dt A prepared DE data.table from prep_de/read_de_results.
#' @param cell_type_a First cell type.
#' @param cell_type_b Second cell type.
#' @param region_a Region for A, or NA for region-combined results.
#' @param region_b Region for B, or NA for region-combined results.
#' @param fdr_cutoff Adjusted p-value threshold used to define the subset for
#'   rho^2 calculation and plot limits.
#' @param xlab_prefix Optional prefix string added to the x-axis label.
#'
#' @return A ggplot object, or NULL if no genes are significant at the specified
#'   adjusted p-value cutoff in either comparison.
#' @export
plot_de_scatter_gg <- function(de_dt,
                               cell_type_a,
                               cell_type_b,
                               region_a = NA,
                               region_b = NA,
                               fdr_cutoff = 0.05,
                               xlab_prefix = NULL) {

  cell_type <- region <- NULL
  adj_p_val.x <- adj_p_val.y <- log_fc.x <- log_fc.y <- NULL

  if (is.na(region_a)) {
    x <- de_dt[cell_type == cell_type_a & is.na(region), ]
  } else {
    x <- de_dt[cell_type == cell_type_a & region == region_a, ]
  }

  if (is.na(region_b)) {
    y <- de_dt[cell_type == cell_type_b & is.na(region), ]
  } else {
    y <- de_dt[cell_type == cell_type_b & region == region_b, ]
  }

  name_a <- paste0(toupper(substr(cell_type_a, 1, 1)),
                   substr(cell_type_a, 2, nchar(cell_type_a)))
  name_b <- paste0(toupper(substr(cell_type_b, 1, 1)),
                   substr(cell_type_b, 2, nchar(cell_type_b)))

  m <- merge(x, y, by = c("chr", "gene"))

  sig_idx <- (m$adj_p_val.x < fdr_cutoff) | (m$adj_p_val.y < fdr_cutoff)

  if (!any(sig_idx)) {
    return(NULL)
  }

  rng <- m[sig_idx, max(abs(c(log_fc.x, log_fc.y)), na.rm = TRUE)]
  rng <- c(-rng, rng)

  xlab_string <- paste("Effect size, log2",
                       paste(name_a, region_a, sep = ", "),
                       sep = "\n")
  if (!is.null(xlab_prefix)) {
    xlab_string <- paste0(xlab_prefix, xlab_string)
  }

  ylab_string <- paste(paste(name_b, region_b, sep = ", "),
                       "Effect size, log2",
                       sep = "\n")

  m_plot <- data.table::as.data.table(m)

  ct <- m_plot[sig_idx,
               stats::cor.test(log_fc.x, log_fc.y, method = "spearman")]

  rho_sqrd <- round(as.numeric(ct$estimate)^2, 2)
  rho_label <- paste0("rho^2 == ", rho_sqrd)

  p <- ggplot2::ggplot(
    m_plot,
    ggplot2::aes(x = log_fc.x, y = log_fc.y)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = 2) +
    ggplot2::scale_x_continuous(limits = rng) +
    ggplot2::scale_y_continuous(limits = rng) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      x = xlab_string,
      y = ylab_string
    ) +
    ggplot2::annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = rho_label,
      parse = TRUE,
      hjust = -0.05,
      vjust = 1.1,
      size = 8
    ) +
    ggplot2::geom_point(
      color = "black",
      size = 1.2
    )

  p
}

#' Compute correlation matrix across cell_type x region groups
#'
#' Correlations are computed on log_fc among genes passing an FDR filter in either
#' of the compared groups. The returned value is signed rho^2.
#'
#' @param de_dt A prepared DE data.table from prep_de/read_de_results.
#' @param cell_types_use Character vector of cell types to include.
#' @param regions_use Character vector of regions to include.
#' @param non_neuron_types Character vector of non-neuronal cell types.
#' @param fdr_cutoff Adjusted p-value threshold.
#' @return A square numeric matrix with dimnames "cell_type__region".
#' @export
compute_de_cor_mat <- function(de_dt,
                               cell_types_use,
                               regions_use,
                               non_neuron_types,
                               fdr_cutoff = 0.05) {

    cell_type <- region <- cr <- gene <- adj_p_val <- log_fc <- NULL

    dt <- data.table::copy(de_dt)

    dt <- dt[cell_type %in% cell_types_use]
    dt <- dt[region %in% regions_use]

    # Exclude ic for neurons (preserved from original)
    dt <- dt[region != "ic" | cell_type %in% non_neuron_types]

    dt[, cell_type := factor(cell_type, levels = cell_types_use)]
    dt[, region := factor(region, levels = regions_use)]
    dt[, cr := paste(cell_type, region, sep = "__")]

    data.table::setorderv(dt, c("cell_type", "region"))

    keys <- unique(dt$cr)
    n <- length(keys)

    out_mat <- matrix(NA_real_,
                            nrow = n,
                            ncol = n,
                            dimnames = list(keys, keys))

    for (i in keys) {
        message(i)
        for (j in keys) {
            if (identical(i, j)) {
                out_mat[i, j] <- 1
                next
            }

            a <- dt[cr == i]
            b <- dt[cr == j]
            m <- merge(a, b, by = "gene")

            # Make R CMD CHECK Happy
            adj_p_val.x <- adj_p_val.y <- log_fc.x <- log_fc.y <- NULL

            ctest <- m[adj_p_val.x < fdr_cutoff | adj_p_val.y < fdr_cutoff,
                       stats::cor.test(log_fc.x, log_fc.y, method = "spearman")]

            r <- sign(ctest$estimate) * ctest$estimate^2
            out_mat[i, j] <- r
        }
    }

    out_mat
}

#' Plot a correlation heatmap with pheatmap
#'
#' @param cor_mat A numeric matrix produced by compute_de_cor_mat.
#' @param clustering_method Clustering method passed to pheatmap.
#' @param breaks Numeric breakpoints for the color scale.
#' @param palette_colors Vector of colors used for the palette.
#' @param show_dendrograms Logical; if FALSE, dendrograms are hidden but clustering order is preserved.
#' @return The value returned by pheatmap::pheatmap.
#' @export
plot_de_cor_heatmap <- function(cor_mat,
                                clustering_method = "complete",
                                breaks = seq(-1, 1, length.out = 101),
                                palette_colors = c("steelblue", "white", "darkorange"),
                                show_dendrograms = TRUE) {

  pal_fun <- grDevices::colorRampPalette(palette_colors)

  treeheight_row <- if (show_dendrograms) 50 else 0
  treeheight_col <- if (show_dendrograms) 50 else 0

  pheatmap::pheatmap(
    cor_mat,
    breaks = breaks,
    color = pal_fun(length(breaks) - 1),
    clustering_method = clustering_method,
    treeheight_row = treeheight_row,
    treeheight_col = treeheight_col
  )
}

#' Plot a correlation heatmap
#'
#' @param cor_mat A numeric matrix produced by compute_de_cor_mat.
#' @param clustering_method Clustering method passed to hclust.
#' @param breaks Numeric breakpoints for the color scale.
#' @param palette_colors Vector of colors used for the palette.
#' @param legend_title Title for the color legend.  If set to NULL, suppress the legend entirely.
#' @param show_dendrograms Logical; if FALSE, dendrograms are hidden but clustering order is preserved.
#' @return A ComplexHeatmap heatmap object.
#' @export
plot_de_cor_heatmap_complex <- function(cor_mat,
                                        clustering_method = "complete",
                                        breaks = seq(-1, 1, length.out = 101),
                                        palette_colors = c("steelblue", "white", "darkorange"),
                                        legend_title = "Correlation",
                                        show_dendrograms = TRUE) {

  col_fun <- circlize::colorRamp2(
    breaks = seq(min(breaks), max(breaks), length.out = length(palette_colors)),
    palette_colors
  )

  show_legend <- !is.null(legend_title)

  ComplexHeatmap::Heatmap(
    cor_mat,
    name = "cor_metric",
    col = col_fun,

    clustering_method_rows = clustering_method,
    clustering_method_columns = clustering_method,
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",

    row_dend_reorder = FALSE,
    column_dend_reorder = FALSE,
    show_row_dend = show_dendrograms,
    show_column_dend = show_dendrograms,

    column_names_rot = 90,
    column_names_gp = grid::gpar(fontsize = 10),
    row_names_gp = grid::gpar(fontsize = 10),

    rect_gp = grid::gpar(col = "grey85", lwd = 1),

    show_heatmap_legend = show_legend,

    heatmap_legend_param = list(
      title = legend_title
    )
  )
}

###############################################################################
# META CELL PARSING AND KMEANS CLUSTERING
###############################################################################

#' Read cell metadata
#'
#' @param cell_metadata_file Path to annotated cell metadata file.
#' @return A data.table of cell metadata.
#' @export
read_cell_metadata <- function(cell_metadata_file) {
  data.table::fread(cell_metadata_file)
}

#' Extract donor ages from cell metadata
#'
#' @param cell_metadata A data.table containing donor_external_id and age.
#' @return A named numeric vector of donor ages.
#' @export
extract_donor_ages <- function(cell_metadata) {
  donor_external_id <- age <- NULL

  donor_ages_dt <- unique(cell_metadata[, list(donor_external_id, age)])

  donor_ages <- as.numeric(donor_ages_dt$age)
  names(donor_ages) <- donor_ages_dt$donor_external_id

  sort(donor_ages)
}

#' Read and aggregate metacells
#' (originally load_metacells)
#'
#' @param path Path to metacell count matrix file.
#' @param cell_types_use Character vector of cell types to retain.
#' @param regions_use Character vector of regions to retain.
#' @return A list with elements:
#'   - metacells: dgCMatrix of aggregated TPM
#'   - col_metadata: data.table of donor/cell_type/region
#' @export
read_metacells <- function(path,
                           cell_types_use,
                           regions_use = c("CaH", "Pu", "NAC", "ic", "DFC")) {

  gene <- donor <- village <- cell_type <- region <- single_cell_assay <- group <- NULL

  dt <- data.table::fread(path)
  names(dt)[1] <- "gene"

  mat <- as.matrix(dt, rownames = "gene")

  stopifnot(requireNamespace("Matrix", quietly = TRUE))
  mat <- methods::as(mat, "dgCMatrix")

  col_metadata <- data.table::as.data.table(
    do.call(rbind, strsplit(colnames(mat), "__"))
  )
  data.table::setnames(
    col_metadata,
    c("donor", "village", "cell_type", "region", "single_cell_assay")
  )

  keep <- col_metadata[, cell_type %in% cell_types_use & region %in% regions_use]
  col_metadata <- col_metadata[keep]
  mat <- mat[, keep]

  col_metadata[, group := paste(donor, cell_type, region, sep = "__")]
  f <- factor(col_metadata$group)

  dm <- Matrix::sparseMatrix(
    i = seq_along(f),
    j = as.integer(f),
    x = 1,
    dims = c(length(f), nlevels(f)),
    dimnames = list(NULL, levels(f))
  )

  mat <- mat %*% dm

  cs <- Matrix::colSums(mat)
  d <- Matrix::Diagonal(x = 1e6 / cs)
  mat_tpm <- mat %*% d
  colnames(mat_tpm) <- colnames(mat)

  col_metadata <- data.table::as.data.table(
    do.call(rbind, strsplit(colnames(mat_tpm), "__"))
  )
  data.table::setnames(col_metadata, c("donor", "cell_type", "region"))

  list(metacells = mat_tpm, col_metadata = col_metadata)
}


#' Fast row-wise summary statistics for dense matrices
#'
#' Computes median, MAD, selected quantiles, and non-missing counts
#' for each row of a numeric matrix.
#'
#' @param mat Numeric matrix (genes x samples).
#' @return A list with elements:
#'   - median
#'   - mad
#'   - q_10
#'   - q_25
#'   - q_75
#'   - q_90
#'   - n (number of non-NA values per row)
#'
#' @keywords internal
row_stats_block_fast <- function(mat) {

  med <- matrixStats::rowMedians(mat, na.rm = TRUE)

  mad <- matrixStats::rowMads(
    mat,
    na.rm = TRUE,
    constant = 1.4826
  )

  q <- matrixStats::rowQuantiles(
    mat,
    probs = c(0.10, 0.25, 0.75, 0.90),
    na.rm = TRUE,
    type = 7
  )

  n <- rowSums(!is.na(mat))

  list(
    median = med,
    mad = mad,
    q_10 = q[, 1],
    q_25 = q[, 2],
    q_75 = q[, 3],
    q_90 = q[, 4],
    n = n
  )
}

#' Summarize metacells by gene and cell_type-region
#'
#' @param metacells dgCMatrix of TPM values.
#' @param col_metadata data.table returned by read_metacells.
#' @param donor_ages Named numeric vector of donor ages.
#' @return A data.table summary.
#' @export
summarize_metacells <- function(metacells,
                                col_metadata,
                                donor_ages) {

  donor <- age <- age_bin <- cell_type <- region <- cr <- genes <- NULL

  col_metadata <- data.table::copy(col_metadata)

  col_metadata[, age := unname(donor_ages[donor])]

  col_metadata[, age_bin := cut(
    age,
    breaks = c(-Inf, 39, 49, 59, 69, 79, 89, Inf),
    labels = c("30", "40", "50", "60", "70", "80", "90"),
    right = TRUE
  )]

  col_metadata[, cr := paste(cell_type, region, sep = "__")]

  groups <- split(seq_len(nrow(col_metadata)), col_metadata$cr)

  genes <- rownames(metacells)
  age_bins <- c("30", "40", "50", "60", "70", "80", "90")

  res_list <- vector("list", length(groups))
  names(res_list) <- names(groups)

  i <- 0L

  for (gn in names(groups)) {

    message(gn)

    i <- i + 1L
    cols <- groups[[gn]]

    sub_sp <- metacells[, cols, drop = FALSE]
    sub_dense <- as.matrix(sub_sp)

    st <- row_stats_block_fast(sub_dense)

    cr_parts <- data.table::tstrsplit(gn, "__", fixed = TRUE)
    ct <- cr_parts[[1]][1]
    rg <- cr_parts[[2]][1]

    dt_out <- data.table::data.table(
      gene = genes,
      cell_type = ct,
      region = rg,
      median = as.numeric(st$median),
      mad = as.numeric(st$mad),
      q_10 = as.numeric(st$q_10),
      q_25 = as.numeric(st$q_25),
      q_75 = as.numeric(st$q_75),
      q_90 = as.numeric(st$q_90),
      n_donors = as.integer(st$n)
    )

    bins_for_cols <- col_metadata$age_bin[cols]

    for (b in age_bins) {

      idx <- which(bins_for_cols == b)
      colname <- paste0("median_", b)

      if (length(idx) == 0L) {
        dt_out[[colname]] <- NA_real_
      } else if (length(idx) == 1L) {
        dt_out[[colname]] <- as.numeric(sub_dense[, idx])
      } else {
        dt_out[[colname]] <- matrixStats::rowMedians(
          sub_dense[, idx, drop = FALSE]
        )
      }
    }

    res_list[[i]] <- dt_out
  }

  #Make R CMD CHECK Happy
  gene <- cell_type <- region <- NULL

  out <- data.table::rbindlist(res_list, use.names = TRUE, fill = TRUE)
  data.table::setorder(out, gene, cell_type, region)

  out
}

#' Split metacells by cell_type and region
#'
#' @param metacells dgCMatrix of TPM values.
#' @param col_metadata data.table returned by read_metacells.
#' @param donor_ages Named numeric vector of donor ages.
#' @return A named list of dgCMatrix objects, sorted by donor age.
#' @export
split_metacells_by_cell_type_region <- function(metacells,
                                                col_metadata,
                                                donor_ages) {

  # Make R CMD CHECK Happy
  donor <- cell_type <- region <- cr <- age<- NULL

  col_metadata <- data.table::copy(col_metadata)

  col_metadata[, age := unname(donor_ages[donor])]
  col_metadata[, cr := paste(cell_type, region, sep = "__")]

  groups <- split(seq_len(nrow(col_metadata)), col_metadata$cr)

  out <- vector("list", length(groups))
  names(out) <- names(groups)

  i <- 0L

  for (gn in names(groups)) {

    i <- i + 1L
    cols <- groups[[gn]]

    sub_sp <- metacells[, cols, drop = FALSE]

    donors <- col_metadata$donor[cols]
    ages <- unname(donor_ages[donors])

    ord <- order(ages, donors)

    sub_sp <- sub_sp[, ord, drop = FALSE]
    donors <- donors[ord]

    colnames(sub_sp) <- donors

    out[[i]] <- sub_sp
  }

  out
}

#' Prepare DE matrices for clustering and heatmaps
#' (originally prepare_de_matrices)
#'
#' @param de_dt Prepared DE data.table.
#' @param metacell_summary Output from summarize_metacells.
#' @param cell_types_use Character vector of cell types.
#' @param fdr_cutoff Adjusted p-value threshold.
#' @param abs_lfc_cutoff Absolute log2 fold-change threshold.
#' @param min_tpm Minimum median TPM in either CaH or DFC.
#' @param regions_use Regions to consider (currently CaH and DFC logic retained).
#' @return A list with lfc_mat, fdr_mat, lfc_mat_z.
#' @export
prep_de_matrices <- function(de_dt,
                             metacell_summary,
                             cell_types_use,
                             fdr_cutoff = 0.01,
                             abs_lfc_cutoff = log2(1.05),
                             min_tpm = 10,
                             regions_use = c("CaH", "DFC")) {

  gene <- cell_type <- region <- median <- log_fc <- adj_p_val <- NULL

  dt <- data.table::copy(de_dt)
  dt <- dt[cell_type %in% cell_types_use]

  dt <- merge(
    dt,
    metacell_summary[region == "CaH", list(gene, cell_type, median_tpm_ca = median)],
    by = c("gene", "cell_type"),
    all.x = TRUE
  )

  dt <- merge(
    dt,
    metacell_summary[region == "DFC", list(gene, cell_type, median_tpm_dfc = median)],
    by = c("gene", "cell_type"),
    all.x = TRUE
  )

  # Make R CMD CHECK Happy
  median_tpm_ca <- median_tpm_dfc <- NULL

  dt <- dt[median_tpm_ca > min_tpm | median_tpm_dfc > min_tpm]

  lfc_mat <- as.matrix(
    data.table::dcast(dt, gene ~ cell_type, value.var = "log_fc"),
    rownames = "gene"
  )[, cell_types_use, drop = FALSE]

  fdr_mat <- as.matrix(
    data.table::dcast(dt, gene ~ cell_type, value.var = "adj_p_val"),
    rownames = "gene"
  )[, cell_types_use, drop = FALSE]

  n_sig <- rowSums(
    (fdr_mat < fdr_cutoff) & (abs(lfc_mat) > abs_lfc_cutoff),
    na.rm = TRUE
  )

  lfc_mat <- lfc_mat[n_sig > 0, , drop = FALSE]
  fdr_mat <- fdr_mat[n_sig > 0, , drop = FALSE]

  lfc_mat[is.na(lfc_mat)] <- 0
  fdr_mat[is.na(fdr_mat)] <- 1

  lfc_mat_z <- scale(lfc_mat, scale = TRUE, center = TRUE)

  list(lfc_mat = lfc_mat, fdr_mat = fdr_mat, lfc_mat_z = lfc_mat_z)
}


#' Prepare region-specific LFC matrix
#' (originally prepare_region_lfc_matrix)
#'
#' @param de_dt Prepared DE data.table.
#' @param genes_use Genes to retain.
#' @param cell_types_use Character vector of cell types.
#' @param regions_use Character vector of regions.
#' @return A numeric matrix of log fold-changes.
#' @export
prep_region_lfc_matrix <- function(de_dt,
                                   genes_use,
                                   cell_types_use,
                                   regions_use) {

  gene <- cell_type <- region <- log_fc <- NULL

  dt <- data.table::copy(de_dt)

  dt <- dt[
    cell_type %in% cell_types_use &
      region %in% regions_use
  ]

  lfc_dt <- data.table::dcast(
    dt,
    gene ~ cell_type + region,
    value.var = "log_fc",
    sep = "__"
  )

  ## IMPORTANT: avoid coercion to character by using gene as rownames
  lfc_mat <- as.matrix(lfc_dt, rownames = "gene")

  ## optional: stop early if genes are missing
  missing_genes <- setdiff(genes_use, rownames(lfc_mat))
  if (length(missing_genes) > 0L) {
    stop(
      "prep_region_lfc_matrix: ", length(missing_genes),
      " genes in genes_use are missing from the cast matrix. Example: ",
      missing_genes[[1]]
    )
  }

  lfc_mat <- lfc_mat[genes_use, , drop = FALSE]

  tmp <- data.table::as.data.table(
    do.call(rbind, strsplit(colnames(lfc_mat), "__"))
  )

  # Make R CMD CHECK Happy
  V1 <- V2 <- NULL

  tmp[, V1 := factor(V1, levels = cell_types_use)]
  tmp[, V2 := factor(V2, levels = regions_use)]

  col_order <- order(tmp$V1, tmp$V2)
  lfc_mat <- lfc_mat[, col_order, drop = FALSE]

  lfc_mat
}
#' Plot average silhouette width across k for k-means clustering
#' (originally plot_k_means_silhouette)
#'
#' Computes k-means clustering across a range of cluster numbers and
#' visualizes the average silhouette width for each k.
#'
#' @param mat Numeric matrix used for clustering.
#'   Rows represent the units being clustered (typically genes),
#'   and columns represent features (e.g., scaled log fold-changes).
#'   The matrix should be numeric and is typically centered and scaled
#'   prior to calling this function.
#'
#' @param ks Integer vector of cluster numbers (k) to evaluate. Default is 10:30.
#'
#' @return Invisibly returns NULL. Produces a base R plot as a side effect.
#'
#' @details
#' - Euclidean distance is computed on rows of `mat`.
#' - K-means clustering is performed for k = 10:30.
#' - For each k, clustering is run with nstart = 200 and iter.max = 20.
#' - Average silhouette width is computed using cluster::silhouette.
#' - A line plot of k versus average silhouette width is produced.
#'
#' @export
plot_kmeans_silhouette <- function(mat, ks = 10:30) {

  d <- stats::dist(mat)

  avg_sil <- sapply(ks, function(k) {

    set.seed(42)
    km <- stats::kmeans(mat, centers = k, nstart = 200, iter.max = 20)

    sil <- cluster::silhouette(km$cluster, d)

    mean(sil[, "sil_width"])
  })

  graphics::plot(
    ks,
    avg_sil,
    type = "b",
    xlab = "Number of clusters",
    ylab = "Average silhouette width"
  )

  invisible(NULL)
}

#' Plot k-means clustered heatmap of log fold-changes
#' (originally plot_k_means_heatmap)
#'
#' Performs k-means clustering on a gene-by-feature matrix and visualizes
#' the corresponding log fold-change matrix as a heatmap, ordered by cluster.
#'
#' @param k_means_mat Numeric matrix used for clustering.
#'   Rows represent genes (or other units being clustered). Rownames must be
#'   non-NULL and must correspond to gene identifiers. Typically this matrix is
#'   a scaled (e.g., z-scored) version of `lfc_mat`, since k-means is
#'   scale-sensitive.
#'
#' @param lfc_mat Numeric matrix of log fold-changes used for visualization.
#'   Rows represent genes and must include all genes present in `k_means_mat`
#'   (at least those that remain after optional cluster dropping/releveling).
#'   Columns typically represent cell types or cell_type-region combinations.
#'
#' @param scaling_factor Numeric multiplier applied to `lfc_mat` prior to
#'   plotting. This is used purely for visualization scaling and does not
#'   affect clustering.
#'
#' @param k Integer number of clusters (centers) for k-means.
#'
#' @param cluster_level_order Optional integer vector defining the desired
#'   ordering of k-means cluster labels. If NULL, cluster labels are used as
#'   returned by k-means (no reordering and no cluster dropping). If non-NULL,
#'   clusters not present in this vector will be dropped unless
#'   `allow_drop_clusters = FALSE`.
#'
#' @param allow_drop_clusters Logical; if FALSE, stop when `cluster_level_order`
#'   would drop any cluster labels present in the k-means solution.
#'
#' @return A named integer vector of cluster assignments for genes, where names
#'   correspond to gene identifiers (after any optional dropping implied by
#'   `cluster_level_order`).
#'
#' @details
#' - Euclidean k-means clustering is performed with `stats::kmeans()` using
#'   `centers = k`, `nstart = 200`, and `iter.max = 20`.
#' - If `cluster_level_order` is provided, cluster labels are reordered
#'   according to that order; optionally, clusters not included may be dropped.
#' - The heatmap is generated using `pheatmap::pheatmap()` without column or row
#'   clustering. Columns are shown in the existing order of `lfc_mat`.
#'
#' @export
plot_kmeans_heatmap <- function(k_means_mat,
                                lfc_mat,
                                scaling_factor,
                                k = 19,
                                cluster_level_order = c(2, 10, 6, 3, 5, 14, 9, 13, 4, 1, 15, 19, 18, 7, 17, 8, 11, 12),
                                allow_drop_clusters = TRUE) {

  ## -----------------------
  ## Assertions / validation
  ## -----------------------

  if (!is.matrix(k_means_mat) || !is.numeric(k_means_mat)) {
    stop("k_means_mat must be a numeric matrix.")
  }

  if (is.null(rownames(k_means_mat))) {
    stop("k_means_mat must have non-NULL rownames (gene identifiers).")
  }

  if (!is.matrix(lfc_mat) || !is.numeric(lfc_mat)) {
    stop("lfc_mat must be a numeric matrix.")
  }

  if (is.null(rownames(lfc_mat))) {
    stop("lfc_mat must have non-NULL rownames (gene identifiers).")
  }

  if (!is.numeric(scaling_factor) || length(scaling_factor) != 1L || is.na(scaling_factor)) {
    stop("scaling_factor must be a single non-NA numeric value.")
  }

  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 2) {
    stop("k must be a single numeric/integer value >= 2.")
  }
  k <- as.integer(k)

  if (!is.null(cluster_level_order)) {

    if (!is.numeric(cluster_level_order) || length(cluster_level_order) < 1L) {
      stop("cluster_level_order must be NULL or a non-empty numeric/integer vector.")
    }

    if (anyDuplicated(cluster_level_order)) {
      stop("cluster_level_order must not contain duplicate entries.")
    }
  }

  ## -----------------------
  ## K-means
  ## -----------------------

  set.seed(42)
  km <- stats::kmeans(
    k_means_mat,
    centers = k,
    nstart = 200,
    iter.max = 20
  )

  gn <- rownames(k_means_mat)

  missing_in_lfc <- setdiff(gn, rownames(lfc_mat))
  if (length(missing_in_lfc) > 0L) {
    stop(
      "lfc_mat is missing ", length(missing_in_lfc),
      " genes present in k_means_mat. Example: ",
      missing_in_lfc[[1]]
    )
  }

  present_clusters <- sort(unique(km$cluster))

  ## -----------------------
  ## Cluster labeling / ordering
  ## -----------------------

  if (is.null(cluster_level_order)) {

    ## No releveling; keep k-means labels as-is
    k_use <- as.integer(km$cluster)
    names(k_use) <- gn

  } else {

    dropped_clusters <- setdiff(present_clusters, cluster_level_order)

    if (length(dropped_clusters) > 0L && !isTRUE(allow_drop_clusters)) {
      stop(
        "cluster_level_order would drop these cluster labels present in the k-means solution: ",
        paste(dropped_clusters, collapse = ", "),
        ". Set allow_drop_clusters=TRUE to allow this."
      )
    }

    ## Relevel and optionally drop clusters not in the ordering vector
    k_use <- factor(km$cluster, levels = cluster_level_order)
    k_use <- as.numeric(k_use)
    names(k_use) <- gn

    keep <- !is.na(k_use)
    k_use <- k_use[keep]

    if (length(k_use) == 0L) {
      stop("After applying cluster_level_order, no genes remain (all clusters were dropped).")
    }
  }

  gene_order <- names(k_use)[order(k_use)]
  boundaries <- cumsum(table(k_use))

  ## -----------------------
  ## Plot
  ## -----------------------

  pheatmap::pheatmap(
    t(lfc_mat[gene_order, , drop = FALSE] * scaling_factor),
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    breaks = seq(-1, 1, length.out = 101),
    color = grDevices::colorRampPalette(c("steelblue", "white", "darkorange"))(100),
    show_colnames = FALSE,
    gaps_col = boundaries
  )

  k_use
}

#' Plot k-means clustered heatmap of log fold-changes
#' (originally plot_k_means_heatmap)
#'
#' Performs k-means clustering on a gene-by-feature matrix and visualizes
#' the corresponding log fold-change matrix as a heatmap, ordered by cluster.
#'
#' @param k_means_mat Numeric matrix used for clustering.
#'   Rows represent genes (or other units being clustered). Rownames must be
#'   non-NULL and must correspond to gene identifiers. Typically this matrix is
#'   a scaled (e.g., z-scored) version of `lfc_mat`, since k-means is
#'   scale-sensitive.
#'
#' @param lfc_mat Numeric matrix of log fold-changes used for visualization.
#'   Rows represent genes and must include all genes present in `k_means_mat`
#'   (at least those that remain after optional cluster dropping/releveling).
#'   Columns typically represent cell types or cell_type-region combinations.
#'
#' @param scaling_factor Numeric multiplier applied to `lfc_mat` prior to
#'   plotting. This is used purely for visualization scaling and does not
#'   affect clustering.
#'
#' @param k Integer number of clusters (centers) for k-means.
#'
#' @param cluster_level_order Optional integer vector defining the desired
#'   ordering of k-means cluster labels. If NULL, cluster labels are used as
#'   returned by k-means (no reordering and no cluster dropping). If non-NULL,
#'   clusters not present in this vector will be dropped unless
#'   `allow_drop_clusters = FALSE`.
#'
#' @param allow_drop_clusters Logical; if FALSE, stop when `cluster_level_order`
#'   would drop any cluster labels present in the k-means solution.
#'
#' @param fontsize_col Integer font size for column labels in the heatmap.
#' @param fontsize_row Integer font size for row labels in the heatmap.
#'
#' @return A named integer vector of cluster assignments for genes, where names
#'   correspond to gene identifiers (after any optional dropping implied by
#'   `cluster_level_order`).
#'
#' @details
#' - Euclidean k-means clustering is performed with `stats::kmeans()` using
#'   `centers = k`, `nstart = 200`, and `iter.max = 20`.
#' - If `cluster_level_order` is provided, cluster labels are reordered
#'   according to that order; optionally, clusters not included may be dropped.
#' - The heatmap is generated using `pheatmap::pheatmap()` without column or row
#'   clustering. Columns are shown in the existing order of `lfc_mat`.
#' - Vertical separators between clusters are implemented by inserting
#'   explicit `NA` columns into the plotted matrix. These columns are rendered
#'   as solid black bars via `na_col = "black"`. This approach avoids fragile
#'   grid coordinate overlays and ensures that separators scale correctly and
#'   remain device-independent, since they are part of the matrix itself rather
#'   than graphical annotations layered on top.
#'
#' @export
plot_kmeans_heatmap_with_cluster_labels <- function(k_means_mat,
                                lfc_mat,
                                scaling_factor,
                                k = 19,
                                cluster_level_order = c(2, 10, 6, 3, 5, 14, 9, 13, 4, 1, 15, 19, 18, 7, 17, 8, 11, 12),
                                allow_drop_clusters = TRUE,
                                fontsize_col=10,
                                fontsize_row=10) {

    ## -----------------------
    ## Assertions / validation
    ## -----------------------

    if (!is.matrix(k_means_mat) || !is.numeric(k_means_mat)) {
      stop("k_means_mat must be a numeric matrix.")
    }

    if (is.null(rownames(k_means_mat))) {
      stop("k_means_mat must have non-NULL rownames (gene identifiers).")
    }

    if (!is.matrix(lfc_mat) || !is.numeric(lfc_mat)) {
      stop("lfc_mat must be a numeric matrix.")
    }

    if (is.null(rownames(lfc_mat))) {
      stop("lfc_mat must have non-NULL rownames (gene identifiers).")
    }

    if (!is.numeric(scaling_factor) || length(scaling_factor) != 1L || is.na(scaling_factor)) {
      stop("scaling_factor must be a single non-NA numeric value.")
    }

    if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 2) {
      stop("k must be >= 2.")
    }
    k <- as.integer(k)

    if (!is.null(cluster_level_order)) {

      if (!is.numeric(cluster_level_order) || length(cluster_level_order) < 1L) {
        stop("cluster_level_order must be NULL or numeric.")
      }

      if (anyDuplicated(cluster_level_order)) {
        stop("cluster_level_order must not contain duplicates.")
      }
    }

    ## -----------------------
    ## K-means
    ## -----------------------
    num_genes=dim (lfc_mat)[1]

    set.seed(42)
    km <- stats::kmeans(
      k_means_mat,
      centers = k,
      nstart = 200,
      iter.max = 20
    )

    gn <- rownames(k_means_mat)

    missing_in_lfc <- setdiff(gn, rownames(lfc_mat))
    if (length(missing_in_lfc) > 0L) {
      stop("lfc_mat missing genes. Example: ", missing_in_lfc[[1]])
    }

    present_clusters <- sort(unique(km$cluster))

    ## -----------------------
    ## Cluster labeling / ordering
    ## -----------------------

    if (is.null(cluster_level_order)) {

      k_use <- as.integer(km$cluster)
      names(k_use) <- gn

    } else {

      dropped_clusters <- setdiff(present_clusters, cluster_level_order)

      if (length(dropped_clusters) > 0L && !isTRUE(allow_drop_clusters)) {
        stop("cluster_level_order would drop clusters: ",
             paste(dropped_clusters, collapse = ", "))
      }

      k_use <- factor(km$cluster, levels = cluster_level_order)
      k_use <- as.numeric(k_use)
      names(k_use) <- gn

      keep <- !is.na(k_use)
      k_use <- k_use[keep]

      if (length(k_use) == 0L) {
        stop("All clusters dropped.")
      }
    }

    gene_order <- names(k_use)[order(k_use)]
    boundaries <- cumsum(table(k_use))

    ## -----------------------
    ## Plot with thick separators
    ## -----------------------

    mat_plot <- t(lfc_mat[gene_order, , drop = FALSE] * scaling_factor)

    ## Replace pre-existing NA values (from lfc_mat) so na_col can be reserved
    ## exclusively for separator columns.
    na_before <- sum(is.na(mat_plot))
    if (na_before > 0L) {
      message("plot_kmeans_heatmap: replacing ", na_before, " NA values in mat_plot with 0 for visualization.")
      mat_plot[is.na(mat_plot)] <- 0
    }

    n_cols <- ncol(mat_plot)

    gap_after <- boundaries
    if (length(gap_after) >= 1L) {
      gap_after <- gap_after[-length(gap_after)]
    }
    gap_after <- as.integer(gap_after)

    sep_width <- 6L

    message("plot_kmeans_heatmap: n_cols = ", n_cols,
            ", n_clusters = ", length(boundaries),
            ", n_separators = ", length(gap_after),
            ", sep_width = ", sep_width)

    n_sep_total <- length(gap_after) * sep_width
    n_cols_exp <- n_cols + n_sep_total

    mat_exp <- matrix(NA_real_, nrow = nrow(mat_plot), ncol = n_cols_exp)
    rownames(mat_exp) <- rownames(mat_plot)

    orig_to_exp <- integer(n_cols)

    exp_col <- 1L
    sep_set <- rep(FALSE, n_cols)
    if (length(gap_after) > 0L) {
      sep_set[gap_after] <- TRUE
    }

    for (j in seq_len(n_cols)) {

      mat_exp[, exp_col] <- mat_plot[, j]
      orig_to_exp[j] <- exp_col
      exp_col <- exp_col + 1L

      if (sep_set[j]) {
        exp_col <- exp_col + sep_width
      }
    }

    ## Compute expanded cluster midpoints
    starts <- c(1L, boundaries[-length(boundaries)] + 1L)
    ends <- boundaries

    start_exp <- orig_to_exp[starts]
    end_exp <- orig_to_exp[ends]
    mid_exp <- as.integer(round((start_exp + end_exp) / 2))

    labels_col <- rep("", n_cols_exp)
    labels_col[mid_exp] <- as.character(seq_along(boundaries))

    #replace the fixed breaks definition with one that scales
    #with the scaling value passed in.
    #this will preserve the same results as scaling_factor=5 and breaks = seq(-1, 1, length.out = 101)
    #for the original unscaled lfc_mat.

    lfc_limit <- 0.2

    breaks_vec <- seq(
      -lfc_limit * scaling_factor,
      lfc_limit * scaling_factor,
      length.out = 101
    )

    #add some extra margin for text.
    ph <- pheatmap::pheatmap(
      mat_exp,
      cluster_cols = FALSE,
      cluster_rows = FALSE,
      breaks = breaks_vec,
      color = grDevices::colorRampPalette(c("steelblue", "white", "darkorange"))(100),
      border_color = NA,
      na_col = "black",
      show_colnames = TRUE,
      labels_col = labels_col,
      angle_col = 0,
      silent = TRUE,
      fontsize_col = fontsize_col,
      fontsize_row = fontsize_row
    )

    gt <- ph$gtable

    ## Add top margin
    gt <- gtable::gtable_add_rows(
      gt,
      heights = grid::unit(6, "mm"),
      pos = 0
    )

    ## Add bottom margin (this is where the label will go)
    gt <- gtable::gtable_add_rows(
      gt,
      heights = grid::unit(8, "mm"),
      pos = nrow(gt)
    )

    ## Shrink the original pheatmap content so margins fit on the device
    shrink_factor <- 0.95  # try 0.85–0.95 as needed
    if (nrow(gt) > 2L) {
      gt$heights[2:(nrow(gt) - 1L)] <- gt$heights[2:(nrow(gt) - 1L)] * shrink_factor
    }

    labelStr <- paste0("Genes (n=", dim(lfc_mat)[1], ")")

    gt <- gtable::gtable_add_grob(
      gt,
      grobs = grid::textGrob(
        labelStr,
        y = grid::unit(0.8, "npc"),
        gp = grid::gpar(fontsize = 12)
      ),
      t = nrow(gt),
      l = 1,
      r = ncol(gt)
    )

    grid::grid.newpage()
    grid::grid.draw(gt)

    invisible(k_use)
}

#' Write lightweight DE outputs (summary + top up/down gene tables)
#'
#' Writes three files:
#' 1) A plain-text summary of inputs and DE results
#' 2) A TSV of top upregulated genes (by t-statistic)
#' 3) A TSV of top downregulated genes (by t-statistic)
#'
#' The DE table is subset to `cell_type_use`. Expression summaries are taken
#' from `metacell_summary` for the specified `cell_type_use` and `region_use`.
#'
#' @param de_dt data.table containing DE results with at least columns:
#'   gene, cell_type, log_fc, log_fc_se, t, adj_p_val.
#' @param metacell_summary data.table containing per-gene expression summaries
#'   with at least columns: gene, cell_type, region, n_donors, and
#'   median_30 ... median_90.
#' @param donor_ages Named numeric vector of donor ages (used for summary text).
#' @param cell_type_use Character scalar; cell type to extract.
#' @param region_use Character scalar; region to extract expression medians from.
#' @param out_name Basename for output files (no extension).
#' @param out_dir Output directory path.
#' @param n_top Maximum number of top up and top down genes to write.
#' @param fdr_thresh Adjusted p-value threshold.
#'
#' @return Invisibly returns a list with file paths and the two output tables.
#' @export
write_de_lite <- function(de_dt,
                          metacell_summary,
                          donor_ages,
                          cell_type_use,
                          region_use,
                          out_name,
                          out_dir,
                          n_top = 200,
                          fdr_thresh = 0.05) {

  gene <- cell_type <- region <- adj_p_val <- log_fc <- log_fc_se <- t <- median <- n_donors <- NULL
  median_30 <- median_40 <- median_50 <- median_60 <- median_70 <- median_80 <- median_90 <- NULL

  if (!dir.exists(out_dir)) {
    stop("out_dir does not exist: ", out_dir)
  }

  if (!is.numeric(n_top) || length(n_top) != 1L || is.na(n_top) || n_top < 1) {
    stop("n_top must be a positive scalar.")
  }
  n_top <- as.integer(n_top)

  if (!is.numeric(fdr_thresh) || length(fdr_thresh) != 1L || is.na(fdr_thresh) || fdr_thresh <= 0 || fdr_thresh > 1) {
    stop("fdr_thresh must be a scalar in (0, 1].")
  }

  de_dt <- data.table::copy(de_dt)
  metacell_summary <- data.table::copy(metacell_summary)

  ## subset DE to requested cell type
  de_dt <- de_dt[cell_type == cell_type_use]

  if (nrow(de_dt) == 0L) {
    stop("No DE rows found for cell_type_use='", cell_type_use, "'.")
  }

  ## medians for requested cell type + region
  medians <- metacell_summary[
    cell_type == cell_type_use & region == region_use,
    list(gene, median_30, median_40, median_50, median_60, median_70, median_80, median_90)
  ]

  if (nrow(medians) == 0L) {
    stop(
      "No metacell_summary rows found for cell_type_use='",
      cell_type_use,
      "' and region_use='",
      region_use,
      "'."
    )
  }

  ## rename medians columns:
  ## median_30 -> median_30s, etc. (preserving the collaborator’s intent)
  idx_dec <- which(!(names(medians) %in% c("gene")))
  data.table::setnames(medians, idx_dec, paste0(names(medians)[idx_dec], "s"))

  ## build up table (up)
  up_genes <- de_dt[
    adj_p_val < fdr_thresh & log_fc > 0,
    list(gene, log_fc, log_fc_se, t_stat = t, fdr = adj_p_val)
  ]

  up_genes <- merge(up_genes, medians, by = "gene", all.x = TRUE)
  data.table::setorderv(up_genes, "t_stat", -1)
  if (nrow(up_genes) > 0L) {
    up_genes <- up_genes[seq_len(min(n_top, nrow(up_genes)))]
  }

  ## build down table (down)
  down_genes <- de_dt[
    adj_p_val < fdr_thresh & log_fc < 0,
    list(gene, log_fc, log_fc_se, t_stat = t, fdr = adj_p_val)
  ]

  down_genes <- merge(down_genes, medians, by = "gene", all.x = TRUE)
  data.table::setorderv(down_genes, "t_stat", 1)
  if (nrow(down_genes) > 0L) {
    down_genes <- down_genes[seq_len(min(n_top, nrow(down_genes)))]
  }

  ## capture summary messages as text
  summary_lines <- character(0)

  add_line <- function(...) {
    summary_lines <<- c(summary_lines, paste0(...))
  }

  add_line("Input summary")
  add_line("cell type: ", cell_type_use)
  add_line("brain region: ", region_use)

  n_donors_val <- metacell_summary[
    cell_type == cell_type_use & region == region_use,
    n_donors
  ][1]

  add_line("N donors: ", n_donors_val)

  donor_ages_num <- as.numeric(donor_ages)
  donor_ages_num <- donor_ages_num[!is.na(donor_ages_num)]

  if (length(donor_ages_num) == 0L) {
    add_line("donor age range: NA - NA")
    add_line("donors median age: NA")
    add_line("donors per decade bin: NA")
  } else {
    add_line("donor age range: ", min(donor_ages_num), " - ", max(donor_ages_num))
    add_line("donors median age: ", stats::median(donor_ages_num))
    str <- convert_ages_to_decade_string(donor_ages_num)
    add_line("donors per decade bin: ", str)
  }

  # Make R CMD CHECK Happy
  .N <- NULL

  add_line("")
  add_line("Results summary")
  add_line("genes tested: ", de_dt[, .N])

  n_sig <- de_dt[adj_p_val < fdr_thresh, .N]
  n_up <- de_dt[adj_p_val < fdr_thresh & log_fc > 0, .N]
  n_dn <- de_dt[adj_p_val < fdr_thresh & log_fc < 0, .N]

  add_line("FDR < ", fdr_thresh, " total: ", n_sig)
  add_line("FDR < ", fdr_thresh, " up: ", n_up)
  add_line("FDR < ", fdr_thresh, " down: ", n_dn)

  pct_down <- (n_dn / n_sig) * 100
  add_line(
    "percent down FDR < ",
    fdr_thresh,
    ": ",
    ifelse(is.na(pct_down), "NA", round(pct_down, 3))
  )

  ## out paths
  summary_path <- file.path(out_dir, paste0(out_name, "_summary.txt"))
  up_path <- file.path(out_dir, paste0(out_name, "_up_genes.txt"))
  down_path <- file.path(out_dir, paste0(out_name, "_down_genes.txt"))

  summary_text <- paste(
    paste(summary_lines, collapse = "\n"),
    "",
    "Additional files",
    basename(up_path),
    basename(down_path),
    sep = "\n"
  )

  writeLines(summary_text, con = summary_path)

  data.table::fwrite(up_genes, file = up_path, sep = "\t", quote = FALSE, na = "NA")
  data.table::fwrite(down_genes, file = down_path, sep = "\t", quote = FALSE, na = "NA")

  invisible(list(
    summary_path = summary_path,
    up_path = up_path,
    down_path = down_path,
    up_genes = up_genes,
    down_genes = down_genes,
    summary_lines = summary_lines
  ))
}


#' Convert ages to a compact decade-bin summary string
#' (originally convert_ages_to_decade_string)
#'
#' @param donor_ages Numeric vector of donor ages.
#' @return Character scalar like "30s N = 12, 40s N = 20, ..."
#' @keywords internal
convert_ages_to_decade_string <- function(donor_ages) {

  donor_ages_chr <- as.character(donor_ages)
  dd <- substr(donor_ages_chr, 1, 1)

  ## preserve original behavior: treat 20s as 30s
  dd[dd == "2"] <- "3"

  dd <- as.numeric(dd)
  dd <- dd * 10

  tbl <- table(dd)
  paste0(names(tbl), "s N = ", tbl, collapse = ", ")
}


#' Run GSEA across cell types and GMT files
#'
#'
#' @param de_dt data.table with at least columns gene, cell_type, t.
#' @param fgsea_cell_types Character vector of cell types to run.
#' @param gmt_files Character vector of GMT file paths.
#' @param seed Random seed for reproducibility.
#' @return data.table of fgsea results with added columns cell_type, gmt, and id.
#' @export
run_gsea <- function(de_dt, fgsea_cell_types, gmt_files, seed=42) {
  set.seed(seed)
  gene <- cell_type <- t <- NULL

  gsea_results <- list()
  for (i in fgsea_cell_types) {

    message(i)

    ranks <- de_dt[cell_type == i, t] # using t-statistic to rank genes
    names(ranks) <- de_dt[cell_type == i, gene]
    ranks <- sort(ranks)

    for (j in gmt_files) {

      message(j)

      pathways <- fgsea::gmtPathways(j) # EDIT: would be faster to load all of these once; but fgsea is rate limiting

      new_results <- fgsea::fgsea(pathways = pathways,
                           stats = ranks,
                           minSize = 15,
                           maxSize = 500)

      # colnames(new_results) = tolower(colnames(new_results)) # leaving these alone for now; would be nice if they were snake_case like everything else

      new_results$cell_type <- i
      new_results$gmt <- basename(j)

      data.table::setorderv(new_results, "pval")

      gsea_results[[length(gsea_results) + 1]] <- new_results

    }

  }

  gsea_results <- data.table::rbindlist(gsea_results)

  # Make R CMD CHECK Happy
  id <- NULL
  gsea_results[, id := 1:nrow(gsea_results)]
}

#' Write lightweight GSEA outputs per cell type
#' (originally write_gsea_lite)
#'
#' Writes two TSV files per cell type:
#' - "<cell_type>_pos_gsea.txt": positive NES, padj < threshold
#' - "<cell_type>_neg_gsea.txt": negative NES, padj < threshold
#'
#' Currently filters to GMTs whose basename contains "c5.go".
#'
#' @param gsea_results data.table returned by run_gsea.
#' @param out_dir Output directory.
#' @param padj_thresh Adjusted p-value threshold.
#' @param gmt_pattern Fixed-pattern string applied to the `gmt` column.
#' @export
write_gsea_lite <- function(gsea_results,
                            out_dir,
                            padj_thresh = 0.05,
                            gmt_pattern = "c5.go") {

  cell_type <- gmt <- NES <- padj <- pathway <- size <- NULL

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed.")
  }

  if (!data.table::is.data.table(gsea_results)) {
    gsea_results <- data.table::as.data.table(gsea_results)
  }

  req_cols <- c("cell_type", "gmt", "NES", "padj", "pathway", "size")
  missing_cols <- setdiff(req_cols, names(gsea_results))
  if (length(missing_cols) > 0L) {
    stop("gsea_results is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!dir.exists(out_dir)) {
    stop("out_dir does not exist: ", out_dir)
  }

  if (!is.numeric(padj_thresh) || length(padj_thresh) != 1L || is.na(padj_thresh) || padj_thresh <= 0 || padj_thresh > 1) {
    stop("padj_thresh must be a scalar in (0, 1].")
  }

  for (ct in sort(unique(gsea_results$cell_type))) {

    message(ct)

    pos <- gsea_results[
      cell_type == ct &
        grepl(gmt_pattern, gmt, fixed = TRUE) &
        NES > 0 &
        padj < padj_thresh
    ]

    pos <- pos[, list(pathway, p_adj = padj, nes = NES, size, gmt)]

    utils::write.table(
      pos,
      file = paste0(out_dir, "/", ct, "_pos_gsea.txt"),
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE,
      sep = "\t"
    )

    neg <- gsea_results[
      cell_type == ct &
        grepl(gmt_pattern, gmt, fixed = TRUE) &
        NES < 0 &
        padj < padj_thresh
    ]

    neg <- neg[, list(pathway, p_adj = padj, nes = NES, size, gmt)]

    utils::write.table(
      neg,
      file = paste0(out_dir, "/", ct, "_neg_gsea.txt"),
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE,
      sep = "\t"
    )
  }

  invisible(NULL)
}

#' Donor-level GEX heatmap ordered by donor age
#' (originally plot_donor_gex_age_heatmap)
#'
#' Plots a heatmap of donor-level expression (TPM) for a gene set, with donors
#' implicitly ordered by column order in `metacells` and decade boundaries
#' indicated by gaps.
#'
#' The expression is scaled per gene using the 10th and 90th percentiles across donors:
#' values are shifted by p10, divided by (p90 - p10), and clipped to [0, 1].
#'
#' @param metacells Numeric matrix (or Matrix) with rows = genes and columns = donors.
#'   Column names must be donor IDs that match names(donor_ages).
#' @param gs Character vector of gene symbols to plot.
#' @param donor_ages Named numeric vector of donor ages (names are donor IDs).
#' @param gs_gaps Optional integer vector giving gap positions between gene groups.
#' @param cluster_gs Logical; if TRUE, cluster genes using Spearman correlation.
#' @param transpose Logical; if TRUE, plot donors as rows and genes as columns.
#'
#' @return Invisibly returns NULL.
#' @export
plot_donor_gex_age_heatmap <- function(metacells,
                                       gs,
                                       donor_ages,
                                       gs_gaps = NULL,
                                       cluster_gs = FALSE,
                                       transpose = FALSE) {

  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' is required but not installed.")
  }
  if (!requireNamespace("viridisLite", quietly = TRUE)) {
    stop("Package 'viridisLite' is required but not installed.")
  }

  if (is.null(rownames(metacells))) {
    stop("metacells must have non-NULL rownames (gene identifiers).")
  }
  if (is.null(colnames(metacells))) {
    stop("metacells must have non-NULL colnames (donor identifiers).")
  }
  if (is.null(names(donor_ages))) {
    stop("donor_ages must be a named numeric vector (names are donor IDs).")
  }

  missing_genes <- gs[!(gs %in% rownames(metacells))]
  if (length(missing_genes) > 0L) {
    message("WARNING: missing ", paste(missing_genes, collapse = ", "))
  }

  gs <- gs[gs %in% rownames(metacells)]
  if (length(gs) == 0L) {
    stop("No genes from `gs` were found in rownames(metacells).")
  }

  metacells <- metacells[gs, , drop = FALSE]

  donors <- colnames(metacells)
  if (!all(donors %in% names(donor_ages))) {
    bad <- donors[!(donors %in% names(donor_ages))]
    stop("Some metacells colnames are missing from donor_ages. Example: ", bad[[1]])
  }

  ## Scale each gene: (x - p10) / (p90 - p10), clamp to [0,1]
  ## Use apply to preserve original behavior exactly.
  p_10 <- apply(metacells, 1, stats::quantile, probs = 0.10, na.rm = TRUE)
  p_90 <- apply(metacells, 1, stats::quantile, probs = 0.90, na.rm = TRUE)

  range <- p_90 - p_10

  ## Avoid divide-by-zero: if range==0, set to 1 so the row becomes all 0 after subtraction.
  range[range == 0] <- 1

  metacells <- sweep(metacells, 1, p_10, "-")
  metacells <- sweep(metacells, 1, range, "/")

  metacells[metacells < 0] <- 0
  metacells[metacells > 1] <- 1

  ## Decade transition gaps
  donor_decades <- donor_ages[donors] %/% 10
  idx <- which(diff(donor_decades) != 0)

  if (isTRUE(cluster_gs)) {
    ## cor() needs dense; keep conversion local to this branch
    dist_mat <- stats::as.dist(1 - stats::cor(t(as.matrix(metacells)), method = "spearman"))
    cluster_arg <- stats::hclust(dist_mat, method = "average")
  } else {
    cluster_arg <- FALSE
  }

  if (isTRUE(transpose)) {
    pheatmap::pheatmap(
      t(metacells),
      cluster_rows = FALSE,
      cluster_cols = cluster_arg,
      color = viridisLite::viridis(100),
      show_rownames = FALSE,
      gaps_col = gs_gaps,
      gaps_row = idx
    )
  } else {
    pheatmap::pheatmap(
      metacells,
      cluster_rows = cluster_arg,
      cluster_cols = FALSE,
      color = viridisLite::viridis(100),
      show_colnames = FALSE,
      gaps_col = idx,
      gaps_row = gs_gaps
    )
  }

  invisible(NULL)
}


#' Donor-level GEX scatterplot vs age (ggplot2)
#' (originally plot_donor_gex_age_scatterplot)
#'
#' @param exp_vector Named numeric vector of expression values (names = donor IDs).
#' @param donor_ages Named numeric vector of donor ages (names = donor IDs).
#' @param main Character plot title.
#' @param show_spearman Logical; if TRUE, compute Spearman correlation and display it above the panel.
#' @param size If show_spearman is TRUE, dictates the size of the correlation text.
#' @param rho_threshold Numeric threshold for |rho| above which points are black; otherwise light grey.
#'   Default is 0.2.
#' @param y_axis_floor Numeric minimum upper Y-axis limit. If NULL, no Y limits are imposed.
#'   Default is 10.
#' @return A ggplot object.
#' @export
plot_donor_gex_age_scatterplot <- function(exp_vector,
                                           donor_ages,
                                           main = "",
                                           show_spearman = FALSE,
                                           size = 6,
                                           rho_threshold = 0.2,
                                           y_axis_floor = 10) {

  age <- expression <- NULL

  if (is.null(names(exp_vector))) {
    stop("exp_vector must be a named numeric vector (names are donor IDs).")
  }
  if (is.null(names(donor_ages))) {
    stop("donor_ages must be a named numeric vector (names are donor IDs).")
  }

  donors <- names(exp_vector)

  if (!all(donors %in% names(donor_ages))) {
    bad <- donors[!(donors %in% names(donor_ages))]
    stop("Some exp_vector names are missing from donor_ages. Example: ", bad[[1]])
  }

  df <- data.frame(
    age = as.numeric(donor_ages[donors]),
    expression = as.numeric(exp_vector),
    stringsAsFactors = FALSE
  )

  point_color <- "black"
  rho <- NA_real_
  subtitle_txt <- NULL

  if (isTRUE(show_spearman)) {

    ok <- stats::complete.cases(df$age, df$expression)

    if (sum(ok) >= 3L) {

      ct <- stats::cor.test(
        df$age[ok],
        df$expression[ok],
        method = "spearman",
        exact = FALSE
      )

      rho <- unname(ct$estimate)
      subtitle_txt <- paste0("rho = ", formatC(rho, format = "f", digits = 2))

      if (is.finite(rho) && abs(rho) < rho_threshold) {
        point_color <- "lightgrey"
      }

    } else {
      warning("Not enough non-missing points to compute Spearman correlation.", call. = FALSE)
    }
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = age, y = expression)) +
    ggplot2::geom_point(size = 1.0, color = point_color) +
    ggplot2::labs(
      x = "Age",
      y = "Expression, TPM",
      title = main,
      subtitle = subtitle_txt
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.subtitle = ggplot2::element_text(size = size, hjust = 0)
    )

  if (!is.null(y_axis_floor)) {
    y_max <- max(y_axis_floor, max(df$expression, na.rm = TRUE))
    p <- p + ggplot2::scale_y_continuous(limits = c(0, y_max))
  }

  p
}
