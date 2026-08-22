# in_dir1 <- "differential_expression/results/LEVEL_6/age_bmi_complete_cases/age_bmi"
# in_dir2 <- "differential_expression/results/LEVEL_6/age_smoking_complete_cases/age_smoking"
# ct_file <- "differential_expression/metadata/cell_types_for_de_filtering_plot.txt"
#
# de_counts <- compute_de_region_contrast_counts(
#   paths = c(in_dir2, in_dir1, in_dir2),
#   contrasts = c("age", "bmi", "yes_vs_no_smoking"),
#   cellTypeListFile = ct_file
# )
#
# plot_de_region_contrast_heatmap(de_counts, label_threshold = 10000)

#' Compute per-cell-type DE gene counts across many contrasts and regions
#'
#' Assembles a long table of total significant differential expression gene
#' counts per cell type, spanning many \code{(path, contrast)} pairs at once.
#' Each pair is read via \code{\link{compute_de_result_counts}} (which in
#' turn uses \code{\link{parse_de_inputs}}), so a single pair already picks
#' up every region present in that directory for that contrast. The same
#' region can be contributed by multiple different paths/contrasts (e.g. one
#' directory supplying a "bmi" contrast and another supplying "age" and
#' smoking contrasts for the same set of regions) — the results are merged
#' by region, not siloed per source path.
#'
#' \code{region} and \code{contrast} are returned as ordered factors —
#' \code{region} by first appearance across \code{paths}, \code{contrast} by
#' \code{unique(contrasts)} (i.e. \code{contrasts}' own order is the display
#' order; list it in the order you want, no separate ordering parameter
#' needed) — this ordering is what lets
#' \code{\link{plot_de_region_contrast_heatmap}} produce a hierarchical
#' "region, then its contrasts" axis with no further bookkeeping.
#'
#' @param paths Character vector of DE result directories, one per
#'   \code{(path, contrast)} pair. The same directory may repeat multiple
#'   times, once per contrast found in it.
#' @param contrasts Character vector of plain contrast names (e.g.
#'   \code{"age"}, \code{"yes_vs_no_smoking"}), the same length as
#'   \code{paths}. Each is turned into an anchored file-matching pattern
#'   internally (\code{"__<contrast>_DE_results\\.txt$"}), so directories
#'   containing stray non-DE files (e.g. \code{*_volcano_plots.pdf}) are
#'   handled safely. \code{unique(contrasts)} also becomes the \code{contrast}
#'   column's factor level order, so list contrasts in your desired display
#'   order.
#' @param cellTypeListFile Optional path to a file specifying the cell types
#'   to include, applied to every pair. When supplied, its order is also
#'   preserved as the \code{cell_type} column's factor level order (see
#'   \code{\link{plot_de_region_contrast_heatmap}}'s \code{cell_type_order}).
#' @param fdr_cutoff Adjusted P-value threshold used to call a gene
#'   significant. Defaults to 0.05.
#' @param min_num_genes Minimum total DE hits (\code{n_total}, summed across
#'   all regions and contrasts) a cell type must have to be kept. Cell types
#'   below this are dropped entirely. Defaults to 20; set to \code{NULL} to
#'   disable.
#' @param filter_ic_to_non_neurons If \code{TRUE} (default), rows for the
#'   \code{"ic"} region are dropped except for non-neuronal cell types
#'   (\code{astrocyte}, \code{OPC}, \code{oligodendrocyte}, \code{microglia})
#'   — i.e. neuronal cell types are set to missing for \code{ic} before
#'   plotting. Mirrors the same manuscript-specific convention used
#'   elsewhere (e.g. \code{bican.mccarroll.de.analysis::compute_de_cor_mat()}
#'   and the TRADE-across-regions analysis in
#'   \code{bican.mccarroll.figures::trade_plots.R}).
#'
#' @return A data.table with columns \code{region}, \code{contrast}
#'   (ordered factors, see above), \code{cell_type}, \code{n_up},
#'   \code{n_down}, and \code{n_total} (\code{n_up + n_down}).
#'
#' @export
compute_de_region_contrast_counts <- function(paths,
                                              contrasts,
                                              cellTypeListFile = NULL,
                                              fdr_cutoff = 0.05,
                                              min_num_genes = 20,
                                              filter_ic_to_non_neurons = TRUE) {
  # Make R CMD CHECK happy
  interaction <- region <- contrast <- n_up <- n_down <- n_total <- cell_type <- NULL

  if (length(paths) != length(contrasts)) {
    stop("paths and contrasts must be the same length.")
  }

  region_levels <- character(0)
  seen_keys <- character(0)

  results <- vector("list", length(paths))

  for (i in seq_along(paths)) {
    file_pattern <- sprintf("__%s_DE_results\\.txt$", contrasts[i])

    counts_i <- compute_de_result_counts(
      in_dir = paths[i],
      file_pattern = file_pattern,
      cellTypeListFile = cellTypeListFile,
      fdr_cutoff = fdr_cutoff
    )

    data.table::setnames(counts_i, "interaction", "region")

    region_values <- unique(as.character(counts_i$region))
    region_levels <- union(region_levels, region_values)

    key <- paste(region_values, contrasts[i], sep = "__")
    dup_keys <- intersect(key, seen_keys)
    if (length(dup_keys) > 0) {
      warning(
        "Duplicate (region, contrast) combination(s) from paths/contrasts input, keeping first occurrence: ",
        paste(dup_keys, collapse = ", ")
      )
      counts_i <- counts_i[!(paste(region, contrasts[i], sep = "__") %in% dup_keys)]
    }
    seen_keys <- c(seen_keys, setdiff(key, dup_keys))

    results[[i]] <- counts_i
  }

  de_counts <- data.table::rbindlist(results, use.names = TRUE)
  de_counts[, n_total := n_up + n_down]

  contrast_levels <- unique(contrasts)

  de_counts[, region := factor(region, levels = region_levels)]
  de_counts[, contrast := factor(contrast, levels = contrast_levels)]

  if (filter_ic_to_non_neurons) {
    non_neuron_types <- c("astrocyte", "OPC", "oligodendrocyte", "microglia")
    de_counts <- de_counts[!(region == "ic" & !(cell_type %in% non_neuron_types))]
  }

  if (!is.null(min_num_genes)) {
    keep_cell_types <- de_counts[, .(total = sum(n_total)), by = cell_type][total >= min_num_genes, cell_type]
    de_counts <- de_counts[cell_type %in% keep_cell_types]
  }

  # If a cell type list was supplied, preserve its order as the cell_type
  # factor levels, so plot_de_region_contrast_heatmap() can default the Y
  # axis to that order (the manuscript-standard cell type ordering) instead
  # of alphabetical.
  if (!is.null(cellTypeListFile)) {
    ct_order <- scan(cellTypeListFile, what = character(), quiet = TRUE)
    present <- unique(as.character(de_counts$cell_type))
    de_counts[, cell_type := factor(cell_type, levels = intersect(ct_order, present))]
  }

  data.table::setorder(de_counts, region, contrast, cell_type)

  de_counts[]
}

#' Plot a cell type x region x contrast DE gene-count heatmap
#'
#' Exploratory heatmap/table view of the counts produced by
#' \code{\link{compute_de_region_contrast_counts}}: cell types on the Y
#' axis, and a hierarchical X axis grouped first by region, then by contrast
#' within that region (via \code{region}/\code{contrast} factor level
#' order). Returns a \code{ggplot} object only — no file output, not wired
#' into any figure pipeline.
#'
#' @param de_counts A data.table produced by
#'   \code{\link{compute_de_region_contrast_counts}}.
#' @param metric Which column to use as the tile fill value. One of
#'   \code{"n_total"} (default), \code{"n_up"}, or \code{"n_down"}.
#' @param cell_type_order Optional character vector fixing the Y-axis row
#'   order. If \code{NULL} (default), uses \code{de_counts$cell_type}'s
#'   factor level order when it is already a factor (i.e. when
#'   \code{cellTypeListFile} was passed to
#'   \code{\link{compute_de_region_contrast_counts}}), otherwise alphabetical.
#' @param fill_transform Transform applied to the fill color scale. Defaults
#'   to \code{"log1p"} (handles zero counts), which spreads out the
#'   near-zero-to-low-hundreds range that a linear scale washes out when some
#'   contrasts have thousands of DE genes and others have almost none. Pass
#'   \code{"identity"} for a plain linear scale.
#' @param fill_breaks Numeric vector of legend break points for the fill
#'   scale. Defaults to \code{c(0, 10, 100, 1000, 10000)}, giving a clean
#'   log10-style legend (values outside the data range are simply not shown).
#' @param label_threshold Tiles whose \code{metric} value is strictly less
#'   than this are annotated with the count as text (useful for
#'   distinguishing small counts that a compressed color scale can't show
#'   precisely). Defaults to 100; set to a high value (e.g. 10000) to label
#'   every tile. Set to \code{NULL} to disable labels entirely.
#' @param label_size Font size (in mm, per \code{ggplot2::geom_text}) for the
#'   cell labels. Defaults to 2.2, small enough to fit a 4-digit count.
#' @param label_color_threshold Tiles whose \code{metric} value is at least
#'   this are labeled in black instead of white, so labels stay legible over
#'   both the dark (low-count) and bright (high-count) ends of the viridis
#'   fill scale. Defaults to 1000.
#' @param legend_title Title for the fill legend. Defaults to a friendly
#'   label based on \code{metric} (\code{"DE genes"} for \code{"n_total"},
#'   \code{"Upregulated DE genes"} for \code{"n_up"}, \code{"Downregulated DE
#'   genes"} for \code{"n_down"}).
#' @param contrast_labels Optional named character vector remapping
#'   \code{contrast} values to display labels on the X axis, e.g.
#'   \code{c(age = "Age", bmi = "BMI", yes_vs_no_smoking = "Smoking")}.
#'   Contrasts not named in the vector keep their original value. Defaults
#'   to \code{NULL} (show the raw contrast values).
#' @param no_color If \code{TRUE}, skip color entirely: every tile is white
#'   with a black border, and every cell is labeled in black (regardless of
#'   \code{label_threshold}/\code{label_color_threshold}, which are ignored
#'   in this mode) — the numbers are the only way to read the values. No
#'   fill legend is drawn. Defaults to \code{FALSE}.
#'
#' @return A \code{ggplot} object.
#'
#' @export
plot_de_region_contrast_heatmap <- function(de_counts,
                                            metric = c("n_total", "n_up", "n_down"),
                                            cell_type_order = NULL,
                                            fill_transform = "log1p",
                                            fill_breaks = c(0, 10, 100, 1000, 10000),
                                            label_threshold = 100,
                                            label_size = 2.2,
                                            label_color_threshold = 1000,
                                            legend_title = NULL,
                                            contrast_labels = NULL,
                                            no_color = FALSE) {
  metric <- match.arg(metric)

  if (is.null(legend_title)) {
    legend_title <- c(
      n_total = "DE genes",
      n_up = "Upregulated DE genes",
      n_down = "Downregulated DE genes"
    )[[metric]]
  }

  # Make R CMD CHECK happy
  cell_type <- contrast <- region <- label_color <- NULL

  d <- data.table::copy(de_counts)

  if (is.null(cell_type_order)) {
    if (is.factor(d$cell_type)) {
      # compute_de_region_contrast_counts() sets this when cellTypeListFile
      # was supplied, preserving its order (e.g. the manuscript-standard
      # cell type ordering) rather than falling back to alphabetical.
      cell_type_order <- levels(d$cell_type)
    } else {
      cell_type_order <- sort(unique(as.character(d$cell_type)))
    }
  }
  d[, cell_type := factor(cell_type, levels = rev(cell_type_order))]

  if (no_color) {
    p <- ggplot2::ggplot(d, ggplot2::aes(x = contrast, y = cell_type)) +
      ggplot2::geom_tile(fill = "white", color = "black") +
      ggplot2::facet_grid(~region, scales = "free_x", space = "free_x") +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        strip.background = ggplot2::element_rect(fill = "white", color = "black")
      )

    if (!is.null(contrast_labels)) {
      p <- p + ggplot2::scale_x_discrete(
        labels = function(x) ifelse(x %in% names(contrast_labels), contrast_labels[x], x)
      )
    }

    if (!is.null(label_threshold)) {
      p <- p + ggplot2::geom_text(
        data = d,
        ggplot2::aes(label = .data[[metric]]),
        color = "black",
        size = label_size
      )
    }

    return(p)
  }

  p <- ggplot2::ggplot(d, ggplot2::aes(x = contrast, y = cell_type, fill = .data[[metric]])) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::facet_grid(~region, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_viridis_c(name = legend_title, transform = fill_transform, breaks = fill_breaks) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.background = ggplot2::element_rect(fill = "white", color = "black")
    )

  if (!is.null(contrast_labels)) {
    p <- p + ggplot2::scale_x_discrete(
      labels = function(x) ifelse(x %in% names(contrast_labels), contrast_labels[x], x)
    )
  }

  if (!is.null(label_threshold)) {
    d_label <- d[d[[metric]] < label_threshold]
    d_label$label_color <- ifelse(d_label[[metric]] >= label_color_threshold, "black", "white")

    p <- p + ggplot2::geom_text(
      data = d_label,
      ggplot2::aes(label = .data[[metric]], color = label_color),
      size = label_size
    ) +
      ggplot2::scale_color_identity()
  }

  p
}
