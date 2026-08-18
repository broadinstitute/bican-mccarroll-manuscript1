## ------------------------------------------------------------------
## Example usage
## ------------------------------------------------------------------

# plot_de_primary_secondary_manifest(
#     manifest_file = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/external_comparison_snap200/metadata/bican_dfc_snap_age_de_overlap_manifest.tsv",
#     out_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/external_comparison_snap200/results",
#     primary_dataset = "bican",
#     secondary_dataset = "snap",
#     primary_label = "BICAN",
#     secondary_label = "SNAP",
#     effect_name = "age effects"
# )

# plot_de_primary_secondary_manifest(
#     manifest_file = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/external_comparison_PMID_40903571/metadata/PMID_40903571_bican_dfc_age_de_overlap_manifest.tsv",
#     out_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/external_comparison_PMID_40903571/results",
#     primary_dataset = "PMID_40903571",
#     secondary_dataset = "bican",
#     primary_label = "PMID 40903571",
#     secondary_label = "BICAN",
#     effect_name = "age effects"
# )

# plot_de_primary_secondary_manifest(
#     manifest_file = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/external_comparison_PMID_39227716/metadata/PMID_39227716_bican_dfc_age_de_overlap_manifest.tsv",
#     out_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/external_comparison_PMID_39227716/results",
#     primary_dataset = "bican",
#     secondary_dataset = "PMID_39227716",
#     primary_label = "BICAN",
#     secondary_label = "PMID 39227716",
#     effect_name = "age effects",
#     min_num_genes=1 # to include the one type that is filtered.
# )

# plot_de_primary_secondary_manifest(
#     manifest_file = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/external_comparison_snap200/metadata/bican_dfc_snap_sex_de_overlap_manifest.tsv",
#     out_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/external_comparison_snap200/results",
#     primary_dataset = "bican",
#     secondary_dataset = "snap",
#     primary_label = "BICAN",
#     secondary_label = "Ling et al.",
#     effect_name = "sex effects",
#     gene_set = "autosome",
#     contig_yaml_file = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/metadata/GRCh38_ensembl_v43.contig_groups.yaml",
#     reduced_gtf_file = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/metadata/GRCh38_ensembl_v43.reduced.gtf.gz"
# )

# plot_de_primary_secondary_manifest(
#     manifest_file = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/external_comparison_snap200/metadata/bican_dfc_snap_sex_de_overlap_manifest.tsv",
#     out_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/external_comparison_snap200/results",
#     primary_dataset = "bican",
#     secondary_dataset = "snap",
#     primary_label = "BICAN",
#     secondary_label = "Ling et al.",
#     effect_name = "sex effects",
#     gene_set = "both",
#     contig_yaml_file = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/metadata/GRCh38_ensembl_v43.contig_groups.yaml",
#     reduced_gtf_file = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/metadata/GRCh38_ensembl_v43.reduced.gtf.gz"
# )

# plot_de_primary_secondary_manifest(
#     manifest_file = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/external_comparison_snap200/metadata/bican_dfc_snap_sex_de_overlap_manifest.tsv",
#     out_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/external_comparison_snap200/results",
#     primary_dataset = "bican",
#     secondary_dataset = "snap",
#     primary_label = "BICAN",
#     secondary_label = "Ling et al.",
#     effect_name = "sex effects",
#     gene_set = "xy",
#     contig_yaml_file = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/metadata/GRCh38_ensembl_v43.contig_groups.yaml",
#     reduced_gtf_file = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/metadata/GRCh38_ensembl_v43.reduced.gtf.gz"
# )

## ------------------------------------------------------------------
## Top-level functions
## ------------------------------------------------------------------

#' Compare a primary dataset's significant DE genes against a secondary dataset
#'
#' For each cell type listed in \code{manifest_file}, reads a primary and
#' secondary differential expression result file, restricts to genes that are
#' significant in the primary dataset (optionally further restricted to
#' \code{genes_to_keep}), and plots secondary logFC against primary logFC with
#' a sign-concordance test. Produces one scatter plot per cell type plus a
#' summary bar plot of sign concordance across cell types.
#'
#' @param manifest_file Path to a tab-delimited manifest with at least the
#'   columns \code{cell_type}, \code{primary_dataset}, and
#'   \code{secondary_dataset}, where the latter two name columns containing
#'   per-cell-type DE result file paths.
#' @param out_dir Directory to write per-cell-type scatter SVGs, the summary
#'   SVG, and the summary TSV. If \code{NULL}, no files are written and only
#'   the returned list is available.
#' @param primary_dataset,secondary_dataset Column names in \code{manifest_file}
#'   identifying the primary and secondary dataset file-path columns.
#' @param primary_label,secondary_label Display labels for the two datasets.
#'   If \code{NULL}, derived from \code{primary_dataset}/\code{secondary_dataset}.
#' @param effect_name Human-readable name of the effect being compared (e.g.
#'   \code{"age effects"}, \code{"sex effects"}), used in plot titles and
#'   output file names.
#' @param gene_set One of \code{"both"}, \code{"autosome"}, or \code{"xy"}.
#'   Restricts comparisons to genes on autosomes or sex chromosomes. Requires
#'   \code{contig_yaml_file} and \code{reduced_gtf_file} unless \code{"both"}.
#' @param contig_yaml_file Path to the contig-groups YAML file used to resolve
#'   \code{gene_set}. Required (together with \code{reduced_gtf_file}) unless
#'   \code{gene_set = "both"}.
#' @param reduced_gtf_file Path to a tab-delimited reduced GTF-like file with
#'   at least the columns \code{chr}, \code{annotationType}, and
#'   \code{gene_name}. Required (together with \code{contig_yaml_file}) unless
#'   \code{gene_set = "both"}.
#' @param min_num_genes Minimum number of overlapping genes required for a
#'   cell type to be included in the sign-concordance summary.
#' @param alpha Adjusted P-value threshold used to select primary-significant
#'   genes.
#' @param width,height Dimensions (inches) for the per-cell-type scatter SVGs.
#' @param scale_effects Logical scalar. If \code{TRUE}, each cell type's
#'   primary and secondary logFC values are independently rescaled to
#'   \code{[-1, 1]} before plotting (see
#'   \code{\link{plot_de_primary_secondary_from_data}}). Defaults to
#'   \code{FALSE}.
#'
#' @return Invisibly returns a list with elements \code{per_cell_type} (a list
#'   of per-cell-type results), \code{sign_results} (a data frame of
#'   sign-concordance statistics per cell type), \code{sign_summary_plot}
#'   (the sign-concordance summary bar plot), and \code{cor_summary_plot}
#'   (the analogous summary bar plot of per-cell-type correlation). Both
#'   summary plots are \code{NULL} if no cell type passed \code{min_num_genes}.
#'
#' @export
plot_de_primary_secondary_manifest <- function(manifest_file,
                                               out_dir = NULL,
                                               primary_dataset = "bican",
                                               secondary_dataset = "snap",
                                               primary_label = NULL,
                                               secondary_label = NULL,
                                               effect_name = "age effects",
                                               gene_set = c("both", "autosome", "xy"),
                                               contig_yaml_file = NULL,
                                               reduced_gtf_file = NULL,
                                               min_num_genes = 20,
                                               alpha = 0.05,
                                               width = 7,
                                               height = 7,
                                               scale_effects = FALSE) {
  gene_set <- match.arg(gene_set)

  de_data <- compute_de_primary_secondary_data(
    manifest_file = manifest_file,
    primary_dataset = primary_dataset,
    secondary_dataset = secondary_dataset,
    gene_set = gene_set,
    contig_yaml_file = contig_yaml_file,
    reduced_gtf_file = reduced_gtf_file,
    alpha = alpha
  )

  plot_de_primary_secondary_from_data(
    de_data,
    out_dir = out_dir,
    primary_dataset = primary_dataset,
    secondary_dataset = secondary_dataset,
    primary_label = primary_label,
    secondary_label = secondary_label,
    effect_name = effect_name,
    gene_set = gene_set,
    min_num_genes = min_num_genes,
    width = width,
    height = height,
    scale_effects = scale_effects
  )
}


#' Gather per-gene primary/secondary DE comparison data for a manifest
#'
#' For each cell type listed in \code{manifest_file}, reads the primary and
#' secondary differential expression result files and joins them on shared,
#' primary-significant genes. This is the data-gathering half of
#' \code{plot_de_primary_secondary_manifest()}, split out so callers can
#' cache the joined gene-level data (as one table spanning all cell types)
#' and skip re-parsing the raw DE result files on subsequent plot renders.
#'
#' @inheritParams plot_de_primary_secondary_manifest
#'
#' @return A data.frame with columns \code{cell_type}, \code{gene},
#'   \code{secondary_logFC}, \code{secondary_adjP}, \code{primary_logFC}, and
#'   \code{primary_adjP} (one row per gene per cell type). Cell types with no
#'   qualifying genes are omitted.
#'
#' @export
compute_de_primary_secondary_data <- function(manifest_file,
                                              primary_dataset = "bican",
                                              secondary_dataset = "snap",
                                              gene_set = c("both", "autosome", "xy"),
                                              contig_yaml_file = NULL,
                                              reduced_gtf_file = NULL,
                                              alpha = 0.05) {
  gene_set <- match.arg(gene_set)

  validate_gene_set_metadata(
    gene_set = gene_set,
    contig_yaml_file = contig_yaml_file,
    reduced_gtf_file = reduced_gtf_file
  )

  genes_to_keep <- get_de_genes_for_gene_set(
    gene_set = gene_set,
    contig_yaml_file = contig_yaml_file,
    reduced_gtf_file = reduced_gtf_file
  )

  manifest <- utils::read.table(
    manifest_file,
    header = TRUE,
    stringsAsFactors = FALSE,
    sep = "\t",
    check.names = FALSE
  )

  required_cols <- c("cell_type", primary_dataset, secondary_dataset)
  missing_cols <- setdiff(required_cols, colnames(manifest))

  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Manifest is missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  validate_manifest_files(
    manifest = manifest,
    file_cols = c(primary_dataset, secondary_dataset)
  )

  per_cell_type <- lapply(seq_len(nrow(manifest)), function(i) {
    cell_type <- manifest$cell_type[i]

    logger::log_info("Processing cell type ", cell_type)

    df <- make_de_primary_secondary_df(
      primary_file = manifest[[primary_dataset]][i],
      secondary_file = manifest[[secondary_dataset]][i],
      genes_to_keep = genes_to_keep,
      alpha = alpha
    )

    if (nrow(df) == 0) {
      return(NULL)
    }

    df$cell_type <- cell_type
    df
  })

  de_data <- do.call(rbind, per_cell_type)

  if (is.null(de_data)) {
    de_data <- data.frame(
      cell_type = character(0),
      gene = character(0),
      secondary_logFC = numeric(0),
      secondary_adjP = numeric(0),
      primary_logFC = numeric(0),
      primary_adjP = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  de_data[, c(
    "cell_type", "gene", "secondary_logFC", "secondary_adjP",
    "primary_logFC", "primary_adjP"
  )]
}


#' Plot a primary/secondary DE comparison from precomputed per-gene data
#'
#' Builds per-cell-type scatter plots and a sign-concordance summary from
#' gene-level data produced by \code{compute_de_primary_secondary_data()}.
#' This is the plotting half of \code{plot_de_primary_secondary_manifest()}.
#'
#' @param de_data A data.frame produced by
#'   \code{compute_de_primary_secondary_data()}, with columns
#'   \code{cell_type}, \code{gene}, \code{secondary_logFC},
#'   \code{secondary_adjP}, \code{primary_logFC}, and \code{primary_adjP}.
#' @param out_dir Directory to write per-cell-type scatter SVGs, the summary
#'   SVG, and the summary TSV. If \code{NULL}, no files are written and only
#'   the returned list is available.
#' @param primary_dataset,secondary_dataset Dataset identifiers recorded in
#'   the summary TSV and used to derive default labels and output filenames;
#'   the join itself was already performed when \code{de_data} was gathered.
#' @param primary_label,secondary_label Display labels for the two datasets.
#'   If \code{NULL}, derived from \code{primary_dataset}/\code{secondary_dataset}.
#' @param effect_name Human-readable name of the effect being compared (e.g.
#'   \code{"age effects"}, \code{"sex effects"}), used in plot titles and
#'   output file names.
#' @param gene_set Recorded in the summary TSV and output filenames
#'   (informational only at this stage; restricting to autosomal/sex-linked
#'   genes already happened when \code{de_data} was gathered).
#' @param min_num_genes Minimum number of overlapping genes required for a
#'   cell type to be included in the sign-concordance summary.
#' @param width,height Dimensions (inches) for the per-cell-type scatter SVGs.
#' @param scale_effects Logical scalar. If \code{TRUE}, each cell type's
#'   primary and secondary logFC values are independently rescaled to
#'   \code{[-1, 1]} (dividing by their own maximum absolute value) before
#'   plotting, so effect sizes are visually comparable across cell types and
#'   datasets on very different native scales. Does not affect the
#'   sign-concordance statistics, which are scale-invariant. Defaults to
#'   \code{FALSE}.
#'
#' @return Invisibly returns a list with elements \code{per_cell_type} (a list
#'   of per-cell-type results), \code{sign_results} (a data frame of
#'   sign-concordance statistics per cell type), \code{sign_summary_plot}
#'   (the sign-concordance summary bar plot), and \code{cor_summary_plot}
#'   (the analogous summary bar plot of per-cell-type correlation). Both
#'   summary plots are \code{NULL} if no cell type passed \code{min_num_genes}.
#'
#' @export
plot_de_primary_secondary_from_data <- function(de_data,
                                                out_dir = NULL,
                                                primary_dataset = "bican",
                                                secondary_dataset = "snap",
                                                primary_label = NULL,
                                                secondary_label = NULL,
                                                effect_name = "age effects",
                                                gene_set = "both",
                                                min_num_genes = 20,
                                                width = 7,
                                                height = 7,
                                                scale_effects = FALSE) {
  if (is.null(primary_label)) {
    primary_label <- make_dataset_label(primary_dataset)
  }

  if (is.null(secondary_label)) {
    secondary_label <- make_dataset_label(secondary_dataset)
  }

  de_data <- as.data.frame(de_data)

  cell_types <- unique(de_data$cell_type)

  results <- vector("list", length(cell_types))
  names(results) <- cell_types

  sign_results <- vector("list", length(cell_types))
  names(sign_results) <- cell_types

  if (!is.null(out_dir) && !dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  output_prefix <- make_output_prefix(
    primary_dataset,
    secondary_dataset,
    gene_set = gene_set
  )

  for (cell_type in cell_types) {
    df <- de_data[de_data$cell_type == cell_type, , drop = FALSE]

    res <- plot_de_primary_secondary_from_df(
      cell_type = cell_type,
      df = df,
      primary_dataset = primary_dataset,
      secondary_dataset = secondary_dataset,
      primary_label = primary_label,
      secondary_label = secondary_label,
      effect_name = effect_name,
      gene_set = gene_set,
      scale_effects = scale_effects
    )

    results[[cell_type]] <- res
    sign_results[[cell_type]] <- res$sign_test

    if (!is.null(out_dir)) {
      out_file <- file.path(
        out_dir,
        sprintf(
          "%s_%s_scatter_%s.svg",
          output_prefix,
          .sanitize_filename(effect_name),
          .sanitize_filename(cell_type)
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
  sign_results <- sign_results[sign_results$n >= min_num_genes, ]

  if (is.null(sign_results) || nrow(sign_results) == 0) {
    warning("No cell types had overlapping primary-significant genes.")
    return(invisible(list(
      per_cell_type = results,
      sign_results = sign_results,
      sign_summary_plot = NULL,
      cor_summary_plot = NULL
    )))
  }

  sign_summary_plot <- plot_sign_test_summary(
    sign_results = sign_results,
    primary_label = primary_label,
    secondary_label = secondary_label,
    effect_name = effect_name
  )

  cor_summary_plot <- plot_correlation_summary(
    sign_results = sign_results,
    primary_label = primary_label,
    secondary_label = secondary_label,
    effect_name = effect_name
  )

  if (!is.null(out_dir)) {
    ggplot2::ggsave(
      filename = file.path(
        out_dir,
        sprintf(
          "%s_%s_sign_test_summary.svg",
          output_prefix,
          .sanitize_filename(effect_name)
        )
      ),
      plot = sign_summary_plot,
      device = "svg",
      width = 8,
      height = max(4, 0.35 * nrow(sign_results) + 2)
    )

    ggplot2::ggsave(
      filename = file.path(
        out_dir,
        sprintf(
          "%s_%s_correlation_summary.svg",
          output_prefix,
          .sanitize_filename(effect_name)
        )
      ),
      plot = cor_summary_plot,
      device = "svg",
      width = 8,
      height = max(4, 0.35 * nrow(sign_results) + 2)
    )

    utils::write.table(
      sign_results,
      file = file.path(
        out_dir,
        sprintf(
          "%s_%s_sign_test_summary.tsv",
          output_prefix,
          .sanitize_filename(effect_name)
        )
      ),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  invisible(list(
    per_cell_type = results,
    sign_results = sign_results,
    sign_summary_plot = sign_summary_plot,
    cor_summary_plot = cor_summary_plot
  ))
}


plot_de_primary_secondary_cell_type <- function(cell_type,
                                                primary_file,
                                                secondary_file,
                                                primary_dataset = "bican",
                                                secondary_dataset = "snap",
                                                primary_label = NULL,
                                                secondary_label = NULL,
                                                effect_name = "age effects",
                                                gene_set = "both",
                                                genes_to_keep = NULL,
                                                alpha = 0.05,
                                                scale_effects = FALSE) {
  df <- make_de_primary_secondary_df(
    primary_file = primary_file,
    secondary_file = secondary_file,
    genes_to_keep = genes_to_keep,
    alpha = alpha
  )

  if (nrow(df) == 0) {
    return(NULL)
  }

  plot_de_primary_secondary_from_df(
    cell_type = cell_type,
    df = df,
    primary_dataset = primary_dataset,
    secondary_dataset = secondary_dataset,
    primary_label = primary_label,
    secondary_label = secondary_label,
    effect_name = effect_name,
    gene_set = gene_set,
    scale_effects = scale_effects
  )
}


plot_de_primary_secondary_from_df <- function(cell_type,
                                              df,
                                              primary_dataset = "bican",
                                              secondary_dataset = "snap",
                                              primary_label = NULL,
                                              secondary_label = NULL,
                                              effect_name = "age effects",
                                              gene_set = "both",
                                              scale_effects = FALSE) {
  # Make R CMD CHECK Happy
  secondary_logFC <- primary_logFC <- NULL

  if (is.null(primary_label)) {
    primary_label <- make_dataset_label(primary_dataset)
  }

  if (is.null(secondary_label)) {
    secondary_label <- make_dataset_label(secondary_dataset)
  }

  # Sign concordance and correlation are scale-invariant (dividing by a
  # positive constant changes neither sign nor Pearson correlation), so these
  # are always computed on the raw, unscaled effect sizes.
  sign_res <- sign_test_fraction(
    x = df$secondary_logFC,
    y = df$primary_logFC
  )

  cor_val <- stats::cor(
    df$secondary_logFC,
    df$primary_logFC,
    use = "complete.obs"
  )

  plot_df <- df
  axis_suffix <- ""

  if (scale_effects) {
    scale_to_unit <- function(x) {
      m <- max(abs(x), na.rm = TRUE)
      if (!is.finite(m) || m == 0) m <- 1
      x / m
    }

    plot_df$secondary_logFC <- scale_to_unit(plot_df$secondary_logFC)
    plot_df$primary_logFC <- scale_to_unit(plot_df$primary_logFC)
    axis_suffix <- " (scaled to [-1, 1])"
  }

  lims <- if (scale_effects) {
    c(-1, 1)
  } else {
    range(c(plot_df$secondary_logFC, plot_df$primary_logFC), na.rm = TRUE)
  }

  subtitle <- sprintf(
    "%s-significant genes: %.1f%% same sign in %s (%d/%d); cor=%.2f",
    primary_label,
    100 * sign_res$frac,
    secondary_label,
    sign_res$k,
    sign_res$n,
    cor_val
  )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = secondary_logFC, y = primary_logFC)) +
    ggplot2::geom_point(na.rm = TRUE, size = 1.2, alpha = 0.6) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "black") +
    ggplot2::coord_cartesian(xlim = lims, ylim = lims) +
    ggplot2::labs(
      title = cell_type,
      subtitle = subtitle,
      x = sprintf("%s logFC%s", secondary_label, axis_suffix),
      y = sprintf("%s logFC%s", primary_label, axis_suffix)
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 13, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 11),
      axis.title = ggplot2::element_text(size = 13),
      axis.text = ggplot2::element_text(size = 11)
    )

  list(
    plot = p,
    sign_test = data.frame(
      cell_type = cell_type,
      primary_dataset = primary_dataset,
      secondary_dataset = secondary_dataset,
      primary_label = primary_label,
      secondary_label = secondary_label,
      effect_name = effect_name,
      gene_set = gene_set,
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

plot_sign_test_summary <- function(sign_results,
                                   primary_label = "BICAN",
                                   secondary_label = "SNAP",
                                   effect_name = "age effects") {
  # Make R CMD CHECK Happy
  cell_type <- percent <- k <- n <- NULL

  sign_results <- sign_results[order(sign_results$percent), ]

  sign_results$cell_type <- factor(
    sign_results$cell_type,
    levels = sign_results$cell_type
  )

  title <- sprintf(
    "%s-significant %s: sign concordance in %s",
    primary_label,
    effect_name,
    secondary_label
  )

  ggplot2::ggplot(sign_results, ggplot2::aes(x = cell_type, y = percent)) +
    ggplot2::geom_col(color = "black") +
    ggplot2::geom_text(
      ggplot2::aes(label = paste0(k, "/", n)),
      hjust = -0.1,
      size = 4
    ) +
    ggplot2::coord_flip(ylim = c(0, 105)) +
    ggplot2::labs(
      x = NULL,
      y = "% same sign",
      title = title
    ) +
    ggplot2::theme_bw(base_size = 12)
}


plot_correlation_summary <- function(sign_results,
                                     primary_label = "BICAN",
                                     secondary_label = "SNAP",
                                     effect_name = "age effects") {
  # Make R CMD CHECK Happy
  cell_type <- cor <- NULL

  sign_results <- sign_results[order(sign_results$cor), ]

  sign_results$cell_type <- factor(
    sign_results$cell_type,
    levels = sign_results$cell_type
  )

  title <- sprintf(
    "%s-significant %s: correlation with %s",
    primary_label,
    effect_name,
    secondary_label
  )

  ggplot2::ggplot(sign_results, ggplot2::aes(x = cell_type, y = cor)) +
    ggplot2::geom_col(color = "black") +
    ggplot2::geom_text(
      ggplot2::aes(
        label = sprintf("%.2f", cor),
        hjust = ifelse(cor >= 0, -0.1, 1.1)
      ),
      size = 4
    ) +
    ggplot2::coord_flip(ylim = c(-1.15, 1.15)) +
    ggplot2::labs(
      x = NULL,
      y = "Correlation (r)",
      title = title
    ) +
    ggplot2::theme_bw(base_size = 12)
}


## ------------------------------------------------------------------
## Data helpers
## ------------------------------------------------------------------

make_de_primary_secondary_df <- function(primary_file,
                                         secondary_file,
                                         genes_to_keep = NULL,
                                         alpha = 0.05) {
  primary <- read_de_result(primary_file)
  secondary <- read_de_result(secondary_file)

  logger::log_info("Primary result number of signifcant genes [", length(which(primary$adj.P.Val <= 0.05)), "]")

  common_genes <- intersect(rownames(primary), rownames(secondary))

  genes <- rownames(primary)[primary$adj.P.Val < alpha]
  genes <- intersect(genes, common_genes)

  if (!is.null(genes_to_keep)) {
    genes <- intersect(genes, genes_to_keep)
  }

  data.frame(
    gene = genes,
    secondary_logFC = secondary[genes, "logFC"],
    secondary_adjP = secondary[genes, "adj.P.Val"],
    primary_logFC = primary[genes, "logFC"],
    primary_adjP = primary[genes, "adj.P.Val"],
    stringsAsFactors = FALSE
  )
}


read_de_result <- function(file) {
  x <- utils::read.table(
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

  if (is.null(rownames(x)) || any(rownames(x) == seq_len(nrow(x)))) {
    warning(
      "Input file may not have gene IDs as row names: ",
      file,
      call. = FALSE
    )
  }

  x
}


validate_manifest_files <- function(manifest, file_cols) {
  missing_cols <- setdiff(file_cols, colnames(manifest))

  if (length(missing_cols) > 0) {
    stop(
      "Missing expected manifest columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  missing_files <- lapply(file_cols, function(col) {
    files <- manifest[[col]]
    missing_idx <- which(is.na(files) | !file.exists(files))

    if (length(missing_idx) == 0) {
      return(NULL)
    }

    data.frame(
      row = missing_idx,
      cell_type = manifest$cell_type[missing_idx],
      column = col,
      file = files[missing_idx],
      stringsAsFactors = FALSE
    )
  })

  missing_files <- do.call(rbind, missing_files)

  if (!is.null(missing_files) && nrow(missing_files) > 0) {
    print(missing_files, row.names = FALSE)

    stop(
      "Missing files found in manifest: ",
      nrow(missing_files),
      call. = FALSE
    )
  }

  logger::log_info("All manifest files found")
  invisible(TRUE)
}


## ------------------------------------------------------------------
## Utility helpers
## ------------------------------------------------------------------

make_dataset_label <- function(x) {
  toupper(gsub("_", " ", x))
}


make_output_prefix <- function(primary_dataset,
                               secondary_dataset,
                               gene_set = "both") {
  prefix <- sprintf(
    "%s_significant_validated_in_%s",
    .sanitize_filename(primary_dataset),
    .sanitize_filename(secondary_dataset)
  )

  if (gene_set != "both") {
    prefix <- sprintf(
      "%s_%s",
      prefix,
      .sanitize_filename(gene_set)
    )
  }

  prefix
}


validate_gene_set_metadata <- function(gene_set,
                                       contig_yaml_file,
                                       reduced_gtf_file) {
  one_missing <- xor(
    is.null(contig_yaml_file),
    is.null(reduced_gtf_file)
  )

  if (one_missing) {
    stop(
      "contig_yaml_file and reduced_gtf_file must either both be supplied or both be NULL.",
      call. = FALSE
    )
  }

  if (gene_set != "both" && is.null(contig_yaml_file)) {
    stop(
      sprintf(
        "gene_set = '%s' requires contig_yaml_file and reduced_gtf_file.",
        gene_set
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


get_de_genes_for_gene_set <- function(gene_set = c("both", "autosome", "xy"),
                                      contig_yaml_file = NULL,
                                      reduced_gtf_file = NULL) {
  gene_set <- match.arg(gene_set)

  if (gene_set == "both") {
    return(NULL)
  }

  contig_groups <- yaml::yaml.load_file(contig_yaml_file)

  if (gene_set == "autosome") {
    contigs <- names(Filter(
      function(x) identical(unname(x), "autosome"),
      contig_groups
    ))
  } else {
    contigs <- names(Filter(
      function(x) any(unlist(x, use.names = FALSE) %in% c("X", "Y")),
      contig_groups
    ))
  }

  gtf <- data.table::fread(
    reduced_gtf_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )

  unique(
    gtf[
      gtf$chr %in% contigs & gtf$annotationType == "gene",
      "gene_name"
    ][[1L]]
  )
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
