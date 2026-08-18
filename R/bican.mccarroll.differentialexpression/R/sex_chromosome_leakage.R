# source ("~/bican-mccarroll-manuscript1/R/bican.mccarroll.differentialexpression/R/differential_expression.R")
# source ("~/bican-mccarroll-manuscript1/R/bican.mccarroll.differentialexpression/R/select_marker_genes.R")
# source ("~/bican-mccarroll-manuscript1/R/bican.mccarroll.differentialexpression/R/createDGEList.R")
# source ("~/bican-mccarroll-manuscript1/R/bican.mccarroll.differentialexpression/R/donor_age_prediction.R")
# data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/metacells/LEVEL_6"; data_name="donor_rxn_DGEList"
# cellTypeListFile="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/metadata/cell_types_for_de_filtering_plot.txt"
# outFile="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/sex_expression_leakage/results_average.txt"
# outFileDonorLevel="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/sex_expression_leakage/results_donor.txt"
# outSVG="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/sex_expression_leakage/results.svg"
# xist_gene = "XIST";y_genes = c("DDX3Y", "RPS4Y1", "USP9Y")


# Sex-chromosome expression "leakage" analysis for village pseudobulk data.
#
# Villages pool many donors into shared physical/sequencing space. If reads
# leak across donors (e.g. via ambient RNA or demultiplexing error), a
# female-specific gene like XIST will show spurious nonzero expression in
# male donors, and male-specific Y-chromosome genes will show spurious
# nonzero expression in female donors. This file quantifies that leakage,
# per village and per cell type, as a ratio of median CPM between sexes.

# data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/metacells/LEVEL_6"; data_name="donor_rxn_DGEList"
# cellTypeListFile="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/metadata/cell_types_for_de_filtering_plot.txt"

#' Compute sex-chromosome expression leakage ratios (file interface)
#'
#'
#' @param data_dir Directory containing the serialized DGEList (passed to
#'   \code{\link{loadDGEList}} as \code{dir}).
#' @param data_name Prefix used by \code{\link{loadDGEList}} to locate the
#'   \code{*_counts.tsv.gz} / \code{*_samples.tsv.gz} files (passed as
#'   \code{prefix}).
#' @param cellTypeListFile Optional file path containing a list of cell
#'   types to include. If \code{NULL}, uses all cell types present.
#' @param outFile Optional tab-delimited output path for the village-level
#'   leakage results table. If \code{NULL}, no file is written.
#' @param outFileDonorLevel Optional tab-delimited output path for the
#'   donor-level leakage results table. If \code{NULL}, no file is written.
#' @param outSVG Optional SVG output path for the leakage summary boxplot. If
#'   \code{NULL}, no file is written.
#' @inheritParams compute_sex_chromosome_leakage
#'
#' @return Invisibly returns a list with elements:
#'   \itemize{
#'     \item \code{leakage}: the village-level data frame produced by
#'       \code{\link{compute_sex_chromosome_leakage}}.
#'     \item \code{donor_leakage}: the donor-level data frame produced by
#'       \code{\link{compute_sex_chromosome_leakage_donor_level}}.
#'     \item \code{plot}: the \code{ggplot2} object produced by
#'       \code{\link{plot_sex_chromosome_leakage_boxplot}}.
#'   }
#'
#' @export
sex_chromosome_leakage_from_files <- function(data_dir, data_name,
                                              cellTypeListFile = NULL,
                                              outFile = NULL,
                                              outFileDonorLevel = NULL,
                                              outSVG = NULL,
                                              xist_gene = "XIST",
                                              y_genes = c("DDX3Y", "RPS4Y1", "USP9Y")) {
  validate_writable_file(outFile)
  validate_writable_file(outFileDonorLevel)
  validate_writable_file(outSVG)

  dge <- loadDGEList(dir = data_dir, prefix = data_name)
  dge <- filter_dgelist_by_celltype_list(dge, cellTypeListFile)

  leakage_df <- compute_sex_chromosome_leakage(dge, xist_gene = xist_gene, y_genes = y_genes)
  donor_leakage_df <- compute_sex_chromosome_leakage_donor_level(dge, xist_gene = xist_gene, y_genes = y_genes)

  if (!is.null(outFile)) {
    write_sex_chromosome_leakage(leakage_df, outFile)
  }
  if (!is.null(outFileDonorLevel)) {
    write_sex_chromosome_leakage(donor_leakage_df, outFileDonorLevel)
  }

  p <- plot_sex_chromosome_leakage_boxplot(leakage_df, outSVG = outSVG)


  plots <- plot_leakage_numerator_denominator(leakage_df)

  plots$expression
  plots$expression_by_cell_type

  invisible(list(leakage = leakage_df, donor_leakage = donor_leakage_df, plot = p))
}

# Collapses dge to one pseudobulk observation per donor x village x cell type
# (one cell type at a time, see compute_sex_chromosome_leakage() for why) and
# returns a donor-level data frame with columns donor, village, cell_type,
# imputed_sex, xist_cpm, y_cpm. Shared by compute_sex_chromosome_leakage()
# and compute_sex_chromosome_leakage_donor_level() so the (expensive)
# collapsing step is implemented in exactly one place.
.donor_level_sex_chromosome_cpm <- function(dge, xist_gene, y_genes) {
  stopifnot("DGEList" %in% class(dge))

  required_cols <- c("donor", "village", "cell_type", "imputed_sex")
  missing_cols <- setdiff(required_cols, colnames(dge$samples))
  if (length(missing_cols) > 0) {
    stop("dge$samples is missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  genes_needed <- c(xist_gene, y_genes)
  missing_genes <- setdiff(genes_needed, rownames(dge$counts))
  if (length(missing_genes) > 0) {
    stop("The following gene(s) were not found in rownames(dge$counts): ",
      paste(missing_genes, collapse = ", "),
      call. = FALSE
    )
  }

  # Subset to one cell type at a time before collapsing donors, mirroring
  # simplify_dge_for_marker_genes() / predict_age_celltype(). collapse_by_donor()
  # is O(n_groups x n_samples); running it once across every cell type at once
  # makes that cost scale with the total sample count across the whole DGEList
  # instead of being bounded per cell type, which is prohibitively slow.
  cell_type_list <- unique(as.character(dge$samples$cell_type))

  dge_list <- vector("list", length(cell_type_list))
  names(dge_list) <- cell_type_list

  for (ct in cell_type_list) {
    logger::log_info(paste("Collapsing donors for cell type:", ct))
    dge_cell <- dge[, dge$samples$cell_type == ct, keep.lib.sizes = TRUE]

    group_key <- paste(as.character(dge_cell$samples$donor), as.character(dge_cell$samples$village), sep = "__")
    dge_cell$samples$.leakage_group_key <- group_key

    dge_list[[ct]] <- collapse_by_donor(
      dge_cell,
      donor_col = ".leakage_group_key",
      keep_cols = c("donor", "village", "cell_type", "imputed_sex")
    )
  }

  dge_collapsed <- do.call(cbind, dge_list)

  cpm_mat <- edgeR::cpm(dge_collapsed, log = FALSE)

  donor_df <- data.frame(
    donor = as.character(dge_collapsed$samples$donor),
    village = as.character(dge_collapsed$samples$village),
    cell_type = as.character(dge_collapsed$samples$cell_type),
    imputed_sex = as.character(dge_collapsed$samples$imputed_sex),
    xist_cpm = as.numeric(cpm_mat[xist_gene, ]),
    y_cpm = as.numeric(colSums(cpm_mat[y_genes, , drop = FALSE])),
    stringsAsFactors = FALSE
  )

  unrecognized_sex <- setdiff(unique(donor_df$imputed_sex), c("1", "2"))
  if (length(unrecognized_sex) > 0) {
    warning("Dropping donors with unrecognized imputed_sex value(s): ",
      paste(unrecognized_sex, collapse = ", "),
      call. = FALSE
    )
    donor_df <- donor_df[donor_df$imputed_sex %in% c("1", "2"), , drop = FALSE]
  }

  donor_df
}

#' Compute sex-chromosome expression leakage ratios per village and cell type
#'
#' For each village x cell type combination, computes counts-per-million
#' (CPM) of a female-specific marker gene (default \code{XIST}) and a set of
#' male-specific Y-chromosome genes (default \code{DDX3Y}, \code{RPS4Y1},
#' \code{USP9Y}) for every donor, then summarizes across donors by sex.
#'
#' Donors within the same village and cell type are first collapsed to a
#' single pseudobulk observation per donor (summing raw counts across any
#' remaining samples, e.g. distinct chemistries or regions) using
#' \code{\link{collapse_by_donor}}, applied one cell type at a time (as in
#' \code{simplify_dge_for_marker_genes()}) to keep each call's cost bounded
#' to that cell type's samples rather than the whole \code{DGEList}. CPM is
#' then computed once across the recombined collapsed \code{DGEList} via
#' \code{edgeR::cpm(log = FALSE)} (library-size normalization only;
#' \code{norm.factors} are left at 1, matching
#' \code{\link{collapse_by_donor}}'s output).
#'
#' A donor's Y-chromosome expression is the \strong{sum} of CPM across all of
#' \code{y_genes}.
#'
#' Two leakage ratios are computed per village x cell type:
#' \itemize{
#'   \item \code{xist_leakage_ratio} = median(XIST CPM in males) / median(XIST
#'     CPM in females). Values near 0 indicate little leakage of the female
#'     signal into male donors.
#'   \item \code{y_leakage_ratio} = median(Y CPM in females) / median(Y CPM
#'     in males). Values near 0 indicate little leakage of the male signal
#'     into female donors.
#' }
#'
#' Sex is read from \code{dge$samples$imputed_sex}, encoded \code{1 = male},
#' \code{2 = female}; any other value is dropped with a warning. If a village
#' x cell type group has zero donors of a given sex, the corresponding
#' median(s) and ratio(s) are \code{NA}.
#'
#' @param dge An \code{edgeR::DGEList} with raw counts (rows = genes,
#'   identified by gene symbol in \code{rownames(dge$counts)}) and
#'   \code{dge$samples} containing columns \code{donor}, \code{village},
#'   \code{cell_type}, and \code{imputed_sex}.
#' @param xist_gene Character scalar naming the female-specific marker gene.
#'   Default \code{"XIST"}.
#' @param y_genes Character vector naming the male-specific Y-chromosome
#'   genes to sum. Default \code{c("DDX3Y", "RPS4Y1", "USP9Y")}.
#'
#' @return A data frame with one row per village x cell type combination and
#'   columns:
#'   \itemize{
#'     \item \code{village}, \code{cell_type}
#'     \item \code{n_male}, \code{n_female}: number of donors of each sex
#'       contributing to this group.
#'     \item \code{median_xist_male}, \code{median_xist_female}
#'     \item \code{median_y_male}, \code{median_y_female}
#'     \item \code{xist_leakage_ratio}, \code{y_leakage_ratio}
#'   }
#'
#' @importFrom edgeR cpm
#' @importFrom stats median
#' @export
compute_sex_chromosome_leakage <- function(dge,
                                           xist_gene = "XIST",
                                           y_genes = c("DDX3Y", "RPS4Y1", "USP9Y")) {
  logger::log_info("Computing sex-chromosome leakage ratios per village and cell type.")

  donor_df <- .donor_level_sex_chromosome_cpm(dge, xist_gene = xist_gene, y_genes = y_genes)

  group_ids <- unique(donor_df[, c("village", "cell_type")])
  group_ids <- group_ids[order(group_ids$village, group_ids$cell_type), , drop = FALSE]

  results <- vector("list", nrow(group_ids))
  for (i in seq_len(nrow(group_ids))) {
    v <- group_ids$village[i]
    ct <- group_ids$cell_type[i]

    grp <- donor_df[donor_df$village == v & donor_df$cell_type == ct, , drop = FALSE]
    male <- grp[grp$imputed_sex == "1", , drop = FALSE]
    female <- grp[grp$imputed_sex == "2", , drop = FALSE]

    median_xist_male <- if (nrow(male) > 0) stats::median(male$xist_cpm, na.rm = TRUE) else NA_real_
    median_xist_female <- if (nrow(female) > 0) stats::median(female$xist_cpm, na.rm = TRUE) else NA_real_
    median_y_male <- if (nrow(male) > 0) stats::median(male$y_cpm, na.rm = TRUE) else NA_real_
    median_y_female <- if (nrow(female) > 0) stats::median(female$y_cpm, na.rm = TRUE) else NA_real_

    results[[i]] <- data.frame(
      village = v,
      cell_type = ct,
      n_male = nrow(male),
      n_female = nrow(female),
      median_xist_male = median_xist_male,
      median_xist_female = median_xist_female,
      median_y_male = median_y_male,
      median_y_female = median_y_female,
      xist_leakage_ratio = median_xist_male / median_xist_female,
      y_leakage_ratio = median_y_female / median_y_male,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results)
}

#' Compute donor-level sex-chromosome expression leakage
#'
#' For each donor x village x cell type observation, reports a single
#' leakage measurement using whichever marker is biologically meaningful for
#' that donor's sex: male donors (\code{imputed_sex == 1}) are scored on
#' their own \code{XIST} CPM (should be near 0), and female donors
#' (\code{imputed_sex == 2}) are scored on their own summed Y-chromosome CPM
#' (should be near 0). The other marker is not meaningful for that donor and
#' is not reported.
#'
#' Because donors vary in library size, sequencing depth, and baseline
#' expression, a donor's raw leakage CPM is normalized against a
#' \strong{same-sex peer median}: the median raw leakage CPM among the
#' \emph{other} donors of the same sex, in the same village and cell type
#' (i.e. a leave-one-out median, excluding the donor being scored). This
#' flags donors that are outliers relative to their own peer group -- e.g. a
#' single male donor with unusually high \code{XIST} compared to other male
#' donors in the same village and cell type -- which is a different signal
#' from the population-level \code{xist_leakage_ratio}/\code{y_leakage_ratio}
#' computed by \code{\link{compute_sex_chromosome_leakage}} (which compares
#' against the opposite sex's median to estimate what fraction of the "true"
#' signal leaked in).
#'
#' If fewer than 2 donors of a given sex are present in a village x cell
#' type group, there are no peers to compute a median from, and
#' \code{peer_median_leakage_cpm} / \code{donor_leakage_ratio} are \code{NA}
#' for that group.
#'
#' @inheritParams compute_sex_chromosome_leakage
#'
#' @return A data frame with one row per donor x village x cell type
#'   observation, containing:
#'   \itemize{
#'     \item \code{donor}, \code{village}, \code{cell_type}, \code{imputed_sex}
#'     \item \code{leakage_gene}: \code{"XIST"} for male donors, \code{"Y"}
#'       for female donors.
#'     \item \code{raw_leakage_cpm}: the donor's own CPM for
#'       \code{leakage_gene} (summed across \code{y_genes} for \code{"Y"}).
#'     \item \code{n_peer_donors}: number of other same-sex donors in the
#'       same village x cell type used to compute the peer median.
#'     \item \code{peer_median_leakage_cpm}: leave-one-out median of
#'       \code{raw_leakage_cpm} among those peer donors.
#'     \item \code{donor_leakage_ratio}: \code{raw_leakage_cpm /
#'       peer_median_leakage_cpm}.
#'   }
#'
#' @importFrom stats median
#' @export
compute_sex_chromosome_leakage_donor_level <- function(dge,
                                                       xist_gene = "XIST",
                                                       y_genes = c("DDX3Y", "RPS4Y1", "USP9Y")) {
  logger::log_info("Computing donor-level sex-chromosome leakage ratios.")

  donor_df <- .donor_level_sex_chromosome_cpm(dge, xist_gene = xist_gene, y_genes = y_genes)

  donor_df$leakage_gene <- ifelse(donor_df$imputed_sex == "1", "XIST", "Y")
  donor_df$raw_leakage_cpm <- ifelse(donor_df$imputed_sex == "1", donor_df$xist_cpm, donor_df$y_cpm)

  group_ids <- unique(donor_df[, c("village", "cell_type", "imputed_sex")])
  group_ids <- group_ids[order(group_ids$village, group_ids$cell_type, group_ids$imputed_sex), , drop = FALSE]

  donor_df$n_peer_donors <- NA_integer_
  donor_df$peer_median_leakage_cpm <- NA_real_

  for (i in seq_len(nrow(group_ids))) {
    v <- group_ids$village[i]
    ct <- group_ids$cell_type[i]
    sx <- group_ids$imputed_sex[i]

    idx <- which(donor_df$village == v & donor_df$cell_type == ct & donor_df$imputed_sex == sx)
    vals <- donor_df$raw_leakage_cpm[idx]

    for (j in seq_along(idx)) {
      peer_vals <- vals[-j]
      donor_df$n_peer_donors[idx[j]] <- length(peer_vals)
      donor_df$peer_median_leakage_cpm[idx[j]] <- if (length(peer_vals) > 0) {
        stats::median(peer_vals, na.rm = TRUE)
      } else {
        NA_real_
      }
    }
  }

  donor_df$donor_leakage_ratio <- donor_df$raw_leakage_cpm / donor_df$peer_median_leakage_cpm

  out <- donor_df[, c(
    "donor", "village", "cell_type", "imputed_sex", "leakage_gene",
    "raw_leakage_cpm", "n_peer_donors", "peer_median_leakage_cpm",
    "donor_leakage_ratio"
  )]
  out <- out[order(out$village, out$cell_type, out$donor), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Plot sex-chromosome expression leakage ratios
#'
#' Visualizes the output of \code{\link{compute_sex_chromosome_leakage}} as a
#' barplot faceted by cell type, with XIST and Y-chromosome leakage ratios
#' side by side on a log10 y-axis and one bar per village.
#'
#' @param leakage_df A data frame produced by
#'   \code{\link{compute_sex_chromosome_leakage}} (or
#'   \code{\link{sex_chromosome_leakage_from_files}}), containing at least
#'   \code{village}, \code{cell_type}, \code{xist_leakage_ratio}, and
#'   \code{y_leakage_ratio}.
#' @param outSVG Optional SVG output path. If \code{NULL}, no file is
#'   written.
#'
#' @return A \code{ggplot2} object.
#'
#' @importFrom grDevices svg dev.off
#' @import ggplot2
#' @export
plot_sex_chromosome_leakage <- function(leakage_df, outSVG = NULL) {
  # Make R CMD CHECK happy
  ratio_type <- ratio_value <- village <- cell_type <- NULL

  xist_long <- data.frame(
    village = leakage_df$village,
    cell_type = leakage_df$cell_type,
    ratio_type = "XIST",
    ratio_value = leakage_df$xist_leakage_ratio,
    stringsAsFactors = FALSE
  )
  y_long <- data.frame(
    village = leakage_df$village,
    cell_type = leakage_df$cell_type,
    ratio_type = "Y",
    ratio_value = leakage_df$y_leakage_ratio,
    stringsAsFactors = FALSE
  )
  plot_df <- rbind(xist_long, y_long)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = ratio_type, y = ratio_value, fill = village)) +
    ggplot2::geom_col(position = ggplot2::position_dodge()) +
    ggplot2::facet_wrap(~cell_type) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = NULL, y = "Leakage ratio (log10 scale)", fill = "Village") +
    ggplot2::theme_classic()

  if (!is.null(outSVG)) {
    logger::log_info(paste("Saving sex-chromosome leakage plot to SVG:", outSVG))
    grDevices::svg(outSVG, width = 11, height = 8)
    on.exit(grDevices::dev.off(), add = TRUE)
    print(p)
  }

  p
}

#' Plot sex-chromosome expression leakage ratios as a boxplot
#'
#' Visualizes the output of \code{\link{compute_sex_chromosome_leakage}} as a
#' boxplot faceted by cell type, with XIST and Y-chromosome leakage ratios
#' side by side, one jittered point per village overlaid on the boxplot, and
#' ratios expressed as a percentage on a linear y-axis. This is the plot
#' produced by \code{\link{sex_chromosome_leakage_from_files}}.
#'
#' @param leakage_df A data frame produced by
#'   \code{\link{compute_sex_chromosome_leakage}} (or
#'   \code{\link{sex_chromosome_leakage_from_files}}), containing at least
#'   \code{village}, \code{cell_type}, \code{xist_leakage_ratio}, and
#'   \code{y_leakage_ratio}. Non-finite ratio values are dropped before
#'   plotting.
#' @param outSVG Optional SVG output path. If \code{NULL}, no file is
#'   written.
#' @param rare_cluster_threshold Optional numeric threshold. If provided, a
#'   dashed horizontal reference line is drawn at this value and the y-axis
#'   is expanded to include it.
#'
#' @return A \code{ggplot2} object.
#'
#' @importFrom grDevices svg dev.off
#' @import ggplot2
#' @export
plot_sex_chromosome_leakage_boxplot <- function(leakage_df,
                                                outSVG = NULL,
                                                rare_cluster_threshold = NULL) {
  ratio_type <- ratio_value <- village <- cell_type <- NULL

  xist_long <- data.frame(
    village = leakage_df$village,
    cell_type = leakage_df$cell_type,
    ratio_type = "XIST",
    ratio_value = leakage_df$xist_leakage_ratio,
    stringsAsFactors = FALSE
  )

  y_long <- data.frame(
    village = leakage_df$village,
    cell_type = leakage_df$cell_type,
    ratio_type = "Y",
    ratio_value = leakage_df$y_leakage_ratio,
    stringsAsFactors = FALSE
  )

  plot_df <- rbind(xist_long, y_long)
  plot_df <- plot_df[is.finite(plot_df$ratio_value), ]
  plot_df$ratio_type <- factor(plot_df$ratio_type, levels = c("XIST", "Y"))

  ymax <- max(plot_df$ratio_value, rare_cluster_threshold, na.rm = TRUE)
  if (!is.finite(ymax)) {
    ymax <- max(plot_df$ratio_value, na.rm = TRUE)
  }

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = ratio_type, y = ratio_value)
  ) +
    ggplot2::geom_boxplot(
      outlier.shape = NA,
      width = 0.55
    ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = village),
      position = ggplot2::position_jitter(width = 0.12, height = 0),
      shape = 21,
      size = 1.8,
      alpha = 0.8
    ) +
    ggplot2::facet_wrap(~cell_type) +
    ggplot2::scale_y_continuous(
      labels = function(x) paste0(round(100 * x, 1), "%")
    ) +
    ggplot2::coord_cartesian(ylim = c(0, ymax * 1.1)) +
    ggplot2::labs(
      x = NULL,
      y = "Leakage ratio (%)",
      fill = "Village"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = ggplot2::rel(1.1)),
      strip.background = ggplot2::element_rect(fill = "white", color = "black"),
      strip.text = ggplot2::element_text(size = ggplot2::rel(1.0)),
      legend.position = "right"
    )

  if (!is.null(rare_cluster_threshold)) {
    p <- p +
      ggplot2::geom_hline(
        yintercept = rare_cluster_threshold,
        linetype = "dashed",
        linewidth = 0.4
      )
  }

  if (!is.null(outSVG)) {
    logger::log_info(paste("Saving sex-chromosome leakage boxplot to SVG:", outSVG))
    grDevices::svg(outSVG, width = 11, height = 8)
    on.exit(grDevices::dev.off(), add = TRUE)
    print(p)
  }

  p
}

#' Write sex-chromosome leakage results to a tab-delimited file
#'
#' @param df A data frame produced by
#'   \code{\link{compute_sex_chromosome_leakage}}.
#' @param file Output file path.
#'
#' @return Invisibly returns \code{TRUE}.
#'
#' @importFrom utils write.table
#' @keywords internal
write_sex_chromosome_leakage <- function(df, file) {
  logger::log_info(paste("Writing sex-chromosome leakage results to:", file))
  utils::write.table(df, file = file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  invisible(TRUE)
}

###################################################################
# Explore the asymetry between XIST and the Y chromosome leakage
###################################################################

plot_leakage_numerator_denominator <- function(leakage_df) {
  # Make R CMD CHECK happy
  role <- value <- marker <- leakage_ratio <- contrast <- NULL

  stopifnot(is.data.frame(leakage_df))

  required_cols <- c(
    "village", "cell_type", "n_male", "n_female",
    "median_xist_male", "median_xist_female",
    "median_y_male", "median_y_female",
    "xist_leakage_ratio", "y_leakage_ratio"
  )

  missing_cols <- setdiff(required_cols, colnames(leakage_df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  df <- as.data.frame(leakage_df)

  numeric_cols <- c(
    "n_male", "n_female",
    "median_xist_male", "median_xist_female",
    "median_y_male", "median_y_female",
    "xist_leakage_ratio", "y_leakage_ratio"
  )

  for (col in numeric_cols) {
    df[[col]] <- as.numeric(df[[col]])
  }

  keep <- complete.cases(df[, numeric_cols]) &
    df$n_male > 0 &
    df$n_female > 0 &
    df$median_xist_male > 0 &
    df$median_xist_female > 0 &
    df$median_y_male > 0 &
    df$median_y_female > 0 &
    df$xist_leakage_ratio > 0 &
    df$y_leakage_ratio > 0

  df <- df[keep, , drop = FALSE]

  if (nrow(df) == 0) {
    stop("No complete positive rows remain.")
  }

  expr_long <- rbind(
    data.frame(
      village = df$village,
      cell_type = df$cell_type,
      marker = "XIST",
      role = "Wrong-sex numerator",
      quantity = "male XIST",
      value = df$median_xist_male,
      stringsAsFactors = FALSE
    ),
    data.frame(
      village = df$village,
      cell_type = df$cell_type,
      marker = "XIST",
      role = "True-sex denominator",
      quantity = "female XIST",
      value = df$median_xist_female,
      stringsAsFactors = FALSE
    ),
    data.frame(
      village = df$village,
      cell_type = df$cell_type,
      marker = "Y genes",
      role = "Wrong-sex numerator",
      quantity = "female Y",
      value = df$median_y_female,
      stringsAsFactors = FALSE
    ),
    data.frame(
      village = df$village,
      cell_type = df$cell_type,
      marker = "Y genes",
      role = "True-sex denominator",
      quantity = "male Y",
      value = df$median_y_male,
      stringsAsFactors = FALSE
    )
  )

  expr_long$marker <- factor(expr_long$marker, levels = c("XIST", "Y genes"))
  expr_long$role <- factor(
    expr_long$role,
    levels = c("Wrong-sex numerator", "True-sex denominator")
  )

  ratio_long <- rbind(
    data.frame(
      village = df$village,
      cell_type = df$cell_type,
      marker = "XIST",
      leakage_ratio = df$xist_leakage_ratio,
      stringsAsFactors = FALSE
    ),
    data.frame(
      village = df$village,
      cell_type = df$cell_type,
      marker = "Y genes",
      leakage_ratio = df$y_leakage_ratio,
      stringsAsFactors = FALSE
    )
  )

  ratio_long$marker <- factor(ratio_long$marker, levels = c("XIST", "Y genes"))

  contrast_long <- rbind(
    data.frame(
      village = df$village,
      cell_type = df$cell_type,
      contrast = "Wrong-sex signal: female Y / male XIST",
      value = df$median_y_female / df$median_xist_male,
      stringsAsFactors = FALSE
    ),
    data.frame(
      village = df$village,
      cell_type = df$cell_type,
      contrast = "True-sex signal: female XIST / male Y",
      value = df$median_xist_female / df$median_y_male,
      stringsAsFactors = FALSE
    ),
    data.frame(
      village = df$village,
      cell_type = df$cell_type,
      contrast = "Leakage ratio: Y / XIST",
      value = df$y_leakage_ratio / df$xist_leakage_ratio,
      stringsAsFactors = FALSE
    )
  )

  contrast_long$contrast <- factor(
    contrast_long$contrast,
    levels = c(
      "Wrong-sex signal: female Y / male XIST",
      "True-sex signal: female XIST / male Y",
      "Leakage ratio: Y / XIST"
    )
  )

  p_expr <- ggplot2::ggplot(
    expr_long,
    ggplot2::aes(x = role, y = value)
  ) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.15, height = 0, alpha = 0.35, size = 0.9) +
    ggplot2::facet_wrap(~marker) +
    ggplot2::scale_y_continuous(trans = "log10") +
    ggplot2::labs(
      x = NULL,
      y = "Median expression, log10 scale",
      title = "Wrong-sex numerator versus true-sex denominator"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      panel.grid.minor = ggplot2::element_blank()
    )

  p_expr_by_cell_type <- ggplot2::ggplot(
    expr_long,
    ggplot2::aes(x = marker, y = value)
  ) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.15, height = 0, alpha = 0.3, size = 0.7) +
    ggplot2::facet_grid(role ~ cell_type, scales = "free_y") +
    ggplot2::scale_y_continuous(trans = "log10") +
    ggplot2::labs(
      x = NULL,
      y = "Median expression, log10 scale",
      title = "Numerator and denominator values by cell type"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
      strip.text.x = ggplot2::element_text(size = 8),
      panel.grid.minor = ggplot2::element_blank()
    )

  p_ratios <- ggplot2::ggplot(
    ratio_long,
    ggplot2::aes(x = marker, y = leakage_ratio)
  ) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.15, height = 0, alpha = 0.4, size = 0.9) +
    ggplot2::scale_y_continuous(trans = "log10") +
    ggplot2::labs(
      x = NULL,
      y = "Leakage ratio, log10 scale",
      title = "Observed leakage ratios"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  p_contrasts <- ggplot2::ggplot(
    contrast_long,
    ggplot2::aes(x = contrast, y = value)
  ) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.15, height = 0, alpha = 0.4, size = 0.9) +
    ggplot2::scale_y_continuous(trans = "log10") +
    ggplot2::labs(
      x = NULL,
      y = "Fold difference, log10 scale",
      title = "Decomposing the Y/XIST leakage-ratio asymmetry"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      panel.grid.minor = ggplot2::element_blank()
    )

  list(
    expression = p_expr,
    expression_by_cell_type = p_expr_by_cell_type,
    ratios = p_ratios,
    contrasts = p_contrasts,
    plot_data = list(
      expression = expr_long,
      ratios = ratio_long,
      contrasts = contrast_long
    )
  )
}
