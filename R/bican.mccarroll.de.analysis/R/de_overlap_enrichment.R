#########################
# Directional DE overlap enrichment between two datasets (e.g. BICAN age
# effects vs. an external AD dataset). Complements the continuous-effect
# views in de_age_sex_plots.R (compute_de_cor_mat_datasets(),
# read_de_results_datasets()) with a discrete gene-set view: are the genes
# that are significantly up/down in one dataset enriched among the genes
# that are significantly up/down (same direction) in the other, on the
# intersected tested-gene universe.
#
# Setting up a 2x2 fisher's test this way is essentially the same
# as using a hypergeometric test (see gene set enrichment), but also
# generates an odds ratio, which is nice to give us an effect size.
#
# See original code:
# https://github.com/AnnaSophieFroehlich/single_cell_aging/blob/main/07_Overlap
# _Age_DE_Alzheimer_DE/01_Enrichment_age-genes_in_AD_v2.ipynb
######################

## ------------------------------------------------------------------
## Top-level functions
## ------------------------------------------------------------------

#' Build the joined per-gene table used for directional overlap enrichment
#'
#' Reshapes a combined long DE table (as produced by
#' \code{\link{read_de_results_datasets}}) into one row per
#' \code{(cell_type, gene)}, restricted to genes with finite effect size and
#' finite adjusted P value in \emph{both} datasets. This is the
#' threshold-independent join step of the overlap enrichment analysis, split
#' out from \code{\link{compute_de_overlap_enrichment}} so callers can cache
#' it (one table spanning all cell types) and skip re-parsing the raw DE
#' result files when only significance thresholds change.
#'
#' @param de_dt A combined, prepared DE data.table (rbind of two
#'   \code{prep_de()}/\code{read_de_results()} outputs) with an added
#'   \code{dataset} column, e.g. the return value of
#'   \code{read_de_results_datasets()}.
#' @param cell_types_use Character vector of cell types to include.
#' @param dataset_age Value of the \code{dataset} column identifying the
#'   first ("age") dataset.
#' @param dataset_ad Value of the \code{dataset} column identifying the
#'   second ("AD"/comparison) dataset.
#'
#' @return A data.table with columns \code{cell_type}, \code{gene},
#'   \code{beta_age}, \code{adj_p_age}, \code{beta_ad}, \code{adj_p_ad}, one
#'   row per gene tested in both datasets for that cell type.
#'
#' @export
build_de_overlap_gene_table <- function(de_dt,
                                        cell_types_use,
                                        dataset_age = "BICAN",
                                        dataset_ad = "PMID_39402379") {
  # Make R CMD CHECK happy
  cell_type <- dataset <- gene <- log_fc <- adj_p_val <- NULL
  beta_age <- adj_p_age <- beta_ad <- adj_p_ad <- N <- NULL

  dt <- data.table::as.data.table(de_dt)
  dt <- dt[cell_type %in% cell_types_use]

  age_dt <- dt[dataset == dataset_age, list(cell_type, gene, beta_age = log_fc, adj_p_age = adj_p_val)]
  ad_dt <- dt[dataset == dataset_ad, list(cell_type, gene, beta_ad = log_fc, adj_p_ad = adj_p_val)]

  dup_age <- age_dt[, .N, by = list(cell_type, gene)][N > 1]
  if (nrow(dup_age) > 0) {
    stop(
      "Duplicate (cell_type, gene) rows in dataset_age (", dataset_age, "): ",
      nrow(dup_age)
    )
  }

  dup_ad <- ad_dt[, .N, by = list(cell_type, gene)][N > 1]
  if (nrow(dup_ad) > 0) {
    stop(
      "Duplicate (cell_type, gene) rows in dataset_ad (", dataset_ad, "): ",
      nrow(dup_ad)
    )
  }

  merged <- merge(age_dt, ad_dt, by = c("cell_type", "gene"))

  merged <- merged[
    is.finite(beta_age) & is.finite(adj_p_age) &
      is.finite(beta_ad) & is.finite(adj_p_ad)
  ]

  merged[]
}


#' Directional DE overlap enrichment between two datasets
#'
#' For each cell type, tests whether genes significantly upregulated in the
#' first ("age") dataset are enriched among genes significantly upregulated
#' in the second ("AD"/comparison) dataset, and likewise for downregulated
#' genes, using one-sided Fisher exact tests on the intersected tested-gene
#' universe. Also reports the Pearson/Spearman correlation of the continuous
#' effect sizes as a separate descriptive statistic, and lists genes that are
#' significant in both datasets with opposite effect directions.
#'
#' A cell type is included in the \code{fisher} result only if both datasets
#' have at least \code{min_de_genes} significant genes (summed across both
#' directions) in the shared universe for that cell type; excluded cell types
#' are still reported in \code{counts} with \code{included = FALSE}, so
#' nothing tested is silently dropped from the record. Benjamini-Hochberg FDR
#' correction (\code{p_adj}) is applied across all \code{2 * K} Fisher tests
#' for the \code{K} included cell types, controlling the false discovery rate
#' over the full set of directional overlap tests reported at once (not per
#' cell type) -- some individual tests (e.g. cell types with few overlapping
#' genes) may look marginal on the raw \code{p_value} but should be judged on
#' \code{p_adj}.
#'
#' \code{pearson_r}/\code{spearman_rho} are correlations of \code{beta_age}
#' vs. \code{beta_ad} among genes significant in \emph{both} datasets
#' (age-significant and AD-significant, either direction) within the shared
#' tested universe. \code{spearman_rho} additionally gets a
#' \code{stats::cor.test()} P value (\code{spearman_p_value}), BH-corrected
#' across the \code{K} included cell types as its own family
#' (\code{spearman_p_adj}), separate from the Fisher-test correction. This
#' matches the reference implementation for this comparison (Frohlich et al.
#' 2024, Nat Neurosci; github.com/AnnaSophieFroehlich/single_cell_aging,
#' \code{07_Overlap_Age_DE_Alzheimer_DE}), which merges each cell type's
#' significant comparison-dataset genes against its significant age genes
#' before running \code{cor.test(method = "spearman")} and BH-correcting the
#' resulting P values across cell types.
#'
#' @param gene_dt A data.table produced by \code{build_de_overlap_gene_table()},
#'   with columns \code{cell_type}, \code{gene}, \code{beta_age},
#'   \code{adj_p_age}, \code{beta_ad}, \code{adj_p_ad}.
#' @param fdr_cutoff_age Adjusted P value threshold used to call a gene
#'   significant in the age dataset. Defaults to 0.05, matching this
#'   manuscript's main age DE calling; Frohlich et al. 2024 uses 0.1 for the
#'   age side specifically in their AD/SCZ overlap comparisons, so pass
#'   \code{fdr_cutoff_age = 0.1} to replicate that exactly.
#' @param fdr_cutoff_ad Adjusted P value threshold used to call a gene
#'   significant in the AD/comparison dataset. Defaults to 0.05.
#' @param min_de_genes Minimum number of significant genes (summed across up
#'   and down) required in \emph{both} datasets for a cell type to be
#'   included in the \code{fisher} result and the FDR correction. Defaults
#'   to 10.
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{fisher}}{One row per included cell type per direction
#'     (\code{"up"}, \code{"down"}): \code{cell_type}, \code{direction},
#'     \code{universe_size}, \code{n_age_de}, \code{n_ad_de},
#'     \code{n_overlap}, \code{odds_ratio}, \code{p_value}, \code{p_adj},
#'     \code{pearson_r}, \code{spearman_rho}, \code{spearman_p_value},
#'     \code{spearman_p_adj} (the four correlation columns are computed on
#'     the age-and-AD-significant subset, see Details, and repeated on both
#'     direction rows for a given cell type).}
#'   \item{\code{counts}}{One row per cell type present in \code{gene_dt}
#'     (including excluded ones): \code{cell_type}, \code{included},
#'     \code{universe_size}, \code{n_age_up}, \code{n_age_down},
#'     \code{n_ad_up}, \code{n_ad_down}, \code{n_up_up}, \code{n_down_down},
#'     \code{n_age_up_ad_down}, \code{n_age_down_ad_up}, \code{pearson_r},
#'     \code{spearman_rho}, \code{spearman_p_value} (\code{NA} if a cell
#'     type has fewer than 2 genes in the age-and-AD-significant subset),
#'     \code{spearman_p_adj} (\code{NA} for excluded cell types, which are
#'     not part of that correction family).}
#'   \item{\code{discordant}}{One row per gene significant in both datasets
#'     with opposite effect signs, for included cell types only:
#'     \code{cell_type}, \code{gene}, \code{beta_age}, \code{adj_p_age},
#'     \code{beta_ad}, \code{adj_p_ad}, \code{pattern} (one of
#'     \code{"age_up_ad_down"}, \code{"age_down_ad_up"}),
#'     \code{n_cell_types_same_pattern} (how many included cell types show
#'     that same gene/pattern combination).}
#' }
#'
#' @export
compute_de_overlap_enrichment <- function(gene_dt,
                                          fdr_cutoff_age = 0.05,
                                          fdr_cutoff_ad = 0.05,
                                          min_de_genes = 10) {
  cell_type <- NULL

  dt <- data.table::as.data.table(gene_dt)

  cell_types_present <- unique(as.character(dt$cell_type))

  per_cell_type <- lapply(cell_types_present, function(ct) {
    .compute_de_overlap_enrichment_one_cell_type(
      dt[cell_type == ct],
      cell_type = ct,
      fdr_cutoff_age = fdr_cutoff_age,
      fdr_cutoff_ad = fdr_cutoff_ad,
      min_de_genes = min_de_genes
    )
  })

  fisher_dt <- data.table::rbindlist(lapply(per_cell_type, `[[`, "fisher"), use.names = TRUE)
  counts_dt <- data.table::rbindlist(lapply(per_cell_type, `[[`, "counts"), use.names = TRUE)
  discordant_dt <- data.table::rbindlist(lapply(per_cell_type, `[[`, "discordant"), use.names = TRUE)

  # Make R CMD CHECK happy
  p_value <- p_adj <- spearman_p_value <- spearman_p_adj <- NULL

  # Correction denominator is the number of tests actually run (the
  # included cell types), not the full matched-cell-type list -- raising
  # min_de_genes to run fewer, stronger tests should correct for exactly
  # those fewer tests.
  fisher_dt[, p_adj := stats::p.adjust(p_value, method = "BH")]

  # Spearman correlation FDR is its own family (one value per cell type, not
  # duplicated by direction), separate from the Fisher-test family above.
  cell_type_rho <- unique(fisher_dt[, list(cell_type, spearman_p_value)])
  cell_type_rho[, spearman_p_adj := stats::p.adjust(spearman_p_value, method = "BH")]

  fisher_dt <- merge(fisher_dt, cell_type_rho[, list(cell_type, spearman_p_adj)], by = "cell_type", sort = FALSE)
  counts_dt <- merge(counts_dt, cell_type_rho[, list(cell_type, spearman_p_adj)], by = "cell_type", all.x = TRUE, sort = FALSE)

  discordant_dt <- .annotate_discordant_recurrence(discordant_dt)

  list(
    fisher = fisher_dt[],
    counts = counts_dt[],
    discordant = discordant_dt[]
  )
}


## ------------------------------------------------------------------
## Plotting
## ------------------------------------------------------------------

#' Plot a directional DE overlap enrichment heatmap
#'
#' Three-panel figure modeled on Frohlich et al. 2024's Fig. 6a: one
#' single-column heatmap each for the up-direction odds ratio, the
#' down-direction odds ratio, and the Spearman correlation, each on its own
#' linear, sequential color scale (odds ratios are never log-transformed or
#' capped). Separate scales per panel matter here because the Fisher test is
#' one-sided (\code{alternative = "greater"}): it only ever tests for
#' enrichment, so a single shared/diverging scale would visually imply a
#' tested depletion signal that doesn't exist, and up- and down-direction
#' odds ratios in this data run on very different typical magnitudes.
#'
#' Odds-ratio cells are labeled \code{"<odds ratio> [p=1e<N>]"}, where
#' \code{N = floor(log10(p_adj))} is the order of magnitude of the
#' BH-corrected P value (e.g. \code{"4.8 [p=1e-25]"}), rather than a
#' significant/not-significant boolean, so the reader can judge the
#' strength of the evidence directly. \code{p_adj >= 1} labels as
#' \code{"[p=1]"}; \code{p_adj} underflowing to \code{0} labels as
#' \code{"[p<1e-300]"}. The bracketed P value is colored black when
#' \code{p_adj < p_adj_threshold} and grey otherwise, to visually emphasize
#' the FDR-significant cells; the odds-ratio number itself is unstyled. The
#' correlation cell is still marked with a plain asterisk using
#' \code{spearman_p_adj < p_adj_threshold}.
#'
#' @param fisher_dt The \code{fisher} element of
#'   \code{\link{compute_de_overlap_enrichment}}'s return value.
#' @param cell_type_order Optional character vector giving the top-to-bottom
#'   display order of cell types (first entry plotted at the top). If
#'   \code{NULL}, uses the order cell types first appear in \code{fisher_dt}.
#'   Cell types in \code{fisher_dt} but not in this vector are dropped.
#' @param p_adj_threshold Adjusted P value threshold used to mark the
#'   correlation cell with an asterisk. Defaults to 0.05.
#' @param up_colors,down_colors,cor_colors Two-color \code{low, high} vectors
#'   for the up odds-ratio, down odds-ratio, and correlation panels
#'   respectively. Default to the same palette as the reference
#'   implementation's figure.
#' @param label_size Font size for the printed cell labels.
#' @param panel_rel_widths Relative widths of the three panels, passed to
#'   \code{cowplot::plot_grid()}. Defaults to \code{c(2, 1, 1)} rather than
#'   equal widths, because only the "Up" panel carries the cell type axis
#'   labels (\code{show_y_axis = TRUE}) - at equal widths, those labels eat
#'   into the "Up" panel's tile/legend space and, for longer cell type names
#'   (e.g. BICAN's \verb{SPN_*}/\verb{GABA_*} names), can visually overlap the
#'   in-tile odds-ratio/P-value text.
#'
#' @return A \code{cowplot} plot grid object combining all three panels.
#'
#' @export
plot_de_overlap_enrichment_heatmap <- function(fisher_dt,
                                               cell_type_order = NULL,
                                               p_adj_threshold = 0.05,
                                               up_colors = c("#FFFFE6", "#FDE725"),
                                               down_colors = c("#D7F5B7", "#4AC16D"),
                                               cor_colors = c("#ACDDF4", "#365C8D"),
                                               label_size = 3,
                                               panel_rel_widths = c(2, 1, 1)) {
  # Make R CMD CHECK happy
  cell_type <- direction <- odds_ratio <- p_adj <- label_or_val <- label_p <- p_color <- NULL
  spearman_rho <- spearman_p_adj <- label_rho <- NULL

  d <- data.table::copy(data.table::as.data.table(fisher_dt))

  if (is.null(cell_type_order)) {
    cell_type_order <- unique(as.character(d$cell_type))
  }

  d[, cell_type := factor(cell_type, levels = rev(cell_type_order))]
  d <- d[!is.na(cell_type)]

  d[, label_or_val := ifelse(
    odds_ratio == 0, "0",
    ifelse(is.infinite(odds_ratio), "Inf", sprintf("%.1f", odds_ratio))
  )]
  d[, label_p := paste0("[", .format_p_order_of_magnitude(p_adj), "]")]
  d[, p_color := ifelse(p_adj < p_adj_threshold, "black", "grey35")]

  .or_panel <- function(dir_value, colors, title, show_y_axis) {
    dd <- d[direction == dir_value]

    p <- ggplot2::ggplot(dd, ggplot2::aes(x = 1, y = cell_type, fill = odds_ratio)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::geom_text(ggplot2::aes(x = 0.58, label = label_or_val), hjust = 0, size = label_size) +
      ggplot2::geom_text(ggplot2::aes(x = 1.42, label = label_p, color = p_color), hjust = 1, size = label_size) +
      ggplot2::scale_color_identity(guide = "none") +
      ggplot2::scale_fill_gradient(low = colors[1], high = colors[2], name = "Odds ratio") +
      ggplot2::labs(x = NULL, y = NULL, title = title) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())

    if (!show_y_axis) {
      p <- p + ggplot2::theme(axis.text.y = ggplot2::element_blank())
    }

    p
  }

  p_up <- .or_panel("up", up_colors, "Up", show_y_axis = TRUE)
  p_down <- .or_panel("down", down_colors, "Down", show_y_axis = FALSE)

  cor_dt <- unique(d[, list(cell_type, spearman_rho, spearman_p_adj)])
  cor_dt[, label_rho := sprintf("%.2f", spearman_rho)]
  cor_dt[, label_rho := ifelse(spearman_p_adj < p_adj_threshold, paste0(label_rho, "*"), label_rho)]

  p_cor <- ggplot2::ggplot(cor_dt, ggplot2::aes(x = 1, y = cell_type, fill = spearman_rho)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = label_rho), size = label_size) +
    ggplot2::scale_fill_gradient(low = cor_colors[1], high = cor_colors[2], name = "Spearman rho") +
    ggplot2::labs(x = NULL, y = NULL, title = "rho") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank()
    )

  cowplot::plot_grid(p_up, p_down, p_cor, nrow = 1, align = "h", axis = "tb", rel_widths = panel_rel_widths)
}


## ------------------------------------------------------------------
## Helpers
## ------------------------------------------------------------------

.compute_de_overlap_enrichment_one_cell_type <- function(g,
                                                         cell_type,
                                                         fdr_cutoff_age,
                                                         fdr_cutoff_ad,
                                                         min_de_genes) {
  universe_size <- nrow(g)

  age_up <- g$adj_p_age <= fdr_cutoff_age & g$beta_age > 0
  age_down <- g$adj_p_age <= fdr_cutoff_age & g$beta_age < 0
  ad_up <- g$adj_p_ad <= fdr_cutoff_ad & g$beta_ad > 0
  ad_down <- g$adj_p_ad <= fdr_cutoff_ad & g$beta_ad < 0

  n_age_de <- sum(age_up) + sum(age_down)
  n_ad_de <- sum(ad_up) + sum(ad_down)
  included <- (n_age_de >= min_de_genes) && (n_ad_de >= min_de_genes)

  # Correlation is computed on the age-and-AD-significant subset; see
  # Details in compute_de_overlap_enrichment()'s documentation for why.
  both_sig <- (age_up | age_down) & (ad_up | ad_down)
  pearson_r <- .cor_on_subset(g$beta_age, g$beta_ad, both_sig, method = "pearson")
  spearman_test <- .cor_test_on_subset(g$beta_age, g$beta_ad, both_sig, method = "spearman")
  spearman_rho <- spearman_test$estimate
  spearman_p_value <- spearman_test$p_value

  counts <- data.table::data.table(
    cell_type = cell_type,
    included = included,
    universe_size = universe_size,
    n_age_up = sum(age_up),
    n_age_down = sum(age_down),
    n_ad_up = sum(ad_up),
    n_ad_down = sum(ad_down),
    n_up_up = sum(age_up & ad_up),
    n_down_down = sum(age_down & ad_down),
    n_age_up_ad_down = sum(age_up & ad_down),
    n_age_down_ad_up = sum(age_down & ad_up),
    pearson_r = pearson_r,
    spearman_rho = spearman_rho,
    spearman_p_value = spearman_p_value
  )

  if (!included) {
    return(list(
      fisher = data.table::data.table(
        cell_type = character(0), direction = character(0),
        universe_size = integer(0), n_age_de = integer(0), n_ad_de = integer(0),
        n_overlap = integer(0), odds_ratio = numeric(0), p_value = numeric(0),
        pearson_r = numeric(0), spearman_rho = numeric(0), spearman_p_value = numeric(0)
      ),
      counts = counts,
      discordant = .empty_discordant_dt()
    ))
  }

  fisher_up <- .fisher_2x2_greater(age_up, ad_up)
  fisher_down <- .fisher_2x2_greater(age_down, ad_down)

  fisher <- data.table::data.table(
    cell_type = cell_type,
    direction = c("up", "down"),
    universe_size = universe_size,
    n_age_de = c(sum(age_up), sum(age_down)),
    n_ad_de = c(sum(ad_up), sum(ad_down)),
    n_overlap = c(sum(age_up & ad_up), sum(age_down & ad_down)),
    odds_ratio = c(fisher_up$odds_ratio, fisher_down$odds_ratio),
    p_value = c(fisher_up$p_value, fisher_down$p_value),
    pearson_r = pearson_r,
    spearman_rho = spearman_rho,
    spearman_p_value = spearman_p_value
  )

  age_up_ad_down <- age_up & ad_down
  age_down_ad_up <- age_down & ad_up

  discordant <- data.table::rbindlist(list(
    data.table::data.table(
      cell_type = cell_type,
      gene = g$gene[age_up_ad_down],
      beta_age = g$beta_age[age_up_ad_down],
      adj_p_age = g$adj_p_age[age_up_ad_down],
      beta_ad = g$beta_ad[age_up_ad_down],
      adj_p_ad = g$adj_p_ad[age_up_ad_down],
      pattern = "age_up_ad_down"
    ),
    data.table::data.table(
      cell_type = cell_type,
      gene = g$gene[age_down_ad_up],
      beta_age = g$beta_age[age_down_ad_up],
      adj_p_age = g$adj_p_age[age_down_ad_up],
      beta_ad = g$beta_ad[age_down_ad_up],
      adj_p_ad = g$adj_p_ad[age_down_ad_up],
      pattern = "age_down_ad_up"
    )
  ), use.names = TRUE)

  list(fisher = fisher, counts = counts, discordant = discordant)
}


.format_p_order_of_magnitude <- function(p) {
  ifelse(
    p >= 1, "p=1",
    ifelse(
      p <= 0, "p<1e-300",
      sprintf("p=1e%d", floor(log10(p)))
    )
  )
}


.fisher_2x2_greater <- function(set_a, set_b) {
  a <- sum(set_a & set_b)
  b <- sum(set_a & !set_b)
  c <- sum(!set_a & set_b)
  d <- sum(!set_a & !set_b)

  tab <- matrix(c(a, b, c, d), nrow = 2, ncol = 2, byrow = TRUE)

  ft <- stats::fisher.test(tab, alternative = "greater")

  list(odds_ratio = unname(ft$estimate), p_value = ft$p.value)
}


.cor_on_subset <- function(beta_x, beta_y, subset, method) {
  if (sum(subset) < 2) {
    return(NA_real_)
  }

  stats::cor(beta_x[subset], beta_y[subset], method = method)
}


.cor_test_on_subset <- function(beta_x, beta_y, subset, method) {
  if (sum(subset) < 2) {
    return(list(estimate = NA_real_, p_value = NA_real_))
  }

  ct <- stats::cor.test(beta_x[subset], beta_y[subset], method = method)

  list(estimate = unname(ct$estimate), p_value = ct$p.value)
}


.empty_discordant_dt <- function() {
  data.table::data.table(
    cell_type = character(0),
    gene = character(0),
    beta_age = numeric(0),
    adj_p_age = numeric(0),
    beta_ad = numeric(0),
    adj_p_ad = numeric(0),
    pattern = character(0)
  )
}


.annotate_discordant_recurrence <- function(discordant_dt) {
  gene <- pattern <- n_cell_types_same_pattern <- NULL

  if (nrow(discordant_dt) == 0) {
    discordant_dt[, n_cell_types_same_pattern := integer(0)]
    return(discordant_dt[])
  }

  recurrence <- discordant_dt[, list(n_cell_types_same_pattern = .N), by = list(gene, pattern)]
  merge(discordant_dt, recurrence, by = c("gene", "pattern"))
}
