# Example usage:
#   devtools::load_all()
#   plot_fisher_exact(
#       fisher_table_path = "/broad/.../manuscript_data/SCZ_eur_fisher_contingency_counts_gene_clusters.tsv",
#       plot_disease_label = "schizophrenia",
#       cluster_order = c(5, 0, 6, 2, 7, 8, 10, 1, 9, 4, 3),
#       output_path = "/broad/.../manuscript_data/SCZ_eur_cluster_enrichment.svg"
#   )

#' Plot Fisher's exact test enrichment of colocalized genes by cluster
#'
#' For each K-means cluster, runs Fisher's exact test (one-sided enrichment)
#' on the contingency counts produced by
#' \code{build_fisher_contingency_table()} and plots odds ratios with 95\%
#' confidence intervals as a forest-style dot plot.
#'
#' @param fisher_table_path Path to the contingency counts TSV (columns:
#'   \code{cluster}, \code{coloc_in_cluster}, \code{noncoloc_in_cluster},
#'   \code{coloc_not_in_cluster}, \code{noncoloc_not_in_cluster}).
#' @param plot_disease_label Display label for the GWAS trait
#'   (e.g., \code{"schizophrenia"}, \code{"Alzheimer's disease"}).
#' @param cluster_order Character vector giving the desired display order
#'   of clusters (e.g., \code{c("5","0","6","2","7","8","10","1","9","4","3")}).
#' @param output_path If not \code{NULL}, saves the plot to this path.
#' @param width Plot width in inches (default 4).
#' @param height Plot height in inches (default 7).
#' @param dpi Plot resolution (default 150).
#'
#' @return A \code{data.table} with Fisher test results (odds ratio,
#'   p-value, adjusted p-value, confidence intervals) invisibly.
#'
#' @importFrom data.table fread
#' @importFrom stats fisher.test p.adjust
#' @importFrom ggplot2 ggplot aes geom_vline geom_errorbarh geom_point
#'   geom_text scale_color_manual labs coord_cartesian theme_classic theme
#'   element_text
#' @noRd

# plot_fisher_exact <- function(fisher_table_path,
#                               plot_disease_label,
#                               cluster_order,
#                               output_path = NULL,
#                               width = 4,
#                               height = 7,
#                               dpi = 150) {
#
#     df <- fread(fisher_table_path)
#
#     ## Fisher's exact test per cluster
#     oddsratio <- numeric(nrow(df))
#     p_value   <- numeric(nrow(df))
#     ci_low    <- numeric(nrow(df))
#     ci_high   <- numeric(nrow(df))
#
#     for (i in seq_len(nrow(df))) {
#         m <- matrix(
#             c(df$coloc_in_cluster[i],     df$noncoloc_in_cluster[i],
#               df$coloc_not_in_cluster[i], df$noncoloc_not_in_cluster[i]),
#             nrow = 2, byrow = TRUE
#         )
#         ft_p  <- fisher.test(m, alternative = "greater")
#         ft_ci <- fisher.test(m, alternative = "two.sided")
#
#         oddsratio[i] <- unname(ft_p$estimate)
#         p_value[i]   <- ft_p$p.value
#         ci_low[i]    <- ft_ci$conf.int[1]
#         ci_high[i]   <- ft_ci$conf.int[2]
#     }
#
#     res <- df
#     res$oddsratio        <- oddsratio
#     res$p_value          <- p_value
#     res$ci_low           <- ci_low
#     res$ci_high          <- ci_high
#     res$adjusted_p_value <- p.adjust(res$p_value, method = "BH")
#     res$is_highlight     <- res$adjusted_p_value < 0.05
#     res$p_label          <- paste0("p=", signif(res$adjusted_p_value, 2))
#
#     ## Sort: highlight first, then by p-value
#     res <- res[order(!res$is_highlight, res$adjusted_p_value), ]
#
#     ## Cluster ordering
#     order_vec <- as.character(cluster_order)
#     res$cluster <- as.character(res$cluster)
#     res$cluster_label <- paste0("Cluster ", res$cluster)
#     label_order <- paste0("Cluster ", order_vec)
#     res$cluster_label <- factor(res$cluster_label, levels = rev(label_order))
#     res <- res[order(res$cluster_label), ]
#
#     ## Cap CI for plotting
#     x_cap  <- 9
#     x_text <- 9.5
#     res$ci_high_plot <- pmin(res$ci_high, x_cap)
#
#     ## Plot (clusters on y-axis, odds ratio on x-axis)
#     p <- ggplot(res, aes(y = cluster_label)) +
#         ggplot2::geom_vline(xintercept = 1, linetype = "dashed",
#                    linewidth = 0.4, color = "grey40") +
#         ggplot2::geom_errorbarh(
#             aes(xmin = ci_low, xmax = ci_high_plot, color = is_highlight),
#             height = 0, linewidth = 0.9
#         ) +
#         geom_point(
#             aes(x = oddsratio, color = is_highlight),
#             size = 2.6
#         ) +
#         scale_color_manual(
#             values = c(`TRUE` = "#d62728", `FALSE` = "grey70"),
#             guide  = "none"
#         ) +
#         labs(
#             x     = "Odds ratio with 95% CI",
#             y     = NULL,
#             title = paste0("Enrichment of ", plot_disease_label,
#                            "\ncolocalized genes by cluster")
#         ) +
#         coord_cartesian(xlim = c(0, x_cap), clip = "off") +
#         theme_classic(base_size = 12) +
#         theme(
#             axis.text.y  = element_text(size = 10),
#             plot.title   = element_text(size = 12, face = "bold",
#                                         hjust = 0.5,
#                                         margin = ggplot2::margin(b = 15)),
#             plot.margin  = ggplot2::margin(5.5, 60, 5.5, 5.5)
#         ) +
#         geom_text(
#             data     = subset(res, as.numeric(adjusted_p_value) < 0.05),
#             aes(x = x_text, label = p_label, color = is_highlight),
#             hjust = 0, size = 3.4, fontface = "bold"
#         )
#
#     if (!is.null(output_path)) {
#         ggplot2::ggsave(output_path, p, width = width, height = height, dpi = dpi)
#     }
#
#     invisible(res)
# }
