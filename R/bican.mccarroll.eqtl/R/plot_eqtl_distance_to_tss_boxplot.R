# Example usage:
#   devtools::load_all()
#   plot_eqtl_distance_to_tss_boxplot(
#       index_snp_matrix_path = "/broad/.../index_snp_slope_matrix_with_zero_impute_qval_0.01.tsv",
#       cluster_assignments_path = "/broad/.../cluster_assignments_qval_0.01_k11.tsv",
#       start_distance_path = "/broad/.../index_snp_start_distance_qval_0.01.tsv",
#       output_path = "/broad/.../eqtl_distance_to_tss_boxplot.svg",
#       orientation = "horizontal"
#   )

#' Plot eQTL variant-to-TSS distance by K-means cluster
#'
#' Boxplot + beeswarm of absolute distance from variant to TSS for each
#' K-means cluster.  Shows full distribution per cluster; boxplot center
#' line = median.  Sorted by increasing median distance, 200 kb cap,
#' linear scale.
#'
#' @param index_snp_matrix_path Character scalar.  Path to the index SNP slope
#'   matrix TSV (must have \code{phenotype_id} and \code{variant_id} columns).
#' @param cluster_assignments_path Character scalar.  Path to cluster
#'   assignments TSV.  Supports 2-column (\code{gene, cluster}) or 3-column
#'   (\code{gene, variant_id, cluster}) formats.
#' @param start_distance_path Character scalar.  Path to the start distance TSV
#'   (output of \code{\link{get_index_snp_start_distance}}).
#' @param output_path Character scalar or \code{NULL}.  If non-NULL, saves the
#'   plot to this path.  Format is inferred from the file extension.
#' @param orientation Character scalar.  \code{"vertical"} (default) or
#'   \code{"horizontal"}.
#'
#' @return The \code{ggplot} object (invisibly).
#'
#' @importFrom data.table fread setnames uniqueN data.table
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_text scale_x_continuous
#'   scale_y_continuous expansion labs theme_classic theme element_text margin
#'   ggsave
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom logger log_info
#' @noRd

###################################################################################
# THIS CONTAINS HARD CODED CLUSTER LABELS AND NEEDS TO BE REFINED BEFORE SAFE USE.
###################################################################################
# plot_eqtl_distance_to_tss_boxplot <- function(index_snp_matrix_path,
#                                                cluster_assignments_path,
#                                                start_distance_path,
#                                                output_path = NULL,
#                                                orientation = "vertical") {
#
#     # ---- 1. Load index SNP list ----
#     index_snp_df <- data.table::fread(index_snp_matrix_path,
#                                       select = c("phenotype_id", "variant_id"))
#
#     # ---- 2. Load cluster assignments ----
#     cluster_df <- data.table::fread(cluster_assignments_path)
#     if (ncol(cluster_df) == 2) {
#         data.table::setnames(cluster_df, c("phenotype_id", "cluster"))
#     } else {
#         data.table::setnames(cluster_df, names(cluster_df)[1], "phenotype_id")
#         cluster_df <- cluster_df[, .(phenotype_id, cluster)]
#     }
#
#     # ---- 3. Get TSS distance ----
#     dist_dt <- data.table::fread(start_distance_path)
#
#     # ---- 4. Merge: index SNPs + start_distance + clusters ----
#     merged_dt <- merge(index_snp_df, dist_dt,
#                        by = c("phenotype_id", "variant_id"), all.x = TRUE)
#     merged_dt[, abs_distance_kb := abs(start_distance) / 1000]
#     merged_dt <- merge(merged_dt, cluster_df, by = "phenotype_id", all.x = TRUE)
#     merged_dt <- merged_dt[!is.na(cluster) & !is.na(abs_distance_kb)]
#
#     logger::log_info("Plotting {nrow(merged_dt)} genes across {data.table::uniqueN(merged_dt$cluster)} clusters")
#
#     # ---- 4b. Cap distance at 200 kb ----
#     distance_cap_kb <- 200
#     n_beyond_cap <- merged_dt[abs_distance_kb > distance_cap_kb, .N]
#     logger::log_info("{n_beyond_cap} / {nrow(merged_dt)} eQTLs ({round(100 * n_beyond_cap / nrow(merged_dt), 1)}%) beyond {distance_cap_kb} kb cap")
#     merged_dt[, abs_distance_kb_capped := pmin(abs_distance_kb, distance_cap_kb)]
#
#     # ---- 5. Cluster labels ----
#     cluster_annotation <- c(
#         "0"  = "Pan-cell-type (0)",
#         "1"  = "Astrocyte",
#         "2"  = "Neuron (2)",
#         "3"  = "Microglia",
#         "4"  = "OPC",
#         "5"  = "Pan-cell-type (5)",
#         "6"  = "Neuron (6)",
#         "7"  = "MSN (7)",
#         "8"  = "MSN (8)",
#         "9"  = "Oligodendrocyte",
#         "10" = "Cortical neuron (10)"
#     )
#
#     merged_dt[, cluster_label := cluster_annotation[as.character(cluster)]]
#
#     # ---- 5b. Sort clusters by increasing median distance ----
#     median_by_cluster <- merged_dt[, .(median_dist = median(abs_distance_kb)),
#                                    by = cluster_label]
#     median_by_cluster <- median_by_cluster[order(median_dist)]
#     sorted_cluster_labels <- median_by_cluster$cluster_label
#
#     counts <- merged_dt[, .N, by = cluster_label]
#     cluster_order <- data.table::data.table(cluster_label = sorted_cluster_labels)
#     cluster_order <- merge(cluster_order, counts,
#                            by = "cluster_label", sort = FALSE)
#     cluster_order[, display_label := paste0(cluster_label, " (n=", N, ")")]
#
#     merged_dt <- merge(merged_dt,
#                        cluster_order[, .(cluster_label, display_label)],
#                        by = "cluster_label")
#     merged_dt[, display_label := factor(display_label,
#                                         levels = cluster_order$display_label)]
#
#     # ---- 6. Compute median per cluster (on original, uncapped distances) ----
#     median_df <- merged_dt[, .(median_dist = median(abs_distance_kb)),
#                            by = display_label]
#
#     # ---- 7. Plot ----
#     if (orientation == "horizontal") {
#         merged_dt[, display_label := factor(display_label,
#                                             levels = rev(levels(display_label)))]
#         median_df[, display_label := factor(display_label,
#                                             levels = levels(merged_dt$display_label))]
#
#         p <- ggplot2::ggplot(merged_dt,
#                              ggplot2::aes(x = abs_distance_kb_capped,
#                                           y = display_label)) +
#             ggplot2::geom_boxplot(
#                 fill = "steelblue", outlier.shape = NA,
#                 alpha = 0.5, width = 0.5, linewidth = 0.4
#             ) +
#             ggbeeswarm::geom_quasirandom(
#                 groupOnX = FALSE, size = 0.05, alpha = 0.15, bandwidth = 0.8
#             ) +
#             ggplot2::geom_text(
#                 data = median_df,
#                 ggplot2::aes(x = distance_cap_kb * 1.02, y = display_label,
#                              label = paste0(round(median_dist), " kb")),
#                 size = 3.2, fontface = "bold", hjust = 0
#             ) +
#             ggplot2::scale_x_continuous(
#                 limits = c(0, distance_cap_kb * 1.12),
#                 expand = ggplot2::expansion(mult = c(0, 0))
#             ) +
#             ggplot2::labs(x = "Absolute distance from TSS (kb)", y = NULL) +
#             ggplot2::theme_classic(base_size = 14) +
#             ggplot2::theme(
#                 axis.text.y = ggplot2::element_text(size = 10),
#                 axis.title.x = ggplot2::element_text(size = 14),
#                 axis.text.x = ggplot2::element_text(size = 11),
#                 plot.margin = ggplot2::margin(10, 30, 10, 10)
#             )
#
#         if (!is.null(output_path)) {
#             ggplot2::ggsave(output_path, plot = p,
#                             width = 10, height = 7, dpi = 300)
#         }
#
#     } else {
#         p <- ggplot2::ggplot(merged_dt,
#                              ggplot2::aes(x = display_label,
#                                           y = abs_distance_kb_capped)) +
#             ggplot2::geom_boxplot(
#                 fill = "steelblue", outlier.shape = NA,
#                 alpha = 0.5, width = 0.4, linewidth = 0.4
#             ) +
#             ggbeeswarm::geom_quasirandom(
#                 size = 0.05, alpha = 0.15, bandwidth = 0.8
#             ) +
#             ggplot2::geom_text(
#                 data = median_df,
#                 ggplot2::aes(x = display_label,
#                              y = distance_cap_kb * 1.02,
#                              label = paste0(round(median_dist), " kb")),
#                 size = 3.2, fontface = "bold", vjust = 0
#             ) +
#             ggplot2::scale_y_continuous(
#                 limits = c(0, distance_cap_kb * 1.08),
#                 expand = ggplot2::expansion(mult = c(0, 0))
#             ) +
#             ggplot2::labs(x = NULL, y = "Absolute distance from TSS (kb)") +
#             ggplot2::theme_classic(base_size = 14) +
#             ggplot2::theme(
#                 axis.text.x = ggplot2::element_text(size = 10, angle = 45,
#                                                      hjust = 1),
#                 axis.title.y = ggplot2::element_text(size = 14),
#                 axis.text.y = ggplot2::element_text(size = 11),
#                 plot.margin = ggplot2::margin(20, 10, 10, 10)
#             )
#
#         if (!is.null(output_path)) {
#             ggplot2::ggsave(output_path, plot = p,
#                             width = 9, height = 7, dpi = 300)
#         }
#     }
#
#     logger::log_info("Done. Plot saved to: {output_path}")
#     invisible(p)
# }
