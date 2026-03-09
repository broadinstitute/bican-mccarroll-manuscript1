# Boxplot + beeswarm of eQTL variant-to-TSS absolute distance by K-means cluster.
# Shows full distribution per cluster; boxplot center line = median.
# Sorted by increasing median distance, 200 kb cap, linear scale.
#
# Usage:
#   Rscript plot_eqtl_distance_to_tss_boxplot.R <index_snp_matrix> <cluster_assignments> <start_distance> <output_path> [orientation]
#
# orientation: "vertical" (default) or "horizontal"

library(ggplot2)
library(ggbeeswarm)
library(data.table)

# ---- Inputs ----
args <- commandArgs(trailingOnly = TRUE)
index_snp_matrix_path    <- args[1]  # index_snp_slope_matrix (has phenotype_id, variant_id)
cluster_assignments_path <- args[2]  # cluster_assignments TSV (gene, cluster)
start_distance_path      <- args[3]  # index_snp_start_distance TSV (phenotype_id, variant_id, start_distance)
output_path              <- args[4]  # output path with extension (e.g. .png or .svg)
orientation              <- ifelse(length(args) >= 5, args[5], "vertical")

# ---- 1. Load index SNP list ----
index_snp_df <- fread(index_snp_matrix_path, select = c("phenotype_id", "variant_id"))

# ---- 2. Load cluster assignments ----
cluster_df <- fread(cluster_assignments_path)
# Handle both 2-col (gene, cluster) and 3-col (gene, variant_id, cluster) formats
if (ncol(cluster_df) == 2) {
    setnames(cluster_df, c("phenotype_id", "cluster"))
} else {
    setnames(cluster_df, names(cluster_df)[1], "phenotype_id")
    cluster_df <- cluster_df[, .(phenotype_id, cluster)]
}

# ---- 3. Get TSS distance ----
dist_dt <- fread(start_distance_path)

# ---- 4. Merge: index SNPs + start_distance + clusters ----
merged_dt <- merge(index_snp_df, dist_dt, by = c("phenotype_id", "variant_id"), all.x = TRUE)
merged_dt[, abs_distance_kb := abs(start_distance) / 1000]
merged_dt <- merge(merged_dt, cluster_df, by = "phenotype_id", all.x = TRUE)
merged_dt <- merged_dt[!is.na(cluster) & !is.na(abs_distance_kb)]

message(sprintf("Plotting %d genes across %d clusters", nrow(merged_dt), uniqueN(merged_dt$cluster)))

# ---- 4b. Cap distance at 200 kb ----
distance_cap_kb <- 200
n_beyond_cap <- merged_dt[abs_distance_kb > distance_cap_kb, .N]
message(sprintf("%d / %d eQTLs (%.1f%%) beyond %d kb cap",
                n_beyond_cap, nrow(merged_dt),
                100 * n_beyond_cap / nrow(merged_dt), distance_cap_kb))
merged_dt[, abs_distance_kb_capped := pmin(abs_distance_kb, distance_cap_kb)]

# ---- 5. Cluster labels ----
cluster_annotation <- c(
  "0"  = "Pan-cell-type (0)",
  "1"  = "Astrocyte",
  "2"  = "Neuron (2)",
  "3"  = "Microglia",
  "4"  = "OPC",
  "5"  = "Pan-cell-type (5)",
  "6"  = "Neuron (6)",
  "7"  = "MSN (7)",
  "8"  = "MSN (8)",
  "9"  = "Oligodendrocyte",
  "10" = "Cortical neuron (10)"
)

merged_dt[, cluster_label := cluster_annotation[as.character(cluster)]]

# ---- 5b. Sort clusters by increasing median distance ----
median_by_cluster <- merged_dt[, .(median_dist = median(abs_distance_kb)), by = cluster_label]
median_by_cluster <- median_by_cluster[order(median_dist)]
sorted_cluster_labels <- median_by_cluster$cluster_label

message("Cluster order by increasing median distance (kb):")
for (i in seq_len(nrow(median_by_cluster))) {
  message(sprintf("  %s: %.1f kb", median_by_cluster$cluster_label[i],
                  median_by_cluster$median_dist[i]))
}

counts <- merged_dt[, .N, by = cluster_label]
cluster_order <- data.table(cluster_label = sorted_cluster_labels)
cluster_order <- merge(cluster_order, counts, by = "cluster_label", sort = FALSE)
cluster_order[, display_label := paste0(cluster_label, " (n=", N, ")")]

merged_dt <- merge(merged_dt, cluster_order[, .(cluster_label, display_label)], by = "cluster_label")
merged_dt[, display_label := factor(display_label, levels = cluster_order$display_label)]

# ---- 6. Compute median per cluster (on original, uncapped distances) ----
median_df <- merged_dt[, .(median_dist = median(abs_distance_kb)), by = display_label]

# ---- 7. Plot ----
if (orientation == "horizontal") {
  # Horizontal: distance on x-axis, clusters on y-axis (reversed so smallest median at top)
  merged_dt[, display_label := factor(display_label, levels = rev(levels(display_label)))]
  median_df[, display_label := factor(display_label, levels = levels(merged_dt$display_label))]

  p <- ggplot(merged_dt, aes(x = abs_distance_kb_capped, y = display_label)) +
    geom_boxplot(
      fill = "steelblue",
      outlier.shape = NA,
      alpha = 0.5,
      width = 0.5,
      linewidth = 0.4
    ) +
    geom_quasirandom(
      groupOnX = FALSE,
      size = 0.05,
      alpha = 0.15,
      bandwidth = 0.8
    ) +
    geom_text(
      data = median_df,
      aes(x = distance_cap_kb * 1.02, y = display_label,
          label = paste0(round(median_dist), " kb")),
      size = 3.2, fontface = "bold", hjust = 0
    ) +
    scale_x_continuous(
      limits = c(0, distance_cap_kb * 1.12),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      x = "Absolute distance from TSS (kb)",
      y = NULL
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 14),
      axis.text.x = element_text(size = 11),
      plot.margin = margin(10, 30, 10, 10)
    )

  ggsave(output_path, plot = p, width = 10, height = 7, dpi = 300)

} else {
  # Vertical (default): clusters on x-axis, distance on y-axis
  p <- ggplot(merged_dt, aes(x = display_label, y = abs_distance_kb_capped)) +
    geom_boxplot(
      fill = "steelblue",
      outlier.shape = NA,
      alpha = 0.5,
      width = 0.4,
      linewidth = 0.4
    ) +
    geom_quasirandom(
      size = 0.05,
      alpha = 0.15,
      bandwidth = 0.8
    ) +
    geom_text(
      data = median_df,
      aes(x = display_label, y = distance_cap_kb * 1.02,
          label = paste0(round(median_dist), " kb")),
      size = 3.2, fontface = "bold", vjust = 0
    ) +
    scale_y_continuous(
      limits = c(0, distance_cap_kb * 1.08),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      x = NULL,
      y = "Absolute distance from TSS (kb)"
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.title.y = element_text(size = 14),
      axis.text.y = element_text(size = 11),
      plot.margin = margin(20, 10, 10, 10)
    )

  ggsave(output_path, plot = p, width = 9, height = 7, dpi = 300)
}

message("Done. Plot saved to: ", output_path)
