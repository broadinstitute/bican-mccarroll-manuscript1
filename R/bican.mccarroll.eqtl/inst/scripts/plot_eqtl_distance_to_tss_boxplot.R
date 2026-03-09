# Boxplot + beeswarm of eQTL variant-to-TSS absolute distance by K-means cluster.
# Shows full distribution per cluster; boxplot center line = median.
# Sorted by increasing median distance, 200 kb cap, linear scale.

library(ggplot2)
library(ggbeeswarm)
library(data.table)

# ---- Inputs ----
args <- commandArgs(trailingOnly = TRUE)
index_snp_matrix_path    <- args[1]  # index_snp_slope_matrix (has phenotype_id, variant_id)
cluster_assignments_path <- args[2]  # cluster_assignments TSV (gene, cluster)
start_distance_path      <- args[3]  # index_snp_start_distance TSV (phenotype_id, variant_id, start_distance)
output_path              <- args[4]  # output path without extension

# ---- 1. Load index SNP list ----
index_snp_df <- fread(index_snp_matrix_path, select = c("phenotype_id", "variant_id"))

# ---- 2. Load cluster assignments ----
cluster_df <- fread(cluster_assignments_path)
setnames(cluster_df, c("phenotype_id", "cluster"))

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
  "6"  = "Neuronal (6)",
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
cluster_order[, display_label := paste0(cluster_label, "\n(n=", N, ")")]

merged_dt <- merge(merged_dt, cluster_order[, .(cluster_label, display_label)], by = "cluster_label")
merged_dt[, display_label := factor(display_label, levels = cluster_order$display_label)]

# ---- 6. Compute median per cluster (on original, uncapped distances) ----
median_df <- merged_dt[, .(median_dist = median(abs_distance_kb)), by = display_label]

# ---- 7. Plot (200 kb cap, linear scale) ----
p <- ggplot(merged_dt, aes(x = display_label, y = abs_distance_kb_capped)) +
  geom_boxplot(
    fill = "steelblue",
    outlier.shape = NA,
    alpha = 0.5,
    width = 0.6,
    linewidth = 0.4
  ) +
  geom_quasirandom(
    size = 0.15,
    alpha = 0.3,
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
    axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 11),
    plot.margin = margin(20, 10, 10, 10)
  )

ggsave(paste0(output_path, ".png"), plot = p, width = 14, height = 7, dpi = 300)

message("Done. Plot saved to: ", output_path)
