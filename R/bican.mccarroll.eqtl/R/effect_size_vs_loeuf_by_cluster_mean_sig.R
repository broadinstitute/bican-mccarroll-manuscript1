#' Effect size vs LOEUF by K-means cluster (mean of significant effects)
#'
#' For each gene-SNP pair, computes the mean absolute effect size across cell
#' types where the nominal p-value is below the per-gene empirical threshold
#' from tensorQTL.  Merges with cluster assignments, filters to protein-coding
#' genes with gnomAD LOEUF scores, then computes per-cluster Spearman
#' correlations between LOEUF and mean effect size.  Also fits an overall
#' linear model regressing effect size on LOEUF + cluster.  Produces a
#' faceted scatterplot and TSV outputs.
#'
#' @param slope_matrix_path Character scalar.  Path to the index-SNP slope
#'   matrix TSV.
#' @param pval_nominal_path Character scalar.  Path to the nominal p-value
#'   matrix TSV (same row/column layout as slope matrix).
#' @param pval_threshold_path Character scalar.  Path to the per-gene
#'   empirical threshold matrix TSV (same layout).
#' @param index_snp_path Character scalar.  Path to the index-SNP matrix TSV
#'   with columns \code{phenotype_id} and \code{variant_id}.
#' @param cluster_path Character scalar.  Path to the cluster assignments TSV
#'   with columns \code{gene} (or \code{phenotype_id}) and \code{cluster}.
#' @param gnomad_path Character scalar.  Path to gnomAD constraint metrics TSV.
#' @param gene_biotype_path Character scalar.  Path to gene biotype TSV (output of
#'   \code{\link{get_gene_biotype_from_gtf}}) with columns
#'   \code{phenotype_id} and \code{gene_type}.
#' @param regression_output_path Character scalar.  Path for output TSV.
#' @param plot_output_path Character scalar.  Path for output plot.
#' @param cluster_order Integer vector.  Order of K-means cluster IDs for
#'   display.  Defaults to \code{c(11, 0, 5, 4, 2, 12, 9, 1, 6, 3, 7, 10, 8)}.
#' @param plot_width Numeric.  Plot width in inches.  Default 14.
#' @param plot_height Numeric.  Plot height in inches.  Default 8.5.
#' @param ncol Integer.  Number of facet columns.  Default 5.
#'
#' @return A list with components \code{per_cluster} (per-cluster Spearman
#'   results) and \code{overall_regression} (LOEUF coefficient from linear
#'   model), both as \code{data.table}s.  Called for side effects (writes
#'   TSVs and SVG plot).
#'
#' @export
#' @importFrom data.table fread fwrite data.table rbindlist setnames
#' @importFrom stats cor.test lm p.adjust
#' @importFrom logger log_info
effect_size_vs_loeuf_by_cluster_mean_sig <- function(
    slope_matrix_path,
    pval_nominal_path,
    pval_threshold_path,
    index_snp_path,
    cluster_path,
    gnomad_path,
    gene_biotype_path,
    regression_output_path,
    plot_output_path,
    cluster_order = c(11, 0, 5, 4, 2, 12, 9, 1, 6, 3, 7, 10, 8),
    plot_width = 14,
    plot_height = 8.5,
    ncol = 5
) {

    # ---- 1. Read matrices ----
    logger::log_info("Reading slope matrix...")
    slope_dt <- data.table::fread(slope_matrix_path)

    logger::log_info("Reading p-value nominal matrix...")
    pval_dt <- data.table::fread(pval_nominal_path)

    logger::log_info("Reading p-value threshold matrix...")
    thresh_dt <- data.table::fread(pval_threshold_path)

    # Subset to index SNP gene-SNP pairs
    logger::log_info("Reading index SNP list...")
    index_dt <- data.table::fread(index_snp_path, select = c("phenotype_id", "variant_id"))
    logger::log_info("{nrow(index_dt)} index SNP gene-SNP pairs")

    index_key <- paste(index_dt$phenotype_id, index_dt$variant_id, sep = "|")
    slope_key <- paste(slope_dt$phenotype_id, slope_dt$variant_id, sep = "|")

    keep <- slope_key %in% index_key
    slope_dt  <- slope_dt[keep]
    pval_dt   <- pval_dt[keep]
    thresh_dt <- thresh_dt[keep]
    logger::log_info("{nrow(slope_dt)} rows after subsetting to index SNPs")

    ct_cols <- setdiff(names(slope_dt), c("phenotype_id", "variant_id"))

    # ---- 2. Mask non-significant slopes and compute mean |slope| ----
    logger::log_info("Masking non-significant slopes and computing mean |slope| per gene...")

    slope_m  <- as.matrix(slope_dt[, ct_cols, with = FALSE])
    pval_m   <- as.matrix(pval_dt[, ct_cols, with = FALSE])
    thresh_m <- as.matrix(thresh_dt[, ct_cols, with = FALSE])

    sig_mask <- !is.na(pval_m) & !is.na(thresh_m) & (pval_m < thresh_m)

    slope_sig <- abs(slope_m)
    slope_sig[!sig_mask] <- NA

    mean_abs_slope <- rowMeans(slope_sig, na.rm = TRUE)
    mean_abs_slope[is.nan(mean_abs_slope)] <- NA

    n_sig <- rowSums(sig_mask, na.rm = TRUE)

    effect_dt <- data.table::data.table(
        phenotype_id    = slope_dt$phenotype_id,
        mean_abs_slope  = mean_abs_slope,
        n_sig_celltypes = n_sig
    )

    # ---- 3. Read cluster assignments ----
    logger::log_info("Reading cluster assignments...")
    cluster_dt <- data.table::fread(cluster_path)
    if ("gene" %in% names(cluster_dt)) data.table::setnames(cluster_dt, "gene", "phenotype_id")

    # ---- 4. Filter to protein-coding genes ----
    logger::log_info("Reading gene biotype (filtering to protein-coding)...")
    map_dt <- data.table::fread(gene_biotype_path)
    pc_genes <- map_dt[gene_type == "protein_coding", .(phenotype_id)]
    logger::log_info("{nrow(pc_genes)} protein-coding genes")

    # ---- 5. Read gnomAD LOEUF ----
    logger::log_info("Reading gnomAD constraint metrics...")
    gnomad_dt <- data.table::fread(gnomad_path)
    loeuf_dt <- gnomad_dt[canonical == TRUE & mane_select == TRUE &
                           !is.na(gene) & !is.na(lof.oe_ci.upper),
                          .(loeuf = unique(lof.oe_ci.upper)), by = gene]
    data.table::setnames(loeuf_dt, "gene", "phenotype_id")
    logger::log_info("{nrow(loeuf_dt)} genes with LOEUF (MANE Select canonical)")

    # ---- 6. Merge all ----
    dt <- merge(effect_dt, cluster_dt[, .(phenotype_id, cluster)], by = "phenotype_id")
    dt <- merge(dt, pc_genes, by = "phenotype_id")
    dt <- merge(dt, loeuf_dt, by = "phenotype_id")
    logger::log_info("{nrow(dt)} protein-coding genes with LOEUF + cluster + mean sig effect size")

    # ---- 7. Per-cluster Spearman correlations ----
    clusters <- sort(unique(dt$cluster))
    results <- list()

    for (cl in clusters) {
        sub <- dt[cluster == cl]
        sp <- stats::cor.test(sub$loeuf, sub$mean_abs_slope, method = "spearman", exact = FALSE)

        results[[length(results) + 1]] <- data.table::data.table(
            cluster      = cl,
            n_total      = nrow(sub),
            spearman_rho = round(sp$estimate, 4),
            spearman_p   = sp$p.value
        )
    }

    res_dt <- data.table::rbindlist(results)
    res_dt[, spearman_p_BH := stats::p.adjust(spearman_p, method = "BH")]

    logger::log_info("Per-cluster results:\n{paste(capture.output(print(res_dt)), collapse = '\n')}")

    data.table::fwrite(res_dt, regression_output_path, sep = "\t")
    logger::log_info("Saved to: {regression_output_path}")

    # ---- 8. Overall regression: effect size ~ LOEUF + cluster ----
    logger::log_info("Fitting overall linear model: mean_abs_slope ~ loeuf + factor(cluster)...")
    overall_lm <- stats::lm(mean_abs_slope ~ loeuf + factor(cluster), data = dt)
    loeuf_coef <- summary(overall_lm)$coefficients["loeuf", ]
    logger::log_info("LOEUF coefficient = {round(loeuf_coef['Estimate'], 4)}, p = {format(loeuf_coef['Pr(>|t|)'], digits = 3)}")

    overall_reg_dt <- data.table::data.table(
        term       = "loeuf",
        estimate   = loeuf_coef["Estimate"],
        se         = loeuf_coef["Std. Error"],
        t_value    = loeuf_coef["t value"],
        p_value    = loeuf_coef["Pr(>|t|)"],
        n_genes    = nrow(dt),
        n_clusters = length(unique(dt$cluster))
    )
    data.table::fwrite(overall_reg_dt, regression_output_path, sep = "\t")
    logger::log_info("Saved to: {regression_output_path}")

    # ---- 9. Scatterplot ----
    logger::log_info("Generating scatterplot...")

    dt[, cluster_label := paste0("Cluster ", cluster)]
    dt[, cluster_label := factor(cluster_label, levels = paste0("Cluster ", cluster_order))]

    cluster_display <- setNames(paste0("Cluster ", seq_along(cluster_order)),
                                paste0("Cluster ", cluster_order))
    dt[, cluster_label := factor(cluster_display[as.character(cluster_label)],
                                 levels = paste0("Cluster ", seq_along(cluster_order)))]

    anno_dt <- res_dt[, .(cluster,
        label = sprintf("rho == %.3f * ',' ~ italic(p) == '%.2g'", spearman_rho, spearman_p_BH))]
    anno_dt[, cluster_label := paste0("Cluster ", cluster)]
    anno_dt[, cluster_label := factor(cluster_display[as.character(cluster_label)],
        levels = levels(dt$cluster_label))]

    p <- ggplot2::ggplot(dt, ggplot2::aes(x = loeuf, y = mean_abs_slope)) +
        ggplot2::geom_point(size = 0.5, alpha = 0.3, color = "black") +
        ggplot2::geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.6) +
        ggplot2::geom_text(data = anno_dt, ggplot2::aes(label = label),
                  x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3,
                  size = 4, color = "red", parse = TRUE) +
        ggplot2::facet_wrap(~ cluster_label, ncol = ncol, scales = "fixed") +
        ggplot2::labs(x = "LOEUF", y = "Mean |eQTL effect size|") +
        ggplot2::theme_classic(base_size = 11, base_family = "Arial") +
        ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", size = 11))

    ggplot2::ggsave(plot_output_path, plot = p, width = plot_width, height = plot_height, device = "svg")
    logger::log_info("Plot saved to: {plot_output_path}")

    invisible(list(per_cluster = res_dt, overall_regression = overall_reg_dt))
}
