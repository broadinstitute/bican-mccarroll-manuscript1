# We can generate differential expression with the data merged across regions and include a covariate for region,
# or directly extract effect sizes as the interaction of a tested variable with region.
# How do those results compare to each other?  It would be nice if they were similar, and the per-region
# estimates were noisy versions of the merged estimates - but it's possible there are region specific
# effects that are masked by merging.

library (data.table)
library(ComplexHeatmap)
library (ggplot2)
library(circlize)
#library (pals)
library(Polychrome)



merged_input_DE="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age/cell_type/aggregated_results/age_DE_logFC_uncorrected_matrix_simple2.txt"
merged_input_clusters="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age/cell_type/aggregated_results/age_DE_logFC_uncorrected_matrix_simple2.gene_clusters.txt"

region_input_DE="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age/cell_type_region_interaction_absolute_effects/aggregated_results/age_DE_logFC_uncorrected_matrix_simple2.txt"
region_input_clusters="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age/cell_type_region_interaction_absolute_effects/aggregated_results/age_DE_logFC_uncorrected_matrix_simple2.gene_clusters.txt"
k=14 #the number of clusters we used in the original analysis.


# Raw vs Mash
experiment_1_name="Raw"
input_1_DE="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age/cell_type/aggregated_results/age_DE_logFC_uncorrected_matrix_simple2.txt"
input_1_clusters="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age/cell_type/aggregated_results/age_DE_logFC_uncorrected_matrix_simple2.gene_clusters.txt"

experiment_2_name="mash-corrected"
input_2_DE="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age/cell_type/aggregated_results/age_DE_logFC_corrected_matrix_simple2.txt"
input_2_clusters="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age/cell_type/aggregated_results/age_DE_logFC_corrected_matrix_simple2.gene_clusters.txt"
k=16


# experiment_1_name="mash-corrected"
# input_1_DE="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age/cell_type/aggregated_results/age_DE_logFC_corrected_matrix_simple2.txt"
# input_1_clusters="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age/cell_type/aggregated_results/age_DE_logFC_corrected_matrix_simple2.gene_clusters.txt"
#
# experiment_2_name="mash-corrected with null correction"
# input_2_DE="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age/cell_type/aggregated_results/age_DE_logFC_corrected_Vhat_matrix_simple2.txt"
# input_2_clusters="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age/cell_type/aggregated_results/age_DE_logFC_corrected_Vhat_matrix_simple2.gene_clusters.txt"
# k=16

logFCRaw="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age/cell_type/aggregated_results/age_DE_logFC_allgenes_uncorrected_matrix_simple2.txt"
tstatRaw="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age/cell_type/aggregated_results/age_DE_tstat_allgenes_uncorrected_matrix_simple2.txt"
logFCMash="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age/cell_type/aggregated_results/age_DE_logFC_allgenes_corrected_matrix_simple2.txt"
tstatMash="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age/cell_type/aggregated_results/age_DE_tstat_allgenes_corrected_matrix_simple2.txt"



compare_merged_vs_region_de<-function (merged_input_clusters, region_input_clusters, k=14) {
    m=read.table(merged_input_clusters, header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE)
    r=read.table(region_input_clusters, header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE)

    # plot the intersection and disjoint sets of genes
    p=plot_gene_set_overlaps(rownames(m), rownames(r), title="Gene Set Overlaps: Merged vs Region-Specific DE")

    geneClustersMerged=m[[paste0("k", k)]]
    names (geneClustersMerged)=rownames(m)
    geneClustersRegion=r[[paste0("k", k)]]
    names (geneClustersRegion)=rownames(r)

    cluster_mapping=map_clusters_by_gene_overlap(geneClustersMerged, geneClustersRegion)

    plot_clustter_overlap(cluster_mapping$similarity_matrix, cluster_mapping$optimal_mapping)

    # Try plotting the effect sizes ordered so that the clusters line up.
    mDE=read.table(merged_input_DE, header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE, check.names = FALSE)
    rDE=read.table(region_input_DE, header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE, check.names = FALSE)
    commonGenes=intersect(rownames(mDE), rownames(rDE))
    mDE=mDE[commonGenes, ]
    rDE=rDE[commonGenes, ]

    # Reorder columns to match
    rDE=reorder_by_base_conditions(mDE, rDE)

    cluster_order_merged=as.character(cluster_mapping$optimal_mapping$row_factor)
    plot_mash_heatmap(effect_size_matrix=mDE, clusters=geneClustersMerged[commonGenes], cluster_order=cluster_order_merged, column_title="Merged DE: Age Effects")

    cluster_order_region=as.character(cluster_mapping$optimal_mapping$best_match)
    plot_mash_heatmap(effect_size_matrix=rDE, clusters=geneClustersRegion[commonGenes], cluster_order=cluster_order_region, column_title="Region-Specific DE: Age Effects")
}

compare_two_de_clustering_results<-function (experiment_1_name, input_1_DE, input_1_clusters,
                                  experiment_2_name, input_2_DE, input_2_clusters,
                                  k=14) {

    d=read_paired_cluster_data(input_1_DE, input_1_clusters,
                             input_2_DE, input_2_clusters,
                             k)

    # plot the intersection and disjoint sets of genes
    p=plot_gene_set_overlaps(rownames(d$clusters1), rownames(d$clusters2), title="Gene Set Overlap", experiment_1_name, experiment_2_name)

    cluster_mapping=map_clusters_by_gene_overlap(d$geneClusters1, d$geneClusters2)

    plot_clustter_overlap(cluster_mapping$similarity_matrix, cluster_mapping$optimal_mapping, experiment_1_name, experiment_2_name)

    #this is generally cluster 0...end.  It would be interesting to order by specificity instead.
    #for each cluster, compute average DE per condition, then (somehow) count the number of similar conditions.

    cluster_order_1=as.character(cluster_mapping$optimal_mapping$row_factor)
    plot_mash_heatmap(effect_size_matrix=d$de1, clusters=d$geneClusters1[d$commonGenes], cluster_order=cluster_order_1, column_title=paste(experiment_1_name, "DE: Age Effects"))

    cluster_order_2=as.character(cluster_mapping$optimal_mapping$best_match)
    plot_mash_heatmap(effect_size_matrix=d$de2, clusters=d$geneClusters2[d$commonGenes], cluster_order=cluster_order_2, column_title=paste(experiment_2_name, "DE: Age Effects"))


}


#explore how genes shift from one cluster to another after mash is applied.
kMeansclusterShifts<-function (experiment_1_name, input_1_DE, input_1_clusters,
                 experiment_2_name, input_2_DE, input_2_clusters,
                 k=14) {

    d=read_paired_cluster_data(input_1_DE, input_1_clusters,
                               input_2_DE, input_2_clusters,
                               k)

    source=d$geneClusters1
    destination=d$geneClusters2

    clusterSource=0
    condition="microglia"
    genesInCluster=names(source[source==clusterSource])
    dest=destination[genesInCluster]

    barplot(table(dest), xlab=paste("cluster in ", experiment_2_name), main=paste("Genes in cluster ", clusterSource))
    #let's look at the raw and mash effect sizes for these genes.
    df=data.frame(rawLogFC=d$de1[genesInCluster, condition], mashLogFC=d$de2[genesInCluster, condition], destination=factor(dest))

    limx=range (df$rawLogFC, df$mashLogFC, na.rm=TRUE)
    corVal=cor(df$rawLogFC, df$mashLogFC, use="complete.obs")

    #load up some high-contrast colors.
    #nColsNeeded=length(unique (df$destination))
    #data(glasbey)
    #cols=glasbey[1:nColsNeeded]
    #names (cols)=levels(df$destination)

    #get the list of colors for all labels
    cluster_mapping=map_clusters_by_gene_overlap(d$geneClusters1, d$geneClusters2)
    clusterLabels=cluster_mapping$optimal_mapping$best_match
    cols=category_palette(clusterLabels)

    ggplot(df, aes(x = rawLogFC, y = mashLogFC, color = destination)) +
        geom_point(alpha = 1, size = 2.5) +
        scale_color_manual(values = cols) +
        geom_abline(slope = 1, intercept = 0, color = "gray50", linetype = "dotted") +
        geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
        labs(
            x = "Raw effect size",
            y = "Mash posterior effect size",
            color = "Destination",
            title = paste0(condition, "  (r = ", round(corVal, 3), ")")
        ) +
        theme_minimal(base_size = 14) +
        theme(
            plot.title = element_text(hjust = 0.5),
            panel.grid.minor = element_blank()
        ) +
        xlim(limx) + ylim(limx)




}

adhocGeneOrder<-function () {
    rtstat=read.table(tstatRaw, header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE, check.names = FALSE)
    mtstat=read.table(tstatMash, header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE, check.names = FALSE)

    all(rownames (rtstat)==rownames(mtstat))
    all (colnames(rtstat)==colnames(mtstat))

    rLogFC=read.table(logFCRaw, header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE, check.names = FALSE)
    mLogFC=read.table(logFCMash, header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE, check.names = FALSE)

    all(rownames (rLogFC)==rownames(mLogFC))
    all (colnames(rLogFC)==colnames(mLogFC))
    all (rownames(rLogFC)==rownames(rtstat))

    condition="microglia"

    df=data.frame(rawLogFC=rLogFC[[condition]], mashLogFC=mLogFC[[condition]], rawTstat=rtstat[[condition]], mashTstat=mtstat[[condition]])
    rownames(df)=rownames(rLogFC)
    df$rawSe=abs(df$rawLogFC  / df$rawTstat)
    df$mashSe=abs(df$mashLogFC / df$mashTstat)


    corVal=cor(df$rawLogFC, df$mashLogFC, use="complete.obs")
    ggplot(df, aes(x = rawLogFC, y = mashLogFC)) +
        geom_point(alpha = 0.4, size = 1.2) +
        geom_abline(slope = 1, intercept = 0, color = "gray50", linetype = "dotted") +
        geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
        labs(
            x = "Raw effect size",
            y = "Mash posterior effect size",
            title = paste0(condition, "  (r = ", round(corVal, 3), ")")
        ) +
        theme_minimal(base_size = 14) +
        theme(
            plot.title = element_text(hjust = 0.5),
            panel.grid.minor = element_blank()
        )

    corVal=cor(df$rawTstat, df$mashTstat, use="complete.obs")
    ggplot(df, aes(x = rawTstat, y = mashTstat)) +
        geom_point(alpha = 0.4, size = 1.2) +
        geom_abline(slope = 1, intercept = 0, color = "gray50", linetype = "dotted") +
        geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
        labs(
            x = "Raw t-stat",
            y = "Mash t-stat",
            title = paste0(condition, "  (r = ", round(corVal, 3), ")")
        ) +
        theme_minimal(base_size = 14) +
        theme(
            plot.title = element_text(hjust = 0.5),
            panel.grid.minor = element_blank()
        )

    #pick out a few exemplar genes that have high mash t-stat and low raw t-stat.
    round (df[df$rawTstat<3 & df$mashTstat>5,],4)

    #plot_es_compare(df, condition, corVal, cap_dt = 6)

}


plot_es_compare <- function(df, condition, corVal, cap_dt = 6) {
    # ---- derived quantities ----
    df$se_raw  <- abs(df$rawLogFC  / df$rawTstat)
    df$se_mash <- abs(df$mashLogFC / df$mashTstat)
    df$d_t     <- df$mashTstat - df$rawTstat
    df$d_se_l2 <- log2(df$se_mash / df$se_raw)   # <0 => SE shrank
    df$size_dt <- pmin(abs(df$d_t), cap_dt)

    # ---- plot ----
    ggplot(df, aes(x = rawLogFC, y = mashLogFC)) +
        geom_point(
            aes(color = d_se_l2,
                size  = size_dt),
            alpha = 0.75
        ) +
        geom_abline(slope = 1, intercept = 0,
                    color = "gray50", linetype = "dotted", linewidth = 0.6) +
        geom_smooth(method = "lm", se = FALSE,
                    color = "red", linetype = "dashed", linewidth = 0.7,
                    inherit.aes = FALSE,
                    mapping = aes(x = df$rawLogFC, y = df$mashLogFC)) +
        scale_color_gradient2(
            low = "#3b4cc0", mid = "white", high = "#b40426",
            midpoint = 0, name = "log2(SE_mash/SE_raw)",
            oob = squish
        ) +
        scale_size(range = c(0.6, 2.2), name = "|Δt| (capped)") +
        labs(
            x = "Raw effect size",
            y = "Mash posterior effect size",
            title = paste0(condition, "  (r = ", round(corVal, 3), ")")
        ) +
        theme_minimal(base_size = 14) +
        theme(plot.title = element_text(hjust = 0.5),
              panel.grid.minor = element_blank())
}


read_paired_cluster_data<-function (input_1_DE, input_1_clusters,
                                    input_2_DE, input_2_clusters,
                                    k=14) {

    clusters1=read.table(input_1_clusters, header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE)
    clusters2=read.table(input_2_clusters, header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE)

    geneClusters1=clusters1[[paste0("k", k)]]
    names (geneClusters1)=rownames(clusters1)
    geneClusters2=clusters2[[paste0("k", k)]]
    names (geneClusters2)=rownames(clusters2)

    # Try plotting the effect sizes ordered so that the clusters line up.
    de1=read.table(input_1_DE, header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE, check.names = FALSE)
    de2=read.table(input_2_DE, header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE, check.names = FALSE)
    commonGenes=intersect(rownames(de1), rownames(de2))
    de1=de1[commonGenes, ]
    de2=de2[commonGenes, ]

    # Reorder columns to match
    commonConditions=intersect(colnames(de1), colnames(de2))
    de1=de1[, commonConditions, drop=FALSE]
    de2=de2[, commonConditions, drop=FALSE]

    result=list(commonGenes=commonGenes, clusters1=clusters1, clusters2=clusters2,
                geneClusters1=geneClusters1, geneClusters2=geneClusters2,
                de1=de1, de2=de2)

    return (result)

}

#effect_size_matrix=mDE; clusters=geneClustersMerged[commonGenes]; column_title="Merged DE: Age Effects"
plot_mash_heatmap <- function(effect_size_matrix,
                              clusters,
                              cluster_order,
                              column_title = NULL,
                              sort_within_cluster = FALSE) {
    if (is.null(colnames(effect_size_matrix)) || is.null(rownames(effect_size_matrix)))
        stop("effect_size_matrix needs rownames (genes) and colnames (cell types).")
    if (is.null(names(clusters)))
        stop("'clusters' must be named by gene IDs.")

    mat <- t(effect_size_matrix)

    gene_ids <- colnames(mat)
    gene_clusters <- as.character(clusters[gene_ids])
    if (anyNA(gene_clusters)) stop("Some genes in effect_size_matrix are missing from 'clusters'.")

    # enforce split order using provided cluster_order
    present_levels <- intersect(cluster_order, unique(gene_clusters))
    if (!length(present_levels)) stop("No overlap between cluster_order and clusters.")
    col_split <- factor(gene_clusters, levels = present_levels)

    # optional: reorder columns within each split by effect strength
    if (isTRUE(sort_within_cluster)) {
        strength <- colMeans(abs(mat), na.rm = TRUE)
        ord <- order(col_split, -strength)
        mat <- mat[, ord, drop = FALSE]
        col_split <- col_split[ord]
    } else {
        # keep current column order from effect_size_matrix, only set splits
        mat <- mat[, seq_len(ncol(mat)), drop = FALSE]
    }

    # shared categorical colors for cluster labels
    ann_cols <- category_palette(levels(col_split))

    # bottom labels per cluster on X axis
    bottom_anno <- HeatmapAnnotation(
        cluster = col_split,
        cluster_label = anno_block(
            labels = levels(col_split),
            labels_gp = gpar(fontsize = 10),
            gp = gpar(fill = NA, col = NA)
        ),
        col = list(cluster = ann_cols),
        show_annotation_name = FALSE,
        show_legend = FALSE,
        which = "column"
    )

    # blue–white–red scale
    rng <- range(mat, na.rm = TRUE)
    brks <- c(seq(rng[1], 0, length.out = 5), seq(0, rng[2], length.out = 5)[-1])
    col_fun <- circlize::colorRamp2(
        brks,
        c("#00008B","#4040FF","#8080FF","#C0C0FF","#FFFFFF",
          "#FFC0C0","#FF8080","#FF4040","#8B0000")
    )

    ht <- Heatmap(
        mat,
        name = "Effect",
        column_split = col_split,
        bottom_annotation = bottom_anno,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        show_column_dend = FALSE,
        show_column_names = FALSE,
        col = col_fun,
        column_title = column_title,
        row_title = "Cell types",
        gap = unit(1.5, "mm"),
        use_raster = TRUE,
        heatmap_legend_param = list(color_bar = "continuous",
                                    legend_height = unit(5, "cm"))
    )

    draw(ht)
}

# ---- helper: default categorical palette (edit here if you want different defaults) ----
category_palette <- function(levels_vec) {
    # load curated Polychrome Glasbey
    data("glasbey", package = "Polychrome", envir = environment())
    base_cols <- get("glasbey", envir = environment())

    # drop unreadable colors (whites / pale yellows / very light pastels)
    bad <- c("#FFFF00", "#FFD700", "#FFFFFF", "#F0F0F0")
    base_cols <- base_cols[!toupper(base_cols) %in% bad]

    if (length(base_cols) < length(levels_vec)) {
        stop("Not enough distinct colors after filtering")
    }

    cols <- base_cols[seq_along(levels_vec)]
    names(cols) <- as.character(levels_vec)
    cols
}

#sim=cluster_mapping$similarity_matrix; mapping=cluster_mapping$optimal_mapping
plot_clustter_overlap<- function(sim, mapping, experiment_1_name="Merged", experiment_2_name="Region-Specific") {

    # sim: square similarity matrix with row/colnames as string cluster labels
    # mapping: data.frame with columns row_factor, best_match

    rf <- as.character(mapping$row_factor)
    bm <- as.character(mapping$best_match)

    ro <- match(rf, rownames(sim))
    co <- match(bm, colnames(sim))
    if (any(is.na(ro)) || any(is.na(co))) {
        stop("Mapping labels not found in matrix dimnames.")
    }

    rng <- range(sim, na.rm = TRUE)
    col_fun <- colorRamp2(
        seq(rng[1], rng[2], length.out = 3),
        c("white", "steelblue3", "navy")
    )

    Heatmap(
        sim[ro, co, drop = FALSE],
        name = "similarity",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_labels = rf,
        column_labels = bm,
        row_names_side = "left",
        column_names_side = "top",
        col = col_fun,
        row_title = experiment_1_name,
        column_title = experiment_2_name,
    )
}

#genesOne=rownames(m); genesTwo=rownames(r)
plot_gene_set_overlaps<-function (genesOne, genesTwo, title="Gene Set Overlaps",
                                  experiment_1_name="Merged",
                                  experiment_2_name="Region-Specific") {
    g1=setdiff(genesOne, genesTwo)
    g2=setdiff(genesTwo, genesOne)
    gboth=intersect(genesOne, genesTwo)

    df=data.frame(
        category=c(paste("Only", experiment_1_name), paste("Only", experiment_2_name), "Both"),
        count=c(length(g1), length(g2), length(gboth))
    )
    p<-ggplot2::ggplot(df, ggplot2::aes(x=category, y=count, fill=category)) +
        ggplot2::geom_bar(stat="identity") +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none") +
        ggplot2::labs(title=title, x="", y="Number of Genes")
    return(p)
}

#clusters1=geneClustersMerged; clusters2=geneClustersRegion
map_clusters_by_gene_overlap <- function(clusters1, clusters2) {
    # Check inputs are named vectors
    if (is.null(names(clusters1)) || is.null(names(clusters2))) {
        stop("Both vectors must be named with gene symbols.")
    }

    # Find genes in common
    shared_genes <- intersect(names(clusters1), names(clusters2))
    if (length(shared_genes) == 0) {
        stop("No overlapping genes found.")
    }

    # Subset to shared genes
    clusters1 <- clusters1[shared_genes]
    clusters2 <- clusters2[shared_genes]

    # Create cluster → gene sets
    cluster_genes_1 <- split(names(clusters1), clusters1)
    cluster_genes_2 <- split(names(clusters2), clusters2)

    labels1 <- names(cluster_genes_1)
    labels2 <- names(cluster_genes_2)

    # Initialize similarity matrix
    sim_mat <- matrix(0, nrow = length(labels1), ncol = length(labels2),
                      dimnames = list(labels1, labels2))

    # Compute Jaccard index for each cluster pair
    for (i in labels1) {
        for (j in labels2) {
            g1 <- cluster_genes_1[[i]]
            g2 <- cluster_genes_2[[j]]
            sim_mat[i, j] <- length(intersect(g1, g2)) / length(union(g1, g2))
        }
    }

    # Solve assignment
    optimal_mapping=get_optimal_factor_mapping(sim_mat)
    return (list(similarity_matrix=sim_mat, optimal_mapping=optimal_mapping))
}

get_optimal_factor_mapping <- function(cor_mat) {
    if (!requireNamespace("clue", quietly = TRUE)) {
        stop("Install 'clue' with install.packages('clue')")
    }

    row_names <- rownames(cor_mat)
    col_names <- colnames(cor_mat)
    n_rows <- length(row_names)
    n_cols <- length(col_names)
    n <- max(n_rows, n_cols)

    # Convert correlation to cost
    max_corr <- max(cor_mat, na.rm = TRUE)
    cost_mat <- max_corr - cor_mat

    # Pad to square matrix
    padded <- matrix(max(cost_mat, na.rm = TRUE) + 1, n, n)
    padded[1:n_rows, 1:n_cols] <- cost_mat

    # Solve
    assignment <- clue::solve_LSAP(padded)

    # Filter: only real assignments
    valid_rows <- seq_len(n_rows)
    valid_cols <- assignment[valid_rows]
    valid <- valid_cols <= n_cols

    data.frame(
        row_factor = row_names[valid],
        best_match = col_names[valid_cols[valid]],
        correlation = mapply(function(i, j) cor_mat[i, j], row_names[valid], col_names[valid_cols[valid]]),
        stringsAsFactors = FALSE
    )
}

# mDE: matrix with "base" condition names
# rDE: matrix with suffixed condition names
# return: rDE reordered with grouping that matches mDE columns
reorder_by_base_conditions <- function(mDE, rDE) {
    # base names from mDE
    base_names <- colnames(mDE)

    # strip suffixes from rDE
    rDE_base <- sub("_[^_]+$", "", colnames(rDE))

    # build ordering
    ord <- unlist(lapply(base_names, function(b) {
        which(rDE_base == b)
    }))

    rDE[, ord, drop = FALSE]
}

# Example:
# rDE_aligned <- reorder_by_base_conditions(mDE, rDE)

