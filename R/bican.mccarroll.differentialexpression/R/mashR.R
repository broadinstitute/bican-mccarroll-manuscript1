# Can we use mashR to get better joint estimates of age and sex effects by stabilizing results across
# multiple cell types that have different levels of power?

# packages
# library(mashr)
# library (data.table)
# library (profmem) #if fit_mash_with_V add.mem.profile=TRUE
# library(ComplexHeatmap)
# library(circlize)
# library(cluster)
# library (ggplot2)
#
# in_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/cell_type_results_sex_age"
# sample_metadata_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metacells/donor_rxn_DGEList_samples.tsv.gz"
# file_pattern="age"
#
# gene_cluster_file_raw="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/mash/metadata/age_DE_logFC_K16_gene_clusters.csv"
# gene_cluster_file_post_mash="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/mash/metadata/age_DE_mash_corrected_effects_7539genes_K17_gene_clusters.csv"
# gene_cluster_labels_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/mash/metadata/age_DE_logFC_K16_gene_clusters_labels.txt"


run_mashr<-function (in_dir, file_pattern="age") {

    #get the list of files from a directory that match some pattern.
    f=list.files(in_dir, pattern = paste0(file_pattern), full.names = TRUE)
    d=lapply(f, read.table, sep="\t", header=TRUE,
             stringsAsFactors = FALSE, check.names = FALSE)
    #need to infer the names for each data.frame in the list
    names(d)=drop_common_suffix(basename(f))

    #make mash inputs
    mash_inputs_union<-make_mash_inputs(d, coef_col = "logFC", t_col = "t", fdr_col = "adj.P.Val", gene_mode="union")
    #how many genes pass an FDR < 0.05 in any condition?
    nGenesFDR=length(which(apply(mash_inputs_union$FDR, 1, function(x) any(x<0.05))))
    logger::log_info("Number of genes with at least one condition that passes FDR [", nGenesFDR, "]")

    idxPassFDRAny<-which(apply(mash_inputs_union$FDR, 1, function (x) any(x<0.01)))
    genesPassFDR<-rownames(mash_inputs_union$Bhat[idxPassFDRAny,])

    #TODO: Maybe what I want to do is fit the data-determined covariance matrixes to the intersection
    # data

    # Fit mash to the data.
    mash_fit<-fit_mash(Bhat=mash_inputs_union$Bhat, Shat=mash_inputs_union$Shat, FDR=mash_inputs_union$FDR, missing_mask=mash_inputs_union$missing_mask, npc=5)
    m<-mash_fit$mash
    covariance_matrix_list<-mash_fit$covariance_matrix_list

    #the learned gene covariance matrixes are the interesting part.
    learned_matrix_names=names (covariance_matrix_list)[grep ("ED_", names (covariance_matrix_list))]

    mixing_proportions=get_estimated_pi(m)
    barplot(mixing_proportions, las=2, cex.names=0.7, main="Mixing proportions for mash components")

    for (component in learned_matrix_names) {
        p1=plot_mash_component(component = component, covariance_matrix_list, mixing_proportions=mixing_proportions, cex=1)
        p2=plot_top_eigenvectors(covariance_matrix_list[[component]], k=1)[[1]]
        p=cowplot::plot_grid(p1, p2, ncol=1)
        if (mixing_proportions[component]>0)
            print (p)
    }


    ################
    # RESULTS
    ################
    #local false sign rate - how confident mash is in the direction of effect for each celltype x gene.
    #https://academic.oup.com/biostatistics/article-abstract/18/2/275/2557030?redirectedFrom=fulltext
    lfsr<-get_lfsr(m)
    pm<-get_pm(m)
    sds<-get_psd(m)

    #Filter results to genes that pass an FDR threshold in at least one condition
    m_fdr=filter_mash_result_by_genes(m, genesPassFDR)

    result <- cluster_mash_effects(m_fdr, lfsr_threshold = 1, k = 16, method="kmeans", dist_metric = "euclidean", nstart=200, seed=2)
    plot_mash_heatmap(result$pm_filtered, result$clusters, column_title="distance [euclidean] clustering [kmeans] (lfsr=1)")

    result <- cluster_mash_effects(m_fdr, lfsr_threshold = 0.05, k = 16, method="kmeans", dist_metric = "euclidean", nstart=200, seed=2)
    plot_mash_heatmap(result$pm_filtered, result$clusters, column_title="distance [euclidean] clustering [kmeans] (lfsr<0.05")
    #plot_mash_heatmap(pm[rownames (result$pm_filtered),], result$clusters, column_title="distance [euclidean] clustering [kmeans] unfiltered data")












    #TODO: Plot the minimum local false sign rate per gene against the minimum FDR from DE.
    # min_lfsr<-apply(lfsr, 1, min, na.rm=T)
    # min_fdr<-apply(FDR, 1, min, na.rm=T)
    # plot (-log10(min_fdr),-log10(min_lfsr))
    # plot (-log10(min_fdr), -log10(min_q))


    #ESTIMATE THE Directional FDR OF THE DATA SET
    # Called sets
    pos_idx <- which(lfsr < 0.05 & pm > 0, arr.ind = TRUE)
    neg_idx <- which(lfsr < 0.05 & pm < 0, arr.ind = TRUE)

    # Expected errors and directional FDR among calls
    exp_errors <- sum(lfsr[pos_idx], na.rm = TRUE) + sum(lfsr[neg_idx], na.rm = TRUE)
    n_calls    <- nrow(pos_idx) + nrow(neg_idx)
    dir_fdr    <- if (n_calls > 0) exp_errors / n_calls else NA_real_

    #write the outputs
    # idx=sort(unique(c(pos_idx,neg_idx)))
    # write.table(pm, "/broad/mccarroll/haley/BICAN/DE/age_DE_logFC_mash_corrected_effects_matrix.txt", row.names=T, col.names = T, quote=F, sep="\t")
    # write.table(pm[idx,], "/broad/mccarroll/haley/BICAN/DE/age_DE_logFC_mash_corrected_effects_filtered_matrix.txt", row.names=T, col.names = T, quote=F, sep="\t")
    #
    # pos_idx <- which(lfsr < 0.01 & pm > 0, arr.ind = TRUE)
    # neg_idx <- which(lfsr < 0.01 & pm < 0, arr.ind = TRUE)
    # idx=sort(unique(c(pos_idx,neg_idx)))
    # write.table(pm[idx,], "/broad/mccarroll/haley/BICAN/DE/age_DE_logFC_mash_corrected_effects_filtered_strict_matrix.txt", row.names=T, col.names = T, quote=F, sep="\t")
    # write.table(pm[clusters$gene,], "/broad/mccarroll/haley/BICAN/DE/age_DE_logFC_mash_corrected_effects_filtered_Haley_matrix.txt", row.names=T, col.names = T, quote=F, sep="\t")

    #How much sharing is there between cell types?
    #the factor is how close the results are , where factor=0 indicates they share the same direction of effect.
    #This is a fast way to get an estimate of how much sharing there is by cell type to see structure
    sharing=get_pairwise_sharing(m, factor=0.5)
    plot_sharing (sharing, strTitle=paste("Pairwise sharing by magnitude (factor=0.5) for", file_pattern, "effects"))

    #TODO: can I look at the eigenvectors of each covariance matrix


    #What covariance structures are being used?
    #This is the proportions of each covariance structure that are learned BEFORE the individual gene
    #posterior means/variances are computed.  These are the priors used during that fitting stage.
    mixing_proportions=get_estimated_pi(m)
    barplot(mixing_proportions, las=2, cex.names=0.7, main="Mixing proportions for mash components")
    heatmap(covariance_matrix_list[["ED_PCA_1"]])
    #plot_cov_image(covariance_matrix_list[["ED_PCA_1"]])
    plot_mash_component(component = "ED_PCA_1", covariance_matrix_list, mixing_proportions=mixing_proportions)

    #lapply(names(covariance_matrix_list), function (x) plot_mash_component(x, covariance_matrix_list, mixing_proportions=mixing_proportions))
    # Partition mash calls by component.
    # only keep genes significant in at least one condition
    res <- partition_mash_calls_by_component(m, thresh = 0.05, keep_null = T)
    heatmap(t(res$weights_collapsed), scale="column")

    z=table (res$assignment$label)
    z=z/sum(z)
    barplot(z, las=2, cex.names=0.7, main="Fraction of significant genes by mash component")

    #are there any genes that are shared?  Look for a lower max weight for a gene
    gene_weight_max= apply(res$weights_collapsed, 1, max, na.rm = TRUE)
    which.min(gene_weight_max)
    mash_forest_with_labels(m,  gene="ENSG00000254733", order_by = "effect", call_rule="both")
    mash_forest_with_labels(m,  gene="A2M", order_by = "effect", call_rule="metric", plot_raw=T)
    mash_forest_with_labels(m,  gene="MYRIP", order_by = "effect", call_rule="metric", plot_raw=T)

    #how many genes are assigned to each mash component?
    table(res$assignment$label) # also : res$counts

    #APOE is assigned to both SPNs and microglia and astocytes.  Complicated.
    par(mar=c(10,4,3,1))
    barplot(res$weights_collapsed["APOE",], las=2, main="APOE weights by mash component")
    plot_mash_component(component = "ED_PCA_1", covariance_matrix_list, mixing_proportions=mixing_proportions)
    plot_top_eigenvectors(covariance_matrix_list[["ED_PCA_1"]], k=1)
    mash_forest_with_labels(m,  gene="APOE", order_by = "effect", call_rule="metric", plot_raw=T)


    plot_mash_component(component = "ED_PCA_2", covariance_matrix_list, mixing_proportions=mixing_proportions)
    plot_top_eigenvectors(covariance_matrix_list[["ED_PCA_2"]], k=1)
    plot_mash_component(component = "ED_PCA_3", covariance_matrix_list, mixing_proportions=mixing_proportions)
    plot_top_eigenvectors(covariance_matrix_list[["ED_PCA_3"]], k=3)

    plot_mash_component(component = "ED_tPCA", covariance_matrix_list, mixing_proportions=mixing_proportions)
    plot_top_eigenvectors(covariance_matrix_list[["ED_tPCA"]], k=5)



    #interesting, ED_PCA_3 looks like microglia and astrocytes, but independent?
    plot_cov_image(covariance_matrix_list[["ED_PCA_3"]])
    head (res$assignment[res$assignment$label=="ED_PCA_3",])
    mash_forest_with_labels(m,  gene="ASRGL1", order_by = "effect", call_rule="metric", plot_raw=T)
    mash_forest_with_labels(m,  gene="ANOS1", order_by = "effect", call_rule="metric", plot_raw=T)
    mash_forest_with_labels(m,  gene="ARRB2", order_by = "effect", call_rule="metric", plot_raw=T)
    mash_forest_with_labels(m,  gene="ASRGL1", order_by = "effect", call_rule="metric", plot_raw=T)

    #TODO: For genes that were originally not ascertained, map their NA results back into
    #the raw data so they aren't plotted.
    # this gene is flagged as from PCA_ED_3, but only the microglia effect is significant.
    mash_forest_with_labels(m,  gene="TGFBI", order_by = "effect", call_rule="metric", plot_raw=T)
    res$weights_collapsed["TGFBI",]
    pm["TGFBI",]
    lfsr["TGFBI",]
    mash_forest_with_labels(m,  gene="TGFBI", order_by = "effect", call_rule="metric", plot_raw=F)

    # Per-gene table (hard assignment, with posterior weight)
    head(res$assignment)

    # If you want a barplot of fractions (same labels as get_estimated_pi):
    op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
    par(mar = c(8, 4, 3, 1))
    barplot(res$fractions,
            las = 2, ylab = "% of significant genes",
            main = "Significant genes partitioned by mash component")


    #try sampling the posterior
    #https://stephenslab.github.io/mashr/articles/mash_sampling.html
    #I will need to fit the model a different way to do this!
    #x = get_pairwise_sharing_from_samples(m, factor=0.5, lfsr_thresh = 1)
    #corrplot(x, method='color', col.lim=c(0,1), type='upper', addCoef.col = "black", tl.col="black", tl.srt=45, title = 'Pairwise Sharing by Magnitude', mar = c(4,0,4,0))


    sig_results<-get_significant_results(m, thresh = 0.05, conditions = NULL, sig_fn = get_lfsr)

    mash_forest_with_labels(m,  gene=names(sig_results[1]), order_by = "effect", call_rule="both")

    #how are significant effects shared across cell types?
    # Draw many samples for each gene using the original data the all the patterns
    # from the covariance matrixes that are blended to get new means / sd many times
    # These can be used to compute cell type correlations of effects that respects
    # the mash covariance structure detected.
    logger::log_info("Computing posterior matrices for mashr samples...")
    m$result <- mashr::mash_compute_posterior_matrices(
        m, data,
        algorithm.version   = "R",
        posterior_samples   = 100,     # increase for smoother estimates
        output_posterior_cov = FALSE
    )
    logger::log_info("Done computing posterior matrices.")


    ## 3) compute pairwise sharing (by magnitude)
    M <- mashr::get_pairwise_sharing_from_samples(
        m,
        factor = 0.5,        # within 2× in magnitude
        lfsr_thresh = 0.05   # include genes active in at least one condition
    )

    ## 4) plot (Fig. 6 style)
    corrplot::corrplot(
        M, method = "color", type = "upper",
        col.lim = c(0.3, 1), addCoef.col = "black",
        tl.col = "black", tl.srt = 45,
        title = "Pairwise sharing by magnitude", mar = c(4, 0, 4, 0)
    )


    #TODO: If I were to just fit the data to the cannoical covariances, what would happen?
    #for genes that were fit to astrocytes and microglia (ED_PCA_3), would they be assigned equally to glia and astrocytes.

    mash_fit_can<-fit_mash(Bhat, Shat, FDR=FDR, missing_mask=missing_mask, only_cannonical_matrixes=T)
    m2<-mash_fit_can$mash
    covariance_matrix_list2<-mash_fit$covariance_matrix_list
    res2 <- partition_mash_calls_by_component(m2, thresh = 0.05, keep_null = T)

    plot_cov_image(covariance_matrix_list[["ED_PCA_3"]])
    head (res$assignment[res$assignment$label=="ED_PCA_3",])
    mash_forest_with_labels(m,  gene="ASRGL1", order_by = "effect", call_rule="metric", plot_raw=T)
    mash_forest_with_labels(m,  gene="ANOS1", order_by = "effect", call_rule="metric", plot_raw=T)
    mash_forest_with_labels(m,  gene="ARRB2", order_by = "effect", call_rule="metric", plot_raw=T)
    mash_forest_with_labels(m,  gene="ASRGL1", order_by = "effect", call_rule="metric", plot_raw=T)

    mash_forest_with_labels(m2,  gene="ASRGL1", order_by = "effect", call_rule="metric", plot_raw=T)
    mash_forest_with_labels(m,  gene="ANOS1", order_by = "effect", call_rule="metric", plot_raw=T)
    mash_forest_with_labels(m,  gene="ARRB2", order_by = "effect", call_rule="metric", plot_raw=T)
    mash_forest_with_labels(m,  gene="ASRGL1", order_by = "effect", call_rule="metric", plot_raw=T)

    barplot(res$weights_collapsed["ASRGL1",], las=2)
    barplot(res2$weights_collapsed["ASRGL1",], las=2)



    #TODO: look at haley's original microglia gene cluster pre-mash, and
    #compare to the mash result genes assigned to microglia.
    #need to more carefully filter the mash results.
    #res$assignment$label
    # maybe a confusion matrix of haley cluster vs mash component?
}




##########################
# CLUSTER GENES POST-MASH
##########################


cluster_mash_effects <- function(m,
                                 lfsr_threshold = 0.05,
                                 method = c("kmeans", "hclust"),
                                 dist_metric = c("euclidean", "correlation"),
                                 k = 10,
                                 n_pcs_plot = 4,
                                 plot = TRUE,
                                 nstart = 100,
                                 seed = 1,
                                 test_signs_only = FALSE,
                                 use_raw_values = FALSE) {

    # Reproducibility
    set.seed(seed)

    method <- match.arg(method)
    dist_metric <- match.arg(dist_metric)

    if (use_raw_values==FALSE) {
        lfsr <- ashr::get_lfsr(m)
        pm <- ashr::get_pm(m)

        # Zero-out non-significant effects
        pm[lfsr > lfsr_threshold] <- 0

        # Drop genes with no effects
        keep <- rowSums(pm != 0) > 0
        pm_filtered <- pm[keep, , drop = FALSE]
        logger::log_info("Filtered to {nrow(pm_filtered)} genes out of {nrow(pm)} with lfsr < {lfsr_threshold} in at least one condition.")
    } else {
        pm_filtered=m$data$Bhat
        logger::log_info("Using raw effect estimates {nrow(pm_filtered)} genes without any filtering.")
    }

    gene_names <- rownames(pm_filtered)



    if (nrow(pm_filtered) < k) {
        stop("Fewer genes remaining after filtering than number of clusters requested.")
    }

    # Optionally trinarize the data
    if (test_signs_only) {
        tri_mat <- matrix(0, nrow = nrow(pm_filtered), ncol = ncol(pm_filtered),
                          dimnames = dimnames(pm_filtered))
        tri_mat[pm_filtered > 0] <- 1
        tri_mat[pm_filtered < 0] <- -1
        data_for_clustering <- tri_mat
    } else {
        data_for_clustering <- pm_filtered
    }

    # Distance matrix
    if (dist_metric == "euclidean") {
        d <- dist(data_for_clustering)
    } else {
        cor_matrix <- cor(t(data_for_clustering), use = "pairwise.complete.obs")
        d <- as.dist(1 - cor_matrix)
    }

    # Clustering
    if (method == "kmeans") {
        km <- kmeans(data_for_clustering, centers = k, nstart = nstart, iter.max = 100)
        clusters <- km$cluster
    } else {
        hc <- hclust(d)
        clusters <- cutree(hc, k = k)
    }

    # PCA (always on the real values, not trinarized matrix)
    pca <- prcomp(pm_filtered, scale. = TRUE)
    pca_df <- as.data.frame(pca$x)
    pca_df$cluster <- factor(clusters)
    pca_df$gene <- gene_names

    # Optional PCA plots
    plots <- NULL
    if (plot) {
        if (!requireNamespace("ggplot2", quietly = TRUE)) {
            warning("ggplot2 not available for plotting.")
        } else {
            plots <- list()
            np <- min(n_pcs_plot, ncol(pca$x) %/% 2 * 2)
            for (i in seq(1, np, by = 2)) {
                x <- paste0("PC", i)
                y <- paste0("PC", i + 1)
                p <- ggplot2::ggplot(pca_df, ggplot2::aes_string(x = x, y = y, color = "cluster")) +
                    ggplot2::geom_point(alpha = 0.8, size = 2) +
                    ggplot2::labs(title = paste("PCA:", x, "vs", y)) +
                    ggplot2::theme_classic()
                plots[[paste0(x, "_vs_", y)]] <- p
            }
        }
    }

    return(list(
        clusters = structure(clusters, names = gene_names),
        pca = pca,
        pca_df = pca_df,
        plots = plots,
        pm_filtered = pm_filtered,
        dist = d
    ))
}


find_optimal_k_for_mash_effects <- function(m,
                                            start_num_clusters = 2,
                                            end_num_clusters = 20,
                                            lfsr_threshold = 0.05,
                                            dist_metric = c("euclidean", "correlation"),
                                            method = "kmeans",
                                            nstart = 100,
                                            seed = 1,
                                            verbose = TRUE) {
    if (!requireNamespace("cluster", quietly = TRUE)) {
        stop("Install the 'cluster' package to compute silhouette scores.")
    }

    dist_metric <- match.arg(dist_metric)

    results <- list()

    for (k in start_num_clusters:end_num_clusters) {
        res <- tryCatch({
            cluster_result <- cluster_mash_effects(
                m = m,
                lfsr_threshold = lfsr_threshold,
                method = method,
                dist_metric = dist_metric,
                k = k,
                plot = FALSE,
                nstart = nstart,
                seed = seed
            )

            # Validate distance matrix exists
            if (is.null(cluster_result$dist)) {
                stop("Distance matrix not returned by cluster_mash_effects().")
            }

            sil <- cluster::silhouette(cluster_result$clusters, cluster_result$dist)
            mean_sil <- mean(sil[, "sil_width"])

            if (verbose) {
                message("k = ", k, ", silhouette = ", round(mean_sil, 4))
            }

            list(k = k, silhouette = mean_sil, result = cluster_result)

        }, error = function(e) {
            if (verbose) message("k = ", k, " failed: ", e$message)
            NULL
        })

        if (!is.null(res)) {
            results[[as.character(k)]] <- res
        }
    }

    if (length(results) == 0) {
        stop("No clustering results succeeded for the specified k range.")
    }

    # Select best based on silhouette score
    best_k <- names(which.max(sapply(results, function(x) x$silhouette)))
    best_result <- results[[best_k]]

    if (verbose) {
        message("Best k = ", best_result$k, " with silhouette = ", round(best_result$silhouette, 4))
    }

    return(list(
        best_k = best_result$k,
        best_silhouette = best_result$silhouette,
        best_result = best_result$result,
        all_results = results
    ))
}


plot_mash_heatmap <- function(pm_filtered,
                              clusters,
                              column_title = NULL) {

    #to get rid of annoying warning message about rasterization
    ht_opt$message = FALSE

    # Transpose to: cell types = rows, genes = columns
    mat <- t(pm_filtered)

    # Step 1: Cluster the cluster labels
    cluster_ids <- unique(clusters)
    cluster_centroids <- sapply(cluster_ids, function(cid) {
        genes_in_cluster <- names(clusters)[clusters == cid]
        if (length(genes_in_cluster) == 1) {
            pm_filtered[genes_in_cluster, ]
        } else {
            colMeans(pm_filtered[genes_in_cluster, , drop = FALSE])
        }
    })
    cluster_centroids <- t(cluster_centroids)

    # Cluster cluster centroids
    dist_mat <- dist(cluster_centroids)
    cluster_order <- hclust(dist_mat)$order
    reordered_cluster_ids <- cluster_ids[cluster_order]

    # Set col_split with reordered levels
    col_split <- factor(clusters[colnames(mat)], levels = reordered_cluster_ids)

    # Colors for each cluster
    n_clusters <- length(reordered_cluster_ids)
    cluster_colors <- setNames(
        colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n_clusters),
        reordered_cluster_ids
    )

    # Annotation bar at bottom with cluster labels
    bottom_anno <- ComplexHeatmap::HeatmapAnnotation(
        gene_cluster = col_split,
        cluster_label = ComplexHeatmap::anno_block(
            gp = grid::gpar(fill = NA, col = NA),
            labels = reordered_cluster_ids,
            labels_gp = grid::gpar(fontsize = 10)
        ),
        col = list(gene_cluster = cluster_colors),
        show_annotation_name = FALSE,
        show_legend = FALSE,
        which = "column"
    )

    # colors and mapping (as you have)
    seismic_colors <- c(
        "#00008B", "#4040FF", "#8080FF", "#C0C0FF",
        "#FFFFFF",
        "#FFC0C0", "#FF8080", "#FF4040", "#8B0000"
    )

    rng <- range(mat, na.rm = TRUE)
    vmin <- rng[1]; vmax <- rng[2]

    neg_breaks <- seq(vmin, 0, length.out = 5)
    pos_breaks <- seq(0, vmax, length.out = 5)
    brks <- c(neg_breaks, pos_breaks[-1])

    col_fun <- circlize::colorRamp2(brks, seismic_colors, space = "RGB")

    # legend ticks: exact ends plus neat 0.1 steps inside
    inner <- seq(ceiling(vmin*2)/2, floor(vmax*2)/2, by = 0.25)
    ats <- unique(c(vmin, inner[inner > vmin & inner < vmax], vmax))
    labs <- sprintf("%.2f", ats)

    ht <- ComplexHeatmap::Heatmap(
        mat,
        name = "Effect",
        column_split = col_split,
        cluster_columns = FALSE,
        cluster_rows = TRUE,
        show_column_dend = FALSE,
        show_column_names = FALSE,
        bottom_annotation = bottom_anno,
        column_title = column_title,
        column_title_gp = grid::gpar(fontsize = 20, fontface = "bold"),
        row_title = "Cell types",
        col = col_fun,
        column_names_gp = grid::gpar(fontsize = 6),
        row_names_gp = grid::gpar(fontsize = 10),
        gap = grid::unit(1.5, "mm"),
        use_raster = TRUE,
        heatmap_legend_param = list(
            color_bar = "continuous",
            legend_height = grid::unit(5, "cm")
        )
    )

    ComplexHeatmap::draw(ht)
}


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



##############
# ADHOC
#############

simpliestKmeansTest<-function () {

    d=read.table("/downloads/effect_sizes.txt", header=T, sep="\t", row.names=1, check.names=F)
    clusters_python=read.table("/downloads/cluster_labels.txt", header=T, sep="\t", stringsAsFactors = F, col.names=c("gene", "cluster"))
    cluster_labels_python<-clusters_python$cluster
    names (cluster_labels_python)<-clusters_python$gene

    centroids_python=

    k=length(unique(cluster_labels_python))

    km <- kmeans(
        d,
        centers = k,
        nstart = 100,        # sklearn default (or auto, which is 10 in most versions)
        iter.max = 300,     # sklearn default
        algorithm = "Lloyd" # closest to sklearn
    )

    cluster_labels_R<-km$cluster
    map_clusters_by_gene_overlap(cluster_labels_python, cluster_labels_R)
    confusion_matrix=table(cluster_labels_python, cluster_labels_R, dnn=list("python", "R"))

    library(ClusterR)
    set.seed(12)
    x=as.matrix(d)
    km <- KMeans_rcpp(
        x,
        clusters = 17,
        num_init = 100,              # = sklearn n_init
        max_iters = 300,             # = sklearn max_iter
        initializer = "kmeans++",    # sklearn default
        tol = 1e-4,                  # sklearn default
        verbose = FALSE
    )

    cluster_labels_R2<-km$cluster
    names(cluster_labels_R2)<-rownames(d)

    #almost exactly the same results.
    map_clusters_by_gene_overlap(cluster_labels_R, cluster_labels_R2)
    map_clusters_by_gene_overlap(cluster_labels_python, cluster_labels_R)



}


#for adhoc testing of clustering
cluster_playground<-function (m, genesPassFDR) {

    #load up Haley's results for comparison.
    getClusterLabels<-function (inFile) {
        clusters<-read.table(inFile, header=T, stringsAsFactors = F, sep=",")
        colnames (clusters)=c("gene", "cluster")
        cluster_labels<-clusters$cluster
        names (cluster_labels)<-clusters$gene
        return (cluster_labels)
    }

    haley_cluster_labels_raw<-getClusterLabels(gene_cluster_file_raw)
    haley_cluster_labels_mash<-getClusterLabels(gene_cluster_file_post_mash)

    #surprisingly, clusters aren't very similar between raw and mash corrected data.
    map_clusters_by_gene_overlap(haley_cluster_labels_raw, haley_cluster_labels_mash)

    #cluster the mash results by correlation or euclidean distance.
    #turns out these are approximately the same!
    #this uses the filtered results only, then sets effect sizes to 0 if the sign test isn't significant.
    m_fdr=filter_mash_result_by_genes(m, genesPassFDR)
    result <- cluster_mash_effects(m_fdr, lfsr_threshold = 0.05, k = 16, method="kmeans", dist_metric = "euclidean", nstart=200, seed=2)
    plot_mash_heatmap(result2$pm_filtered, result2$clusters, column_title="distance [euclidean] clustering [kmeans]")

    result2 <- cluster_mash_effects(m_fdr, lfsr_threshold = 0.05, k = 16, method="kmeans", dist_metric = "correlation", nstart=200, seed=1)
    plot_mash_heatmap(result$pm_filtered, result$clusters, column_title="distance [correlation] clustering [kmeans]")

    map_clusters_by_gene_overlap(result$clusters, result2$clusters)

    #write the effect size matrix out for Haley
    write.table(result2$pm_filtered, "/broad/mccarroll/haley/BICAN/DE/age_DE_logFC_mash_corrected_effects_filtered_sig_effects_matrix.txt", row.names=T, col.names = T, quote=F, sep="\t")


    #######################
    #compare to Haley's results
    ##########################
    #subset to haley's explicit list of genes.
    m_genes=filter_mash_result_by_genes(m, names (haley_cluster_labels_raw))
    plot_mash_heatmap(m_genes$data$Bhat, haley_cluster_labels_raw, column_title="Haley Clustering of raw effect sizes")


    #try to replicate?
    k=length(unique(haley_cluster_labels_raw))
    #km is almost exactly the same result as cluster_mash_effects.
    km <- kmeans(
        m_genes$data$Bhat,
        centers = k,
        nstart = 100,        # sklearn default (or auto, which is 10 in most versions)
        iter.max = 300,     # sklearn default
        algorithm = "Lloyd" # closest to sklearn
    )

    map_clusters_by_gene_overlap(haley_cluster_labels_raw, km$cluster)

    result <- cluster_mash_effects(m_genes, lfsr_threshold = 1, k = k, method="kmeans", dist_metric = "euclidean", nstart=200, seed=1, use_raw_values=TRUE)
    plot_mash_heatmap(m_genes$data$Bhat, result$clusters, column_title="Jim Clustering of raw effect sizes")
    map_clusters_by_gene_overlap(haley_cluster_labels_raw, result$clusters)
    map_clusters_by_gene_overlap(km$cluster, result$clusters)




    #What if I just cluster on the sign of the effect?  It's yucky.
    result_sign <- cluster_mash_effects(m_fdr, lfsr_threshold = 0.05, k = 15, method="kmeans", dist_metric = "euclidean", nstart=200, seed=1, test_signs_only=TRUE)
    plot_mash_heatmap(pm[rownames (result_sign$pm_filtered),], result_sign$clusters, column_title="distance [correlation] clustering [kmeans]")



    r=find_optimal_k_for_mash_effects(m_fdr, lfsr_threshold = 0.05, start_num_clusters=5, end_num_clusters=30, method="kmeans", dist_metric = "correlation", nstart=200, seed=1)

    result3 <- cluster_mash_effects(m_fdr, lfsr_threshold = 0.05, k = 20, method="kmeans", dist_metric = "correlation", nstart=200, seed=1)
    plot_mash_heatmap(result3$pm_filtered, result3$clusters, column_title="distance [correlation] clustering [kmeans]")
    #plot_mash_heatmap(pm[rownames (result3$pm_filtered),], result3$clusters, column_title="distance [correlation] clustering [kmeans] unfiltered data")

    map_clusters_by_gene_overlap(result$clusters, result3$clusters)

    result4 <- cluster_mash_effects(m_fdr, lfsr_threshold = 0.05, k = 11, method="kmeans", dist_metric = "correlation", nstart=200, seed=1)
    plot_mash_heatmap(result4$pm_filtered, result4$clusters, column_title="distance [correlation] clustering [kmeans]")
    #plot_mash_heatmap(pm[rownames (result3$pm_filtered),], result3$clusters, column_title="distance [correlation] clustering [kmeans] unfiltered data")

    map_clusters_by_gene_overlap(result$clusters, result3$clusters)


    resul4t <- cluster_mash_effects(m_fdr, lfsr_threshold = 0.05, k = 20, method="kmeans", dist_metric = "euclidean")
    plot_mash_heatmap(result$pm_filtered, result$clusters, column_title="distance [euclidean] clustering [kmeans]")
    plot_mash_heatmap(pm[rownames (result$pm_filtered),], result$clusters, column_title="distance [euclidean] clustering [kmeans] unfiltered data")

    result <- cluster_mash_effects(m_fdr, lfsr_threshold = 0.05, k = 30, method="kmeans", dist_metric = "correlation")
    plot_mash_heatmap(result$pm_filtered, result$clusters, column_title="distance [correlation] clustering [kmeans]")
    plot_mash_heatmap(pm[rownames (result$pm_filtered),], result$clusters, column_title="distance [correlation] clustering [kmeans] unfiltered data")

    result <- cluster_mash_effects(m_fdr, lfsr_threshold = 0.05, k = 30, method="kmeans", dist_metric = "euclidean")
    plot_mash_heatmap(result$pm_filtered, result$clusters, column_title="distance [euclidean] clustering [kmeans]")
    plot_mash_heatmap(pm[rownames (result$pm_filtered),], result$clusters, column_title="distance [euclidean] clustering [kmeans] unfiltered data")

    #The hclust results look categorically worse across the board.
    #result <- cluster_mash_effects(m_fdr, lfsr_threshold = 0.05, k = 16, method= "hclust", dist_metric = "correlation")
    #plot_mash_heatmap(result$pm_filtered, result$clusters, column_title="distance [correlation] clustering [hclust]")

}

effect_size_playground<-function () {
    #Visualize the effect size distribution of the raw data
    bh=as.vector(abs(Bhat))
    fh=as.vector(FDR)
    bh_filtered=bh[which(fh<0.05)]
    hist (bh_filtered, breaks=50, main="effect size conditioned on FDR<=0.05")

    #maybe you want the max per gene?
    fdr_threshold=0.01
    idx <- apply(Bhat, 1, function(x) which.max(abs(x)))
    max_effect <- abs(Bhat[cbind(seq_len(nrow(Bhat)), idx)])
    fdr_at_max <- FDR[cbind(seq_len(nrow(FDR)), idx)]
    bh_filtered=max_effect[fdr_at_max<=fdr_threshold]
    hist (bh_filtered, breaks=50, main=paste("max effect size conditioned on FDR<=", fdr_threshold, sep=""))
    abline (v=log2(1.02), col='red')
    length(bh_filtered)
    length(which(bh_filtered>log2(1.02)))


    #Visualize the effect size distribution of the posterior data
    pm_filtered=abs(as.vector(pm))[as.vector(lfsr)<0.05]
    hist (pm_filtered, breaks=50, main="posterior effect size conditioned on lsfr<=0.05")


    lfsr_threshold=0.05
    idx <- apply(pm, 1, function(x) which.max(abs(x)))
    max_effect <- abs(pm[cbind(seq_len(nrow(pm)), idx)])
    fdr_at_max <- lfsr[cbind(seq_len(nrow(lfsr)), idx)]
    bh_filtered=max_effect[fdr_at_max<=lfsr_threshold]
    hist (bh_filtered, breaks=50, main=paste("max effect size conditioned on lfsr<=", fdr_threshold, sep=""))
    abline (v=log2(1.02), col='red')
    length(bh_filtered)
    length(which(bh_filtered>log2(1.02)))

}

confusionMatrixWithKMeansClustering<-function (m, gene_cluster_file, gene_cluster_labels_file) {
    clusters<-read.table(gene_cluster_file, header=T, stringsAsFactors = F, sep=",")
    colnames (clusters)=c("gene", "cluster")
    labels<-read.table(gene_cluster_labels_file, header=T, stringsAsFactors = F, sep="\t")
    colnames(labels)<-c("cluster", "kmeans_label")
    clusters<-merge(clusters, labels, by = "cluster")

    #make a confusion matrix of the mash components vs the haley clusters
    #don't threshold the mash components, match Haley's inputs.
    res <- partition_mash_calls_by_component(m, thresh = 1, keep_null = T)
    assignment<-res$assignment

    idx=match(clusters$gene, assignment$gene)
    stopifnot(length(which(is.na(idx)))==0)

    clusters$mash_component <- assignment$label[idx]
    confusion_matrix=table(clusters$kmeans_label, clusters$mash_component, dnn=list("Kmeans", "Mash"))




    plot_confusion_matrix(confusion_matrix, normalize=T, title="Confusion matrix of Kmeans clusters vs mash components")

    confusion_matrix[1,,drop=F]

    #astrocyte down vs ED_PCA1.
    head (clusters[clusters$kmeans_label=="astrocyte_down" & clusters$mash_component=="ED_PCA_1",])
    mash_forest_with_labels(m,  gene="TGIF2", order_by = "effect", call_rule="both", plot_raw=T)
    barplot (res$weights_collapsed["TGIF2",], las=2)
    mash_forest_with_labels(m,  gene="PHACTR1", order_by = "effect", call_rule="both", plot_raw=T)
    barplot (res$weights_collapsed["PHACTR1",], las=2)



}


find_abs_increases <- function(Bhat, pm) {
    if (!all(dim(Bhat) == dim(pm))) {
        stop("Bhat and pm must have the same dimensions")
    }

    # 1) Identify increases, excluding zero-effect Bhat entries
    inc_mat <- (abs(pm) > abs(Bhat)) & (Bhat != 0)
    idx <- which(inc_mat, arr.ind = TRUE)

    # 2) Extract values
    genes <- if (!is.null(rownames(Bhat))) rownames(Bhat)[idx[,1]] else idx[,1]
    conds <- if (!is.null(colnames(Bhat))) colnames(Bhat)[idx[,2]] else idx[,2]

    bh_vals <- Bhat[idx]
    pm_vals <- pm[idx]

    df_inc <- data.frame(
        gene       = genes,
        condition  = conds,
        Bhat       = bh_vals,
        pm         = pm_vals,
        abs_Bhat   = abs(bh_vals),
        abs_pm     = abs(pm_vals),
        delta_abs  = abs(pm_vals) - abs(bh_vals),
        same_sign  = sign(pm_vals) == sign(bh_vals),
        stringsAsFactors = FALSE
    )

    # 3) Order by size of increase
    df_inc <- df_inc[order(df_inc$delta_abs, decreasing = TRUE), ]

    # 4) Summaries
    summary <- list(
        num_increases   = sum(inc_mat),
        frac_increases  = mean(inc_mat),
        row_increases   = rowSums(inc_mat),
        col_increases   = colSums(inc_mat)
    )

    return(list(
        summary = summary,
        increases = df_inc,
        increases_same_sign = df_inc[df_inc$same_sign, ]
    ))
}

plot_gene_effects <- function(Bhat, pm, gene) {
    if (!(gene %in% rownames(Bhat)) || !(gene %in% rownames(pm))) {
        stop("Gene not found in both matrices")
    }

    df <- data.frame(
        Bhat = Bhat[gene, ],
        pm   = pm[gene, ],
        condition = colnames(Bhat)
    )

    ggplot(df, aes(x = Bhat, y = pm, label = condition)) +
        geom_point(size = 3, color = "steelblue") +
        geom_smooth(method = "lm", se = FALSE, color = "red") +
        geom_text(vjust = -0.8, size = 3) +
        labs(
            title = paste("Effect sizes for", gene),
            x = "Original effect size (Bhat)",
            y = "Posterior mean (pm)"
        ) +
        theme_minimal(base_size = 14)
}



#########################
# FIT / Filter MASH MODEL
###########################


#' Fit mash with EM-estimated residual correlation
#'
#' Fits a \pkg{mashr} model to per-gene effects across cell types while
#' estimating the residual correlation matrix \eqn{V} (e.g., from overlapping
#' donors) via \code{mash_estimate_corr_em}. Data-driven covariance components
#' are learned from "strong" signals and combined with canonical covariances for
#' the final mash fit.
#'
#' @param Bhat Numeric matrix of effect estimates (genes x cell types).
#' @param Shat Numeric matrix of standard errors (same dimensions as \code{Bhat});
#'   must be finite and strictly positive.
#' @param FDR Numeric matrix of FDR estimates (genes x cell types)
#' @param npc Integer. Number of principal components used to seed
#'   data-driven covariances (capped at \code{ncol(Bhat)}). Default: \code{5L}.
#' @param n_null_est Integer. Number of rows randomly sampled to estimate the
#'   residual correlation \eqn{V} via EM. Default: \code{5000L}.
#'
#' @return A fitted \code{mash} object (use \code{mashr} accessors such as
#'   \code{get_pm}, \code{get_psd}, and \code{get_lfsr} for summaries).
#'
#' @export
fit_mash<- function(Bhat, Shat, FDR, npc = 5L, n_null_est = 5000L, missing_mask=NULL, add.mem.profile=TRUE, only_cannonical_matrixes=F) {
    data0 <- mash_set_data(Bhat, Shat)

    # pick "strong" rows to learn data-driven covariances
    Z <- Bhat / Shat
    zmax <- apply(abs(Z), 1, max)
    strong_idx <- choose_strong_rows(Bhat, Shat)
    data_strong <- mash_set_data(Bhat[strong_idx, , drop = FALSE], Shat[strong_idx, , drop = FALSE])

    # canonical + data-driven covariances
    U.c  <- cov_canonical(data0)
    U.p  <- cov_pca(data_strong, npc = min(npc, ncol(Bhat)))
    U.ed <- cov_ed(data_strong, U.p)  # extreme deconvolution refines covariances

    if (only_cannonical_matrixes) {
        covariance_matrix_list=c(U.c)
    } else {
        covariance_matrix_list=c(U.c, U.ed)
    }


    # estimate V with EM on a random subset (accounts for donor overlap)
    # Many attempts to construct V capture biological signal of known related cell subtypes
    # so we're going to skip V.
    #set.seed(1)
    #idx <- sample.int(nrow(Bhat), size = min(n_null_est, nrow(Bhat)))
    #data_tmp <- mash_set_data(Bhat[idx, , drop = FALSE], Shat[idx, , drop = FALSE])
    #Vhat_out <- mash_estimate_corr_em(data_tmp, U.ed)   # returns list with $V and $mash.model
    #Vhat <- Vhat_out$V

    # fit on all genes using V
    Vhat=diag(ncol(Bhat)) #this is the default for mashr
    data_all <- mash_set_data(Bhat, Shat, V = Vhat)
    m <- mash(data_all, covariance_matrix_list, usepointmass = TRUE, outputlevel = 2)


    #remask the raw data.


    #add the raw data directly to the mash object
    m$data<-data_all
    result=list(mash=m, covariance_matrix_list=covariance_matrix_list, missing_mask=missing_mask)
    return (result)

}

#' Two-stage mash fit: learn covariances on informative rows, then apply to all genes
#'
#' @param Bhat  numeric matrix (genes x conditions) of effect estimates
#' @param Shat  numeric matrix (genes x conditions) of standard errors
#' @param npc   integer; number of PCs to pass to cov_pca (capped at ncol(Bhat))
#' @param only_cannonical_matrixes logical; if TRUE use only canonical covariances (no ED_tPCA)
#' @param huge  numeric; SE threshold used to define "well-measured" rows (exclude rows with any SE >= huge)
#' @param grid  numeric vector; scale grid passed to \code{mash} (set to \code{NULL} to use mash defaults)
#' @param autoselect.grid logical; passed to \code{mash}; default FALSE to enforce the provided grid
#' @param usepointmass logical; include a point-mass-at-0 component
#'
#' @return list with elements:
#'   \item{m_learn}{mash fit on informative rows (learned prior)}
#'   \item{m}{final mash fit on all rows using the learned prior}
#'   \item{learn_rows}{integer indices of rows used for learning}
#'   \item{well_measured}{logical vector marking rows with all SE < huge}
#'   \item{strong_idx_w}{logical vector (within well_measured subset) chosen by choose_strong_rows}
#'   \item{Ulist}{covariance library used}
#'   \item{grid}{scale grid actually used}
#'
#' @details
#' Stage 1 learns the mixture on rows that are both well-measured and strong.
#' Stage 2 applies the learned prior to all genes, preventing missing/huge-SE rows
#' from steering the covariance learning. A floor on the scale grid helps avoid
#' soaking up tiny correlated drift with a structured component.
fit_mash_two_stage <- function(Bhat, Shat,
                               npc = 5,
                               only_cannonical_matrixes = FALSE,
                               huge = 1e5,
                               usepointmass = TRUE) {
    # basic checks
    if (!is.matrix(Bhat) || !is.matrix(Shat)) stop("Bhat and Shat must be matrices")
    if (!all(dim(Bhat) == dim(Shat))) stop("Bhat and Shat must have identical dimensions")
    if (ncol(Bhat) < 2L) stop("Need at least 2 conditions (columns)")
    if (npc < 1L) stop("npc must be >= 1")

    # 1) define well-measured rows (exclude any with huge/non-finite SE)
    well_measured <- apply(Shat < huge & is.finite(Shat), 1L, all)

    # 2) among well-measured, pick "strong" rows for learning
    Bhat_w <- Bhat[well_measured, , drop = FALSE]
    Shat_w <- Shat[well_measured, , drop = FALSE]
    if (nrow(Bhat_w) == 0L) stop("No well-measured rows under the given 'huge' threshold")
    strong_idx_w <- choose_strong_rows(Bhat_w, Shat_w)

    if (!any(strong_idx_w)) stop("No strong rows found within the well-measured set")
    learn_rows <- which(well_measured)[strong_idx_w]

    # 3) build mash data objects
    data_all   <- mash_set_data(Bhat, Shat)  # V = I by default
    data_learn <- mash_set_data(Bhat[learn_rows, , drop = FALSE],
                                Shat[learn_rows, , drop = FALSE])

    # 4) covariance library (canonical + ED_tPCA learned on informative rows)
    U.c  <- cov_canonical(data_all)
    U.p  <- cov_pca(data_learn, npc = min(npc, ncol(Bhat)))
    U.ed <- cov_ed(data_learn, U.p)
    Ulist <- if (isTRUE(only_cannonical_matrixes)) U.c else c(U.c, U.ed)

    # 5) stage 1: learn prior on informative rows
    m_learn <- mash(data_learn,
                    Ulist = Ulist,
                    usepointmass = usepointmass,
                    outputlevel = 2)

    # 6) stage 2: apply learned prior to all rows
    m <- mash(data_all, g = get_fitted_g(m_learn), outputlevel = 2)

    list(m_learn = m_learn,
         m = m,
         learn_rows = learn_rows,
         well_measured = well_measured,
         strong_idx_w = strong_idx_w,
         Ulist = Ulist,
         grid = grid)
}


#Restore the NA values in the original data to the mash outputs.
mask_mash_inplace <- function(m, missing_mask) {
    if (is.null(missing_mask))
        return (m)

    stopifnot(inherits(m, "mash"), is.matrix(missing_mask))

    # Validate dimensions against one of the result matrices
    if (!is.list(m$result) || is.null(m$result$PosteriorMean))
        stop("m$result$PosteriorMean not found; is this a fitted mash object?")
    if (!all(dim(m$result$PosteriorMean) == dim(missing_mask)))
        stop("missing_mask must have same dim as the posterior matrices.")

    # List of result matrices to mask if present
    fields <- c("PosteriorMean", "PosteriorSD", "lfsr", "NegativeProb")
    for (f in fields) {
        if (!is.null(m$result[[f]])) {
            mat <- m$result[[f]]
            if (!all(dim(mat) == dim(missing_mask)))
                stop("Dimension mismatch for result field: ", f)
            mat[missing_mask] <- NA_real_
            m$result[[f]] <- mat
        }
    }
    return (m)
}


# Choose a "strong" subset adaptively for mash covariance learning
choose_strong_rows <- function(Bhat, Shat,
                               target_n = 5000L,         # desired size of strong set
                               min_per_ct = 300L,        # ensure per-cell-type representation
                               max_n = 20000L) {         # guardrail for very large sets
    stopifnot(is.matrix(Bhat), is.matrix(Shat), all(dim(Bhat) == dim(Shat)))
    Z <- Bhat / Shat
    G <- nrow(Z); C <- ncol(Z)

    # 1) Global top by rowwise max |Z|
    zmax <- apply(abs(Z), 1L, max)
    ord_global <- order(zmax, decreasing = TRUE)
    k_global <- min(target_n, G)
    idx_global <- ord_global[seq_len(k_global)]

    # 2) Per-cell-type top |Z| (representation)
    idx_ct <- integer(0)
    for (j in seq_len(C)) {
        ord_j <- order(abs(Z[, j]), decreasing = TRUE)
        k_j <- min(min_per_ct, G)
        idx_ct <- union(idx_ct, ord_j[seq_len(k_j)])
    }

    # 3) Union and cap
    idx <- union(idx_global, idx_ct)
    if (length(idx) > max_n) idx <- idx[seq_len(max_n)]

    sort(idx)
}


#' Subset a mash fit by per-gene directional significance and effect size
#'
#' @description
#' Keep only genes that have **at least one condition** \(j\) with
#' \code{lfsr[i, j] < lfsr_threshold} **and** \code{abs(pm[i, j]) >= effect_threshold}.
#' Returns a new `mash` object containing just the selected genes (no refit).
#'
#' @param m A fitted `mash` object.
#' @param lfsr_threshold Numeric scalar; directional error cutoff. Default `0.05`.
#' @param effect_threshold Numeric scalar; minimum absolute **posterior mean** on
#'   the same scale as \code{ashr::get_pm(m)}. Default `0.05`.
#'
#' @return A `mash` object subset to the passing genes. If no genes pass, the
#'   original object is returned with a message.
#'
#' @details
#' Posterior summaries are obtained via the `ashr` generics
#' \code{ashr::get_pm} and \code{ashr::get_lfsr}. The function subsets:
#' \itemize{
#'   \item \code{m$result} matrices/vectors that are gene-aligned,
#'   \item \code{m$posterior_weights} (if present),
#'   \item \code{m$data} via \code{mashr::mash_set_data} to keep Bhat/Shat (and V/alpha).
#' }
#' Model-level slots (mixture, \eqn{\pi}, loglik) are kept unchanged.
#'
#' @export
filter_mash_by_lfsr_effect <- function(m,
                                       lfsr_threshold = 0.05,
                                       effect_threshold = 0.05) {
    stopifnot(inherits(m, "mash"))

    # pull posterior summaries (ashr generics dispatch on mash)
    PM   <- ashr::get_pm(m)
    LFSR <- ashr::get_lfsr(m)
    if (is.null(PM) || is.null(LFSR))
        stop("Posterior summaries not found in mash object (pm/lfsr).")

    # per-gene: does ANY condition meet both criteria?
    pass_any <- apply(LFSR < lfsr_threshold & abs(PM) >= effect_threshold,
                      1L, function(x) any(isTRUE(x), na.rm = TRUE))
    pass_any[is.na(pass_any)] <- FALSE

    if (!any(pass_any)) {
        message("No genes met (lfsr < ", lfsr_threshold,
                " & |pm| >= ", effect_threshold, "); returning original object.")
        return(m)
    }
    idx <- which(pass_any)

    # --- build subset without refitting ---
    m2 <- m

    # subset gene-aligned entries in m$result
    if (!is.null(m$result) && length(m$result)) {
        for (nm in names(m$result)) {
            x <- m$result[[nm]]
            if (is.matrix(x) && nrow(x) == nrow(PM)) {
                m2$result[[nm]] <- x[idx, , drop = FALSE]
            } else if (is.array(x) && length(dim(x)) >= 3L && dim(x)[1] == nrow(PM)) {
                # e.g. PosteriorSamples [G x C x D]
                sl <- vector("list", length(dim(x))); sl[] <- quote(expr=)
                sl[[1]] <- idx
                m2$result[[nm]] <- do.call(`[`, c(list(x), sl, list(drop = FALSE)))
            } else if (is.vector(x) && length(x) == nrow(PM)) {
                m2$result[[nm]] <- x[idx]
            } else {
                m2$result[[nm]] <- x  # leave as-is
            }
        }
    }

    # subset posterior weights (if present)
    if (!is.null(m$posterior_weights)) {
        m2$posterior_weights <- m$posterior_weights[idx, , drop = FALSE]
    }

    # rebuild data slot so ashr::get_* continues to work cleanly
    if (!is.null(m$data)) {
        Bhat  <- tryCatch(m$data$Bhat,  error = function(e) NULL)
        Shat  <- tryCatch(m$data$Shat,  error = function(e) NULL)
        V     <- tryCatch(m$data$V,     error = function(e) NULL)
        alpha <- tryCatch(m$data$alpha, error = function(e) NULL)

        if (!is.null(Bhat) && !is.null(Shat)) {
            m2$data <- mashr::mash_set_data(
                Bhat[idx, , drop = FALSE],
                Shat[idx, , drop = FALSE],
                V = V,
                alpha = alpha
            )
        }
    }

    m2
}

filter_mash_result_by_genes <- function(mash_result, genes_to_keep) {
    # Get gene names from mash result
    gene_names <- rownames(ashr::get_pm(mash_result))

    # Create logical vector for subsetting
    keep_logical <- gene_names %in% genes_to_keep

    if (sum(keep_logical) == 0) {
        stop("None of the specified genes were found in the mash result.")
    }

    # Subset result matrices
    new_mash <- mash_result
    new_mash$result$PosteriorMean <- mash_result$result$PosteriorMean[keep_logical, , drop = FALSE]
    new_mash$result$PosteriorSD <- mash_result$result$PosteriorSD[keep_logical, , drop = FALSE]
    new_mash$result$lfsr <- mash_result$result$lfsr[keep_logical, , drop = FALSE]
    new_mash$result$lfdr <- mash_result$result$lfdr[keep_logical, , drop = FALSE]
    new_mash$result$NegativeLog10PosteriorLocalBayesFactor <-
        mash_result$result$NegativeLog10PosteriorLocalBayesFactor[keep_logical, , drop = FALSE]

    # Subset data matrices
    new_mash$data$Bhat <- mash_result$data$Bhat[keep_logical, , drop = FALSE]
    new_mash$data$Shat <- mash_result$data$Shat[keep_logical, , drop = FALSE]

    # Optional fields
    if (!is.null(mash_result$strong_signals)) {
        new_mash$strong_signals <- mash_result$strong_signals[keep_logical]
    }
    if (!is.null(mash_result$posterior_weights)) {
        new_mash$posterior_weights <- mash_result$posterior_weights[keep_logical, , drop = FALSE]
    }
    if (!is.null(mash_result$null_loglik)) {
        new_mash$null_loglik <- mash_result$null_loglik[keep_logical]
    }
    if (!is.null(mash_result$alt_loglik)) {
        new_mash$alt_loglik <- mash_result$alt_loglik[keep_logical]
    }

    return(new_mash)
}


#' Partition significant genes by mash mixture component families
#'
#' @description
#' For genes with at least one condition passing an LFSR threshold, collapse
#' mash posterior weights across scales/variants into **component families**
#' that align with `mashr::get_estimated_pi(m)`, optionally relabel those
#' families, and assign each gene to the family with the largest collapsed
#' posterior weight (also returning the full soft weights).
#'
#' @param m A fitted `mash` object (fit with `outputlevel = 2` so that
#'   `m$posterior_weights` is available).
#' @param thresh Numeric scalar. Genes are considered for partitioning if
#'   they have `ashr::get_lfsr(m) < thresh` in **any** condition.
#'   Default: `0.05`.
#' @param keep_null Logical. If `FALSE` (default), drop the `"null"` family
#'   from the output columns; if `TRUE`, keep it.
#' @param map_labels Character; one of `"requested"` or `"identity"`.
#'   - `"identity"`: keep family names exactly as derived from the mixture
#'     (e.g., `"ED_PCA_1"`, `"equal_effects"`, `"identity"`, cell-type
#'     private components).
#'   - `"requested"`: friendlier relabeling: keep raw cell-type names for
#'     private components, preserve each `ED_*` family as-is, map
#'     `"equal_effects"`/`"EE"` → `"shared"`, and `"identity"` → `"independent"`.
#'
#' @return A list with:
#' \describe{
#'   \item{assignment}{`data.frame(gene, label, post_weight)` — per-gene hard
#'   label (argmax family) and its collapsed posterior weight.}
#'   \item{counts}{`table` of hard-label counts (descending).}
#'   \item{fractions}{Named numeric vector: percentage for each label.}
#'   \item{weights_collapsed}{Matrix [genes × families] of **soft** weights,
#'   obtained by summing posterior weights across all scales/variants that
#'   belong to the same family. Row names are gene IDs.}
#'   \item{labels}{Character vector of family labels (column order of
#'   `weights_collapsed`).}
#' }
#'
partition_mash_calls_by_component <- function(
        m,
        thresh = 0.05,                 # gene is "significant" if any lfsr<thresh across CTs
        keep_null = FALSE,             # keep "null" in the label set
        map_labels = c("requested", "identity")  # relabeling policy
) {
    stopifnot(inherits(m, "mash"))
    map_labels <- match.arg(map_labels)

    # Use ashr generics so S3 dispatch hits mash methods (per earlier note).
    PM   <- ashr::get_pm(m)
    LFSR <- ashr::get_lfsr(m)

    # ---- 1) Select significant genes (any condition LFSR < thresh) ----
    sig_idx <- which(rowSums(is.finite(LFSR) & LFSR < thresh) > 0L)
    if (!length(sig_idx)) {
        return(list(
            assignment = data.frame(gene = character(0), label = character(0),
                                    post_weight = numeric(0), stringsAsFactors = FALSE),
            counts = integer(0), fractions = numeric(0),
            weights_collapsed = matrix(0, 0, 0),
            labels = character(0)
        ))
    }
    gene_names <- if (!is.null(rownames(PM))) rownames(PM)[sig_idx] else as.character(sig_idx)

    # ---- 2) Pull posterior mixture weights for the selected genes ----
    W <- m$posterior_weights
    if (is.null(W)) stop("m$posterior_weights not found; fit mash with outputlevel = 2.")
    W <- W[sig_idx, , drop = FALSE]   # [genes × mixture components]

    # ---- 3) Collapse columns to component 'families' (match get_estimated_pi) ----
    comp_names <- colnames(W)
    if (is.null(comp_names)) comp_names <- names(m$fitted_g$pi)

    # Remove scale/variant suffixes and internal tags to get the base family name
    base_of_comp <- comp_names
    base_of_comp <- sub("\\s*[:_\\.]\\s*\\d+$", "", base_of_comp)   # strip trailing _3/:2/.4
    base_of_comp <- sub("\\s*\\(scale.*\\)$", "", base_of_comp)     # strip "(scale=...)"
    base_of_comp <- sub("\\s*::.*$", "", base_of_comp)              # strip "::internal"

    # Align to the names used by mixing proportions (pi) so columns line up with get_estimated_pi
    pi_names <- names(mashr::get_estimated_pi(m))
    map_to_pi <- vapply(base_of_comp, function(nm) {
        # allow either "family" prefix to match (robust to minor naming differences)
        hit <- which(startsWith(nm, pi_names) | startsWith(pi_names, nm))
        if (length(hit)) pi_names[hit[1]] else nm
    }, character(1))

    # Sum weights over all mixture columns that belong to the same family
    labs_unique <- unique(map_to_pi)
    Wc <- do.call(cbind, lapply(labs_unique, function(lb) {
        cols <- which(map_to_pi == lb)
        rowSums(W[, cols, drop = FALSE])
    }))
    colnames(Wc) <- labs_unique

    # Order families to match get_estimated_pi (with any extras appended)
    ord_cols <- unique(c(pi_names[pi_names %in% colnames(Wc)], setdiff(colnames(Wc), pi_names)))
    Wc <- Wc[, ord_cols, drop = FALSE]

    # Optionally drop the "null" family
    if (!keep_null && "null" %in% colnames(Wc)) {
        Wc <- Wc[, colnames(Wc) != "null", drop = FALSE]
    }
    rownames(Wc) <- gene_names

    # ---- 4) Optional relabeling to friendlier names ----
    final_labels <- colnames(Wc)
    if (map_labels == "requested") {
        celltypes <- colnames(PM)  # recognize private components that equal a CT name
        relabel <- function(lb) {
            if (lb %in% celltypes) lb                   # private to a CT: keep CT name
            else if (grepl("^ED_", lb)) lb              # each ED_* kept separate
            else if (lb %in% c("equal_effects", "EE")) "shared"
            else if (lb == "identity")                 "independent"
            else lb                                    # e.g., simple_het_3 stays as-is
        }
        final_labels <- vapply(colnames(Wc), relabel, character(1))
        colnames(Wc) <- final_labels
    }

    # ---- 5) Hard assignment = argmax family per gene (+ the winning weight) ----
    top_col <- max.col(Wc, ties.method = "first")
    top_lab <- colnames(Wc)[top_col]
    top_w   <- Wc[cbind(seq_len(nrow(Wc)), top_col)]

    assignment <- data.frame(
        gene = gene_names,
        label = top_lab,
        post_weight = top_w,
        stringsAsFactors = FALSE
    )

    # ---- 6) Summary counts and fractions ----
    counts    <- sort(table(assignment$label), decreasing = TRUE)
    fractions <- round(100 * counts / sum(counts), 2)

    # ---- 7) Return both hard calls and soft weights ----
    list(
        assignment = assignment,            # per-gene hard label (+ posterior weight)
        counts = counts,                    # table by label
        fractions = fractions,              # %
        weights_collapsed = Wc,             # genes × labels (soft weights)
        labels = colnames(Wc)               # final label order used
    )
}


###################
# PLOTTING CODE
###################

plot_confusion_matrix <- function(cm, normalize = FALSE, title = "Confusion Matrix") {
    if (!is.matrix(cm)) stop("cm must be a matrix")
    if (is.null(rownames(cm)) || is.null(colnames(cm))) {
        stop("Matrix must have row and column names")
    }

    m <- cm
    if (normalize) {
        m <- sweep(m, 1, rowSums(m), FUN = "/")  # row-wise normalization
    }

    df <- as.data.frame(as.table(m))
    names(df) <- c("True", "Predicted", "Freq")

    ggplot(df, aes(x = Predicted, y = True, fill = Freq)) +
        geom_tile(color = "white") +
        geom_text(aes(label = round(Freq, 2)), size = 3) +
        scale_fill_gradient(low = "white", high = "steelblue") +
        labs(title = title, x = "Predicted (Mash)", y = "True (Kmeans)") +
        theme_minimal(base_size = 12) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = 10)
        )
}


#' Forest plot of mash posteriors (optional raw overlay) with per-condition calls
#'
#' @description
#' For a single gene, plot posterior means and CIs from a fitted `mash` object
#' across conditions (cell types). Conditions can be colored by a significance
#' rule based on a metric (default **lfsr**), CI exclusion of 0, or both.
#' Optionally overlay the **raw** (pre-mash) estimates and CIs taken from the
#' mash input (`m$data$Bhat`, `m$data$Shat`) to visualize shrinkage.
#'
#' @param m A fitted `mash` object.
#' @param gene Character gene ID (matching rownames of `ashr::get_pm(m)`) or
#'   1-based integer row index.
#' @param level Two-sided posterior CI level (default `0.95`).
#' @param metric One of `c("lfsr","custom")`. If `"lfsr"`, uses
#'   `ashr::get_lfsr(m)`; if `"custom"`, supply `metric_matrix` with same
#'   dimensions as `ashr::get_pm(m)`.
#' @param metric_matrix Numeric matrix used when `metric="custom"`.
#' @param thresh Threshold for `metric` (e.g., `0.05` for lfsr).
#' @param call_rule One of `c("metric","ci","both","either")` determining how
#'   per-condition “significant” is defined.
#' @param col_sig,col_ns Colors for significant / not significant **points**.
#' @param ci_sig,ci_ns Colors for significant / not significant **CIs**.
#' @param lwd_sig,lwd_ns Line widths for significant / not significant CIs.
#' @param pch_pts Point shape for posterior means.
#' @param cex_names Size of condition labels on the y-axis.
#' @param order_by One of `c("effect","abs_effect","none")` for row order.
#' @param top_largest Logical. If `TRUE` (default), largest effects appear at top.
#' @param colors Optional list like `meta.colors()` to control bg/axes/text.
#' @param zero_col Color of the vertical zero line.
#' @param plot_raw Logical. If `TRUE`, overlay raw (input) mean ± CI using
#'   `m$data$Bhat` and `m$data$Shat`. Requires that those matrices are present
#'   in `m$data`. Default `FALSE`.
#' @param col_raw Color for raw points (default `"steelblue4"`).
#' @param ci_raw  Color for raw CIs (default `"steelblue3"`).
#' @param lwd_raw Line width for raw CIs (default `1`).
#' @param pch_raw Point shape for raw means (default `16`).
#'
#' @return (Invisibly) a list with the vectors used for plotting:
#'   `sig`, `order`, `lo`, `hi`, `pm`, and if `plot_raw=TRUE`,
#'   `raw_lo`, `raw_hi`, `raw_pm`.
#'
#' @export
mash_forest_with_labels <- function(
        m, gene,
        level = 0.95,
        metric = c("lfsr", "custom"),
        metric_matrix = NULL, thresh = 0.05,
        call_rule = c("metric","ci","both","either"),
        col_sig = "firebrick", col_ns = "black",
        ci_sig = "firebrick", ci_ns = "grey60",
        lwd_sig = 2, lwd_ns = 1, pch_pts = 15,
        cex_names = 0.9, order_by = c("effect","abs_effect","none"),
        top_largest = TRUE,
        colors = NULL, zero_col = "grey70",
        plot_raw = FALSE, col_raw  = "steelblue4", ci_raw = "steelblue3",
        lwd_raw  = 1, pch_raw  = 16,
        raw_y_offset = 0.18, raw_side = c("below","above")
) {
    call_rule <- match.arg(call_rule)
    metric    <- match.arg(metric)
    order_by  <- match.arg(order_by)
    raw_side  <- match.arg(raw_side)

    PM   <- ashr::get_pm(m)
    PSD  <- ashr::get_psd(m)
    LFSR <- ashr::get_lfsr(m)

    gi <- if (is.character(gene)) match(gene, rownames(PM)) else as.integer(gene)
    if (is.na(gi) || gi < 1L || gi > nrow(PM)) stop("Gene not found/invalid.")

    pm  <- PM[gi, ];  psd <- PSD[gi, ];  ct <- names(pm)

    M <- if (metric == "lfsr") LFSR[gi, ] else {
        stopifnot(!is.null(metric_matrix), all(dim(metric_matrix) == dim(PM)))
        metric_matrix[gi, ]
    }
    sig_metric <- is.finite(M) & (M < thresh)

    z  <- stats::qnorm(0.5 + level/2)
    lo <- pm - z*psd
    hi <- pm + z*psd
    sig_ci <- (lo > 0 | hi < 0)

    sig <- switch(call_rule,
                  metric = sig_metric,
                  ci     = sig_ci,
                  both   = sig_metric & sig_ci,
                  either = sig_metric | sig_ci)

    # ----- ordering: ensures most positive at TOP, most negative at BOTTOM -----
    ord <- switch(order_by,
                  effect     = order(pm,      decreasing = TRUE),  # big -> small
                  abs_effect = order(abs(pm), decreasing = TRUE),
                  none       = seq_along(pm))
    pm <- pm[ord]; psd <- psd[ord]
    lo <- lo[ord]; hi <- hi[ord]; ct <- ct[ord]; sig <- sig[ord]

    raw_pm <- raw_lo <- raw_hi <- NULL
    if (isTRUE(plot_raw)) {
        if (is.null(m$data) || is.null(m$data$Bhat) || is.null(m$data$Shat)) {
            stop("plot_raw=TRUE but m$data$Bhat/Shat not found.")
        }
        raw_pm <- m$data$Bhat[gi, ][ord]
        raw_se <- m$data$Shat[gi, ][ord]
        #mask missing data with a standard error of 1e06.
        raw_pm[!is.finite(raw_pm) | raw_se>=1e6 ] <- NA_real_
        raw_se[!is.finite(raw_se) | raw_se <= 0 | raw_se>=1e6] <- NA_real_
        raw_lo <- raw_pm - z*raw_se
        raw_hi <- raw_pm + z*raw_se
    }

    seg_col <- ifelse(sig, ci_sig, ci_ns)
    pt_col  <- ifelse(sig, col_sig, col_ns)
    seg_lwd <- ifelse(sig, lwd_sig, lwd_ns)

    if (!is.null(colors)) {
        bg_col   <- if (is.na(colors$background)) "white" else colors$background
        axes_col <- if (!is.null(colors$axes)) colors$axes else "black"
        text_col <- if (!is.null(colors$text)) colors$text else "black"
    } else {
        bg_col <- "white"; axes_col <- "black"; text_col <- "black"
    }

    op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
    par(mar = c(5, 14, 3, 2), bg = bg_col, col.axis = axes_col, col.lab = axes_col, col.main = axes_col)

    # y positions: using rev(...) puts the first (largest) at the TOP
    y <- if (top_largest) rev(seq_along(pm)) else seq_along(pm)
    dir  <- if (raw_side == "below") -1 else 1
    yraw <- y + dir * raw_y_offset

    x_all <- c(lo, hi, 0)
    if (isTRUE(plot_raw)) x_all <- c(x_all, raw_lo, raw_hi)
    xlim <- range(x_all, finite = TRUE)

    y_all <- if (isTRUE(plot_raw)) c(y, yraw) else y
    ylim  <- range(y_all) + c(-0.5, 0.5)

    main <- if (is.character(gene)) gene else rownames(PM)[gi]

    plot(xlim, ylim, type = "n", bty = "n",
         xlab = if (isTRUE(plot_raw)) "Effect (posterior & raw ± CI)" else "Effect (posterior mean ± CI)",
         ylab = "", yaxt = "n", main = main)
    abline(v = 0, col = zero_col, lty = 2)

    if (isTRUE(plot_raw)) {
        ok <- is.finite(raw_lo) & is.finite(raw_hi) & is.finite(raw_pm)
        segments(raw_lo[ok], yraw[ok], raw_hi[ok], yraw[ok], col = ci_raw, lwd = lwd_raw)
        points(raw_pm[ok],  yraw[ok], pch = pch_raw, col = col_raw)
    }

    segments(lo, y, hi, y, col = seg_col, lwd = seg_lwd)
    points(pm, y, pch = pch_pts, col = pt_col)

    axis(2, at = y, labels = ct, las = 2, cex.axis = cex_names, col.axis = text_col, tick = FALSE)

    leg_txt <- c(paste0("passes '", call_rule, "' rule"), "not significant")
    leg_col <- c(col_sig, col_ns)
    leg_pch <- c(pch_pts, pch_pts)
    if (isTRUE(plot_raw)) {
        leg_txt <- c(leg_txt, sprintf("raw (%s) mean ± CI", raw_side))
        leg_col <- c(leg_col, col_raw)
        leg_pch <- c(leg_pch, pch_raw)
    }

    # moved from "topright" to "bottomright"
    legend("bottomright", legend = leg_txt, col = leg_col, pch = leg_pch, bty = "n", text.col = axes_col)

    invisible(list(sig = sig, order = ord, lo = lo, hi = hi, pm = pm,
                   raw_lo = raw_lo, raw_hi = raw_hi, raw_pm = raw_pm,
                   y = y, yraw = if (isTRUE(plot_raw)) yraw else NULL))
}

#' Pretty heatmap for a mash covariance/correlation component
#'
#' @param U   symmetric matrix (a mash component).
#' @param mode one of "cov", "cor", "trace":
#'   - "cov": plot U directly (magnitude + shape)
#'   - "cor": plot cov2cor(U) (shape only)
#'   - "trace": plot U / trace(U) (magnitude-aware shape; sum diag = 1)
#' @param cluster logical; reorder by hclust on 1 - correlation of the **plotted** matrix.
#' @param main title; if NULL, derived from mode.
#' @param pi optional mixing proportion to annotate (from get_estimated_pi).
#' @param usage optional total posterior weight used by this component (sum over genes).
#' @param cex axis label size.
#' @param palette colors.
plot_cov_image <- function(
        U,
        mode     = c("cov","cor","trace"),
        cluster  = TRUE,
        main     = NULL,
        pi       = NULL,
        usage    = NULL,
        cex      = 0.8,
        palette  = colorRampPalette(c("#ffffe5", "#fdb863", "#b2182b"))(200)
) {
    stopifnot(is.matrix(U), nrow(U) == ncol(U))
    mode <- match.arg(mode)
    labs <- if (is.null(colnames(U))) as.character(seq_len(ncol(U))) else colnames(U)

    # stabilizer for cov2cor
    epsI <- diag(1e-8 * mean(diag(U), na.rm = TRUE), nrow(U))
    # choose what to plot
    M_plot <- switch(mode,
                     cov   = U,
                     cor   = stats::cov2cor(U + epsI),
                     trace = U / sum(diag(U), na.rm = TRUE))
    if (is.null(main)) main <- switch(mode, cov="Covariance", cor="Correlation", trace="Trace-normalized covariance")

    # cluster on the plotted matrix using 1 - cor where appropriate
    if (cluster) {
        # compute correlation of rows of M_plot for ordering
        R_for_dist <- try(stats::cov2cor(M_plot + epsI), silent = TRUE)
        if (inherits(R_for_dist, "try-error") || anyNA(R_for_dist)) {
            D <- stats::dist(M_plot)  # fallback
        } else {
            D <- stats::as.dist(1 - R_for_dist)
        }
        ord <- stats::hclust(D, method = "average")$order
        M_plot <- M_plot[ord, ord, drop = FALSE]
        labs   <- labs[ord]
    }

    # annotate title with pi/usage if provided
    if (!is.null(pi) || !is.null(usage)) {
        bits <- c(if (!is.null(pi)) paste0("π=", signif(pi,3)),
                  if (!is.null(usage)) paste0("usage=", signif(usage,3)))
        main <- paste(main, paste(bits, collapse="  •  "))
    }

    k <- ncol(M_plot)
    op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
    par(mar = c(max(6, 0.6*max(nchar(labs))), max(6, 0.6*max(nchar(labs))), 3, 1))

    graphics::image(x = 1:k, y = 1:k, z = t(M_plot)[, k:1],
                    col = palette, xlab = "", ylab = "", axes = FALSE, main = main)
    axis(1, at = 1:k, labels = labs, las = 2, cex.axis = cex)
    axis(2, at = k:1, labels = labs, las = 2, cex.axis = cex)
    box()
}

#' Plot a mash covariance component by name (no clustering)
#'
#' @title Plot a mash covariance/correlation component
#' @description
#' Look up a component in a list of covariance matrices (e.g., the `Ulist`
#' you used to fit mash) and draw a labeled heatmap. Choose to plot the raw
#' covariance, the correlation (shape only), or a trace-normalized covariance
#' (sum of diagonal = 1). Optionally annotate the title with the component's
#' prior weight from `mixing_proportions` (e.g., `get_estimated_pi(m)`).
#'
#' @param component Character; name of the component to plot. Must match a name
#'   in `covariance_matrix_list` (exact or unique partial match).
#' @param covariance_matrix_list Named list of symmetric matrices (components).
#' @param mixing_proportions Optional named numeric vector (like
#'   `mashr::get_estimated_pi(m)`). If provided and `annotate_pi = TRUE`,
#'   the title includes the prior weight for `component`.
#' @param mode One of `"cov"`, `"cor"`, `"trace"`:
#'   - `"cov"`   plots the covariance (magnitude + shape);
#'   - `"cor"`   plots `cov2cor(U)` (shape only);
#'   - `"trace"` plots `U / trace(U)` (magnitude-aware, comparable scale).
#' @param annotate_pi Logical; if `TRUE` and `mixing_proportions` is given,
#'   append π for this component to the title. Default `FALSE`.
#' @param main Optional title; if `NULL`, a title is auto-generated from the
#'   component name and mode.
#' @param cex Axis label size.
#' @param palette Color vector for `image()`.
#'
#' @return Invisibly returns the matrix that was plotted (after any transform).
plot_mash_component <- function(
        component,
        covariance_matrix_list,
        mixing_proportions = NULL,
        mode = c("cov","cor","trace"),
        main = NULL,
        cex = 0.8,
        palette = grDevices::colorRampPalette(c("#ffffe5", "#fdb863", "#b2182b"))(200)
) {
    stopifnot(is.list(covariance_matrix_list))
    mode <- match.arg(mode)

    # resolve component name (allow unique partial match)
    nm <- names(covariance_matrix_list)
    if (is.null(nm)) stop("covariance_matrix_list must be a *named* list.")
    hit <- which(nm == component)
    if (!length(hit)) hit <- grep(component, nm, fixed = TRUE)
    if (!length(hit)) stop("Component '", component, "' not found.")
    if (length(hit) > 1L) stop("Component name '", component, "' is ambiguous: ",
                               paste(nm[hit], collapse = ", "))
    comp_name <- nm[hit]
    U <- covariance_matrix_list[[hit]]
    if (!is.matrix(U) || nrow(U) != ncol(U)) stop("Component is not a square matrix.")

    labs <- if (is.null(colnames(U))) as.character(seq_len(ncol(U))) else colnames(U)

    # small ridge for numerical stability when converting to correlation
    epsI <- diag(1e-8 * mean(diag(U), na.rm = TRUE), nrow(U))

    # matrix to plot (no clustering)
    M_plot <- switch(mode,
                     cov   = U,
                     cor   = stats::cov2cor(U + epsI),
                     trace = U / sum(diag(U), na.rm = TRUE))

    # title
    if (is.null(main)) {
        mode_label <- switch(mode, cov="Covariance", cor="Correlation", trace="Trace-normalized covariance")
        main <- paste0(comp_name, " — ", mode_label)
    }
    if (!is.null(mixing_proportions)) {
        if (!is.null(names(mixing_proportions)) && comp_name %in% names(mixing_proportions)) {
            main <- paste0(main, "  (π=", signif(mixing_proportions[[comp_name]], 3), ")")
        }
    }

    # plot
    # k <- ncol(M_plot)
    # op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
    # par(mar = c(max(6, 0.6 * max(nchar(labs))),  # bottom
    #             max(6, 0.6 * max(nchar(labs))),  # left
    #             3, 1))                            # top, right
    #
    # graphics::image(x = 1:k, y = 1:k, z = t(M_plot)[, k:1],
    #                 col = palette, xlab = "", ylab = "", axes = FALSE, main = main)
    # axis(1, at = 1:k, labels = labs, las = 2, cex.axis = cex)
    # axis(2, at = k:1, labels = labs, las = 2, cex.axis = cex)
    # box()
    # invisible(M_plot)

    k <- ncol(M_plot)

    # Convert to long format for ggplot2
    df <- as.data.frame(as.table(M_plot))
    names(df) <- c("Row", "Col", "Value")

    # Ensure ordering (so rows plot top→bottom like image())
    df$Row <- factor(df$Row, levels = labs)
    df$Col <- factor(df$Col, levels = labs)

    p=ggplot(df, aes(x = Col, y = Row, fill = Value)) +
        geom_tile() +
        scale_fill_gradientn(colors = palette) +
        labs(title = main, x = "", y = "") +
        theme_minimal(base_size = 14 * cex) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = rel(cex)),
            panel.grid = element_blank()
        )

    return (p)

}


# Core function ---------------------------------------------------------------
plot_top_eigenvectors <- function(S, k = 3, condition_names = NULL, bar_colors = NULL) {

    # Ordinal helper for plot titles
    .ordinal <- function(i) {
        s <- c("st","nd","rd","th")
        idx <- if (i %% 100L %in% 11:13) 4L else min(max(i %% 10L, 1L), 4L)
        paste0(i, s[idx])
    }


    if (!is.matrix(S) || nrow(S) != ncol(S))
        stop("S must be a square covariance matrix")

    # Make sure it's exactly symmetric (tolerant to tiny numeric noise)
    S <- 0.5 * (S + t(S))

    # Eigen-decomposition (symmetric)
    ev <- eigen(S, symmetric = TRUE)
    o  <- order(ev$values, decreasing = TRUE)
    vals <- ev$values[o]
    vecs <- ev$vectors[, o, drop = FALSE]

    # Proportion of variance explained
    # (guard against tiny negative eigenvalues from numerics)
    vals_pos <- pmax(vals, 0)
    pve <- vals_pos / sum(vals_pos)

    # Keep top k
    k <- min(k, ncol(vecs))
    vals_k <- vals[seq_len(k)]
    pve_k  <- pve[seq_len(k)]
    vecs_k <- vecs[, seq_len(k), drop = FALSE]

    # Default names/colors
    if (is.null(condition_names)) condition_names <- rownames(S)
    if (is.null(condition_names)) condition_names <- as.character(seq_len(nrow(S)))
    if (!is.null(bar_colors) && length(bar_colors) != length(condition_names))
        stop("bar_colors must be length nrow(S) or NULL")

    # Flip eigenvector signs so the largest |loading| is positive (for consistency)
    for (j in seq_len(k)) {
        jmax <- which.max(abs(vecs_k[, j]))
        if (vecs_k[jmax, j] < 0) vecs_k[, j] <- -vecs_k[, j]
    }

    # Build ggplots (one per eigenvector)
    plots <- vector("list", k)
    for (j in seq_len(k)) {
        df <- data.frame(
            condition = factor(condition_names, levels = condition_names),
            loading   = vecs_k[, j],
            stringsAsFactors = FALSE
        )

        p <- ggplot(df, aes(x = condition, y = loading)) +
            geom_col(width = 0.8, aes(fill = condition), show.legend = FALSE) +
            geom_hline(yintercept = 0) +
            labs(
                title = sprintf("%s eigenvector (PVE = %.1f%%)", .ordinal(j), 100 * pve_k[j]),
                x = NULL, y = "Loading"
            ) +
            theme_minimal(base_size = 13) +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                plot.title  = element_text(face = "bold")
            )

        if (!is.null(bar_colors)) {
            p <- p + scale_fill_manual(values = setNames(bar_colors, condition_names))
        }

        plots[[j]] <- p
    }

    return (plots)
}




plot_sharing<-function (sharing, strTitle="") {
    # 1) Build a distance from similarity
    d   <- as.dist(1 - sharing)           # dissimilarity

    # 2) Cluster (one tree for both axes)
    hc  <- hclust(d, method = "average")  # or "ward.D2"

    # 3) Reorder rows/cols with that tree
    ord <- hc$order
    S   <- sharing[ord, ord]              # <-- this is the S I used in plots

    # (optional) Upper-triangle only:
    S_up <- S; S_up[lower.tri(S_up)] <- NA

    # choose fixed limits so figures are comparable across runs
    lo  <- 0.20     # lower bound shown in the paper
    mid <- 0.65     # center of the scale (tweak if you like)
    hi  <- 1.00

    # blue – white – red (paper-ish: RdYlBu reversed)
    col_fun <- colorRamp2(
        c(lo,  mid,  hi),
        c("#2c7fb8", "#ffffbf", "#d73027")
    )

    d  <- as.dist(1 - S)
    hc <- hclust(d, "average")

    Heatmap(
        S,
        name = "pairwise sharing",
        col = col_fun,
        na_col = "white",
        cluster_rows = as.dendrogram(hc),
        cluster_columns = as.dendrogram(hc),
        row_names_side = "left",
        show_row_dend=FALSE,
        show_row_names=TRUE,
        row_dend_side = "right",
        column_names_rot = 45,
        show_column_names=FALSE,
        column_title = strTitle,
        column_title_gp=gpar(fontsize=18),
        heatmap_legend_param = list(
            at = seq(lo, hi, by = 0.1),           # tick marks
            labels = sprintf("%.1f", seq(lo, hi, 0.1)),
            title = ""
        ),
        rect_gp = gpar(col = NA)
    )
}


#######################
# READING INPUT DATA
#######################


#' Build Bhat/Shat/FDR matrices from per–cell-type DE tables
#'
#' @description
#' Construct matrices for `mashr` from a list of differential expression tables
#' (one per cell type/condition). You can choose to include the **union** of all
#' genes seen across tables (filling missing entries with beta=0 and a very large
#' SE so they have no influence on mash), or restrict to the **intersection**
#' (genes present in every table).
#'
#' For each table, the function reads a coefficient column (e.g. `logFC`),
#' a t-statistic column (e.g. `t`), and an FDR column (e.g. `adj.P.Val`), and
#' computes `Shat = |beta / t|` (so `t = beta / se` as in limma).
#'
#' @param lst Named `list` of data frames (one per cell type). Each must have
#'   rownames as gene IDs and the columns specified by `coef_col`, `t_col`,
#'   and `fdr_col`.
#' @param coef_col Character scalar. Column name for the effect size (beta),
#'   e.g. `"logFC"`. Default `"logFC"`.
#' @param t_col Character scalar. Column name for the t-statistic, e.g. `"t"`.
#'   Default `"t"`.
#' @param fdr_col Character scalar. Column name for the FDR/q-value,
#'   e.g. `"adj.P.Val"`. Default `"adj.P.Val"`.
#' @param gene_mode One of `c("union","intersect")`. Use `"union"` to include
#'   all genes observed in any table (missing entries are replaced with
#'   `beta=0` and `SE=big_se`). Use `"intersect"` to keep only genes present
#'   in **every** table. Default `"union"`.
#' @param big_se Numeric scalar. Standard error to use when filling missing or
#'   invalid entries. A very large value ensures such entries contribute ~0
#'   information to mash. Default `1e6`.
#' @param fill_missing_fdr Numeric scalar used to fill missing FDR entries in
#'   the returned `FDR` matrix (for plotting/logic only; **not** used by mash).
#'   Default `NA_real_`.
#'
#' @return A list with components:
#' \itemize{
#'   \item `Bhat` — matrix [genes × cell types] of betas.
#'   \item `Shat` — matrix [genes × cell types] of standard errors (`|beta/t|`).
#'   \item `FDR`  — matrix [genes × cell types] of FDR values (filled with
#'         `fill_missing_fdr` where absent).
#'   \item `missing_mask` — logical matrix [genes × cell types] indicating where
#'         input was missing/invalid and got replaced by `beta=0, SE=big_se`.
#' }
#'
#' @details
#' - With `gene_mode = "union"`, the function guarantees rectangular matrices by
#'   inserting neutral placeholders for missing cells (0, `big_se`). This is
#'   a common trick so that mash can ingest more genes without biasing fits.
#' - With `gene_mode = "intersect"`, all retained genes have entries in every
#'   column (aside from pathological `t=0`; those are still treated as missing
#'   and assigned `SE=big_se`).
#'
#' @export
make_mash_inputs<- function(
        lst,
        coef_col = "logFC",
        t_col    = "t",
        fdr_col  = "adj.P.Val",
        gene_mode = c("union", "intersect"),
        big_se   = 1e6,
        fill_missing_fdr = NA_real_
) {
    gene_mode <- match.arg(gene_mode)

    stopifnot(is.list(lst), length(lst) >= 2L)
    # Give unnamed list elements sensible names
    if (is.null(names(lst)) || any(names(lst) == "")) {
        names(lst) <- paste0("CT", seq_along(lst))
    }

    # ---- determine the gene set (union or intersection) ----
    gene_lists <- lapply(lst, rownames)
    if (gene_mode == "union") {
        genes <- sort(unique(unlist(gene_lists, use.names = FALSE)))
    } else { # "intersect"
        genes <- Reduce(intersect, gene_lists)
        genes <- sort(unique(genes))
    }
    G <- length(genes); C <- length(lst)
    if (G == 0L) stop("No gene names found after applying gene_mode = '", gene_mode, "'.")

    # Preallocate outputs
    Bhat <- matrix(NA_real_, G, C, dimnames = list(genes, names(lst)))
    Shat <- matrix(NA_real_, G, C, dimnames = list(genes, names(lst)))
    FDR  <- matrix(NA_real_, G, C, dimnames = list(genes, names(lst)))

    # ---- fill matrices column by column ----
    for (j in seq_along(lst)) {
        df <- lst[[j]]

        # Validate required columns
        required <- c(coef_col, t_col, fdr_col)
        if (!all(required %in% colnames(df))) {
            stop("Missing required columns in lst[['", names(lst)[j], "']]: ",
                 paste(setdiff(required, colnames(df)), collapse = ", "))
        }

        # Which genes appear in this table (restricted to our gene set)
        g <- intersect(genes, rownames(df))
        if (!length(g)) next

        beta <- as.numeric(df[g, coef_col])
        tt   <- as.numeric(df[g, t_col])
        fdr  <- as.numeric(df[g, fdr_col])

        # Convert t-stat to SE (limma: t = beta / se)
        se <- abs(beta / tt)  # Inf/NA handled below

        # Place into the right rows for this column
        Bhat[g, j] <- beta
        Shat[g, j] <- se
        FDR [g, j] <- fdr
    }

    # ---- handle missing/invalid entries ----
    # Missing beta/SE, non-finite SE, or non-positive SE -> neutral placeholders
    missing_mask <- !is.finite(Bhat) | !is.finite(Shat) | (Shat <= 0)
    if (any(missing_mask)) {
        Bhat[missing_mask] <- 0
        Shat[missing_mask] <- big_se
    }

    # Missing FDR gets a chosen fill (for plotting; mash doesn't use FDR)
    nas_fdr <- !is.finite(FDR)
    if (any(nas_fdr)) FDR[nas_fdr] <- fill_missing_fdr

    list(Bhat = Bhat, Shat = Shat, FDR = FDR, missing_mask = missing_mask)
}

#' Drop common suffix from file names
#' Useful to infer the data set names from differential expression files
#' @param files A list of file names to process
#' @return A vector with the same length as the input files containing the name of the file without
#' the common suffix.
drop_common_suffix <- function(files) {
    # Reverse each filename so suffix becomes a prefix
    rev_files <- sapply(files, function(x) paste0(rev(strsplit(x, NULL)[[1]]), collapse = ""))

    # Find common prefix of reversed strings
    common_rev_prefix <- Reduce(function(a, b) {
        i <- 1
        while (i <= nchar(a) && i <= nchar(b) && substr(a, i, i) == substr(b, i, i)) {
            i <- i + 1
        }
        substr(a, 1, i - 1)
    }, rev_files)

    # Reverse back to get the suffix
    common_suffix <- paste0(rev(strsplit(common_rev_prefix, NULL)[[1]]), collapse = "")

    # Drop suffix from each filename
    sub(paste0(common_suffix, "$"), "", files)
}
