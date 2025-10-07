# Input:
#   M: numeric matrix (genes x conditions), log2 effect sizes
# Outputs:
#   list with clusters per gene, Seurat object, UMAP coords
#
# Notes:
# - No normalization; we z-score each gene across conditions.
# - Leiden via Seurat (algorithm = 4).
# - Euclidean on PCs approximates correlation after row z-scoring.


library(Seurat)

in_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/old/sex_age/cell_type_region_interaction_absolute_effects"
file_pattern="age"
cellTypeListFile="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metadata/mash_cell_type_list_simple.txt"

run_me<-function (in_dir, file_pattern="age", cellTypeListFile=NULL) {
    d=parse_de_inputs(in_dir, file_pattern, cellTypeListFile)
    mash_inputs_union<-make_mash_inputs(d, coef_col = "logFC", t_col = "t", fdr_col = "adj.P.Val", gene_mode="union")
    M=mash_inputs_union$Bhat

}


# Row-wise zscore (per gene across conditions)
zscore_rows <- function(mat) {
    m <- rowMeans(mat, na.rm = TRUE)
    s <- sqrt(rowMeans((mat - m)^2, na.rm = TRUE))
    s[s == 0 | !is.finite(s)] <- 1
    sweep(sweep(mat, 1, m, "-"), 1, s, "/")
}

#n_pcs should be lower than the number of conditions, and probably more like the number of cell types
gene_leiden <- function(M, n_pcs = 20, k = 20, resolution = 1.0, umap = TRUE) {
    stopifnot(is.matrix(M))
    stopifnot(!is.null(rownames(M)), !is.null(colnames(M)))  # need gene + condition names

    # 1. Row-zscore genes across conditions
    Mz <- zscore_rows(M)              # genes × conditions

    # 2. Build Seurat object where genes are "cells"
    mm <- t(Mz)                       # conditions × genes
    obj <- CreateSeuratObject(counts = mm, assay = "eff")
    DefaultAssay(obj) <- "eff"

    # 3. Seurat requires a 'data' slot. Copy counts → data so features match
    obj <- SetAssayData(obj, assay = "eff", slot = "data",
                        new.data = GetAssayData(obj, assay = "eff", slot = "counts"))

    # 4. Use all conditions as features
    VariableFeatures(obj) <- rownames(obj)

    # 5. Populate scale.data (no extra scaling since already zscored)
    obj <- ScaleData(obj, features = rownames(obj),
                     do.center = FALSE, do.scale = FALSE, verbose = FALSE)

    # 6. PCA → neighbors → Leiden → UMAP
    obj <- RunPCA(obj, features = rownames(obj), npcs = n_pcs, verbose = FALSE)
    obj <- FindNeighbors(obj, dims = 1:n_pcs, k.param = k, verbose = FALSE)
    obj <- FindClusters(obj, resolution = resolution, algorithm = 4, verbose = FALSE, random.seed = 1)
    obj <- RunUMAP(obj, dims = 1:n_pcs, verbose = FALSE)

    # Return clusters named by gene
    clusters <- obj$seurat_clusters
    names(clusters) <- colnames(obj)  # these are gene names

    list(
        clusters = clusters,
        seurat   = obj,
        umap     = if (umap) Embeddings(obj, "umap") else NULL
    )
}

