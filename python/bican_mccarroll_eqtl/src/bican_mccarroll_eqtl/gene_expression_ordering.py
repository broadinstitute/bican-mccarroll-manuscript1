"""Order genes within K-means clusters by expression correlation."""

import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

from .kmeans_heatmap import DEFAULT_CELLTYPE_ORDER


def order_genes_by_correlation(median_expression_path,
                               cluster_assignments_path,
                               output_path=None,
                               celltype_order=None):
    """Order genes within each K-means cluster by hierarchical clustering
    on expression correlation.

    Within each cluster, genes are reordered so that genes with similar
    expression patterns across cell types are adjacent.  Genes with zero
    variance (no expression variation) are placed at the end of their
    cluster.

    Parameters
    ----------
    median_expression_path : str
        Path to the median expression TSV (output of
        ``get_heatmap_index_snp_median_expression``).  First column
        should be ``Gene``; remaining columns are cell type / region
        groups.
    cluster_assignments_path : str
        Path to the cluster assignments TSV (output of
        ``run_kmeans_heatmap`` with ``cluster_assignments_output_path``).
        Must have columns ``gene`` and ``cluster``.
    output_path : str or None
        If provided, saves the ordered gene list as a single-column TSV
        with header ``gene``.
    celltype_order : list of str or None
        Column order for cell types used when computing correlations.
        Defaults to DEFAULT_CELLTYPE_ORDER.

    Returns
    -------
    ordered_genes : list of str
        Gene identifiers ordered by cluster then by expression
        correlation within each cluster.
    """
    if celltype_order is None:
        celltype_order = DEFAULT_CELLTYPE_ORDER

    # Read inputs
    median_expr_df = pd.read_csv(
        median_expression_path, sep="\t"
    ).set_index("Gene")
    cluster_df = pd.read_csv(cluster_assignments_path, sep="\t")
    cluster_map = pd.Series(
        cluster_df["cluster"].values, index=cluster_df["gene"]
    ).astype(int)

    # Use only columns present in the expression matrix
    cols_in_order = [ct for ct in celltype_order if ct in median_expr_df.columns]

    # Determine cluster order from the categorical ordering in the file
    cluster_order = sorted(cluster_map.unique())

    ordered_genes = []

    for k in cluster_order:
        genes_k = cluster_map[cluster_map == k].index.tolist()

        # Keep only genes present in the expression matrix
        genes_k = [g for g in genes_k if g in median_expr_df.index]
        if len(genes_k) == 0:
            continue

        if len(genes_k) == 1:
            ordered_genes.extend(genes_k)
            continue

        X = median_expr_df.loc[genes_k, cols_in_order].astype(float).values

        # Separate zero-variance rows (can't compute correlation)
        row_var = X.var(axis=1)
        keep = row_var > 0
        genes_keep = [g for g, ok in zip(genes_k, keep) if ok]
        genes_drop = [g for g, ok in zip(genes_k, keep) if not ok]

        if len(genes_keep) <= 1:
            ordered_genes.extend(genes_keep + genes_drop)
            continue

        X_keep = X[keep, :]

        d = pdist(X_keep, metric="correlation")
        Z = linkage(d, method="average")
        leaf_idx = leaves_list(Z)
        ordered_keep = [genes_keep[i] for i in leaf_idx]

        ordered_genes.extend(ordered_keep + genes_drop)

    if output_path is not None:
        out_df = pd.DataFrame({"gene": ordered_genes})
        out_df.to_csv(output_path, sep="\t", index=False)

    return ordered_genes
