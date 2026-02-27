"""Order genes within K-means clusters by expression correlation."""

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

from .kmeans_heatmap import DEFAULT_CELLTYPE_ORDER, DEFAULT_CELLTYPE_LABEL_MAP


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


def plot_expression_heatmap(median_expression_path,
                            ordered_genes_path,
                            cluster_assignments_path,
                            cluster_order,
                            celltype_order=None,
                            celltype_label_map=None,
                            output_path=None,
                            figsize=(12, 9),
                            dpi=400):
    """Plot median expression heatmap with genes ordered by correlation.

    Normalises expression per gene (divides by 95th percentile, clips to
    [0, 1]) and plots a heatmap grouped by K-means cluster using a
    truncated YlGnBu colourmap.

    Parameters
    ----------
    median_expression_path : str
        Path to the median expression TSV (output of
        ``get_heatmap_index_snp_median_expression``).  First column is
        ``Gene``; remaining columns are cell type / region groups.
    ordered_genes_path : str
        Path to the ordered gene list TSV (output of
        ``order_genes_by_correlation``).  Single column ``gene``.
    cluster_assignments_path : str
        Path to the cluster assignments TSV (output of
        ``run_kmeans_heatmap``).  Columns ``gene`` and ``cluster``.
    cluster_order : list of int
        Display order for clusters (e.g., ``[8, 1, 3, 5, ...]``).
    celltype_order : list of str or None
        Column order for cell types.  Defaults to DEFAULT_CELLTYPE_ORDER.
    celltype_label_map : dict or None
        Map from raw cell type names to display labels.
        Defaults to DEFAULT_CELLTYPE_LABEL_MAP.
    output_path : str or None
        If provided, saves the figure to this path.
    figsize : tuple
        Figure size (width, height) in inches.
    dpi : int
        Resolution for saved figure.

    Returns
    -------
    fig : matplotlib Figure
    """
    if celltype_order is None:
        celltype_order = DEFAULT_CELLTYPE_ORDER
    if celltype_label_map is None:
        celltype_label_map = DEFAULT_CELLTYPE_LABEL_MAP

    # Read inputs
    median_expr_df = pd.read_csv(
        median_expression_path, sep="\t"
    ).set_index("Gene")
    ordered_genes_df = pd.read_csv(ordered_genes_path, sep="\t")
    ordered_genes = ordered_genes_df["gene"].tolist()
    cluster_df = pd.read_csv(cluster_assignments_path, sep="\t")
    cluster_map = pd.Series(
        cluster_df["cluster"].values, index=cluster_df["gene"]
    ).astype(int)

    cols_in_order = [ct for ct in celltype_order if ct in median_expr_df.columns]

    # Normalise expression per gene: divide by p95, clip to [0, 1]
    eps = 1e-8

    def _normalize_row(row):
        p95 = np.percentile(row, 95)
        if p95 < eps:
            return np.zeros_like(row, dtype=float)
        return np.clip(row / p95, 0, 1)

    normalized_df = median_expr_df.apply(_normalize_row, axis=1)

    # Build AnnData in gene order
    expr_matrix = normalized_df.loc[ordered_genes, cols_in_order].values
    adata_expr = ad.AnnData(expr_matrix)
    adata_expr.obs_names = pd.Index(ordered_genes)
    adata_expr.var_names = pd.Index(cols_in_order)

    # Assign cluster labels in desired order
    adata_expr.obs["gene_clusters"] = pd.Categorical(
        cluster_map.reindex(ordered_genes).astype(int),
        categories=cluster_order,
        ordered=True,
    )

    # Rename cell type columns to display labels
    adata_expr.var["var_name_original"] = adata_expr.var_names.copy()
    adata_expr.var_names = pd.Index([
        celltype_label_map[v] for v in adata_expr.var["var_name_original"]
    ])

    # Truncated YlGnBu colormap (0 to 0.6 range)
    base_cmap = mpl.colormaps["YlGnBu"]
    cmap = mpl.colors.ListedColormap(base_cmap(np.linspace(0, 0.6, 256)))

    sc.pl.heatmap(
        adata_expr,
        var_names=adata_expr.var_names,
        groupby="gene_clusters",
        dendrogram=False,
        cmap=cmap,
        swap_axes=True,
        figsize=figsize,
        show=False,
    )

    for text in plt.gcf().findobj(match=plt.Text):
        text.set_fontsize(16)

    plt.suptitle(
        "Expression patterns of genes grouped by eQTL effect-size clustering",
        fontsize=20,
        fontweight="bold",
        x=0.55,
        y=0.95,
    )

    fig = plt.gcf()

    if output_path is not None:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")

    return fig
