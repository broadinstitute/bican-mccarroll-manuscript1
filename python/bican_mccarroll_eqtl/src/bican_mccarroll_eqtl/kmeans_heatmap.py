"""K-means clustering and heatmap visualization of eQTL effect-size profiles."""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


# Default cell type display order and label map
DEFAULT_CELLTYPE_ORDER = [
    "MSN_D2_matrix__CaH",
    "MSN_D1_matrix__CaH",
    "MSN_D1_striosome__CaH",
    "MSN_D2_striosome__CaH",
    "GABA_MGE_CAP__CaH",
    "GABA_MGE_DFC__DFC",
    "GABA_CGE_DFC__DFC",
    "glutamatergic_L23IT__DFC",
    "glutamatergic_L5IT__DFC",
    "astrocyte__DFC",
    "astrocyte__CaH",
    "oligodendrocyte__DFC",
    "oligodendrocyte__CaH",
    "OPC__DFC",
    "OPC__CaH",
    "microglia__DFC",
    "microglia__CaH",
]

DEFAULT_CELLTYPE_LABEL_MAP = {
    "MSN_D2_matrix__CaH": "MSN D2 matrix (CaH)",
    "MSN_D1_matrix__CaH": "MSN D1 matrix (CaH)",
    "MSN_D1_striosome__CaH": "MSN D1 striosome (CaH)",
    "MSN_D2_striosome__CaH": "MSN D2 striosome (CaH)",
    "GABA_MGE_CAP__CaH": "MGE-derived GABAergic (CaH)",
    "GABA_MGE_DFC__DFC": "MGE-derived GABAergic (DFC)",
    "GABA_CGE_DFC__DFC": "CGE-derived GABAergic (DFC)",
    "glutamatergic_L23IT__DFC": "Glutamatergic L2/3 IT (DFC)",
    "glutamatergic_L5IT__DFC": "Glutamatergic L5 IT (DFC)",
    "astrocyte__DFC": "Astrocyte (DFC)",
    "astrocyte__CaH": "Astrocyte (CaH)",
    "oligodendrocyte__DFC": "Oligodendrocyte (DFC)",
    "oligodendrocyte__CaH": "Oligodendrocyte (CaH)",
    "OPC__DFC": "OPC (DFC)",
    "OPC__CaH": "OPC (CaH)",
    "microglia__DFC": "Microglia (DFC)",
    "microglia__CaH": "Microglia (CaH)",
}

# 17q21.31 inversion region genes (H1/H2 haplotype)
H1H2_GENES = {
    "LRRC37A", "LRRC37A2", "KANSL1", "ARL17", "ARL17P1", "NSF", "CRHR1",
    "SPPL2C", "STH", "MAPT", "PLEKHM1",
}


def _load_slope_matrix(input_path, celltype_order=None):
    """Read the index-SNP slope matrix and return an AnnData (genes x cell types).

    Parameters
    ----------
    input_path : str
        Path to the TSV output of get_index_snp_slope_matrix_with_median_impute.
    celltype_order : list of str or None
        Column order for cell types. Defaults to DEFAULT_CELLTYPE_ORDER.

    Returns
    -------
    adata : AnnData
        Genes in obs, cell types in var.
    input_matrix : DataFrame
        Original matrix (still has variant_id column).
    """
    if celltype_order is None:
        celltype_order = DEFAULT_CELLTYPE_ORDER

    input_matrix = pd.read_csv(input_path, sep="\t", index_col=0)
    slope_matrix = input_matrix.drop(columns=["variant_id"])
    adata = ad.AnnData(slope_matrix)
    adata = adata[:, celltype_order]
    return adata, input_matrix


def run_k_selection(input_path, k_range=(5, 25), random_state=32,
                    celltype_order=None):
    """Run K-means for a range of K values and compute selection metrics.

    Parameters
    ----------
    input_path : str
        Path to the index-SNP slope matrix TSV.
    k_range : tuple of (int, int)
        (min_k, max_k) range to evaluate (exclusive upper bound).
    random_state : int
        Random seed for KMeans.
    celltype_order : list of str or None
        Column order for cell types.

    Returns
    -------
    silhouette_df : DataFrame
        Columns: k, silhouette_score, wcss.
    """
    adata, _ = _load_slope_matrix(input_path, celltype_order)

    rows = []
    for k in range(k_range[0], k_range[1]):
        km = KMeans(n_clusters=k, random_state=random_state).fit(adata.X)
        score = silhouette_score(adata.X, km.labels_)
        rows.append({"k": k, "silhouette_score": score, "wcss": km.inertia_})

    return pd.DataFrame(rows)


def plot_k_selection(silhouette_df, output_path=None):
    """Plot silhouette score and WCSS for K selection.

    Parameters
    ----------
    silhouette_df : DataFrame
        Output of run_k_selection with columns k, silhouette_score, wcss.
    output_path : str or None
        If provided, saves the figure to this path.

    Returns
    -------
    fig : Figure
    """
    ks = silhouette_df["k"].values
    sil = silhouette_df["silhouette_score"].values
    wcss = silhouette_df["wcss"].values

    fig, ax1 = plt.subplots(figsize=(6.5, 4.5))

    ax1.plot(ks, sil, marker="o", markersize=4, linewidth=1.8, color="black",
             label="Silhouette score")
    ax1.set_xlabel("Number of clusters (K)")
    ax1.set_ylabel("Silhouette score", color="black")
    ax1.tick_params(axis="y", colors="black")
    ax1.spines["top"].set_visible(False)

    ax2 = ax1.twinx()
    ax2.plot(ks, wcss, marker="o", markersize=4, linewidth=1.8, color="tab:blue",
             label="WCSS")
    ax2.set_ylabel("Within-cluster sum of squares (WCSS)", color="tab:blue")
    ax2.tick_params(axis="y", colors="tab:blue")
    ax2.spines["top"].set_visible(False)

    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(handles1 + handles2, labels1 + labels2, frameon=False, loc="best")

    fig.tight_layout()

    if output_path is not None:
        fig.savefig(output_path, dpi=400, bbox_inches="tight")

    return fig


def run_kmeans_heatmap(input_path, K, desired_order=None, random_state=32,
                       celltype_order=None, celltype_label_map=None,
                       heatmap_output_path=None, cluster_counts_output_path=None,
                       cluster_assignments_output_path=None):
    """Run K-means clustering and generate the effect-size heatmap.

    Intended workflow:
      1. Run with ``desired_order=None`` to produce an initial heatmap with
         clusters in numeric order.  Inspect the heatmap and decide on a
         cluster ordering that places the diagonal color blocks in the
         desired sequence.
      2. Re-run with ``desired_order=[6, 8, 1, ...]`` to produce the final
         ordered heatmap.

    Parameters
    ----------
    input_path : str
        Path to the index-SNP slope matrix TSV.
    K : int
        Number of clusters.
    desired_order : list of int or None
        Cluster display order (e.g., [6, 8, 1, 4, 2, 3, 10, 9, 7, 5, 0]).
        If None, clusters are displayed in numeric order (0, 1, 2, ...).
    random_state : int
        Random seed for KMeans.
    celltype_order : list of str or None
        Column order. Defaults to DEFAULT_CELLTYPE_ORDER.
    celltype_label_map : dict or None
        Map from raw cell type names to display labels.
        Defaults to DEFAULT_CELLTYPE_LABEL_MAP.
    heatmap_output_path : str or None
        If provided, saves the heatmap PNG to this path.
    cluster_counts_output_path : str or None
        If provided, saves gene cluster counts TSV to this path.
    cluster_assignments_output_path : str or None
        If provided, saves a TSV with columns ``gene`` and ``cluster``
        mapping each gene to its K-means cluster label.

    Returns
    -------
    adata : AnnData
        With gene_clusters in obs.
    input_matrix : DataFrame
        Original matrix with variant_id (useful for downstream).
    """
    if celltype_label_map is None:
        celltype_label_map = DEFAULT_CELLTYPE_LABEL_MAP

    adata, input_matrix = _load_slope_matrix(input_path, celltype_order)

    # Rename cell type labels for display
    adata.var["var_name_original"] = adata.var_names.copy()
    adata.var_names = [celltype_label_map[v] for v in adata.var_names]

    # K-means clustering
    km = KMeans(n_clusters=K, random_state=random_state).fit(adata.X)
    adata.obs["gene_clusters"] = pd.Categorical(km.labels_.astype(int))

    if desired_order is None:
        desired_order = sorted(adata.obs["gene_clusters"].cat.categories)

    adata.obs["gene_clusters"] = adata.obs["gene_clusters"].cat.reorder_categories(
        desired_order, ordered=True
    )

    # --- Output 1: cluster assignments ---
    if cluster_assignments_output_path is not None:
        assignments_df = pd.DataFrame({
            "gene": adata.obs_names,
            "cluster": adata.obs["gene_clusters"].values,
        })
        assignments_df.sort_values("gene").to_csv(
            cluster_assignments_output_path, sep="\t", index=False
        )

    # --- Output 2: cluster counts ---
    if cluster_counts_output_path is not None:
        counts = (
            adata.obs["gene_clusters"]
            .value_counts()
            .sort_index()
            .reset_index()
        )
        counts.columns = ["gene_cluster", "n_genes"]
        counts.to_csv(cluster_counts_output_path, sep="\t", index=False)

    # --- Output 3: heatmap ---
    sc.pl.heatmap(
        adata,
        var_names=adata.var_names,
        groupby="gene_clusters",
        dendrogram=False,
        cmap="seismic",
        vcenter=0,
        vmin=-2,
        vmax=2,
        swap_axes=True,
        figsize=(12, 9),
        show=False,
    )

    for text in plt.gcf().findobj(match=plt.Text):
        text.set_fontsize(16)

    plt.suptitle(
        "K-means clustering of lead eQTL effect-size profiles across cell types",
        fontsize=20,
        fontweight="bold",
        x=0.55,
        y=0.95,
    )

    if heatmap_output_path is not None:
        plt.savefig(heatmap_output_path, dpi=300, bbox_inches="tight")

    plt.close()

    return adata, input_matrix


def build_fisher_contingency_table(adata, input_matrix, coloc_path,
                                   contingency_output_path=None,
                                   coloc_cluster_output_path=None):
    """Build Fisher's exact test contingency table for coloc gene enrichment.

    For each K-means cluster, counts colocalized vs non-colocalized genes
    inside vs outside the cluster. Collapses 17q21.31 (H1/H2) genes within
    each cluster to avoid inflating enrichment.

    Parameters
    ----------
    adata : AnnData
        Output of run_kmeans_heatmap (must have gene_clusters in obs).
    input_matrix : DataFrame
        Original slope matrix with variant_id column.
    coloc_path : str
        Path to colocalized genes TSV (must have phenotype_id column).
    contingency_output_path : str or None
        If provided, saves the contingency counts TSV.
    coloc_cluster_output_path : str or None
        If provided, saves coloc gene-to-cluster mapping TSV.

    Returns
    -------
    contingency_df : DataFrame
        One row per cluster with contingency counts.
    """
    coloc_df = pd.read_csv(coloc_path, sep="\t")
    coloc_set_raw = coloc_df["phenotype_id"].unique()

    # Collapse H1/H2 genes: keep at most one per cluster
    coloc_set = set(coloc_set_raw)
    clusters = sorted(adata.obs["gene_clusters"].unique())

    for cl in clusters:
        in_cluster = adata.obs["gene_clusters"] == cl
        genes_in_cluster = list(adata.obs_names[in_cluster])

        h1h2_in_cluster = [g for g in genes_in_cluster
                           if (g in H1H2_GENES) and (g in coloc_set)]
        if len(h1h2_in_cluster) > 1:
            for g in h1h2_in_cluster[1:]:
                coloc_set.remove(g)

    heatmap_genes = set(adata.obs_names)
    coloc_set_filtered = coloc_set.intersection(heatmap_genes)

    # Build contingency table
    is_coloc_all = adata.obs_names.isin(coloc_set_filtered)
    coloc_total = int(is_coloc_all.sum())
    noncoloc_total = int((~is_coloc_all).sum())

    rows = []
    for cl in clusters:
        in_cluster = adata.obs["gene_clusters"] == cl
        a = int((in_cluster & is_coloc_all).sum())
        b = int((in_cluster & ~is_coloc_all).sum())
        c = coloc_total - a
        d = noncoloc_total - b

        rows.append({
            "cluster": str(cl),
            "coloc_in_cluster": a,
            "noncoloc_in_cluster": b,
            "coloc_not_in_cluster": c,
            "noncoloc_not_in_cluster": d,
            "coloc_total": coloc_total,
            "noncoloc_total": noncoloc_total,
            "genes_in_cluster": a + b,
            "genes_total": coloc_total + noncoloc_total,
        })

    contingency_df = pd.DataFrame(rows)

    if contingency_output_path is not None:
        contingency_df.to_csv(contingency_output_path, sep="\t", index=False)

    # Coloc gene -> cluster mapping
    if coloc_cluster_output_path is not None:
        coloc_subset_df = adata.obs.loc[
            adata.obs_names.intersection(coloc_set_filtered)
        ].copy()
        gene_cluster_df = coloc_subset_df[["gene_clusters"]].copy()
        gene_cluster_df["phenotype_id"] = coloc_subset_df.index
        gene_cluster_df = gene_cluster_df[["phenotype_id", "gene_clusters"]].rename(
            columns={"gene_clusters": "cluster"}
        )
        gene_cluster_df.to_csv(coloc_cluster_output_path, sep="\t", index=False)

    return contingency_df
