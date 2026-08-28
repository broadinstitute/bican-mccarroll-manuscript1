"""Match K-means gene-cluster IDs across two runs by shared gene membership.

Two runs of `run_kmeans_heatmap` on different input matrices (e.g. the
manuscript's main and supplemental cell-type sets) assign genes to raw
K-means cluster IDs that carry no correspondence to each other on their own
-- the same biological cluster can come out labeled `3` in one run and `7`
in another. This module finds the one-to-one correspondence that maximizes
shared gene membership between a "raw" run's clusters and a "reference"
run's (already-ordered) clusters, so the raw run's `--desired-order` can be
set to reproduce the reference run's visual/biological ordering.
"""

import pandas as pd
from scipy.optimize import linear_sum_assignment


def match_cluster_order(raw_assignments_path, reference_assignments_path):
    """Match raw cluster IDs to a reference run's ordering by gene overlap.

    Parameters
    ----------
    raw_assignments_path : str
        Path to a `cluster_assignments` TSV (columns: gene, variant_id,
        cluster) from a run made without `--desired-order` /
        `--use_sequential_cluster_labels`, so `cluster` holds the raw
        K-means labels (0..K-1).
    reference_assignments_path : str
        Path to a `cluster_assignments` TSV from an already-ordered
        reference run (e.g. the main figure), whose `cluster` values define
        the desired sequence.

    Returns
    -------
    dict with keys:
        n_shared_genes : int
            Number of genes present in both files (matched by `gene`).
        n_matched_genes : int
            Number of shared genes falling in their Hungarian-matched pair.
        overlap_table : DataFrame
            Contingency table of shared gene counts, rows=raw cluster,
            columns=reference cluster.
        raw_to_reference_position : dict
            Raw cluster ID -> its matched reference cluster value.
        desired_order : list
            Raw cluster IDs ordered to reproduce the reference sequence;
            suitable for `--desired-order` on the raw run.
    """
    raw = pd.read_csv(raw_assignments_path, sep="\t")
    reference = pd.read_csv(reference_assignments_path, sep="\t")

    raw_clusters = sorted(raw["cluster"].unique())
    reference_clusters = sorted(reference["cluster"].unique())

    if len(raw_clusters) != len(reference_clusters):
        raise ValueError(
            f"Cluster count mismatch: raw has {len(raw_clusters)} clusters, "
            f"reference has {len(reference_clusters)}"
        )

    merged = raw[["gene", "cluster"]].merge(
        reference[["gene", "cluster"]], on="gene", suffixes=("_raw", "_reference")
    )

    overlap = pd.crosstab(merged["cluster_raw"], merged["cluster_reference"])
    overlap = overlap.reindex(index=raw_clusters, columns=reference_clusters, fill_value=0)

    # Hungarian algorithm: maximize total overlap == minimize negative overlap.
    row_ind, col_ind = linear_sum_assignment(-overlap.to_numpy())

    raw_to_reference_position = {
        raw_clusters[r]: reference_clusters[c] for r, c in zip(row_ind, col_ind)
    }

    desired_order = sorted(
        raw_clusters, key=lambda cluster_id: raw_to_reference_position[cluster_id]
    )

    return {
        "n_shared_genes": len(merged),
        "n_matched_genes": int(overlap.to_numpy()[row_ind, col_ind].sum()),
        "overlap_table": overlap,
        "raw_to_reference_position": raw_to_reference_position,
        "desired_order": desired_order,
    }
