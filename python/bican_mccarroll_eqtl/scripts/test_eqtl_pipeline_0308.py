"""
============================================================
Test pipeline 0308: K-means clustering

Called automatically by test_eqtl_pipeline_0308.R via system2().
============================================================
"""

import os

# -----------------------
# Paths
# -----------------------

base_dir = "/Users/tracyyuan/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls"
out_dir = os.path.join(base_dir, "manuscript_test_0308")

qval = "0.01"
K = 11
random_state = 42
desired_order = [5, 0, 6, 2, 7, 8, 10, 1, 9, 4, 3]

input_path = os.path.join(out_dir, f"index_snp_slope_matrix_with_zero_impute_qval_{qval}.tsv")
heatmap_output = os.path.join(out_dir, f"kmeans_eqtl_heatmap_qval_{qval}_k{K}.svg")
cluster_counts_output = os.path.join(out_dir, f"gene_cluster_counts_qval_{qval}_k{K}.tsv")
cluster_assignments_output = os.path.join(out_dir, f"cluster_assignments_qval_{qval}_k{K}.tsv")

# -----------------------
# Execution
# -----------------------

from bican_mccarroll_eqtl import run_kmeans_heatmap

print("\n===== K-means clustering =====")
adata, input_matrix = run_kmeans_heatmap(
    input_path=input_path,
    K=K,
    desired_order=desired_order,
    random_state=random_state,
    heatmap_output_path=heatmap_output,
    cluster_counts_output_path=cluster_counts_output,
    cluster_assignments_output_path=cluster_assignments_output,
)
print(f"  Heatmap: {heatmap_output}")
print(f"  Cluster counts: {cluster_counts_output}")
print(f"  Cluster assignments: {cluster_assignments_output}")
print(f"  Genes: {adata.n_obs}, Cell types: {adata.n_vars}")

print("\n===== K-means clustering done! =====\n")
