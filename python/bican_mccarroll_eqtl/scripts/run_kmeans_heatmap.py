"""
============================================================
Task: K-means clustering heatmap of eQTL effect-size profiles
============================================================

Workflow:
  Step 1 — Choose K:
      Run k-selection to pick number of clusters.

  Step 2 — Initial heatmap (desired_order = None):
      Run K-means with desired_order=None to see clusters in numeric order.
      Inspect the heatmap and decide on a cluster ordering that places
      the diagonal color blocks in the desired sequence.

  Step 3 — Final heatmap (desired_order = [6, 8, 1, ...]):
      Re-run with the chosen desired_order to produce the final heatmap,
      cluster assignments, and downstream outputs.

  Step 4 — Gene expression ordering:
      Order genes within each cluster by expression correlation for the
      expression heatmap.

Outputs:
  1. K-selection plot (silhouette + WCSS)
  2. K-means heatmap + gene cluster counts + cluster assignments
  3. Fisher contingency table TSV (input for plot_fisher_exact.R)
  4. Ordered gene list (by expression correlation within clusters)
"""

# -----------------------
# Parameters
# -----------------------

input_path = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3/index_snp_slope_matrix_with_median_impute_qval_0.01.tsv"

output_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3"

coloc_path = f"{output_dir}/SCZ_eur_coloc_genes_pp_h4_0.9.tsv"

median_expression_path = f"{output_dir}/heatmap_index_snp_median_expression_qval_0.01.tsv"

K = 11
# Set to None for initial exploration, then fill in after inspecting heatmap.
desired_order = [6, 8, 1, 4, 2, 3, 10, 9, 7, 5, 0]
random_state = 32

k_selection_output = f"{output_dir}/kmeans_cluster_k_selection_qval_0.01.png"
heatmap_output = f"{output_dir}/kmeans_eqtl_heatmap_qval_0.01_k{K}.png"
cluster_counts_output = f"{output_dir}/gene_cluster_counts_qval_0.01_k{K}.tsv"
cluster_assignments_output = f"{output_dir}/cluster_assignments_qval_0.01_k{K}.tsv"
contingency_output = f"{output_dir}/SCZ_eur_fisher_contingency_counts_gene_clusters.tsv"
coloc_cluster_output = f"{output_dir}/SCZ_eur_coloc_genes_with_clusters.tsv"
ordered_genes_output = f"{output_dir}/ordered_genes_by_expression_correlation_k{K}.tsv"

# -----------------------
# Execution
# -----------------------

from bican_mccarroll_eqtl import (
    run_k_selection,
    plot_k_selection,
    run_kmeans_heatmap,
    build_fisher_contingency_table,
    order_genes_by_correlation,
)

# 1. K-selection plot
silhouette_df = run_k_selection(input_path, random_state=random_state)
plot_k_selection(silhouette_df, output_path=k_selection_output)

# 2. K-means heatmap + cluster counts + cluster assignments
adata, input_matrix = run_kmeans_heatmap(
    input_path=input_path,
    K=K,
    desired_order=desired_order,
    random_state=random_state,
    heatmap_output_path=heatmap_output,
    cluster_counts_output_path=cluster_counts_output,
    cluster_assignments_output_path=cluster_assignments_output,
)

# 3. Fisher contingency table
build_fisher_contingency_table(
    adata=adata,
    input_matrix=input_matrix,
    coloc_path=coloc_path,
    contingency_output_path=contingency_output,
    coloc_cluster_output_path=coloc_cluster_output,
)

# 4. Order genes by expression correlation within clusters
order_genes_by_correlation(
    median_expression_path=median_expression_path,
    cluster_assignments_path=cluster_assignments_output,
    output_path=ordered_genes_output,
)
