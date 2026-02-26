"""
============================================================
Task: K-means clustering heatmap of eQTL effect-size profiles
============================================================

Outputs:
  1. K-selection plot (silhouette + WCSS)
  2. K-means heatmap + gene cluster counts TSV
  3. Fisher contingency table TSV (input for plot_fisher_exact.R)
"""

# -----------------------
# Parameters
# -----------------------

input_path = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3/index_snp_slope_matrix_with_median_impute_qval_0.01.tsv"

output_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3"

coloc_path = f"{output_dir}/SCZ_eur_coloc_genes_pp_h4_0.9.tsv"

K = 11
desired_order = [6, 8, 1, 4, 2, 3, 10, 9, 7, 5, 0]
random_state = 32

k_selection_output = f"{output_dir}/kmeans_cluster_k_selection_qval_0.01.png"
heatmap_output = f"{output_dir}/kmeans_eqtl_heatmap_qval_0.01_k{K}.png"
cluster_counts_output = f"{output_dir}/gene_cluster_counts_qval_0.01_k{K}.tsv"
contingency_output = f"{output_dir}/SCZ_eur_fisher_contingency_counts_gene_clusters.tsv"
coloc_cluster_output = f"{output_dir}/SCZ_eur_coloc_genes_with_clusters.tsv"

# -----------------------
# Execution
# -----------------------

from bican_mccarroll_eqtl import (
    run_k_selection,
    plot_k_selection,
    run_kmeans_heatmap,
    build_fisher_contingency_table,
)

# 1. K-selection plot
silhouette_df = run_k_selection(input_path, random_state=random_state)
plot_k_selection(silhouette_df, output_path=k_selection_output)

# 2. K-means heatmap + cluster counts
adata, input_matrix = run_kmeans_heatmap(
    input_path=input_path,
    K=K,
    desired_order=desired_order,
    random_state=random_state,
    heatmap_output_path=heatmap_output,
    cluster_counts_output_path=cluster_counts_output,
)

# 3. Fisher contingency table
build_fisher_contingency_table(
    adata=adata,
    input_matrix=input_matrix,
    coloc_path=coloc_path,
    contingency_output_path=contingency_output,
    coloc_cluster_output_path=coloc_cluster_output,
)
