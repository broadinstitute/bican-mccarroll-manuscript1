"""
============================================================
Test script: Run all bican_mccarroll_eqtl Python functions
============================================================

Usage:
  1. SSH into Broad server
  2. Set up environment and install package:

     python3 -m venv ~/test_eqtl_env
     source ~/test_eqtl_env/bin/activate
     pip install "bican_mccarroll_eqtl @ git+https://github.com/broadinstitute/bican-mccarroll-manuscript1.git@ty_eqtl#subdirectory=python/bican_mccarroll_eqtl"

  3. Make sure you've already run test_eqtl_pipeline.R first
     (the Python steps depend on R outputs from steps 5 and 8)

  4. python3 test_eqtl_pipeline.py
============================================================
"""

import os

# -----------------------
# Paths
# -----------------------

out_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/manuscript_data"

qval = "0.01"

input_path = os.path.join(out_dir, f"index_snp_slope_matrix_with_median_impute_qval_{qval}.tsv")
median_expression_path = os.path.join(out_dir, f"heatmap_index_snp_median_expression_qval_{qval}.tsv")

ad_coloc_path = os.path.join(out_dir, "AD_2022_coloc_genes_pp_h4_0.9.tsv")
scz_coloc_path = os.path.join(out_dir, "SCZ_eur_coloc_genes_pp_h4_0.9.tsv")

K = 9
random_state = 119
desired_order = [2, 1, 7, 4, 8, 5, 6, 3, 0]

k_selection_output = os.path.join(out_dir, f"kmeans_cluster_k_selection_qval_{qval}.png")
heatmap_output = os.path.join(out_dir, f"kmeans_eqtl_heatmap_qval_{qval}_k{K}.png")
cluster_counts_output = os.path.join(out_dir, f"gene_cluster_counts_qval_{qval}_k{K}.tsv")
cluster_assignments_output = os.path.join(out_dir, f"cluster_assignments_qval_{qval}_k{K}.tsv")
ordered_genes_output = os.path.join(out_dir, f"ordered_genes_by_expression_correlation_k{K}.tsv")
expression_heatmap_output = os.path.join(out_dir, "median_expression_heatmap_kmeans_clusters_expr_ordered.png")

# -----------------------
# Execution
# -----------------------

from bican_mccarroll_eqtl import (
    run_k_selection,
    plot_k_selection,
    run_kmeans_heatmap,
    build_fisher_contingency_table,
    order_genes_by_correlation,
    plot_expression_heatmap,
)

# Step 1: K-selection
print("\n===== Step 1: run_k_selection =====")
silhouette_df = run_k_selection(input_path, random_state=random_state)
plot_k_selection(silhouette_df, output_path=k_selection_output)
print(f"  K-selection plot saved to: {k_selection_output}")

# Step 2: K-means heatmap
print("\n===== Step 2: run_kmeans_heatmap =====")
adata, input_matrix = run_kmeans_heatmap(
    input_path=input_path,
    K=K,
    desired_order=desired_order,
    random_state=random_state,
    heatmap_output_path=heatmap_output,
    cluster_counts_output_path=cluster_counts_output,
    cluster_assignments_output_path=cluster_assignments_output,
)
print(f"  Heatmap saved to: {heatmap_output}")
print(f"  Cluster counts saved to: {cluster_counts_output}")
print(f"  Cluster assignments saved to: {cluster_assignments_output}")
print(f"  Genes: {adata.n_obs}, Cell types: {adata.n_vars}")

# Step 3: Fisher contingency tables (AD and SCZ)
print("\n===== Step 3: build_fisher_contingency_table =====")

for label, coloc_path in [("AD_2022", ad_coloc_path), ("SCZ_eur", scz_coloc_path)]:
    if os.path.exists(coloc_path):
        contingency_output = os.path.join(out_dir, f"{label}_fisher_contingency_counts_gene_clusters.tsv")
        coloc_cluster_output = os.path.join(out_dir, f"{label}_coloc_genes_with_clusters.tsv")
        contingency_df = build_fisher_contingency_table(
            adata=adata,
            input_matrix=input_matrix,
            coloc_path=coloc_path,
            contingency_output_path=contingency_output,
            coloc_cluster_output_path=coloc_cluster_output,
        )
        print(f"  {label} contingency table saved to: {contingency_output}")
        print(f"  {label} coloc cluster mapping saved to: {coloc_cluster_output}")
    else:
        print(f"  SKIPPED {label}: coloc file not found at {coloc_path}")

# Step 4: Order genes by expression correlation
print("\n===== Step 4: order_genes_by_correlation =====")
if os.path.exists(median_expression_path):
    ordered_genes = order_genes_by_correlation(
        median_expression_path=median_expression_path,
        cluster_assignments_path=cluster_assignments_output,
        output_path=ordered_genes_output,
    )
    print(f"  Ordered {len(ordered_genes)} genes")
    print(f"  Saved to: {ordered_genes_output}")
else:
    print(f"  SKIPPED: median expression file not found at {median_expression_path}")
    print("  (Run test_eqtl_pipeline.R step 8 first)")

# Step 5: Plot expression heatmap
print("\n===== Step 5: plot_expression_heatmap =====")
if os.path.exists(median_expression_path) and os.path.exists(ordered_genes_output):
    fig = plot_expression_heatmap(
        median_expression_path=median_expression_path,
        ordered_genes_path=ordered_genes_output,
        cluster_assignments_path=cluster_assignments_output,
        cluster_order=desired_order,
        output_path=expression_heatmap_output,
    )
    print(f"  Expression heatmap saved to: {expression_heatmap_output}")
else:
    print("  SKIPPED: requires median expression and ordered genes files from previous steps")

print("\n===== All Python steps completed successfully! =====")
print("\nNow go back to test_eqtl_pipeline.R step 11 to run plot_fisher_exact")
print("(it depends on the contingency tables generated in step 3 above).\n")
