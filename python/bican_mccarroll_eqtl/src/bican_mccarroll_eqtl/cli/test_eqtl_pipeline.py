"""
============================================================
Test pipeline: Run bican_mccarroll_eqtl Python steps
============================================================

Defaults match the original hard-coded script.

Behavior:
  - By default, steps are skipped if their expected output file(s) exist.
  - Use --force to recompute and overwrite outputs.
============================================================
"""

import argparse
import os
import sys

from bican_mccarroll_eqtl import (
    run_k_selection,
    plot_k_selection,
    run_kmeans_heatmap,
    build_fisher_contingency_table,
    order_genes_by_correlation,
    plot_expression_heatmap
)


# -----------------------
# Hard-coded defaults (unchanged)
# -----------------------

OUT_DIR = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/manuscript_data"
QVAL = "0.01"

K = 11
RANDOM_STATE = 51
DESIRED_ORDER = [8, 1, 3, 5, 4, 10, 2, 0, 7, 6, 9]


def _maybe_skip(step_label, outputs, force):
    if force:
        return False
    if outputs and all(os.path.exists(p) for p in outputs):
        print(f"  SKIPPED {step_label}: outputs already exist")
        for p in outputs:
            print(f"    {p}")
        return True
    return False


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="test-eqtl-pipeline",
        description="Run bican_mccarroll_eqtl Python pipeline steps (depends on R outputs).",
    )
    parser.add_argument(
        "--out-dir",
        default=OUT_DIR,
        help="Output directory used by the R pipeline (inputs) and for Python outputs.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Recompute steps even if output files already exist.",
    )

    args = parser.parse_args(argv)
    out_dir = args.out_dir
    force = args.force

    # -----------------------
    # Derived paths (unchanged logic, but rooted at out_dir)
    # -----------------------

    input_path = os.path.join(
        out_dir, f"index_snp_slope_matrix_with_median_impute_qval_{QVAL}.tsv"
    )
    median_expression_path = os.path.join(
        out_dir, f"heatmap_index_snp_median_expression_qval_{QVAL}.tsv"
    )

    ad_coloc_path = os.path.join(out_dir, "AD_2022_coloc_genes_pp_h4_0.9.tsv")
    scz_coloc_path = os.path.join(out_dir, "SCZ_eur_coloc_genes_pp_h4_0.9.tsv")

    k_selection_output = os.path.join(
        out_dir, f"kmeans_cluster_k_selection_qval_{QVAL}.png"
    )
    heatmap_output = os.path.join(
        out_dir, f"kmeans_eqtl_heatmap_qval_{QVAL}_k{K}.png"
    )
    cluster_counts_output = os.path.join(
        out_dir, f"gene_cluster_counts_qval_{QVAL}_k{K}.tsv"
    )
    cluster_assignments_output = os.path.join(
        out_dir, f"cluster_assignments_qval_{QVAL}_k{K}.tsv"
    )
    ordered_genes_output = os.path.join(
        out_dir, f"ordered_genes_by_expression_correlation_k{K}.tsv"
    )
    expression_heatmap_output = os.path.join(
        out_dir, "median_expression_heatmap_kmeans_clusters_expr_ordered.png"
    )

    # -----------------------
    # Step 1: K-selection
    # -----------------------
    print("\n===== Step 1: run_k_selection =====")

    if not os.path.exists(input_path):
        print(f"  ERROR: input file not found: {input_path}")
        sys.exit(2)

    if not _maybe_skip("Step 1", [k_selection_output], force):
        silhouette_df = run_k_selection(input_path, random_state=RANDOM_STATE)
        plot_k_selection(silhouette_df, output_path=k_selection_output)
        print(f"  K-selection plot saved to: {k_selection_output}")

    # -----------------------
    # Step 2: K-means heatmap
    # -----------------------
    print("\n===== Step 2: run_kmeans_heatmap =====")

    step2_outputs = [
        heatmap_output,
        cluster_counts_output,
        cluster_assignments_output,
    ]

    if not _maybe_skip("Step 2", step2_outputs, force):
        adata, input_matrix = run_kmeans_heatmap(
            input_path=input_path,
            K=K,
            desired_order=DESIRED_ORDER,
            random_state=RANDOM_STATE,
            heatmap_output_path=heatmap_output,
            cluster_counts_output_path=cluster_counts_output,
            cluster_assignments_output_path=cluster_assignments_output,
        )
        print(f"  Heatmap saved to: {heatmap_output}")
        print(f"  Cluster counts saved to: {cluster_counts_output}")
        print(f"  Cluster assignments saved to: {cluster_assignments_output}")
        print(f"  Genes: {adata.n_obs}, Cell types: {adata.n_vars}")
    else:
        adata = None
        input_matrix = None

    # -----------------------
    # Step 3: Fisher contingency tables
    # -----------------------
    print("\n===== Step 3: build_fisher_contingency_table =====")

    if adata is None or input_matrix is None:
        print("  SKIPPED Step 3: requires Step 2 to run in this invocation.")
    else:
        for label, coloc_path in [
            ("AD_2022", ad_coloc_path),
            ("SCZ_eur", scz_coloc_path),
        ]:
            if not os.path.exists(coloc_path):
                print(f"  SKIPPED {label}: coloc file not found at {coloc_path}")
                continue

            contingency_output = os.path.join(
                out_dir,
                f"{label}_fisher_contingency_counts_gene_clusters.tsv",
            )
            coloc_cluster_output = os.path.join(
                out_dir,
                f"{label}_coloc_genes_with_clusters.tsv",
            )

            if _maybe_skip(
                f"Step 3 ({label})",
                [contingency_output, coloc_cluster_output],
                force,
            ):
                continue

            build_fisher_contingency_table(
                adata=adata,
                input_matrix=input_matrix,
                coloc_path=coloc_path,
                contingency_output_path=contingency_output,
                coloc_cluster_output_path=coloc_cluster_output,
            )
            print(f"  {label} contingency table saved to: {contingency_output}")
            print(f"  {label} coloc cluster mapping saved to: {coloc_cluster_output}")

    # -----------------------
    # Step 4: Order genes
    # -----------------------
    print("\n===== Step 4: order_genes_by_correlation =====")

    if not os.path.exists(median_expression_path):
        print(f"  SKIPPED: median expression file not found at {median_expression_path}")
    elif not os.path.exists(cluster_assignments_output):
        print(f"  SKIPPED: cluster assignments file not found at {cluster_assignments_output}")
    elif not _maybe_skip("Step 4", [ordered_genes_output], force):
        ordered_genes = order_genes_by_correlation(
            median_expression_path=median_expression_path,
            cluster_assignments_path=cluster_assignments_output,
            output_path=ordered_genes_output,
        )
        print(f"  Ordered {len(ordered_genes)} genes")
        print(f"  Saved to: {ordered_genes_output}")

    # -----------------------
    # Step 5: Expression heatmap
    # -----------------------
    print("\n===== Step 5: plot_expression_heatmap =====")

    prereq_ok = (
        os.path.exists(median_expression_path)
        and os.path.exists(ordered_genes_output)
        and os.path.exists(cluster_assignments_output)
    )

    if not prereq_ok:
        print("  SKIPPED: requires median expression, ordered genes, and cluster assignments.")
    elif not _maybe_skip("Step 5", [expression_heatmap_output], force):
        plot_expression_heatmap(
            median_expression_path=median_expression_path,
            ordered_genes_path=ordered_genes_output,
            cluster_assignments_path=cluster_assignments_output,
            cluster_order=DESIRED_ORDER,
            output_path=expression_heatmap_output,
        )
        print(f"  Expression heatmap saved to: {expression_heatmap_output}")

    print("\n===== All Python steps completed successfully! =====")
    print("\nNow go back to the R pipeline step that runs plot_fisher_exact.\n")


if __name__ == "__main__":
    main()