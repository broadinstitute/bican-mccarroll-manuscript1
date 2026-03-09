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
from bican_mccarroll_eqtl import run_kmeans_heatmap

# -----------------------
# Hard-coded defaults
# -----------------------

DEFAULT_CELLTYPE_ORDER = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/manuscript_data"
QVAL = "0.01"
K_DEFAULT = 11
RANDOM_STATE = 42
DESIRED_ORDER_DEFAULT = [5, 0, 6, 2, 7, 8, 10, 1, 9, 4, 3]


def _parse_desired_order(desired_order_str):
    return [int(x) for x in desired_order_str.split(",")]


def _maybe_skip(step_label, outputs, force):
    if force:
        return False
    if outputs and all(os.path.exists(p) for p in outputs):
        print(f"  SKIPPED {step_label}: outputs already exist")
        for p in outputs:
            print(f"    {p}")
        return True
    return False

def _validate_desired_order(desired_order, K):
    if len(desired_order) != K:
        raise ValueError(
            f"desired_order must have exactly {K} elements, but got {len(desired_order)}: "
            f"{desired_order}"
        )

    if len(set(desired_order)) != len(desired_order):
        raise ValueError(
            f"desired_order must not contain repeated values: {desired_order}"
        )


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="test-eqtl-pipeline",
        description="Run bican_mccarroll_eqtl Python pipeline steps (depends on R outputs).",
    )
    parser.add_argument(
        "--out-dir",
        default=DEFAULT_CELLTYPE_ORDER,
        help="Output directory used by the R pipeline (inputs) and for Python outputs.",
    )
    parser.add_argument(
        "--K",
        type=int,
        default=K_DEFAULT,
        help=f"Number of K-means clusters. Default: {K_DEFAULT}.",
    )
    parser.add_argument(
        "--desired-order",
        default=",".join(str(x) for x in DESIRED_ORDER_DEFAULT),
        help=(
            "Comma-separated cluster order to use in the heatmap. "
            f"Default: {','.join(str(x) for x in DESIRED_ORDER_DEFAULT)}."
        ),
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

    qval = QVAL
    K = args.K
    random_state = RANDOM_STATE
    desired_order = _parse_desired_order(args.desired_order)
    _validate_desired_order(desired_order, K)

    input_path = os.path.join(out_dir, f"index_snp_slope_matrix_with_zero_impute_qval_{qval}.tsv")
    heatmap_output = os.path.join(out_dir, f"kmeans_eqtl_heatmap_qval_{qval}_k{K}.svg")
    cluster_counts_output = os.path.join(out_dir, f"gene_cluster_counts_qval_{qval}_k{K}.tsv")
    cluster_assignments_output = os.path.join(out_dir, f"cluster_assignments_qval_{qval}_k{K}.tsv")

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


if __name__ == "__main__":
    main()