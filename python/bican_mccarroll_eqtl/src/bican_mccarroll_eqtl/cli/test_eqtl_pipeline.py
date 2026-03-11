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
import pandas as pd
from bican_mccarroll_eqtl import run_kmeans_heatmap

# -----------------------
# Hard-coded defaults
# -----------------------

DEFAULT_OUT_DIR = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/eqtl_analysis_pipeline_run_main_figure"
QVAL = "0.01"
K_DEFAULT = 13
RANDOM_STATE = 42
DESIRED_ORDER_DEFAULT = [11,0,5,4,2,12,9,1,6,3,7,10,8]


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

    min_val = min(desired_order)
    max_val = max(desired_order)

    if min_val != 0:
        raise ValueError(
            f"desired_order must contain 0 as the minimum value, but min={min_val}: "
            f"{desired_order}"
        )

    if max_val != K - 1:
        raise ValueError(
            f"desired_order must contain {K - 1} as the maximum value, but max={max_val}: "
            f"{desired_order}"
        )


def _read_celltype_order(celltype_order_file):
    if celltype_order_file is None:
        return None

    if not os.path.exists(celltype_order_file):
        raise FileNotFoundError(
            f"celltype order file does not exist: {celltype_order_file}"
        )

    celltype_order_df = pd.read_csv(
        celltype_order_file,
        sep="\t",
        header=None,
        names=["cell_type_name"],
        dtype=str,
        keep_default_na=False
    )

    celltype_order = [
        value.strip()
        for value in celltype_order_df["cell_type_name"].tolist()
        if value.strip() != ""
    ]

    return celltype_order


def _read_celltype_label_map(celltype_label_map_file):
    if celltype_label_map_file is None:
        return None

    if not os.path.exists(celltype_label_map_file):
        raise FileNotFoundError(
            f"celltype label map file does not exist: {celltype_label_map_file}"
        )

    celltype_label_map_df = pd.read_csv(
        celltype_label_map_file,
        sep="\t",
        dtype=str,
        keep_default_na=False
    )

    expected_columns = ["cell_type_name", "pretty_label"]
    observed_columns = celltype_label_map_df.columns.tolist()

    if observed_columns != expected_columns:
        raise ValueError(
            f"celltype label map file must have header columns "
            f"{expected_columns} in that order, but got {observed_columns}"
        )

    duplicated_celltypes = celltype_label_map_df["cell_type_name"].duplicated()
    if duplicated_celltypes.any():
        duplicated_values = celltype_label_map_df.loc[
            duplicated_celltypes, "cell_type_name"
        ].tolist()
        raise ValueError(
            f"celltype label map file contains duplicated cell_type_name values: "
            f"{duplicated_values}"
        )

    celltype_label_map = dict(
        zip(
            celltype_label_map_df["cell_type_name"],
            celltype_label_map_df["pretty_label"]
        )
    )

    return celltype_label_map


def _validate_celltype_inputs(celltype_order, celltype_label_map):
    if celltype_order is None or celltype_label_map is None:
        return

    missing_labels = [
        cell_type_name
        for cell_type_name in celltype_order
        if cell_type_name not in celltype_label_map
    ]

    if len(missing_labels) > 0:
        raise ValueError(
            "celltype_label_map is missing entries for the following cell types "
            f"from celltype_order: {missing_labels}"
        )


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="test-eqtl-pipeline",
        description="Run bican_mccarroll_eqtl Python pipeline steps (depends on R outputs).",
    )
    parser.add_argument(
        "--out-dir",
        default=DEFAULT_OUT_DIR,
        help="Output directory used by the R pipeline (inputs) and for Python outputs.",
    )

    parser.add_argument(
        "--qval",
        type=float,
        default=0.01,
        help="Q-value threshold used for filtering (default: 0.01)."
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
        "--celltype-order-file",
        default=None,
        help=(
            "Tab-delimited file with no header and one cell type per line. "
            "Blank lines are ignored."
        ),
    )
    parser.add_argument(
        "--celltype-label-map-file",
        default=None,
        help=(
            "Tab-delimited 2-column file with header: "
            "cell_type_name and pretty_label."
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Recompute steps even if output files already exist.",
    )
    parser.add_argument(
        "--use_sequential_cluster_labels",
        action="store_true",
        help="Use sequential cluster labels in the output.",
    )

    args = parser.parse_args(argv)
    out_dir = args.out_dir
    force = args.force
    qval = args.qval
    K = args.K
    random_state = RANDOM_STATE


    desired_order = _parse_desired_order(args.desired_order)
    _validate_desired_order(desired_order, K)
    celltype_order = _read_celltype_order(args.celltype_order_file)
    celltype_label_map = _read_celltype_label_map(args.celltype_label_map_file)
    _validate_celltype_inputs(celltype_order, celltype_label_map)

    # -----------------------
    # Derived paths (rooted at out_dir)
    # -----------------------

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
            celltype_order=celltype_order,
            celltype_label_map=celltype_label_map,
            heatmap_output_path=heatmap_output,
            cluster_counts_output_path=cluster_counts_output,
            cluster_assignments_output_path=cluster_assignments_output,
            use_sequential_cluster_labels=args.use_sequential_cluster_labels,
        )
        print(f"  Heatmap: {heatmap_output}")
        print(f"  Cluster counts: {cluster_counts_output}")
        print(f"  Cluster assignments: {cluster_assignments_output}")
        print(f"  Genes: {adata.n_obs}, Cell types: {adata.n_vars}")

        print("\n===== K-means clustering done! =====\n")


if __name__ == "__main__":
    main()