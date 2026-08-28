"""
============================================================
Match K-means cluster order across two eQTL cluster_assignments files
============================================================

Given a reference cluster_assignments file (already ordered, e.g. produced
with --desired-order and --use_sequential_cluster_labels) and a raw
cluster_assignments file (produced without those flags, so its "cluster"
column holds the raw K-means labels 0..K-1), finds the one-to-one cluster
correspondence that maximizes shared gene membership (Hungarian assignment)
and prints the resulting --desired-order string for the raw run.

Diagnostics (shared/matched gene counts, overlap table) are written to
stderr; only the comma-separated desired-order string is written to stdout,
so this can be captured directly, e.g.:

    ORDER=$(match-cluster-order \\
        --reference-assignments main/cluster_assignments_qval_0.01_k13.tsv \\
        --raw-assignments raw/cluster_assignments_qval_0.01_k13.tsv)
============================================================
"""

import argparse
import sys

from bican_mccarroll_eqtl import match_cluster_order


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="match-cluster-order",
        description=(
            "Match raw K-means cluster IDs to a reference run's cluster "
            "ordering by shared gene membership, and print the resulting "
            "--desired-order string."
        ),
    )
    parser.add_argument(
        "--reference-assignments",
        required=True,
        help=(
            "cluster_assignments TSV from the already-ordered reference run "
            "(e.g. the main figure), with columns gene, variant_id, cluster."
        ),
    )
    parser.add_argument(
        "--raw-assignments",
        required=True,
        help=(
            "cluster_assignments TSV from the new, unordered run (produced "
            "without --desired-order / --use_sequential_cluster_labels), "
            "with columns gene, variant_id, cluster."
        ),
    )

    args = parser.parse_args(argv)

    result = match_cluster_order(
        raw_assignments_path=args.raw_assignments,
        reference_assignments_path=args.reference_assignments,
    )

    print(
        f"Shared genes: {result['n_shared_genes']}; "
        f"matched to best 1:1 pairing: {result['n_matched_genes']}",
        file=sys.stderr,
    )
    print("Overlap table (rows=raw cluster, cols=reference cluster):", file=sys.stderr)
    print(result["overlap_table"], file=sys.stderr)

    print(",".join(str(x) for x in result["desired_order"]))


if __name__ == "__main__":
    main()
