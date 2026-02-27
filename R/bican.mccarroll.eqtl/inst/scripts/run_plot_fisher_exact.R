## ============================================================
## Task: Plot Fisher's exact test enrichment by cluster
## ============================================================

## -----------------------
## Parameters
## -----------------------

output_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/manuscript_data"

cluster_order <- c("6", "8", "1", "4", "2", "3", "10", "9", "7", "5", "0")

## -----------------------
## Execution
## -----------------------

# 1. SCZ_eur
bican.mccarroll.eqtl::plot_fisher_exact(
    fisher_table_path  = file.path(output_dir, "SCZ_eur_fisher_contingency_counts_gene_clusters.tsv"),
    plot_disease_label = "schizophrenia",
    cluster_order      = cluster_order,
    output_path        = file.path(output_dir, "SCZ_eur_cluster_enrichment.png")
)

# 2. AD_2022
bican.mccarroll.eqtl::plot_fisher_exact(
    fisher_table_path  = file.path(output_dir, "AD_2022_fisher_contingency_counts_gene_clusters.tsv"),
    plot_disease_label = "Alzheimer's disease",
    cluster_order      = cluster_order,
    output_path        = file.path(output_dir, "AD_2022_cluster_enrichment.png")
)
