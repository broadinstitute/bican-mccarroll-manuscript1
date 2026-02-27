## ============================================================
## Test script: Run all bican.mccarroll.eqtl R functions in order
## ============================================================
##
## Usage:
##   1. SSH into Broad server
##   2. Start R
##   3. First install the package (only need to do this once per update):
##
##      remotes::install_github(
##        "broadinstitute/bican-mccarroll-manuscript1",
##        subdir = "R/bican.mccarroll.eqtl",
##        ref = "ty_eqtl",
##        dependencies = TRUE
##      )
##
##   4. source("test_eqtl_pipeline.R")
##
## Each step prints a status message. If a step fails, the error
## message will tell you which function broke.
## ============================================================

## -----------------------
## Paths
## -----------------------

eqtl_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/LEVEL_3"
region_cell_type_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/manuscript_data/region_cell_type.tsv"

# Output directory
out_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/manuscript_data"

qval <- 0.01

## VCF and expression paths (for plot_gene_snp)
vcf_path <- "/broad/bican_um1_mccarroll/vcfs/2025-05-05/gvs_concat_outputs_2025-05-05T14-10-02.donors_renamed_filtered_norm.vcf.gz"
combined_expression_path <- file.path(out_dir, "combined_tpm_expression_across_cell_types.tsv")

## Coloc paths (for get_sig_coloc and plot_fisher_exact)
ad_coloc_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/LEVEL_3/coloc/AD_2022"
scz_coloc_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/LEVEL_3/coloc/SCZ_eur"

## -----------------------
## Step 1: get_egene_union_pairs
## -----------------------
cat("\n===== Step 1: get_egene_union_pairs =====\n")

egene_path <- file.path(out_dir, paste0("egene_union_pairs_qval_", qval, ".tsv"))

egene_dt <- bican.mccarroll.eqtl::get_egene_union_pairs(
    eqtl_dir              = eqtl_dir,
    region_cell_type_path = region_cell_type_path,
    qval_threshold        = qval,
    output_path           = egene_path
)
cat("  Rows:", nrow(egene_dt), " Written to:", egene_path, "\n")

## -----------------------
## Step 2: get_slope_matrix
## -----------------------
cat("\n===== Step 2: get_slope_matrix =====\n")

slope_path <- file.path(out_dir, paste0("slope_matrix_qval_", qval, ".tsv"))

slope_dt <- bican.mccarroll.eqtl::get_slope_matrix(
    eqtl_dir               = eqtl_dir,
    region_cell_type_path  = region_cell_type_path,
    egene_union_pairs_path = egene_path,
    output_path            = slope_path
)
cat("  Rows:", nrow(slope_dt), " Cols:", ncol(slope_dt), " Written to:", slope_path, "\n")

## -----------------------
## Step 3: get_pval_nominal_matrix
## -----------------------
cat("\n===== Step 3: get_pval_nominal_matrix =====\n")

pval_path <- file.path(out_dir, paste0("pval_nominal_matrix_qval_", qval, ".tsv"))

pval_dt <- bican.mccarroll.eqtl::get_pval_nominal_matrix(
    eqtl_dir               = eqtl_dir,
    region_cell_type_path  = region_cell_type_path,
    egene_union_pairs_path = egene_path,
    output_path            = pval_path
)
cat("  Rows:", nrow(pval_dt), " Cols:", ncol(pval_dt), " Written to:", pval_path, "\n")

## -----------------------
## Step 4: get_pval_nominal_threshold_matrix
## -----------------------
cat("\n===== Step 4: get_pval_nominal_threshold_matrix =====\n")

pval_thresh_path <- file.path(out_dir, paste0("pval_nominal_threshold_matrix_qval_", qval, ".tsv"))

pval_thresh_dt <- bican.mccarroll.eqtl::get_pval_nominal_threshold_matrix(
    eqtl_dir               = eqtl_dir,
    region_cell_type_path  = region_cell_type_path,
    egene_union_pairs_path = egene_path,
    output_path            = pval_thresh_path
)
cat("  Rows:", nrow(pval_thresh_dt), " Cols:", ncol(pval_thresh_dt), " Written to:", pval_thresh_path, "\n")

## -----------------------
## Step 5: get_index_snp_slope_matrix_with_median_impute
## -----------------------
cat("\n===== Step 5: get_index_snp_slope_matrix_with_median_impute =====\n")

index_snp_path <- file.path(out_dir, paste0("index_snp_slope_matrix_with_median_impute_qval_", qval, ".tsv"))

index_snp_dt <- bican.mccarroll.eqtl::get_index_snp_slope_matrix_with_median_impute(
    slope_matrix_path = slope_path,
    output_path       = index_snp_path
)
cat("  Rows:", nrow(index_snp_dt), " Cols:", ncol(index_snp_dt), " Written to:", index_snp_path, "\n")

## -----------------------
## Step 6: get_cell_type_pairwise_cor_matrix
## -----------------------
cat("\n===== Step 6: get_cell_type_pairwise_cor_matrix =====\n")

r_squared_path <- file.path(out_dir, paste0("cell_type_pairwise_r_squared_qval_", qval, ".tsv"))

r_squared <- bican.mccarroll.eqtl::get_cell_type_pairwise_cor_matrix(
    slope_path            = slope_path,
    pval_path             = pval_path,
    pval_threshold_path   = pval_thresh_path,
    egene_path            = egene_path,
    region_cell_type_path = region_cell_type_path,
    output_path           = r_squared_path
)
cat("  Matrix:", nrow(r_squared), "x", ncol(r_squared), " Written to:", r_squared_path, "\n")

## -----------------------
## Step 7: plot_cell_type_pairwise_cor
## -----------------------
cat("\n===== Step 7: plot_cell_type_pairwise_cor =====\n")

cor_plot_path <- file.path(out_dir, paste0("cell_type_cor_plot_qval_", qval, ".png"))

bican.mccarroll.eqtl::plot_cell_type_pairwise_cor(
    r_squared_path = r_squared_path,
    output_path    = cor_plot_path
)
cat("  Written to:", cor_plot_path, "\n")

## -----------------------
## Step 8: get_heatmap_index_snp_median_expression
## -----------------------
cat("\n===== Step 8: get_heatmap_index_snp_median_expression =====\n")

median_expr_path <- file.path(out_dir, paste0("heatmap_index_snp_median_expression_qval_", qval, ".tsv"))

median_expr_dt <- bican.mccarroll.eqtl::get_heatmap_index_snp_median_expression(
    index_snp_path        = index_snp_path,
    region_cell_type_path = region_cell_type_path,
    expression_path       = combined_expression_path,
    output_path           = median_expr_path
)
cat("  Rows:", nrow(median_expr_dt), " Cols:", ncol(median_expr_dt), " Written to:", median_expr_path, "\n")

## -----------------------
## Step 9: combine_expression_across_cell_types
## -----------------------
cat("\n===== Step 9: combine_expression_across_cell_types =====\n")

if (!file.exists(combined_expression_path)) {
    combined_dt <- bican.mccarroll.eqtl::combine_expression_across_cell_types(
        eqtl_dir              = eqtl_dir,
        region_cell_type_path = region_cell_type_path,
        output_path           = combined_expression_path
    )
    cat("  Rows:", nrow(combined_dt), " Cols:", ncol(combined_dt), " Written to:", combined_expression_path, "\n")
} else {
    cat("  SKIPPED: file already exists at", combined_expression_path, "\n")
}

## -----------------------
## Step 10: get_sig_coloc (AD and SCZ)
## -----------------------
cat("\n===== Step 10: get_sig_coloc =====\n")

ad_coloc_path <- file.path(out_dir, "AD_2022_coloc_genes_pp_h4_0.9.tsv")
scz_coloc_path <- file.path(out_dir, "SCZ_eur_coloc_genes_pp_h4_0.9.tsv")

if (dir.exists(ad_coloc_dir)) {
    ad_coloc_dt <- bican.mccarroll.eqtl::get_sig_coloc(
        coloc_dir   = ad_coloc_dir,
        pp_h4_threshold = 0.9,
        output_path = ad_coloc_path
    )
    cat("  AD coloc genes:", nrow(ad_coloc_dt), " Written to:", ad_coloc_path, "\n")
} else {
    cat("  SKIPPED: AD coloc dir not found at", ad_coloc_dir, "\n")
}

if (dir.exists(scz_coloc_dir)) {
    scz_coloc_dt <- bican.mccarroll.eqtl::get_sig_coloc(
        coloc_dir   = scz_coloc_dir,
        pp_h4_threshold = 0.9,
        output_path = scz_coloc_path
    )
    cat("  SCZ coloc genes:", nrow(scz_coloc_dt), " Written to:", scz_coloc_path, "\n")
} else {
    cat("  SKIPPED: SCZ coloc dir not found at", scz_coloc_dir, "\n")
}

## -----------------------
## Step 11: plot_fisher_exact (AD and SCZ)
## -----------------------
## NOTE: Requires contingency tables from the Python pipeline (build_fisher_contingency_table).
## Run test_eqtl_pipeline.py steps 2-3 first, then come back to this step.
cat("\n===== Step 11: plot_fisher_exact =====\n")

cluster_order <- c(8, 1, 3, 5, 4, 10, 2, 0, 7, 6, 9)

ad_fisher_path <- file.path(out_dir, "AD_2022_fisher_contingency_counts_gene_clusters.tsv")
scz_fisher_path <- file.path(out_dir, "SCZ_eur_fisher_contingency_counts_gene_clusters.tsv")

if (file.exists(ad_fisher_path)) {
    bican.mccarroll.eqtl::plot_fisher_exact(
        fisher_table_path  = ad_fisher_path,
        plot_disease_label = "Alzheimer's disease",
        cluster_order      = cluster_order,
        output_path        = file.path(out_dir, "AD_2022_cluster_enrichment.png")
    )
    cat("  AD enrichment plot saved\n")
} else {
    cat("  SKIPPED: AD contingency table not found. Run Python pipeline first.\n")
}

if (file.exists(scz_fisher_path)) {
    bican.mccarroll.eqtl::plot_fisher_exact(
        fisher_table_path  = scz_fisher_path,
        plot_disease_label = "schizophrenia",
        cluster_order      = cluster_order,
        output_path        = file.path(out_dir, "SCZ_eur_cluster_enrichment.png")
    )
    cat("  SCZ enrichment plot saved\n")
} else {
    cat("  SKIPPED: SCZ contingency table not found. Run Python pipeline first.\n")
}

## -----------------------
## Step 12: plot_gene_snp
## -----------------------
cat("\n===== Step 12: plot_gene_snp =====\n")

gene_snp_cases <- list(
    list(gene = "XRRA1",  chr = "chr11", pos = 74935168),
    list(gene = "NPAS3",  chr = "chr14", pos = 32935820),
    list(gene = "CEP112", chr = "chr17", pos = 66192315)
)

for (case in gene_snp_cases) {
    out_file <- file.path(out_dir, paste0(case$gene, "_", case$chr, "_", case$pos, ".png"))
    bican.mccarroll.eqtl::plot_gene_snp(
        gene            = case$gene,
        chr             = case$chr,
        pos             = case$pos,
        vcf_path        = vcf_path,
        expression_path = combined_expression_path,
        output_path     = out_file
    )
    cat("  ", case$gene, "saved to:", out_file, "\n")
}

## -----------------------
## Done with R steps
## -----------------------
cat("\n===== All R steps completed successfully! =====\n")
cat("\nNext: run the Python pipeline.\n")
cat("Run the following in a terminal:\n\n")
cat("  python3 -m venv ~/test_eqtl_env\n")
cat("  source ~/test_eqtl_env/bin/activate\n")
cat('  pip install "bican_mccarroll_eqtl @ git+https://github.com/broadinstitute/bican-mccarroll-manuscript1.git@ty_eqtl#subdirectory=python/bican_mccarroll_eqtl"\n')
cat("  python3 test_eqtl_pipeline.py\n\n")
