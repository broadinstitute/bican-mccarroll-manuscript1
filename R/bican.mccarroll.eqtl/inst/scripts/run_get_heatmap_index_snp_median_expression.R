## ============================================================
## Task: Get median expression per cell type for heatmap genes
## ============================================================

## -----------------------
## Parameters
## -----------------------

index_snp_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3/index_snp_slope_matrix_with_median_impute_qval_0.01.tsv"

region_cell_type_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/region_cell_type.tsv"

expression_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3/combined_gene_expression_tpm.tsv"

output_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3/heatmap_index_snp_median_expression_qval_0.01.tsv"

## -----------------------
## Execution
## -----------------------

bican.mccarroll.eqtl::get_heatmap_index_snp_median_expression(
    index_snp_path        = index_snp_path,
    region_cell_type_path = region_cell_type_path,
    expression_path       = expression_path,
    output_path           = output_path
)
