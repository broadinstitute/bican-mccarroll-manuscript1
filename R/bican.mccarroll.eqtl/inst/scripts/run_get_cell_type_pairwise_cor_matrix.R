## ============================================================
## Task: Compute pairwise Spearman R-squared matrix across cell types
## ============================================================

## -----------------------
## Parameters
## -----------------------

slope_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3/slope_matrix_qval_0.01.tsv"

pval_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3/pval_nominal_matrix_qval_0.01.tsv"

pval_threshold_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3/pval_nominal_threshold_matrix_qval_0.01.tsv"

egene_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3/egene_union_pairs_qval_0.01.tsv"

region_cell_type_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/region_cell_type.tsv"

output_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3/cell_type_pairwise_r_squared.tsv"

## -----------------------
## Execution
## -----------------------

bican.mccarroll.eqtl::get_cell_type_pairwise_cor_matrix(
    slope_path            = slope_path,
    pval_path             = pval_path,
    pval_threshold_path   = pval_threshold_path,
    egene_path            = egene_path,
    region_cell_type_path = region_cell_type_path,
    output_path           = output_path
)
