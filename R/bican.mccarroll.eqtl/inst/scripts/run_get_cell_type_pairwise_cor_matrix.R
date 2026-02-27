## ============================================================
## Task: Compute pairwise Spearman R-squared matrix across cell types
## ============================================================

## -----------------------
## Parameters
## -----------------------

slope_matrix_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/manuscript_data/slope_matrix_qval_0.01.tsv"

pval_nominal_matrix_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/manuscript_data/pval_nominal_matrix_qval_0.01.tsv"

pval_nominal_threshold_matrix_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/manuscript_data/pval_nominal_threshold_matrix_qval_0.01.tsv"

egene_union_pairs_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/manuscript_data/egene_union_pairs_qval_0.01.tsv"

region_cell_type_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/manuscript_data/region_cell_type.tsv"

output_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/manuscript_data/cell_type_pairwise_r_squared_qval_0.01.tsv"

## -----------------------
## Execution
## -----------------------

bican.mccarroll.eqtl::get_cell_type_pairwise_cor_matrix(
    slope_matrix_path            = slope_matrix_path,
    pval_nominal_matrix_path             = pval_nominal_matrix_path,
    pval_nominal_threshold_matrix_path   = pval_nominal_threshold_matrix_path,
    egene_union_pairs_path            = egene_union_pairs_path,
    region_cell_type_path = region_cell_type_path,
    output_path           = output_path
)
