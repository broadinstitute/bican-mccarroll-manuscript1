## ============================================================
## Task: Select index SNPs and build median-imputed slope matrix
## ============================================================

## -----------------------
## Parameters
## -----------------------

slope_matrix_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/manuscript_data/slope_matrix_qval_0.01.tsv"

output_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/manuscript_data/index_snp_slope_matrix_with_median_impute_qval_0.01.tsv"

## -----------------------
## Execution
## -----------------------

bican.mccarroll.eqtl::get_index_snp_slope_matrix_with_median_impute(
    slope_matrix_path = slope_matrix_path,
    output_path       = output_path
)
