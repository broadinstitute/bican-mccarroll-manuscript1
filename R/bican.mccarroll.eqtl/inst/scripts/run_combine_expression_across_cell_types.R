## ============================================================
## Task: Combine gene expression TPM across cell types
## ============================================================

## -----------------------
## Parameters
## -----------------------

eqtl_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/LEVEL_3"

region_cell_type_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/manuscript_data/region_cell_type.tsv"

output_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/manuscript_data/combined_tpm_expression_across_cell_types.tsv"

## -----------------------
## Execution
## -----------------------

bican.mccarroll.eqtl::combine_expression_across_cell_types(
    eqtl_dir              = eqtl_dir,
    region_cell_type_path = region_cell_type_path,
    output_path           = output_path
)
