## ============================================================
## Task: Build slope matrix for eGene-variant pairs
## ============================================================

## -----------------------
## Parameters
## -----------------------

eqtl_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/LEVEL_3"

region_cell_type_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/manuscript_data/region_cell_type.tsv"

egene_union_pairs_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/manuscript_data/egene_union_pairs_qval_0.01.tsv"

output_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/manuscript_data/slope_matrix_qval_0.01.tsv"

## -----------------------
## Execution
## -----------------------

bican.mccarroll.eqtl::get_slope_matrix(
    eqtl_dir              = eqtl_dir,
    region_cell_type_path = region_cell_type_path,
    egene_union_pairs_path = egene_union_pairs_path,
    output_path           = output_path
)
