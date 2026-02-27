## ============================================================
## Task: Plot heatmap of pairwise R-squared across cell types
## ============================================================

## -----------------------
## Parameters
## -----------------------

r_squared_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/manuscript_data/cell_type_pairwise_r_squared_qval_0.01.tsv"

output_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/manuscript_data/cell_type_cor_plot_qval_0.01.png"

## -----------------------
## Execution
## -----------------------

bican.mccarroll.eqtl::plot_cell_type_pairwise_cor(
    r_squared_path = r_squared_path,
    output_path    = output_path
)
