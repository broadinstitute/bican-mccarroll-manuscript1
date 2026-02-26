## ============================================================
## Task: Get union of eGene-variant pairs across cell types
## ============================================================

## -----------------------
## Parameters
## -----------------------

eqtl_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/LEVEL_3"

region_cell_type_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/region_cell_type_neuron_evo.tsv"

qval_threshold <- 0.05

output_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/script_output/LEVEL_3/egene_union_pairs_neuron_evo_qval_0.05.tsv"

## -----------------------
## Execution
## -----------------------

bican.mccarroll.eqtl::get_egene_union_pairs(
    eqtl_dir              = eqtl_dir,
    region_cell_type_path = region_cell_type_path,
    qval_threshold        = qval_threshold,
    output_path           = output_path
)
