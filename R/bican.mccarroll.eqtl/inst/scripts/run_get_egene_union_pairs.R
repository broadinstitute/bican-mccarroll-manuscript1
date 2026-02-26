## ============================================================
## Task: Get union of eGene-variant pairs across cell types
## ============================================================

## -----------------------
## Parameters
## -----------------------

eqtl_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/LEVEL_3"

region_cell_type_file <- system.file(
    "extdata", "region_cell_type_neuron_evo.tsv",
    package = "bican.mccarroll.eqtl"
)

population_suffix <- "_eur"
qval_threshold    <- 0.01
output_path       <- "/downloads/tmp/egene_union_pairs_neuron_evo_qval_0.01.tsv"

## -----------------------
## Execution
## -----------------------

bican.mccarroll.eqtl::get_egene_union_pairs(
    eqtl_dir             = eqtl_dir,
    region_cell_type_file = region_cell_type_file,
    population_suffix     = population_suffix,
    qval_threshold        = qval_threshold,
    output_path           = output_path
)
