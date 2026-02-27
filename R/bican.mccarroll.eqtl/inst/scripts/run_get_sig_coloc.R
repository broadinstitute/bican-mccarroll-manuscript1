## ============================================================
## Task: Get significant colocalization genes
## ============================================================

## -----------------------
## Parameters
## -----------------------

coloc_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/coloc/results/LEVEL_3_EUR"

pp_h4_threshold <- 0.9

region_cell_type_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/manuscript_data/region_cell_type.tsv"

output_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/manuscript_data"

## -----------------------
## Execution
## -----------------------

# 1. SCZ_eur
bican.mccarroll.eqtl::get_sig_coloc(
    coloc_dir             = coloc_dir,
    gwas_trait            = "SCZ_eur",
    pp_h4_threshold       = pp_h4_threshold,
    region_cell_type_path = region_cell_type_path,
    output_path           = file.path(output_dir, "SCZ_eur_coloc_genes_pp_h4_0.9.tsv")
)

# 2. AD_2022
bican.mccarroll.eqtl::get_sig_coloc(
    coloc_dir             = coloc_dir,
    gwas_trait            = "AD_2022",
    pp_h4_threshold       = pp_h4_threshold,
    region_cell_type_path = region_cell_type_path,
    output_path           = file.path(output_dir, "AD_2022_coloc_genes_pp_h4_0.9.tsv")
)
