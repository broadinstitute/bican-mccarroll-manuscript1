## ============================================================
## Task: Write lightweight DE outputs (selected cell types)
## ============================================================

## -----------------------
## Parameters
## -----------------------

ct_file <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/sburger_tmp/cell_types_use.txt"
gene_to_chr_file <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/sburger_tmp/gene_to_chromosome.txt"

cell_metadata_file <- "/broad/bican_um1_mccarroll/RNAseq/analysis/cellarium_upload/CAP_freeze_3/CAP_cell_metadata.annotated.txt.gz"
metacells_file <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/metacells/LEVEL_3/donor_rxn_DGEList_counts.tsv.gz"

de_results_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results"
de_dir <- base::paste0(de_results_dir, "/LEVEL_3/sex_age/cell_type")
test <- "age"

regions_use_metacells <- c("CaH", "Pu", "NAC", "ic", "DFC")

de_lite_dir <- "/downloads/tmp/de_lite_new"
fdr_thresh <- 0.05
n_top <- 200

tasks <- list(
    list(cell_type_use = "MSN_D1_matrix",        region_use = "CaH", out_name = "MSN_D1_matrix"),
    list(cell_type_use = "glutamatergic_L23IT",  region_use = "DFC", out_name = "glutamatergic_L23IT"),
    list(cell_type_use = "GABA_PTHLH-PVALB",     region_use = "CaH", out_name = "GABA_PTHLH-PVALB"),
    list(cell_type_use = "astrocyte",            region_use = "CaH", out_name = "astrocyte"),
    list(cell_type_use = "OPC",                  region_use = "CaH", out_name = "OPC"),
    list(cell_type_use = "oligodendrocyte",      region_use = "CaH", out_name = "oligodendrocyte"),
    list(cell_type_use = "microglia",            region_use = "CaH", out_name = "microglia")
)

## -----------------------
## Execution
## -----------------------

cell_types_use <- read_cell_types(ct_file)
gene_to_chr <- read_gene_to_chr(gene_to_chr_file)

cell_metadata <- read_cell_metadata(cell_metadata_file)
donor_ages <- extract_donor_ages(cell_metadata)

de_age <- read_de_results(de_dir, test, ct_file, gene_to_chr)

tmp <- read_metacells(
    metacells_file,
    cell_types_use = cell_types_use,
    regions_use = regions_use_metacells
)

metacell_summary <- summarize_metacells(
    tmp$metacells,
    tmp$col_metadata,
    donor_ages
)

for (x in tasks) {
    write_de_lite(
        de_dt = de_age,
        metacell_summary = metacell_summary,
        donor_ages = donor_ages,
        cell_type_use = x$cell_type_use,
        region_use = x$region_use,
        out_name = x$out_name,
        out_dir = de_lite_dir,
        n_top = n_top,
        fdr_thresh = fdr_thresh
    )
}
