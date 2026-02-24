## ============================================================
## Task: DE volcano plots (age, region-combined)
## ============================================================

## -----------------------
## Parameters
## -----------------------

ct_file <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/sburger_tmp/cell_types_use.txt"
gene_to_chr_file <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/sburger_tmp/gene_to_chromosome.txt"

de_results_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results"
de_dir <- base::paste0(de_results_dir, "/LEVEL_3/sex_age/cell_type")

test <- "age"
region_use <- NA

fdr_cutoff <- 0.05
abs_log_fc_cutoff <- base::log2(1.05)

## -----------------------
## Execution
## -----------------------

cell_types_use <- read_cell_types(ct_file)
gene_to_chr <- read_gene_to_chr(gene_to_chr_file)

de_age <- read_de_results(de_dir, test, ct_file, gene_to_chr)

for (ct in cell_types_use) {
    plot_de_volcano(
        de_age,
        cell_type_use = ct,
        region_use = region_use,
        fdr_cutoff = fdr_cutoff,
        abs_log_fc_cutoff = abs_log_fc_cutoff
    )
}
