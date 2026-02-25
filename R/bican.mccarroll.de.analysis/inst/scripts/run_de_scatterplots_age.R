## ============================================================
## Task: DE scatterplots (age)
## ============================================================

## -----------------------
## Parameters
## -----------------------

ct_file <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/sburger_tmp/cell_types_use.txt"
gene_to_chr_file <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/sburger_tmp/gene_to_chromosome.txt"

de_results_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results"
de_region_interaction_dir <- base::paste0(de_results_dir, "/LEVEL_3/sex_age/cell_type_region_interaction_absolute_effects")

test <- "age"

fdr_cutoff <- 0.05
add_fit <- TRUE

out_pdf <- "/downloads/tmp/age_DE_scatterplots.pdf"

## -----------------------
## Execution
## -----------------------

gene_to_chr <- bican.mccarroll.de.analysis::read_gene_to_chr(gene_to_chr_file)
de_ri_age <- bican.mccarroll.de.analysis::read_de_results(de_region_interaction_dir, test, ct_file, gene_to_chr)

grDevices::pdf(out_pdf)

## Within-cell-type, across-region comparisons
bican.mccarroll.de.analysis::plot_de_scatter(de_ri_age, "MSN_D1_matrix", "MSN_D1_matrix", "CaH", "Pu", fdr_cutoff = fdr_cutoff, add_fit = add_fit)
bican.mccarroll.de.analysis::plot_de_scatter(de_ri_age, "MSN_D1_matrix", "MSN_D1_matrix", "CaH", "Pu", fdr_cutoff = fdr_cutoff, add_fit = add_fit, xlab_prefix="DE Age ")

bican.mccarroll.de.analysis::plot_de_scatter(de_ri_age, "MSN_D1_matrix", "MSN_D1_matrix", "CaH", "NAC", fdr_cutoff = fdr_cutoff, add_fit = add_fit)
bican.mccarroll.de.analysis::plot_de_scatter(de_ri_age, "MSN_D1_matrix", "MSN_D1_matrix", "Pu", "NAC", fdr_cutoff = fdr_cutoff, add_fit = add_fit)

bican.mccarroll.de.analysis::plot_de_scatter(de_ri_age, "MSN_D2_matrix", "MSN_D2_matrix", "CaH", "Pu", fdr_cutoff = fdr_cutoff, add_fit = add_fit)
bican.mccarroll.de.analysis::plot_de_scatter(de_ri_age, "MSN_D2_matrix", "MSN_D2_matrix", "CaH", "NAC", fdr_cutoff = fdr_cutoff, add_fit = add_fit)
bican.mccarroll.de.analysis::plot_de_scatter(de_ri_age, "MSN_D2_matrix", "MSN_D2_matrix", "Pu", "NAC", fdr_cutoff = fdr_cutoff, add_fit = add_fit)

bican.mccarroll.de.analysis::plot_de_scatter(de_ri_age, "MSN_D1_striosome", "MSN_D1_striosome", "CaH", "Pu", fdr_cutoff = fdr_cutoff, add_fit = add_fit)
bican.mccarroll.de.analysis::plot_de_scatter(de_ri_age, "MSN_D1_striosome", "MSN_D1_striosome", "CaH", "NAC", fdr_cutoff = fdr_cutoff, add_fit = add_fit)
bican.mccarroll.de.analysis::plot_de_scatter(de_ri_age, "MSN_D1_striosome", "MSN_D1_striosome", "Pu", "NAC", fdr_cutoff = fdr_cutoff, add_fit = add_fit)

bican.mccarroll.de.analysis::plot_de_scatter(de_ri_age, "MSN_D2_striosome", "MSN_D2_striosome", "CaH", "Pu", fdr_cutoff = fdr_cutoff, add_fit = add_fit)
bican.mccarroll.de.analysis::plot_de_scatter(de_ri_age, "MSN_D2_striosome", "MSN_D2_striosome", "CaH", "NAC", fdr_cutoff = fdr_cutoff, add_fit = add_fit)
bican.mccarroll.de.analysis::plot_de_scatter(de_ri_age, "MSN_D2_striosome", "MSN_D2_striosome", "Pu", "NAC", fdr_cutoff = fdr_cutoff, add_fit = add_fit)

bican.mccarroll.de.analysis::plot_de_scatter(de_ri_age, "astrocyte", "astrocyte", "CaH", "Pu", fdr_cutoff = fdr_cutoff, add_fit = add_fit)
bican.mccarroll.de.analysis::plot_de_scatter(de_ri_age, "astrocyte", "astrocyte", "CaH", "NAC", fdr_cutoff = fdr_cutoff, add_fit = add_fit)
bican.mccarroll.de.analysis::plot_de_scatter(de_ri_age, "astrocyte", "astrocyte", "Pu", "NAC", fdr_cutoff = fdr_cutoff, add_fit = add_fit)

bican.mccarroll.de.analysis::plot_de_scatter(de_ri_age, "OPC", "OPC", "CaH", "Pu", fdr_cutoff = fdr_cutoff, add_fit = add_fit)
bican.mccarroll.de.analysis::plot_de_scatter(de_ri_age, "OPC", "OPC", "CaH", "NAC", fdr_cutoff = fdr_cutoff, add_fit = add_fit)
bican.mccarroll.de.analysis::plot_de_scatter(de_ri_age, "OPC", "OPC", "Pu", "NAC", fdr_cutoff = fdr_cutoff, add_fit = add_fit)

bican.mccarroll.de.analysis::plot_de_scatter(de_ri_age, "oligodendrocyte", "oligodendrocyte", "CaH", "Pu", fdr_cutoff = fdr_cutoff, add_fit = add_fit)
bican.mccarroll.de.analysis::plot_de_scatter(de_ri_age, "oligodendrocyte", "oligodendrocyte", "CaH", "NAC", fdr_cutoff = fdr_cutoff, add_fit = add_fit)
bican.mccarroll.de.analysis::plot_de_scatter(de_ri_age, "oligodendrocyte", "oligodendrocyte", "Pu", "NAC", fdr_cutoff = fdr_cutoff, add_fit = add_fit)

bican.mccarroll.de.analysis::plot_de_scatter(de_ri_age, "microglia", "microglia", "CaH", "Pu", fdr_cutoff = fdr_cutoff, add_fit = add_fit)
bican.mccarroll.de.analysis::plot_de_scatter(de_ri_age, "microglia", "microglia", "CaH", "NAC", fdr_cutoff = fdr_cutoff, add_fit = add_fit)
bican.mccarroll.de.analysis::plot_de_scatter(de_ri_age, "microglia", "microglia", "Pu", "NAC", fdr_cutoff = fdr_cutoff, add_fit = add_fit)

bican.mccarroll.de.analysis::plot_de_scatter(de_ri_age, "astrocyte", "astrocyte", "CaH", "DFC", fdr_cutoff = fdr_cutoff, add_fit = add_fit)
bican.mccarroll.de.analysis::plot_de_scatter(de_ri_age, "OPC", "OPC", "CaH", "DFC", fdr_cutoff = fdr_cutoff, add_fit = add_fit)
bican.mccarroll.de.analysis::plot_de_scatter(de_ri_age, "oligodendrocyte", "oligodendrocyte", "CaH", "DFC", fdr_cutoff = fdr_cutoff, add_fit = add_fit)
bican.mccarroll.de.analysis::plot_de_scatter(de_ri_age, "microglia", "microglia", "CaH", "DFC", fdr_cutoff = fdr_cutoff, add_fit = add_fit)

grDevices::dev.off()
