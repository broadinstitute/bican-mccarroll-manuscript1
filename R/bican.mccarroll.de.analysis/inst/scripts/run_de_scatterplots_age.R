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

## -----------------------
## Execution
## -----------------------

gene_to_chr <- bican.mccarroll.de.analysis::read_gene_to_chr(gene_to_chr_file)
de_ri_age <- bican.mccarroll.de.analysis::read_de_results(de_region_interaction_dir, test, ct_file, gene_to_chr)

## Within-cell-type, across-region comparisons
bican.mccarroll.de.analysis::plot_de_scatter(de_ri_age, "MSN_D1_matrix", "MSN_D1_matrix", "CaH", "Pu", fdr_cutoff = fdr_cutoff, add_fit = add_fit)
bican.mccarroll.de.analysis::plot_de_scatter(de_ri_age, "MSN_D1_matrix", "MSN_D1_matrix", "CaH", "Pu", fdr_cutoff = fdr_cutoff, add_fit = add_fit, xlab_prefix="DE Age ")

p1<-bican.mccarroll.de.analysis::plot_de_scatter_gg(de_ri_age, "MSN_D1_matrix", "MSN_D1_matrix", "CaH", "Pu", fdr_cutoff = fdr_cutoff, add_fit = add_fit, xlab_prefix="DE Age ")

print (p1)
p1 <- p1 + ggplot2::theme_classic() +
    ggplot2::theme(
        axis.title.x = ggplot2::element_text(size = ggplot2::rel(1.5)),
        axis.title.y = ggplot2::element_text(size = ggplot2::rel(1.5)),

    )

print (p1)
