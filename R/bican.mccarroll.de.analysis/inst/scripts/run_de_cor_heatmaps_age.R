## ============================================================
## Task: DE correlation matrices + heatmaps (age)
## ============================================================

## -----------------------
## Parameters
## -----------------------

ct_file <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/sburger_tmp/cell_types_use.txt"
gene_to_chr_file <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/sburger_tmp/gene_to_chromosome.txt"

de_results_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results"
de_region_interaction_dir <- base::paste0(de_results_dir, "/LEVEL_3/sex_age/cell_type_region_interaction_absolute_effects")

test <- "age"

non_neuron_types <- c("astrocyte", "OPC", "oligodendrocyte", "microglia")

region_order <- c("CaH", "Pu", "NAC", "ic", "DFC")
regions_main <- c("CaH", "DFC")
regions_supp <- region_order

fdr_cutoff <- 0.05

breaks <- base::seq(-1, 1, length.out = 101)
palette_colors <- c("steelblue", "white", "darkorange")
clustering_method <- "complete"

## -----------------------
## Execution
## -----------------------

cell_types_use <- read_cell_types(ct_file)
gene_to_chr <- read_gene_to_chr(gene_to_chr_file)
de_ri_age <- read_de_results(de_region_interaction_dir, test, ct_file, gene_to_chr)

cor_mat_main <- compute_de_cor_mat(de_ri_age, cell_types_use, regions_main, non_neuron_types, fdr_cutoff = fdr_cutoff)
cor_mat_supp <- compute_de_cor_mat(de_ri_age, cell_types_use, regions_supp, non_neuron_types, fdr_cutoff = fdr_cutoff)

plot_de_cor_heatmap(
    cor_mat_main,
    clustering_method = clustering_method,
    breaks = breaks,
    palette_colors = palette_colors
)

plot_de_cor_heatmap(
    cor_mat_supp,
    clustering_method = clustering_method,
    breaks = breaks,
    palette_colors = palette_colors
)
