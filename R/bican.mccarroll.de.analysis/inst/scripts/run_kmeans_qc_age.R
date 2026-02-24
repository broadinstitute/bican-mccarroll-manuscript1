## ============================================================
## Task: K-means QC plots for DE matrices (age)
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
de_region_interaction_dir <- base::paste0(
    de_results_dir,
    "/LEVEL_3/sex_age/cell_type_region_interaction_absolute_effects"
)

test <- "age"

regions_use_metacells <- c("CaH", "Pu", "NAC", "ic", "DFC")

## Prep matrices parameters
fdr_cutoff <- 0.01
abs_lfc_cutoff <- base::log2(1.05)
min_tpm <- 10
regions_use_de_mats <- c("CaH", "DFC")

## Region LFC matrix regions
regions_use_region_lfc <- c("CaH", "DFC")

## Heatmap scaling
scaling_factor <- 5

## -----------------------
## Execution
## -----------------------

cell_types_use <- read_cell_types(ct_file)

cell_metadata <- read_cell_metadata(cell_metadata_file)
donor_ages <- extract_donor_ages(cell_metadata)

gene_to_chr <- read_gene_to_chr(gene_to_chr_file)

## Load DE (region-combined) for gene selection
de_age <- read_de_results(de_dir, test, ct_file, gene_to_chr)

tmp <- read_metacells(
    metacells_file,
    cell_types_use = cell_types_use,
    regions_use = regions_use_metacells
)

metacells <- tmp$metacells
col_metadata <- tmp$col_metadata

metacell_summary <- summarize_metacells(
    metacells,
    col_metadata,
    donor_ages
)

de_age_mat_list <- prep_de_matrices(
    de_age,
    metacell_summary,
    cell_types_use,
    fdr_cutoff = fdr_cutoff,
    abs_lfc_cutoff = abs_lfc_cutoff,
    min_tpm = min_tpm,
    regions_use = regions_use_de_mats
)


## QC plot: silhouette vs k
plot_kmeans_silhouette(de_age_mat_list$lfc_mat_z)

## QC plot: region-combined heatmap
gene_clusters <- plot_kmeans_heatmap(
    de_age_mat_list$lfc_mat_z,
    de_age_mat_list$lfc_mat,
    scaling_factor = scaling_factor,
    k=19,
    cluster_level_order = c(2, 10, 6, 3, 5, 14, 9, 13, 4, 1, 15, 19, 18, 7, 17, 8, 11, 12)
)

# gene_clusters <- plot_kmeans_heatmap(
#     de_age_mat_list$lfc_mat_z,
#     de_age_mat_list$lfc_mat,
#     scaling_factor = scaling_factor,
#     k=19,
#     cluster_level_order = c(2, 10, 6, 3, 5, 14, 9, 13, 4, 1, 15, 19, 18, 7, 17, 8, 11, 12, 16)
# )

## Load DE (region interaction) for region-specific matrix
de_ri_age <- read_de_results(de_region_interaction_dir, test, ct_file, gene_to_chr)

de_ri_age_flc_mat <- prep_region_lfc_matrix(
    de_dt = de_ri_age,
    genes_use = rownames(de_age_mat_list$lfc_mat),
    cell_types_use = cell_types_use,
    regions_use = regions_use_region_lfc
)

## QC plot: region-specific heatmap
plot_kmeans_heatmap(
    de_age_mat_list$lfc_mat_z,
    de_ri_age_flc_mat,
    scaling_factor = scaling_factor
)
