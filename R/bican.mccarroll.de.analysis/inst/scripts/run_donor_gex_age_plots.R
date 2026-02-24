## ============================================================
## Task: Donor-level GEX plots vs age
## ============================================================

## -----------------------
## Parameters
## -----------------------

ct_file <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/sburger_tmp/cell_types_use.txt"
cell_metadata_file <- "/broad/bican_um1_mccarroll/RNAseq/analysis/cellarium_upload/CAP_freeze_3/CAP_cell_metadata.annotated.txt.gz"
metacells_file <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/metacells/LEVEL_3/donor_rxn_DGEList_counts.tsv.gz"

regions_use_metacells <- c("CaH", "Pu", "NAC", "ic", "DFC")

## Heatmap gene set
gs <- c(
    "HOMER1", "RIMS1", "NLGN4X", "IL1RAPL1", "SLITRK2", "LRRN1", "CSMD2",
    "PPP1R1B", "PPP1R1A", "ADCY1", "ADCY3", "PDE10A", "PDE4B", "PDE7A", "PRKCB",
    "KCNQ3", "KCNAB1", "RELN", "DAB1",
    "PTGDR", "PTGER2", "IFNLR1", "TLR5", "SERPING1", "IL17RB",
    "FKBP5", "BLVRA",
    "ATG9B", "SLC17A5",
    "TGFBR3", "ACVR2B", "NOTCH2", "CHRD",
    "WDR19", "BBS1", "KIZ", "TTC29", "DNAI7"
)

gs_gaps <- c(7, 15, 19, 25, 27, 29, 33)

## Heatmap selection
heatmap_key <- "MSN_D1_matrix__CaH"
heatmap_cluster_gs <- FALSE
heatmap_transpose <- FALSE

## Single-gene scatter
single_gene_key <- "MSN_D1_matrix__CaH"
single_gene <- "LINC01735"
single_gene_main <- ""

## Metagene scatter
microglia_key <- "microglia__CaH"
microglia_priming <- c(
    "MS4A6A", "MS4A4A", "MS4A4E", "CD163", "GPNMB", "IL15", "AIM2", "NOD2",
    "TRIM14", "TNFAIP8", "DOCK5", "ANXA2", "CEACAM1", "PLXNA1", "PLXNC1",
    "GAS7", "ESR1", "IFIH1", "DDX60L", "HERC5", "PARP9", "DTX3L", "HLA-C",
    "ERAP2", "BTN3A2"
)
microglia_main <- ""

## -----------------------
## Execution
## -----------------------

cell_types_use <- read_cell_types(ct_file)

cell_metadata <- read_cell_metadata(cell_metadata_file)
donor_ages <- extract_donor_ages(cell_metadata)

tmp <- read_metacells(
    metacells_file,
    cell_types_use = cell_types_use,
    regions_use = regions_use_metacells
)

metacell_cr_list <- split_metacells_by_cell_type_region(
    tmp$metacells,
    tmp$col_metadata,
    donor_ages
)

## Heatmap
plot_donor_gex_age_heatmap(
    metacell_cr_list[[heatmap_key]],
    gs = gs,
    donor_ages = donor_ages,
    gs_gaps = gs_gaps,
    cluster_gs = heatmap_cluster_gs,
    transpose = heatmap_transpose
)

## Single gene scatter
plot_donor_gex_age_scatterplot(
    metacell_cr_list[[single_gene_key]][single_gene, ],
    donor_ages,
    main = single_gene_main
)

## Metagene scatter (sum)
plot_donor_gex_age_scatterplot(
    Matrix::colSums(metacell_cr_list[[microglia_key]][microglia_priming, , drop = FALSE]),
    donor_ages,
    main = microglia_main
)
