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

cell_types_use <- bican.mccarroll.de.analysis::read_cell_types(ct_file)

cell_metadata <- bican.mccarroll.de.analysis::read_cell_metadata(cell_metadata_file)
donor_ages <- bican.mccarroll.de.analysis::extract_donor_ages(cell_metadata)

tmp <- bican.mccarroll.de.analysis::read_metacells(
    metacells_file,
    cell_types_use = cell_types_use,
    regions_use = regions_use_metacells
)

metacell_cr_list <- bican.mccarroll.de.analysis::split_metacells_by_cell_type_region(
    tmp$metacells,
    tmp$col_metadata,
    donor_ages
)

## Heatmap
bican.mccarroll.de.analysis::plot_donor_gex_age_heatmap(
    metacell_cr_list[[heatmap_key]],
    gs = gs,
    donor_ages = donor_ages,
    gs_gaps = gs_gaps,
    cluster_gs = heatmap_cluster_gs,
    transpose = heatmap_transpose
)

## Single gene scatter
bican.mccarroll.de.analysis::plot_donor_gex_age_scatterplot(
    metacell_cr_list[[single_gene_key]][single_gene, ],
    donor_ages,
    main = single_gene_main,
    show_spearman=TRUE,
    size=8
)

## Metagene scatter (sum)
bican.mccarroll.de.analysis::plot_donor_gex_age_scatterplot(
    Matrix::colSums(metacell_cr_list[[microglia_key]][microglia_priming, , drop = FALSE]),
    donor_ages,
    main = microglia_main
)

#let's test the 5 x 5 layout of single gene scatter plots.
make_panel <- function(metacell_cr_list,
                       key,
                       gene,
                       donor_ages,
                       show_x = FALSE,
                       show_spearman = TRUE,
                       spearman_text_size = 4,
                       rho_threshold = 0.2,
                       y_axis_floor = 10) {

    exp_vec <- metacell_cr_list[[key]][gene, ]

    p <- bican.mccarroll.de.analysis::plot_donor_gex_age_scatterplot(
        exp_vec,
        donor_ages,
        main = "",
        show_spearman = show_spearman,
        size = spearman_text_size,
        rho_threshold = rho_threshold,
        y_axis_floor = y_axis_floor
    )

    p +
        ggplot2::theme_classic() +
        ggplot2::theme(
            axis.title = ggplot2::element_blank(),

            ## X axis: only show on bottom row
            axis.text.x  = if (show_x) ggplot2::element_text(size = 8) else ggplot2::element_blank(),
            axis.ticks.x = if (show_x) ggplot2::element_line() else ggplot2::element_blank(),

            ## Y axis: ALWAYS show for every panel
            axis.text.y  = ggplot2::element_text(size = 8),
            axis.ticks.y = ggplot2::element_line(),

            plot.margin = ggplot2::margin(2, 2, 2, 2)
        )
}

gene_list <- c("FKBP5","RGS9","RYR3","GRIA1","CLEC2B")

keys <- c("MSN_D1_matrix__CaH", "glutamatergic_L23IT__DFC", "astrocyte__CaH",
          "OPC__CaH", "oligodendrocyte__CaH", "microglia__CaH")


nice_names <- c("D1 matrix MSN", "L23IT glutamatergic\nneuron", "Astrocyte",
                "OPC", "Oligodendrocyte", "Microglia")

plot_list <- list()
idx <- 0L
n_row <- length(gene_list)
n_col <- length(keys)

for (r in seq_len(n_row)) {
    for (c in seq_len(n_col)) {
        idx <- idx + 1L
        plot_list[[idx]] <- make_panel(
            metacell_cr_list = metacell_cr_list,
            key = keys[c],
            gene = gene_list[r],
            donor_ages = donor_ages,
            show_x = (r == n_row),
            show_spearman = TRUE,
            spearman_text_size = 4
        )
    }
}

grid <- cowplot::plot_grid(
    plotlist = plot_list,
    nrow = n_row,
    align = "hv"
)


col_label_grobs <- lapply(nice_names, function(x) cowplot::ggdraw() + cowplot::draw_label(x, size = 12))
col_labels <- cowplot::plot_grid(plotlist = col_label_grobs, nrow = 1)

grid_with_cols <- cowplot::plot_grid(grid, col_labels, ncol = 1, rel_heights = c(1, 0.08))

row_label_grobs <- lapply(gene_list, function(g) cowplot::ggdraw() + cowplot::draw_label(g, angle = 90, size = 14))
row_labels <- cowplot::plot_grid(plotlist = row_label_grobs, ncol = 1)

all_plot <- cowplot::plot_grid(row_labels, grid_with_cols, nrow = 1, rel_widths = c(0.06, 1))
