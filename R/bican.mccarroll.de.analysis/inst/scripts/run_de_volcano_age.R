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
ct="microglia"

fdr_cutoff <- 0.05
abs_log_fc_cutoff <- base::log2(1.05)

chr_color_map=c("X"="cornflowerblue", "Y"="tomato", "autosome"="black")

## -----------------------
## Execution
## -----------------------

gene_to_chr <- bican.mccarroll.de.analysis::read_gene_to_chr(gene_to_chr_file)
de_age <- bican.mccarroll.de.analysis::read_de_results(de_dir, test, ct_file, gene_to_chr)

#explicitly remove mitochondria
de_age=de_age[de_age$chr!="M",]
#map remaining genes to "X", "Y", or "autosome"
de_age$chr <- ifelse(de_age$chr %in% c("X", "Y"), de_age$chr, "autosome")

p<-bican.mccarroll.de.analysis::plot_de_volcano_gg(
    de_age,
    cell_type_use = ct,
    region_use = region_use,
    fdr_cutoff = fdr_cutoff,
    abs_log_fc_cutoff = abs_log_fc_cutoff,
    show_title = FALSE,
    chr_color_map=chr_color_map
)

p <- p +
    ggplot2::theme(
        axis.title.x = ggplot2::element_text(size = ggplot2::rel(1.5)),
        axis.title.y = ggplot2::element_text(size = ggplot2::rel(1.5)),
        legend.text = ggplot2::element_text(size = ggplot2::rel(1.25))
    ) +
    ggplot2::guides(
        color = ggplot2::guide_legend(override.aes = list(size = 2))
    )

print (p)

####################
# TEST SEX
####################
de_sex <- bican.mccarroll.de.analysis::read_de_results(de_dir, test="female_vs_male", ct_file, gene_to_chr)

p2<-plot_de_volcano_gg(
    de_sex,
    cell_type_use = ct,
    region_use = region_use,
    fdr_cutoff = fdr_cutoff,
    abs_log_fc_cutoff = abs_log_fc_cutoff,
    show_title = FALSE,
    chr_color_map=chr_color_map
)

print (p2)




