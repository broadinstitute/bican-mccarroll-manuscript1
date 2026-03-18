grey_white_comparison<-function () {
    #bican.mccarroll.differentialexpression::differential_expression_region(data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/emuratog", data_name="OPC_astro_01_clusters_DGEList", randVars=c('donor','village'), fixedVars=c('age', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'pct_intronic', 'frac_contamination', 'imputed_sex', 'single_cell_assay', 'region', 'biobank'), contrast_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/metadata/differential_expression_contrasts_sex_age.txt", cellTypeListFile=NULL, outPDF="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/emuratog/results/subset/volcano_plots.pdf",  result_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/emuratog/results/subset", n_cores = 14)

    gene_to_chr_file ="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/sburger_tmp/gene_to_chromosome.txt"
    de_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/emuratog/results/subset"
    outFile=paste(de_dir, "/grey_vs_white_plots.pdf", sep="")
    test="age"

    gene_to_chr <- bican.mccarroll.de.analysis::read_gene_to_chr(gene_to_chr_file)

    df <- bican.mccarroll.de.analysis::read_de_results(
        de_dir,
        test,
        ct_file=NULL,
        gene_to_chr)

    name_map = c(
        astrocyte_0 = "astrocyte grey",
        astrocyte_1 = "astrocyte white",
        OPC_0        = "OPC grey",
        OPC_1        = "OPC white"
    )

    df$cell_type = ifelse(
        df$cell_type %in% names(name_map),
        name_map[df$cell_type],
        df$cell_type
    )

    # Define the comparisons once, then iterate.
    plot_specs = list(
        # same color between regions
        list(cell_type_a="astrocyte grey",  cell_type_b="astrocyte grey",  region_a="CaH", region_b="Pu"),
        list(cell_type_a="OPC grey",        cell_type_b="OPC grey",        region_a="CaH", region_b="Pu"),

        # white between regions
        list(cell_type_a="astrocyte white", cell_type_b="astrocyte white", region_a="CaH", region_b="Pu"),
        list(cell_type_a="OPC white",       cell_type_b="OPC white",       region_a="CaH", region_b="Pu"),

        # grey vs white within regions
        list(cell_type_a="astrocyte white", cell_type_b="astrocyte grey",  region_a="CaH", region_b="CaH"),
        list(cell_type_a="astrocyte white", cell_type_b="astrocyte grey",  region_a="Pu",  region_b="Pu"),
        list(cell_type_a="OPC white",       cell_type_b="OPC grey",        region_a="CaH", region_b="CaH"),
        list(cell_type_a="OPC white",       cell_type_b="OPC grey",        region_a="Pu",  region_b="Pu")
    )

    plot_specs = list(
        # astrocyte plots
        list(cell_type_a="astrocyte white", cell_type_b="astrocyte grey", region_a="ic", region_b="CaH"),
        list(cell_type_a="astrocyte white", cell_type_b="astrocyte grey", region_a="ic", region_b="Pu"),
        list(cell_type_a="astrocyte white", cell_type_b="astrocyte grey", region_a="ic", region_b="NAC"),
        list(cell_type_a="astrocyte white", cell_type_b="astrocyte grey", region_a="ic", region_b="DFC"),

        # OPC plots
        list(cell_type_a="OPC white", cell_type_b="OPC grey", region_a="ic", region_b="CaH"),
        list(cell_type_a="OPC white", cell_type_b="OPC grey", region_a="ic", region_b="Pu"),
        list(cell_type_a="OPC white", cell_type_b="OPC grey", region_a="ic", region_b="NAC"),
        list(cell_type_a="OPC white", cell_type_b="OPC grey", region_a="ic", region_b="DFC")
    )



    grDevices::pdf(outFile)

    for (spec in plot_specs) {
        bican.mccarroll.de.analysis::plot_de_scatter_gg(
            de_dt = df,
            cell_type_a = spec$cell_type_a,
            cell_type_b = spec$cell_type_b,
            region_a = spec$region_a,
            region_b = spec$region_b,
            fdr_cutoff = 0.05,
            add_fit = TRUE,
            xlab_prefix = "DE age ",
            plot_only_significant=TRUE
        )
    }

    grDevices::dev.off()

    d=bican.mccarroll.differentialexpression::prepare_data_for_differential_expression(data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/emuratog", data_name="OPC_astro_01_clusters_DGEList", randVars=c('donor','village'), fixedVars=c('age', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'pct_intronic', 'frac_contamination', 'imputed_sex', 'single_cell_assay', 'region', 'biobank'))
    dge=d$dge;
    s=dge$samples

    s=s[s$cell_type %in% names(name_map),]

    s$cell_type = ifelse(
        s$cell_type %in% names(name_map),
        name_map[s$cell_type],
        s$cell_type
    )

    cell_type <- region <- lib.size <- NULL

    .get_field <- function(spec, nm) spec[[nm]]

    allowed <- unique(data.frame(
        cell_type = c(
            vapply(plot_specs, .get_field, character(1), nm = "cell_type_a"),
            vapply(plot_specs, .get_field, character(1), nm = "cell_type_b")
        ),
        region = c(
            vapply(plot_specs, .get_field, character(1), nm = "region_a"),
            vapply(plot_specs, .get_field, character(1), nm = "region_b")
        ),
        stringsAsFactors = FALSE
    ))

    # Keep only rows in s that match allowed (cell_type, region) pairs
    key_s <- paste(s$cell_type, s$region, sep = "||")
    key_a <- paste(allowed$cell_type, allowed$region, sep = "||")
    s2 <- s[key_s %in% key_a, , drop = FALSE]

    cell_type <- region <- lib.size <- NULL

    # Split "astrocyte grey" -> base="astrocyte", color="grey"
    s2$color <- sub(".*\\s+", "", s2$cell_type)
    s2$cell_type_base <- sub("\\s+(grey|white)$", "", s2$cell_type)

    # Optional: enforce ordering
    s2$color <- factor(s2$color, levels = c("grey", "white"))
    s2$region <- factor(s2$region, levels = c("CaH", "Pu"))
    s2$cell_type_base <- factor(s2$cell_type_base, levels = c("astrocyte", "OPC"))

    x_rng <- range(log10(s2$lib.size), na.rm = TRUE)

    # Make R CMD CHECK Happy

    p<-ggplot2::ggplot(s2, ggplot2::aes(x = log10(lib.size))) +
        ggplot2::geom_histogram(bins = 50) +
        ggplot2::facet_grid(color ~ region + cell_type_base, switch = "y") +
        ggplot2::coord_cartesian(xlim = x_rng) +
        ggplot2::labs(x = "log10 number of transcripts (lib.size)", y = "Number of samples") +
        ggplot2::theme_bw() +
        ggplot2::theme(
            strip.placement = "outside",
            strip.background = ggplot2::element_rect(fill = "grey90")
        )

    print(p)


}
