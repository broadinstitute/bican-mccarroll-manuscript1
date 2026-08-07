# library (bican.mccarroll.differentialexpression)
# library(data.table)
# library(ggplot2)

# in_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/results/LEVEL_6/pmi_rqs/cell_type_region_interaction_absolute_effects"
# cellTypeListFile="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/metadata/cell_types_for_de_filtering_plot.txt"
# file_pattern="pmi_hr"


# plot_sex_de_by_chromosome(
#     in_dir = paste0(
#         "/broad/bican_um1_mccarroll/RNAseq/analysis/",
#         "CAP_freeze_3.1_analysis/differential_expression/results/",
#         "LEVEL_6/sex_age/cell_type_region_interaction_absolute_effects"
#     ),
#     file_pattern = "female_vs_male",
#     contig_yaml_file = paste0(
#         "/broad/bican_um1_mccarroll/RNAseq/analysis/",
#         "CAP_freeze_3.1_analysis/metadata/",
#         "GRCh38_ensembl_v43.contig_groups.yaml"
#     ),
#     reduced_gtf_file = paste0(
#         "/broad/bican_um1_mccarroll/RNAseq/analysis/",
#         "CAP_freeze_3.1_analysis/metadata/",
#         "GRCh38_ensembl_v43.reduced.gtf.gz"
#     ),
#     cellTypeListFile="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/metadata/cell_types_for_de_filtering_plot.txt",
#     pdf_output_file = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/external_comparison_snap200/results/sex_effects_by_chromosome.pdf"
# )


#' Plot differential expression result counts
#'
#' Parses differential expression result files and creates faceted bar plots
#' showing the number of genes with positive and negative log fold changes
#' passing an adjusted P-value threshold of 0.05.
#'
#' @param in_dir Directory containing the differential expression result files.
#' @param file_pattern Pattern used to identify differential expression result
#'   files.
#' @param cellTypeListFile Optional path to a file specifying the cell types to
#'   include.
#' @param pdf_output_file Optional path for writing the plot to a PDF file.
#'   When `NULL`, no PDF is written.
#'
#' @return A `ggplot` object.
#'
#' @export
barplot_de_results <- function(
        in_dir,
        file_pattern,
        cellTypeListFile = NULL,
        pdf_output_file = NULL) {

    d <- parse_de_inputs(in_dir, file_pattern, cellTypeListFile)

    setDT(d)

    de_counts <- d[
        adj.P.Val < 0.05 & logFC != 0,
        .(
            n_up = sum(logFC > 0),
            n_down = sum(logFC < 0)
        ),
        by = .(cell_type, interaction)
    ]

    # Include cell type/interaction combinations with zero significant genes
    all_groups <- unique(d[, .(cell_type, interaction)])

    de_counts <- merge(
        all_groups,
        de_counts,
        by = c("cell_type", "interaction"),
        all.x = TRUE
    )

    de_counts[is.na(n_up), n_up := 0L]
    de_counts[is.na(n_down), n_down := 0L]

    setorder(de_counts, interaction, cell_type)

    de_counts_long <- melt(
        de_counts,
        id.vars = c("cell_type", "interaction"),
        measure.vars = c("n_up", "n_down"),
        variable.name = "direction",
        value.name = "n_genes"
    )

    de_counts_long[
        ,
        direction := factor(
            direction,
            levels = c("n_up", "n_down"),
            labels = c("Positive logFC", "Negative logFC")
        )
    ]

    contrast_title <- unique(d$contrast)

    if (length(contrast_title) != 1L) {
        stop("Expected exactly one unique contrast value in d.")
    }

    p <- ggplot(
        de_counts_long,
        aes(
            x = cell_type_label,
            y = n_genes,
            fill = direction
        )
    ) +
        geom_col(
            position = position_dodge(width = 0.8),
            width = 0.7
        ) +
        facet_wrap(
            ~ interaction,
            ncol = 2
        ) +
        labs(
            title = contrast_title,
            x = NULL,
            y = "Number of genes with adjusted P-value < 0.05",
            fill = "Direction"
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(
                angle = 45,
                hjust = 1
            ),
            legend.position = "top"
        )

    if (!is.null(pdf_output_file)) {
        ggplot2::ggsave(
            filename = pdf_output_file,
            plot = p
        )
    }

    p
}

#######################################################
# PLOTS FOR SEX EFFECTS PARTITIONED BY X/Y vs Autosome
#######################################################

#' Annotate differential-expression results by chromosome group
#'
#' @param d A data frame or data.table containing a \code{gene} column.
#' @param contig_yaml_file Path to a YAML file mapping contigs to classes.
#' @param reduced_gtf_file Path to the reduced GTF-like annotation file.
#'
#' @return A data.table with an added \code{chromosome_group} column.
#' @export
annotate_de_chromosome_group <- function(
        d,
        contig_yaml_file,
        reduced_gtf_file) {

    gene <- chr <- annotationType <- gene_name <- chromosome_group <- NULL

    data.table::setDT(d)

    contig_groups <- unlist(
        yaml::yaml.load_file(contig_yaml_file)
    )

    autosomes <- names(contig_groups[contig_groups == "autosome"])

    gtf <- data.table::fread(
        reduced_gtf_file,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE
    )

    gene_annotation <- unique(
        gtf[
            annotationType == "gene" &
                chr %in% c(autosomes, "chrX", "chrY"),
            .(
                gene = gene_name,
                chromosome_group = data.table::fifelse(
                    chr %in% autosomes,
                    "Autosome",
                    "X/Y"
                )
            )
        ],
        by = "gene"
    )

    d[
        gene_annotation,
        chromosome_group := i.chromosome_group,
        on = "gene"
    ]

    d[
        ,
        chromosome_group := factor(
            chromosome_group,
            levels = c("Autosome", "X/Y")
        )
    ]

    d
}


#' Count significant genes by chromosome group
#'
#' @param d Annotated differential-expression results.
#'
#' @return A data.table with gene counts by cell type, interaction, and
#'   chromosome group.
make_sex_de_chromosome_counts <- function(d) {

    cell_type <- interaction <- chromosome_group <- n_genes <- NULL

    observed_counts <- d[
        ,
        .(n_genes = .N),
        by = .(
            cell_type,
            interaction,
            chromosome_group
        )
    ]

    all_groups <- unique(
        d[, .(cell_type, interaction)]
    )

    complete_groups <- all_groups[
        ,
        .(
            chromosome_group = factor(
                c("Autosome", "X/Y"),
                levels = c("Autosome", "X/Y")
            )
        ),
        by = .(
            cell_type,
            interaction
        )
    ]

    counts <- merge(
        complete_groups,
        observed_counts,
        by = c(
            "cell_type",
            "interaction",
            "chromosome_group"
        ),
        all.x = TRUE
    )

    counts[
        is.na(n_genes),
        n_genes := 0L
    ]

    data.table::setorder(
        counts,
        interaction,
        cell_type,
        chromosome_group
    )

    counts
}


#' Plot significant-gene counts by chromosome group
#'
#' @param counts Output from \code{make_sex_de_chromosome_counts()}.
#' @param title Plot title.
#' @param alpha Adjusted P-value threshold used to select genes.
#'
#' @return A ggplot object.
plot_sex_de_chromosome_counts <- function(
        counts,
        title,
        alpha = 0.05) {

    #to make the cell type labels look nicer.
    counts[, cell_type_label := gsub("_", " ", cell_type)]

    ggplot2::ggplot(
        counts,
        ggplot2::aes(
            x = cell_type_label,
            y = n_genes,
            fill = chromosome_group
        )
    ) +
        ggplot2::geom_col(
            position = ggplot2::position_dodge(width = 0.8),
            width = 0.7
        ) +
        ggplot2::facet_wrap(
            ~ interaction,
            ncol = 2,
            scales = "free_x"
        ) +
        ggplot2::labs(
            title = title,
            subtitle = sprintf(
                "Significant genes by chromosome group; adjusted P-value < %.2g",
                alpha
            ),
            x = NULL,
            y = "Number of significant genes",
            fill = "Chromosome group"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5),
            plot.subtitle = ggplot2::element_text(hjust = 0.5),
            axis.text.x = ggplot2::element_text(
                angle = 45,
                hjust = 1
            ),
            legend.position = "top"
        )
}


#' Plot signed effect-size densities for one interaction
#'
#' @param d Significant annotated differential-expression results.
#' @param interaction_name Interaction or region to plot.
#' @param title Plot title.
#' @param alpha Adjusted P-value threshold used to select genes.
#'
#' @return A ggplot object.
plot_sex_de_chromosome_density <- function(
        d,
        interaction_name,
        title,
        alpha = 0.05) {

    interaction <- cell_type <- chromosome_group <- logFC <- NULL

    plot_data <- d[
        interaction == interaction_name
    ]

    plot_data[, cell_type_label := gsub("_", " ", cell_type)]

    ggplot2::ggplot(
        plot_data,
        ggplot2::aes(
            x = logFC,
            color = chromosome_group,
            fill = chromosome_group
        )
    ) +
        ggplot2::geom_density(
            alpha = 0.2,
            linewidth = 0.8,
            na.rm = TRUE
        ) +
        ggplot2::geom_vline(
            xintercept = 0,
            linetype = "dashed",
            linewidth = 0.4
        ) +
        ggplot2::facet_wrap(
            ~ cell_type_label,
            scales = "free_y"
        ) +
        ggplot2::labs(
            title = sprintf(
                "%s: %s",
                title,
                interaction_name
            ),
            subtitle = sprintf(
                "Signed logFC distributions; adjusted P-value < %.2g",
                alpha
            ),
            x = "Signed log fold change",
            y = "Density",
            color = "Chromosome group",
            fill = "Chromosome group"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5),
            plot.subtitle = ggplot2::element_text(hjust = 0.5),
            legend.position = "top"
        )
}


#' Plot significant sex effects by chromosome group
#'
#' Parses differential-expression results, annotates genes as autosomal or
#' X/Y, and writes count and signed-effect density plots to a PDF.
#'
#' @param in_dir Directory containing differential-expression result files.
#' @param file_pattern Pattern used by \code{parse_de_inputs()}.
#' @param contig_yaml_file Path to the contig classification YAML file.
#' @param reduced_gtf_file Path to the reduced GTF-like annotation file.
#' @param pdf_output_file Path for the output PDF.
#' @param cellTypeListFile Optional file containing cell types to include.
#' @param alpha Adjusted P-value threshold.
#'
#' @return Invisibly returns the significant results, count data, and plots.
#' @export
plot_sex_de_by_chromosome <- function(
        in_dir,
        file_pattern,
        contig_yaml_file,
        reduced_gtf_file,
        pdf_output_file,
        cellTypeListFile = NULL,
        alpha = 0.05) {

    adj.P.Val <- logFC <- chromosome_group <- contrast <- interaction <- NULL

    d <- parse_de_inputs(
        in_dir,
        file_pattern,
        cellTypeListFile
    )

    data.table::setDT(d)

    d <- annotate_de_chromosome_group(
        d,
        contig_yaml_file,
        reduced_gtf_file
    )

    significant <- d[
        adj.P.Val < alpha &
            is.finite(logFC) &
            logFC != 0 &
            !is.na(chromosome_group)
    ]

    if (nrow(significant) == 0L) {
        stop(
            "No significant autosomal or X/Y genes were found.",
            call. = FALSE
        )
    }

    contrast_title <- "Sex differential expression"

    counts <- make_sex_de_chromosome_counts(significant)

    count_plot <- plot_sex_de_chromosome_counts(
        counts,
        title = contrast_title,
        alpha = alpha
    )

    interaction_names <- sort(
        unique(significant$interaction)
    )

    density_plots <- lapply(
        interaction_names,
        function(interaction_name) {
            plot_sex_de_chromosome_density(
                significant,
                interaction_name = interaction_name,
                title = contrast_title,
                alpha = alpha
            )
        }
    )

    names(density_plots) <- interaction_names

    output_dir <- dirname(pdf_output_file)

    if (!dir.exists(output_dir)) {
        dir.create(
            output_dir,
            recursive = TRUE
        )
    }

    grDevices::pdf(
        pdf_output_file,
        width = 11,
        height = 8.5,
        onefile = TRUE
    )

    on.exit(
        grDevices::dev.off(),
        add = TRUE
    )

    print(count_plot)

    for (p in density_plots) {
        print(p)
    }

    invisible(
        list(
            significant_genes = significant,
            count_summary = counts,
            count_plot = count_plot,
            density_plots = density_plots
        )
    )
}
