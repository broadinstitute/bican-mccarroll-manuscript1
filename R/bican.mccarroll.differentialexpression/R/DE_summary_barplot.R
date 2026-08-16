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
#     pdf_output_file = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/external_comparison_snap200/results/sex_effects_by_chromosome.pdf",
#     svg_barplot_file = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/external_comparison_snap200/results/sex_effects_by_chromosome_barplot.svg",
#     svg_density_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/external_comparison_snap200/results/sex_effects_density"
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
#' @param svg_output_file Optional path for writing the plot to an SVG file.
#'   When `NULL`, no SVG is written.
#'
#' @return A `ggplot` object.
#'
#' @export
barplot_de_results <- function(
        in_dir,
        file_pattern,
        cellTypeListFile = NULL,
        pdf_output_file = NULL,
        svg_output_file = NULL) {

    # Make R CMD CHECK happy
    . <- adj.P.Val <- logFC <- cell_type <- n_up <- n_down <- direction <- n_genes <- NULL

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

    y_max <- max(de_counts_long$n_genes, na.rm = TRUE)

    if (!is.finite(y_max) || y_max <= 0) {
        y_max <- 1
    }

    outlier_region <- find_cell_type_outlier(de_counts_long)

    other_data <- de_counts_long[
        de_counts_long$interaction != outlier_region,
        ,
        drop = FALSE
    ]

    outlier_data <- de_counts_long[
        de_counts_long$interaction == outlier_region,
        ,
        drop = FALSE
    ]

    p_other <- ggplot2::ggplot(
        other_data,
        ggplot2::aes(
            x = cell_type,
            y = n_genes,
            fill = direction
        )
    ) +
        ggplot2::geom_col(
            position = ggplot2::position_dodge(width = 0.8),
            width = 0.7
        ) +
        ggplot2::facet_wrap(
            ~ interaction,
            ncol = 2
        ) +
        ggplot2::scale_x_discrete(
            labels = function(x) gsub("_", " ", x, fixed = TRUE)
        ) +
        ggplot2::scale_y_continuous(
            limits = c(0, y_max)
        ) +
        ggplot2::labs(
            title = contrast_title,
            x = NULL,
            y = NULL,
            fill = "Direction"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5),
            axis.text.x = ggplot2::element_text(
                angle = 45,
                hjust = 1
            ),
            legend.position = "top"
        )

    p_outlier <- ggplot2::ggplot(
        outlier_data,
        ggplot2::aes(
            x = cell_type,
            y = n_genes,
            fill = direction
        )
    ) +
        ggplot2::geom_col(
            position = ggplot2::position_dodge(width = 0.8),
            width = 0.7
        ) +
        ggplot2::facet_wrap(
            ~ interaction
        ) +
        ggplot2::scale_x_discrete(
            labels = function(x) gsub("_", " ", x, fixed = TRUE)
        ) +
        ggplot2::scale_y_continuous(
            limits = c(0, y_max)
        ) +
        ggplot2::labs(
            x = NULL,
            y = NULL
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(
                angle = 45,
                hjust = 1
            ),
            legend.position = "none"
        )

    p_panels <- cowplot::plot_grid(
        p_other,
        p_outlier,
        ncol = 1,
        rel_heights = c(2, 1),
        align = "v",
        axis = "lr"
    )

    .save_de_plot(
        p_panels,
        pdf_output_file = pdf_output_file,
        svg_output_file = svg_output_file
    )

    p_panels
}

#' Write a plot to PDF and/or SVG
#'
#' Shared output helper used across this file's plotting functions. Either
#' path may be \code{NULL}, in which case that format is skipped, so callers
#' can support PDF, SVG, or both from the same pair of parameters.
#'
#' @param plot A `ggplot` object.
#' @param pdf_output_file Optional path for writing the plot to a PDF file.
#' @param svg_output_file Optional path for writing the plot to an SVG file.
#' @param width Optional plot width passed to `ggplot2::ggsave()`.
#' @param height Optional plot height passed to `ggplot2::ggsave()`.
#'
#' @return Invisibly returns `NULL`.
.save_de_plot <- function(
        plot,
        pdf_output_file = NULL,
        svg_output_file = NULL,
        width = NULL,
        height = NULL) {

    save_one <- function(path) {
        output_dir <- dirname(path)
        if (!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE)
        }

        args <- list(filename = path, plot = plot)
        if (!is.null(width)) args$width <- width
        if (!is.null(height)) args$height <- height

        do.call(ggplot2::ggsave, args)
    }

    if (!is.null(pdf_output_file)) save_one(pdf_output_file)
    if (!is.null(svg_output_file)) save_one(svg_output_file)

    invisible(NULL)
}


.sanitize_filename <- function(x) {
    gsub("[^A-Za-z0-9_.-]+", "_", x)
}


#find the region that is the outlier in cell type usage.  Usually DFC.
find_cell_type_outlier <- function(data) {
    cell_types <- split(
        as.character(data$cell_type),
        as.character(data$interaction)
    )
    cell_types <- lapply(cell_types, unique)

    regions <- names(cell_types)
    distances <- matrix(
        0,
        nrow = length(regions),
        ncol = length(regions),
        dimnames = list(regions, regions)
    )

    for (i in seq_along(regions)) {
        for (j in seq_along(regions)) {
            shared <- length(intersect(
                cell_types[[i]],
                cell_types[[j]]
            ))
            total <- length(union(
                cell_types[[i]],
                cell_types[[j]]
            ))

            distances[i, j] <- 1 - shared / total
        }
    }

    regions[which.max(rowMeans(distances))]
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

    # Make R CMD CHECK happy
    . <- gene <- chr <- annotationType <- gene_name <- chromosome_group <- i.chromosome_group <- NULL

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

    # Make R CMD CHECK happy
    . <- cell_type <- interaction <- chromosome_group <- n_genes <- NULL

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

    # Make R CMD CHECK happy
    cell_type_label <- cell_type <- n_genes <- chromosome_group <- NULL

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

    # Make R CMD CHECK happy
    interaction <- cell_type <- chromosome_group <- logFC <- cell_type_label <- NULL

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
#' X/Y, and writes count and signed-effect density plots. Three independent
#' output paths are supported, each defaulting to \code{NULL} (skipped):
#' a multi-page \code{pdf_output_file} (the summary count plot followed by
#' one density plot per region, as before), a single \code{svg_barplot_file}
#' holding just the summary count plot, and \code{svg_density_dir}, a
#' directory that receives one signed-effect density SVG per region (all
#' cell types faceted within a single plot, as in the PDF pages), named
#' \code{<svg_density_prefix><region>.svg}.
#'
#' @param in_dir Directory containing differential-expression result files.
#' @param file_pattern Pattern used by \code{parse_de_inputs()}.
#' @param contig_yaml_file Path to the contig classification YAML file.
#' @param reduced_gtf_file Path to the reduced GTF-like annotation file.
#' @param cellTypeListFile Optional file containing cell types to include.
#' @param alpha Adjusted P-value threshold.
#' @param pdf_output_file Optional path for the output multi-page PDF. When
#'   \code{NULL}, no PDF is written.
#' @param svg_barplot_file Optional path for writing the summary count plot
#'   to an SVG file. When \code{NULL}, no barplot SVG is written.
#' @param svg_density_dir Optional directory for writing one signed-effect
#'   density SVG per region. When \code{NULL}, no density SVGs are
#'   written.
#' @param svg_density_prefix Filename prefix used for each density SVG
#'   written to \code{svg_density_dir}.
#'
#' @return Invisibly returns the significant results, count data, and plots.
#' @export
plot_sex_de_by_chromosome <- function(
        in_dir,
        file_pattern,
        contig_yaml_file,
        reduced_gtf_file,
        cellTypeListFile = NULL,
        alpha = 0.05,
        pdf_output_file = NULL,
        svg_barplot_file = NULL,
        svg_density_dir = NULL,
        svg_density_prefix = "") {

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

    if (!is.null(pdf_output_file)) {
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
    }

    if (!is.null(svg_barplot_file)) {
        .save_de_plot(
            count_plot,
            svg_output_file = svg_barplot_file,
            width = 11,
            height = 8.5
        )
    }

    if (!is.null(svg_density_dir)) {
        if (!dir.exists(svg_density_dir)) {
            dir.create(
                svg_density_dir,
                recursive = TRUE
            )
        }

        for (interaction_name in interaction_names) {
            out_file <- file.path(
                svg_density_dir,
                paste0(
                    svg_density_prefix,
                    .sanitize_filename(interaction_name),
                    ".svg"
                )
            )

            ggplot2::ggsave(
                filename = out_file,
                plot = density_plots[[interaction_name]],
                width = 11,
                height = 8.5
            )
        }
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
