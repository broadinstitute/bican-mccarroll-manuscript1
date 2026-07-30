# library (bican.mccarroll.differentialexpression)
# library(data.table)
# library(ggplot2)

# in_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/results/LEVEL_6/pmi_rqs/cell_type_region_interaction_absolute_effects"
# cellTypeListFile="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/metadata/cell_types_for_de_filtering_plot.txt"
# file_pattern="pmi_hr"



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
            x = cell_type,
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
