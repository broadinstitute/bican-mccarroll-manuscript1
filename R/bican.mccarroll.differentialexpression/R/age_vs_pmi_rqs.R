# library (bican.mccarroll.differentialexpression)
# library(data.table)
# library(ggplot2)
# library(cowplot)
#
# in_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/results/LEVEL_6/age_pmi_rqs_complete_cases/age_pmi_rqs"
# cellTypeListFile="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/metadata/cell_types_for_de_filtering_plot.txt"
# out_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/results/LEVEL_6/age_pmi_rqs_complete_cases/age_pmi_rqs_effect_comparisons.pdf"
#
#
#
#
# in_dir_full_model="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/results/LEVEL_6/age_pmi_rqs_complete_cases/age_pmi_rqs"
# in_dir_partial_model="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/results/LEVEL_6/age_pmi_rqs_complete_cases/age"
# output_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/results/LEVEL_6/age_pmi_rqs_complete_cases/age_effects_with_without_rqs.pdf"

# run age DE with and without RQS, and plot the age coefficients.
# in both cases, we subset to the complete cases of donors with RQS, PMI, age, then test nested models.


#' Compare two nested models with and without the tested trait vs age
#'
#' Does adding the tested trait change the age effect size?
#'
#' The setup for this test is important - it's critical to have the exact same
#' set of data points in both model runs, and only change the model by adding the
#' tested trait to the full model.
#'
#' @param tested_trait The covariate added in the full model
#' @param in_dir_full_model The directory containing the full model results
#' @param in_dir_partial_model The directory containing the partial results
#' @param cellTypeListFile A file containing a list of cell types to generate plots for.
#' If null generate plots for all cell types.
#' @param output_file A PDF report with summary plots and per cell type / region plots of the effect sizes.
#' @export
#' @import cowplot ggplot2 data.table
plot_nested_models<-function (tested_trait="RQS", in_dir_full_model, in_dir_partial_model, cellTypeListFile=NULL, output_file) {

    full_model_de=parse_de_inputs(in_dir_full_model, file_pattern="age_DE", cellTypeListFile)
    partial_model_de=parse_de_inputs(in_dir_partial_model, file_pattern="age_DE", cellTypeListFile)

    write_age_model_comparison_pdf(
        tested_trait=tested_trait,
        full_model_de = full_model_de,
        partial_model_de = partial_model_de,
        output_file = output_file
    )

}


plot_age_vs_pmi_rqs<-function (in_dir, cellTypeListFile, out_file) {

    age=parse_de_inputs(in_dir, file_pattern="age_DE", cellTypeListFile)
    pmi=parse_de_inputs(in_dir, file_pattern="pmi_hr_DE", cellTypeListFile)
    rqs=parse_de_inputs(in_dir, file_pattern="RQS_DE", cellTypeListFile)


    plot_effect_comparison_page(
        age = age,
        pmi = pmi,
        rqs = rqs,
        cell_type_value = "astrocyte",
        region_value = "CaH"
    )

    write_effect_comparison_pdf(
        age = age,
        pmi = pmi,
        rqs = rqs,
        output_file = out_file
    )


}

plot_effect_comparison_page <- function(
        age,
        pmi,
        rqs,
        cell_type_value,
        region_value) {

    age_pmi <- plot_effect_comparison(
        x = age,
        y = pmi,
        x_label = "Age",
        y_label = "PMI",
        cell_type_value = cell_type_value,
        region_value = region_value
    )

    age_rqs <- plot_effect_comparison(
        x = age,
        y = rqs,
        x_label = "Age",
        y_label = "RQS",
        cell_type_value = cell_type_value,
        region_value = region_value
    )

    pmi_rqs <- plot_effect_comparison(
        x = pmi,
        y = rqs,
        x_label = "PMI",
        y_label = "RQS",
        cell_type_value = cell_type_value,
        region_value = region_value
    )

    plot_grid <- cowplot::plot_grid(
        age_pmi,
        age_rqs,
        pmi_rqs,
        NULL,
        ncol = 2,
        labels = c("A", "B", "C", "")
    )

    page_title <- paste(
        region_value,
        cell_type_value,
        sep = " - "
    )

    cowplot::ggdraw() +
        cowplot::draw_plot(
            plot_grid,
            x = 0,
            y = 0,
            width = 1,
            height = 0.94
        ) +
        cowplot::draw_label(
            page_title,
            x = 0.5,
            y = 0.985,
            hjust = 0.5,
            vjust = 1,
            fontface = "bold",
            size = 16
        )
}


plot_effect_comparison <- function(
        x,
        y,
        x_label,
        y_label,
        cell_type_value,
        region_value) {

    # Make R CMD CHECK happy
    cell_type <- gene <- logFC <- x_effect <- y_effect <- NULL

    x <- data.table::as.data.table(x)
    y <- data.table::as.data.table(y)

    x_subset <- x[
        cell_type == cell_type_value &
            interaction == region_value,
        list(
            gene,
            x_effect = logFC
        )
    ]

    y_subset <- y[
        cell_type == cell_type_value &
            interaction == region_value,
        list(
            gene,
            y_effect = logFC
        )
    ]

    plot_data <- merge(
        x_subset,
        y_subset,
        by = "gene"
    )

    comparison_title <- paste(x_label, "vs", y_label)

    if (nrow(plot_data) < 3L) {
        return(
            ggplot() +
                annotate(
                    "text",
                    x = 0,
                    y = 0,
                    label = paste0(
                        comparison_title,
                        "\nFewer than 3 intersecting genes"
                    ),
                    size = 4
                ) +
                xlim(-1, 1) +
                ylim(-1, 1) +
                ggtitle(comparison_title) +
                theme_void() +
                theme(
                    plot.title = element_text(
                        hjust = 0.5,
                        face = "bold"
                    )
                )
        )
    }

    fit <- lm(y_effect ~ x_effect, data = plot_data)

    slope <- unname(coef(fit)[["x_effect"]])

    correlation <- cor(
        plot_data$x_effect,
        plot_data$y_effect,
        method = "pearson"
    )

    annotation <- sprintf(
        "Pearson r = %.3f\nSlope = %.3f\nGenes = %s",
        correlation,
        slope,
        format(nrow(plot_data), big.mark = ",")
    )

    ggplot(
        plot_data,
        aes(x = x_effect, y = y_effect)
    ) +
        geom_hline(
            yintercept = 0,
            linewidth = 0.3,
            color = "grey75"
        ) +
        geom_vline(
            xintercept = 0,
            linewidth = 0.3,
            color = "grey75"
        ) +
        ggrastr::geom_point_rast(
            alpha = 0.5,
            size = 1,
            raster.dpi = 300
        ) +
        geom_smooth(
            method = "lm",
            formula = y ~ x,
            se = FALSE,
            linewidth = 0.7
        ) +
        annotate(
            "text",
            x = Inf,
            y = -Inf,
            label = annotation,
            hjust = 1.05,
            vjust = -0.4,
            size = 3.5
        ) +
        labs(
            title = comparison_title,
            x = paste(x_label, "logFC"),
            y = paste(y_label, "logFC")
        ) +
        theme_cowplot() +
        theme(
            plot.title = element_text(
                hjust = 0.5,
                face = "bold"
            ),
            plot.margin = margin(8, 12, 12, 8)
        )
}




write_effect_comparison_pdf <- function(
        age,
        pmi,
        rqs,
        output_file = "age_pmi_rqs_effect_comparisons.pdf") {

    # Make R CMD CHECK happy
    . <- cell_type <- NULL

    age <- as.data.table(age)
    pmi <- as.data.table(pmi)
    rqs <- as.data.table(rqs)

    required_columns <- c(
        "gene",
        "cell_type",
        "interaction",
        "logFC"
    )

    inputs <- list(
        age = age,
        pmi = pmi,
        rqs = rqs
    )

    missing_columns <- lapply(
        inputs,
        function(x) setdiff(required_columns, names(x))
    )

    has_missing_columns <- lengths(missing_columns) > 0L

    if (any(has_missing_columns)) {
        missing_message <- vapply(
            names(inputs)[has_missing_columns],
            function(input_name) {
                sprintf(
                    "%s: %s",
                    input_name,
                    paste(
                        missing_columns[[input_name]],
                        collapse = ", "
                    )
                )
            },
            character(1)
        )

        stop(
            "Missing required columns:\n",
            paste(missing_message, collapse = "\n"),
            call. = FALSE
        )
    }

    combinations <- Reduce(
        merge,
        list(
            unique(age[, .(cell_type, interaction)]),
            unique(pmi[, .(cell_type, interaction)]),
            unique(rqs[, .(cell_type, interaction)])
        )
    )

    data.table::setorder(
        combinations,
        cell_type,
        interaction
    )

    if (!nrow(combinations)) {
        stop(
            paste(
                "No cell-type and region combinations",
                "are shared across all three data frames."
            ),
            call. = FALSE
        )
    }

    grDevices::pdf(
        output_file,
        width = 12,
        height = 10,
        onefile = TRUE
    )

    on.exit(
        grDevices::dev.off(),
        add = TRUE
    )

    for (i in seq_len(nrow(combinations))) {
        page <- plot_effect_comparison_page(
            age = age,
            pmi = pmi,
            rqs = rqs,
            cell_type_value = combinations$cell_type[[i]],
            region_value = combinations$interaction[[i]]
        )

        print(page)
    }

    invisible(combinations)
}

##########################
# NESTED AGE PLOTS
##########################

plot_age_model_comparison_page <- function(
        full_model_de,
        partial_model_de,
        cell_type_value,
        tested_trait="RQS") {

    # Make R CMD CHECK happy
    cell_type <- NULL

    full_model_de <- data.table::as.data.table(full_model_de)
    partial_model_de <- data.table::as.data.table(partial_model_de)

    full_regions <- unique(
        full_model_de[
            cell_type == cell_type_value,
            interaction
        ]
    )

    partial_regions <- unique(
        partial_model_de[
            cell_type == cell_type_value,
            interaction
        ]
    )

    regions <- sort(intersect(full_regions, partial_regions))

    if (!length(regions)) {
        stop(
            "No shared regions found for cell type: ",
            cell_type_value,
            call. = FALSE
        )
    }

    plots <- lapply(
        regions,
        function(region_value) {
            plot_effect_comparison(
                x = full_model_de,
                y = partial_model_de,
                x_label = "Full model age effect",
                y_label = paste("Model without",tested_trait, "age effect"),
                cell_type_value = cell_type_value,
                region_value = region_value
            ) +
                ggplot2::labs(title = region_value)
        }
    )

    n_plots <- length(plots)

    ncol <- if (n_plots == 1L) {
        1L
    } else if (n_plots <= 4L) {
        2L
    } else {
        3L
    }

    plot_grid <- cowplot::plot_grid(
        plotlist = plots,
        ncol = ncol,
        align = "hv",
        axis = "tblr"
    )

    cowplot::ggdraw() +
        cowplot::draw_plot(
            plot_grid,
            x = 0,
            y = 0,
            width = 1,
            height = 0.94
        ) +
        cowplot::draw_label(
            cell_type_value,
            x = 0.5,
            y = 0.985,
            hjust = 0.5,
            vjust = 1,
            fontface = "bold",
            size = 16
        )
}

write_age_model_comparison_pdf <- function(
        tested_trait="RQS",
        full_model_de,
        partial_model_de,
        output_file = "age_effects_with_without_rqs.pdf") {

    # Make R CMD CHECK happy
    cell_type <- NULL

    full_model_de <- data.table::as.data.table(full_model_de)
    partial_model_de <- data.table::as.data.table(partial_model_de)

    model_summary <- summarize_model_differences(
        full_model_de,
        partial_model_de
    )

    p=plot_summarized_model_dif (model_summary)

    full_combinations <- unique(
        full_model_de[, list(cell_type, interaction)]
    )

    partial_combinations <- unique(
        partial_model_de[, list(cell_type, interaction)]
    )

    combinations <- merge(
        full_combinations,
        partial_combinations,
        by = c("cell_type", "interaction")
    )

    cell_types <- sort(unique(combinations$cell_type))

    if (!length(cell_types)) {
        stop(
            "No cell-type and region combinations are shared by both models.",
            call. = FALSE
        )
    }

    grDevices::pdf(
        output_file,
        width = 12,
        height = 10,
        onefile = TRUE
    )

    on.exit(
        grDevices::dev.off(),
        add = TRUE
    )

    print (p$cor_plot)
    print (p$abs_dif_plot)

    for (cell_type_value in cell_types) {
        page <- plot_age_model_comparison_page(
            full_model_de = full_model_de,
            partial_model_de = partial_model_de,
            cell_type_value = cell_type_value,
            tested_trait=tested_trait
        )

        print(page)
    }

    invisible(combinations)
}

plot_summarized_model_dif<-function (model_summary) {
    # Make R CMD CHECK happy
    comparison <- cell_type <- pearson_r <- median_abs_difference <- NULL

    model_summary[
        ,
        comparison := paste(cell_type, interaction, sep = " - ")
    ]

    data.table::setorder(
        model_summary,
        pearson_r
    )

    model_summary[
        ,
        comparison := factor(
            comparison,
            levels = comparison
        )
    ]

    p_correlation <- ggplot2::ggplot(
        model_summary,
        ggplot2::aes(
            x = 1 - pearson_r,
            y = comparison
        )
    ) +
        ggplot2::geom_point(size = 2) +
        ggplot2::labs(
            title = "Departure from perfect correlation",
            x = "1 - Pearson r",
            y = NULL
        ) +
        cowplot::theme_cowplot()

    p_difference <- ggplot2::ggplot(
        model_summary,
        ggplot2::aes(
            x = median_abs_difference,
            y = comparison
        )
    ) +
        ggplot2::geom_point(size = 2) +
        ggplot2::labs(
            title = "Absolute effect-size change",
            x = "Median absolute logFC difference",
            y = NULL
        ) +
        cowplot::theme_cowplot()

    return (list(cor_plot=p_correlation, abs_dif_plot=p_difference))

}

summarize_model_differences <- function(
        full_model_de,
        partial_model_de) {

    # Make R CMD CHECK happy
    gene <- cell_type <- interaction <- logFC <- adj.P.Val <-
        full_logFC <- partial_logFC <- full_adj_p <- partial_adj_p <- NULL

    full_model_de <- data.table::as.data.table(full_model_de)
    partial_model_de <- data.table::as.data.table(partial_model_de)

    merged <- merge(
        full_model_de[
            ,
            list(
                gene,
                cell_type,
                interaction,
                full_logFC = logFC,
                full_adj_p = adj.P.Val
            )
        ],
        partial_model_de[
            ,
            list(
                gene,
                cell_type,
                interaction,
                partial_logFC = logFC,
                partial_adj_p = adj.P.Val
            )
        ],
        by = c("gene", "cell_type", "interaction")
    )

    merged <- merged[
        is.finite(full_logFC) &
            is.finite(partial_logFC)
    ]

    merged[
        ,
        {
            fit <- stats::lm(partial_logFC ~ full_logFC)

            list(
                n_genes = .N,
                pearson_r = stats::cor(
                    full_logFC,
                    partial_logFC
                ),
                slope = unname(
                    stats::coef(fit)[["full_logFC"]]
                ),
                median_abs_difference = stats::median(
                    abs(partial_logFC - full_logFC)
                ),
                rmse = sqrt(
                    mean(
                        (partial_logFC - full_logFC)^2
                    )
                ),
                direction_change_fraction = mean(
                    sign(full_logFC) != sign(partial_logFC)
                ),
                significance_change_fraction = mean(
                    (full_adj_p < 0.05) !=
                        (partial_adj_p < 0.05)
                )
            )
        },
        by = list(cell_type, interaction)
    ]
}

# are PMI and RQS donor level traits?
# Yes, there's one unique set of values per donor.
tiny_adhoc<-function () {
    # Make R CMD CHECK happy
    pmi_hr <- RQS <- NULL

    a=read.table("/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/metacells/LEVEL_6/donor_rxn_DGEList_samples.tsv.gz", sep="\t", header=T)

    cols=c("pmi_hr", "RQS", "donor")

    aa=unique (a[,cols])
    length(unique(a$donor))

    ggplot2::ggplot(
        aa,
        ggplot2::aes(x = pmi_hr, y = RQS)
    ) +
        ggplot2::geom_point(
            size = 2,
            alpha = 0.6
        ) +
        ggplot2::geom_smooth(
            method = "lm",
            formula = y ~ x,
            se = TRUE,
            linewidth = 0.8
        ) +
        ggplot2::labs(
            title = "Relationship between PMI and RQS",
            x = "Postmortem interval (hours)",
            y = "RQS"
        ) +
        cowplot::theme_cowplot() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(
                hjust = 0.5,
                face = "bold"
            )
        )

}
