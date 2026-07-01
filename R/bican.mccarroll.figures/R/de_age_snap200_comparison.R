library(ggplot2)

# plot_de_primary_secondary_manifest(
#     manifest_file = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/external_comparison_snap200/metadata/bican_dfc_snap_age_de_overlap_manifest.tsv",
#     out_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/external_comparison_snap200/results",
#     primary_dataset = "bican",
#     secondary_dataset = "snap",
#     primary_label = "BICAN",
#     secondary_label = "SNAP",
#     effect_name = "age effects"
# )

# plot_de_primary_secondary_manifest(
#     manifest_file = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/external_comparison_PMID_40903571/metadata/PMID_40903571_bican_dfc_age_de_overlap_manifest.tsv",
#     out_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/external_comparison_PMID_40903571/results",
#     primary_dataset = "PMID_40903571",
#     secondary_dataset = "bican",
#     primary_label = "PMID 40903571",
#     secondary_label = "BICAN",
#     effect_name = "age effects"
# )


## ------------------------------------------------------------------
## Top-level functions
## ------------------------------------------------------------------

plot_de_primary_secondary_manifest <- function(
        manifest_file,
        out_dir = NULL,
        primary_dataset = "bican",
        secondary_dataset = "snap",
        primary_label = NULL,
        secondary_label = NULL,
        effect_name = "age effects",
        alpha = 0.05,
        width = 7,
        height = 7) {

    if (is.null(primary_label)) {
        primary_label <- make_dataset_label(primary_dataset)
    }

    if (is.null(secondary_label)) {
        secondary_label <- make_dataset_label(secondary_dataset)
    }

    manifest <- read.table(
        manifest_file,
        header = TRUE,
        stringsAsFactors = FALSE,
        sep = "\t",
        check.names = FALSE
    )

    required_cols <- c("cell_type", primary_dataset, secondary_dataset)
    missing_cols <- setdiff(required_cols, colnames(manifest))

    if (length(missing_cols) > 0) {
        stop(sprintf(
            "Manifest is missing required columns: %s",
            paste(missing_cols, collapse = ", ")
        ))
    }

    if (!is.null(out_dir) && !dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }

    validate_manifest_files(
        manifest = manifest,
        file_cols = c(primary_dataset, secondary_dataset)
    )

    results <- vector("list", nrow(manifest))
    names(results) <- manifest$cell_type

    sign_results <- vector("list", nrow(manifest))
    names(sign_results) <- manifest$cell_type

    output_prefix <- make_output_prefix(primary_dataset, secondary_dataset)

    for (i in seq_len(nrow(manifest))) {
        cell_type <- manifest$cell_type[i]

        logger::log_info("Processing cell type ", cell_type)

        res <- plot_de_primary_secondary_cell_type(
            cell_type = cell_type,
            primary_file = manifest[[primary_dataset]][i],
            secondary_file = manifest[[secondary_dataset]][i],
            primary_dataset = primary_dataset,
            secondary_dataset = secondary_dataset,
            primary_label = primary_label,
            secondary_label = secondary_label,
            effect_name = effect_name,
            alpha = alpha
        )

        if (is.null(res)) {
            next
        }

        results[[cell_type]] <- res
        sign_results[[cell_type]] <- res$sign_test

        if (!is.null(out_dir)) {
            out_file <- file.path(
                out_dir,
                sprintf(
                    "%s_%s_scatter_%s.svg",
                    output_prefix,
                    sanitize_filename(effect_name),
                    sanitize_filename(cell_type)
                )
            )

            ggplot2::ggsave(
                filename = out_file,
                plot = res$plot,
                device = "svg",
                width = width,
                height = height
            )
        }
    }

    sign_results <- do.call(rbind, sign_results)

    if (is.null(sign_results) || nrow(sign_results) == 0) {
        warning("No cell types had overlapping primary-significant genes.")
        return(invisible(list(
            per_cell_type = results,
            sign_results = sign_results,
            sign_summary_plot = NULL
        )))
    }

    sign_summary_plot <- plot_sign_test_summary(
        sign_results = sign_results,
        primary_label = primary_label,
        secondary_label = secondary_label,
        effect_name = effect_name
    )

    if (!is.null(out_dir)) {
        ggplot2::ggsave(
            filename = file.path(
                out_dir,
                sprintf(
                    "%s_%s_sign_test_summary.svg",
                    output_prefix,
                    sanitize_filename(effect_name)
                )
            ),
            plot = sign_summary_plot,
            device = "svg",
            width = 8,
            height = max(4, 0.35 * nrow(sign_results) + 2)
        )

        write.table(
            sign_results,
            file = file.path(
                out_dir,
                sprintf(
                    "%s_%s_sign_test_summary.tsv",
                    output_prefix,
                    sanitize_filename(effect_name)
                )
            ),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
        )
    }

    invisible(list(
        per_cell_type = results,
        sign_results = sign_results,
        sign_summary_plot = sign_summary_plot
    ))
}


plot_de_primary_secondary_cell_type <- function(
        cell_type,
        primary_file,
        secondary_file,
        primary_dataset = "bican",
        secondary_dataset = "snap",
        primary_label = NULL,
        secondary_label = NULL,
        effect_name = "age effects",
        alpha = 0.05) {

    if (is.null(primary_label)) {
        primary_label <- make_dataset_label(primary_dataset)
    }

    if (is.null(secondary_label)) {
        secondary_label <- make_dataset_label(secondary_dataset)
    }

    df <- make_de_primary_secondary_df(
        primary_file = primary_file,
        secondary_file = secondary_file,
        alpha = alpha
    )

    if (nrow(df) == 0) {
        return(NULL)
    }

    sign_res <- sign_test_fraction(
        x = df$secondary_logFC,
        y = df$primary_logFC
    )

    cor_val <- stats::cor(
        df$secondary_logFC,
        df$primary_logFC,
        use = "complete.obs"
    )

    lims <- range(c(df$secondary_logFC, df$primary_logFC), na.rm = TRUE)

    subtitle <- sprintf(
        "%s-significant genes: %.1f%% same sign in %s (%d/%d); cor=%.2f",
        primary_label,
        100 * sign_res$frac,
        secondary_label,
        sign_res$k,
        sign_res$n,
        cor_val
    )

    p <- ggplot(df, aes(x = secondary_logFC, y = primary_logFC)) +
        geom_point(na.rm = TRUE, size = 1.2, alpha = 0.6) +
        geom_abline(intercept = 0, slope = 1, color = "black") +
        coord_cartesian(xlim = lims, ylim = lims) +
        labs(
            title = cell_type,
            subtitle = subtitle,
            x = sprintf("%s logFC", secondary_label),
            y = sprintf("%s logFC", primary_label)
        ) +
        theme_bw(base_size = 12) +
        theme(
            plot.title = element_text(size = 13, face = "bold"),
            plot.subtitle = element_text(size = 11),
            axis.title = element_text(size = 13),
            axis.text = element_text(size = 11)
        )

    list(
        plot = p,
        sign_test = data.frame(
            cell_type = cell_type,
            primary_dataset = primary_dataset,
            secondary_dataset = secondary_dataset,
            primary_label = primary_label,
            secondary_label = secondary_label,
            effect_name = effect_name,
            k = sign_res$k,
            n = sign_res$n,
            fraction = sign_res$frac,
            percent = 100 * sign_res$frac,
            p_value = sign_res$p_value,
            cor = cor_val,
            stringsAsFactors = FALSE
        ),
        data = df
    )
}


## ------------------------------------------------------------------
## Plotting helpers
## ------------------------------------------------------------------

plot_sign_test_summary <- function(
        sign_results,
        primary_label = "BICAN",
        secondary_label = "SNAP",
        effect_name = "age effects") {

    sign_results <- sign_results[order(sign_results$percent), ]

    sign_results$cell_type <- factor(
        sign_results$cell_type,
        levels = sign_results$cell_type
    )

    title <- sprintf(
        "%s-significant %s: sign concordance in %s",
        primary_label,
        effect_name,
        secondary_label
    )

    ggplot(sign_results, aes(x = cell_type, y = percent)) +
        geom_col(color = "black") +
        geom_text(
            aes(label = paste0(k, "/", n)),
            hjust = -0.1,
            size = 5
        ) +
        coord_flip(ylim = c(0, 105)) +
        labs(
            x = NULL,
            y = "% same sign",
            title = title
        ) +
        theme_bw(base_size = 12)
}


## ------------------------------------------------------------------
## Data helpers
## ------------------------------------------------------------------

make_de_primary_secondary_df <- function(primary_file, secondary_file, alpha = 0.05) {
    primary <- read_de_result(primary_file)
    secondary <- read_de_result(secondary_file)

    common_genes <- intersect(rownames(primary), rownames(secondary))

    genes <- rownames(primary)[primary$adj.P.Val < alpha]
    genes <- intersect(genes, common_genes)

    data.frame(
        gene = genes,
        secondary_logFC = secondary[genes, "logFC"],
        secondary_adjP = secondary[genes, "adj.P.Val"],
        primary_logFC = primary[genes, "logFC"],
        primary_adjP = primary[genes, "adj.P.Val"],
        stringsAsFactors = FALSE
    )
}


read_de_result <- function(file) {
    x <- read.table(
        file,
        header = TRUE,
        stringsAsFactors = FALSE,
        sep = "\t",
        check.names = FALSE
    )

    required_cols <- c("logFC", "adj.P.Val")
    missing_cols <- setdiff(required_cols, colnames(x))

    if (length(missing_cols) > 0) {
        stop(sprintf(
            "Missing required columns in %s: %s",
            file,
            paste(missing_cols, collapse = ", ")
        ))
    }

    if (is.null(rownames(x)) || any(rownames(x) == seq_len(nrow(x)))) {
        warning(
            "Input file may not have gene IDs as row names: ",
            file,
            call. = FALSE
        )
    }

    x
}


validate_manifest_files <- function(manifest, file_cols) {
    missing_cols <- setdiff(file_cols, colnames(manifest))

    if (length(missing_cols) > 0) {
        stop(
            "Missing expected manifest columns: ",
            paste(missing_cols, collapse = ", "),
            call. = FALSE
        )
    }

    missing_files <- lapply(file_cols, function(col) {
        files <- manifest[[col]]
        missing_idx <- which(is.na(files) | !file.exists(files))

        if (length(missing_idx) == 0) {
            return(NULL)
        }

        data.frame(
            row = missing_idx,
            cell_type = manifest$cell_type[missing_idx],
            column = col,
            file = files[missing_idx],
            stringsAsFactors = FALSE
        )
    })

    missing_files <- do.call(rbind, missing_files)

    if (!is.null(missing_files) && nrow(missing_files) > 0) {
        print(missing_files, row.names = FALSE)

        stop(
            "Missing files found in manifest: ",
            nrow(missing_files),
            call. = FALSE
        )
    }

    logger::log_info("All manifest files found")
    invisible(TRUE)
}


## ------------------------------------------------------------------
## Utility helpers
## ------------------------------------------------------------------

sanitize_filename <- function(x) {
    gsub("[^A-Za-z0-9_.-]+", "_", x)
}


make_dataset_label <- function(x) {
    toupper(gsub("_", " ", x))
}


make_output_prefix <- function(primary_dataset, secondary_dataset) {
    sprintf(
        "%s_significant_validated_in_%s",
        sanitize_filename(primary_dataset),
        sanitize_filename(secondary_dataset)
    )
}


## ------------------------------------------------------------------
## Statistics helpers
## ------------------------------------------------------------------

sign_test_fraction <- function(x, y) {
    ok <- is.finite(x) & is.finite(y) & x != 0 & y != 0

    n <- sum(ok)

    if (n == 0) {
        return(list(
            frac = NA_real_,
            n = 0,
            k = 0,
            p_value = NA_real_
        ))
    }

    k <- sum(sign(x[ok]) == sign(y[ok]))

    bt <- stats::binom.test(
        x = k,
        n = n,
        p = 0.5,
        alternative = "greater"
    )

    list(
        frac = k / n,
        n = n,
        k = k,
        p_value = bt$p.value
    )
}
