#select marker genes for each cell type

# source ("R/differential_expression.R")
# source ("R/donor_age_prediction.R")
#
# data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/metacells/LEVEL_3"
# data_name="donor_rxn_DGEList"
# cellTypeListFile="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/metadata/cell_types_for_marker_gene_discovery_level_1.txt"
# markerGeneReportFile="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/cell_type_marker_genes/marker_gene_selection_level_1.txt"
# high_expression_threshold=10;fold_change=10


#' Select marker genes for each cell type (file interface)
#'
#' Wrapper around \code{\link{select_marker_genes}} that loads a DGEList from disk
#' (via \code{prepare_data_for_marker_genes()}) and writes a tab-delimited report
#' (via \code{write_marker_genes()}).
#'
#' Marker genes are identified per target cell type by comparing donor-matched CPM
#' values against every other cell type, then selecting the worst-case comparator
#' per gene (maximizing donor failures; ties broken by lowest median fold change).
#'
#' This implementation is tied closely to the BICAN data preparation pipeline for
#' comparing cell types to find marker genes.
#' @param data_dir Directory containing serialized inputs for
#'   \code{prepare_data_for_marker_genes()}.
#' @param data_name Basename/key used by \code{prepare_data_for_marker_genes()} to
#'   locate inputs.
#' @param cellTypeListFile Optional file path containing a list of cell types to
#'   include (passed to \code{prepare_data_for_marker_genes()}). If \code{NULL},
#'   uses all cell types present.
#' @param markerGeneReportFile Output TSV path. A one-line header beginning with
#'   \code{#} is written first and records \code{high_expression_threshold} and
#'   \code{fold_change}.
#' @param high_expression_threshold CPM threshold used to define "expressed in
#'   target" for donor-level pass/fail metrics.
#' @param fold_change Minimum per-donor fold change (target / other) required to
#'   count a donor as passing the fold-change criterion.
#'
#' @return Invisibly returns the marker gene data frame produced by
#'   \code{\link{select_marker_genes}} (also written to \code{markerGeneReportFile}).
#'
#' @export
select_marker_genes_from_files<-function (data_dir, data_name, cellTypeListFile=NULL, markerGeneReportFile,
                                          high_expression_threshold=10, fold_change=10) {

    #validate the outputs are writable
    validate_writable_file(markerGeneReportFile)
    #validate_writable_file(outPDFFile)

    #read in the data, apply cell type filter if provided
    dge=prepare_data_for_marker_genes (data_dir, data_name, cellTypeListFile, donor_col = "donor")

    #perform the test, write the output.
    marker_gene_df<-select_marker_genes (dge, annotation_column = "cell_type",
                                         high_expression_threshold=high_expression_threshold,
                                         fold_change=fold_change, markerGeneReportFile=markerGeneReportFile)
}

#' Select marker genes by annotation group (DGEList interface)
#'
#' Identifies marker genes for each annotation group in a donor-collapsed
#' \code{edgeR} \code{DGEList}. For a given target annotation value and gene,
#' CPM values are compared against each other annotation value using only donors
#' present in both groups.
#'
#' Donor-level pass/fail metrics are computed using:
#' \itemize{
#'   \item target CPM > \code{high_expression_threshold}
#'   \item (target CPM / other CPM) > \code{fold_change}
#' }
#'
#' For each gene, the worst-case comparator annotation is selected by maximizing
#' the number of failing donors; ties are broken by the lowest median fold change
#' across donors.
#'
#' If \code{markerGeneReportFile} is provided, the results are written as a
#' tab-delimited file with a leading comment line beginning with \code{#} that
#' records \code{high_expression_threshold} and \code{fold_change}.
#'
#' @param dgeList An \code{edgeR::DGEList} with one observation per donor and
#'   annotation group (or an object accepted by
#'   \code{simplify_dge_for_marker_genes()}), and with
#'   \code{dgeList$samples[[annotation_column]]} and
#'   \code{dgeList$samples$donor} present.
#' @param annotation_column Character scalar giving the name of the column in
#'   \code{dgeList$samples} that defines annotation groups (default:
#'   \code{"cell_type"}).
#' @param high_expression_threshold CPM threshold used to define "expressed in
#'   target" for donor-level metrics.
#' @param fold_change Minimum per-donor fold change (target / other) required to
#'   count a donor as passing the fold-change criterion.
#' @param markerGeneReportFile Optional output TSV path. If \code{NULL}, no file
#'   is written.
#'
#' @return A data frame of marker-gene candidates across all target annotation
#'   values. The returned data frame includes at least:
#'   \itemize{
#'     \item \code{gene}
#'     \item \code{annotation_column}
#'     \item \code{target_annotation}
#'     \item \code{closest_annotation} (worst-case comparator)
#'     \item \code{median_cpm_target}
#'     \item \code{median_fold_change_across_donors}
#'     \item \code{n_donors_expr}
#'     \item \code{n_donors_fc_pass}
#'     \item \code{n_donors_pass}
#'     \item \code{n_donors_fail}
#'   }
#'
#' @export
select_marker_genes <- function(dgeList,
                                annotation_column = "cell_type",
                                high_expression_threshold = 10,
                                fold_change = 10,
                                markerGeneReportFile = NULL) {
    dge <- simplify_dge_for_marker_genes(dgeList)

    stopifnot(annotation_column %in% colnames(dge$samples))

    outputs_list <- list()

    annotation_values <- sort(unique(as.character(dge$samples[[annotation_column]])))

    lineStr <- strrep("=", 80)
    logger::log_info(lineStr)
    logger::log_info("Starting marker genes selection")
    logger::log_info(paste("Annotation column:", annotation_column))
    logger::log_info(paste("High expression threshold:", high_expression_threshold))
    logger::log_info(paste("Fold change threshold:", fold_change))
    logger::log_info(lineStr)

    for (targetAnnotation in annotation_values) {
        logger::log_info(
            paste("Learning marker genes for", annotation_column, "=", targetAnnotation)
        )

        marker_genes <- find_marker_genes(
            dge,
            targetAnnotation = targetAnnotation,
            annotation_column = annotation_column,
            high_expression_threshold = high_expression_threshold,
            fold_change = fold_change
        )

        outputs_list[[targetAnnotation]] <- marker_genes
    }

    all_marker_genes <- do.call(rbind, outputs_list)

    if (!is.null(markerGeneReportFile)) {
        write_marker_genes(
            all_marker_genes,
            markerGeneReportFile,
            high_expression_threshold,
            fold_change
        )
    }

    all_marker_genes
}



.validate_find_marker_genes_inputs <- function(dge,
                                               targetAnnotation,
                                               annotation_column,
                                               high_expression_threshold,
                                               fold_change) {
    if (is.null(dge) || is.null(dge$counts) || is.null(dge$samples)) {
        stop("dge must be a DGEList with $counts and $samples.", call. = FALSE)
    }
    if (!("donor" %in% colnames(dge$samples))) {
        stop("dge$samples must contain a 'donor' column.", call. = FALSE)
    }
    if (length(annotation_column) != 1 || is.na(annotation_column) || !nzchar(annotation_column)) {
        stop("annotation_column must be a single, non-empty string.", call. = FALSE)
    }
    if (!(annotation_column %in% colnames(dge$samples))) {
        stop("dge$samples must contain annotation_column='", annotation_column, "'.", call. = FALSE)
    }
    if (length(targetAnnotation) != 1 || is.na(targetAnnotation) || !nzchar(targetAnnotation)) {
        stop("targetAnnotation must be a single, non-empty string.", call. = FALSE)
    }
    if (!is.numeric(high_expression_threshold) || length(high_expression_threshold) != 1 ||
        is.na(high_expression_threshold) || high_expression_threshold < 0) {
        stop("high_expression_threshold must be a single non-negative number.", call. = FALSE)
    }
    if (!is.numeric(fold_change) || length(fold_change) != 1 ||
        is.na(fold_change) || fold_change <= 0) {
        stop("fold_change must be a single positive number.", call. = FALSE)
    }

    ann <- as.character(dge$samples[[annotation_column]])
    idx_target <- which(ann == targetAnnotation)
    if (length(idx_target) == 0) {
        stop("No samples found with ", annotation_column, " == ", targetAnnotation, ".", call. = FALSE)
    }

    invisible(TRUE)
}

# dge is already collapsed to one observation per donor/cell_type.
# This returns a named integer vector mapping donor -> column index.
.donor_index <- function(donors) {
    donors <- as.character(donors)
    idx <- match(unique(donors), donors)
    names(idx) <- donors[idx]
    idx
}

find_marker_genes <- function(dge,
                              targetAnnotation,
                              annotation_column = "cell_type",
                              high_expression_threshold = 10,
                              fold_change = 10) {
    .validate_find_marker_genes_inputs(
        dge,
        targetAnnotation = targetAnnotation,
        annotation_column = annotation_column,
        high_expression_threshold = high_expression_threshold,
        fold_change = fold_change
    )

    samp <- dge$samples
    ann <- as.character(samp[[annotation_column]])

    cpm_mat <- edgeR::cpm(dge, log = FALSE)

    idx_target <- which(ann == targetAnnotation)
    cpm_target <- cpm_mat[, idx_target, drop = FALSE]

    median_cpm_target <- apply(cpm_target, 1, stats::median, na.rm = TRUE)

    other_ann_values <- setdiff(unique(ann), targetAnnotation)
    if (length(other_ann_values) == 0) {
        stop("No other annotation values present besides the target annotation.", call. = FALSE)
    }

    per_ct <- vector("list", length(other_ann_values))
    names(per_ct) <- other_ann_values

    for (ct in other_ann_values) {
        idx_other <- which(ann == ct)
        if (length(idx_other) == 0) next

        # Extract CPM matrix for the other annotation value, keeping one column per donor
        donors_other <- as.character(samp$donor[idx_other])
        idx_other_by_donor <- .donor_index(donors_other)
        cpm_other <- cpm_mat[, idx_other[idx_other_by_donor], drop = FALSE]

        # Restrict to donors present in both the target and comparison groups
        common_donors <- intersect(colnames(cpm_target), colnames(cpm_other))
        if (length(common_donors) == 0) next

        tmat <- cpm_target[, common_donors, drop = FALSE]
        omat <- cpm_other[, common_donors, drop = FALSE]

        # Compute per-donor fold change and pass/fail status for each gene
        ratio <- tmat / omat
        ratio[omat == 0 & tmat > 0] <- Inf
        ratio[omat == 0 & tmat == 0] <- NA_real_

        # Donors with sufficient absolute expression in target
        expr_pass <- tmat > high_expression_threshold

        # Donors with sufficient fold change
        fc_pass <- ratio > fold_change

        # Donors passing both criteria
        pass <- expr_pass & fc_pass
        fail <- !pass

        # Summarize per-gene statistics for this target/other comparison
        per_ct[[ct]] <- data.frame(
            gene = rownames(cpm_mat),
            other_annotation = ct,
            median_fc = apply(ratio, 1, function(x) {
                x <- x[!is.na(x)]
                if (length(x) == 0) return(NA_real_)
                xf <- x[is.finite(x)]
                if (length(xf) > 0) stats::median(xf, na.rm = TRUE) else Inf
            }),
            n_donors_expr = apply(expr_pass, 1, function(x) sum(x, na.rm = TRUE)),
            n_donors_fc_pass = apply(fc_pass, 1, function(x) sum(x, na.rm = TRUE)),
            n_donors_pass = apply(pass, 1, function(x) sum(x, na.rm = TRUE)),
            n_donors_fail = apply(fail, 1, function(x) sum(x, na.rm = TRUE)),
            stringsAsFactors = FALSE
        )
    }

    all_ct <- do.call(rbind, per_ct)
    if (is.null(all_ct) || nrow(all_ct) == 0) {
        stop("No target/other comparisons had overlapping donors.", call. = FALSE)
    }

    # Select the worst-case comparison per gene:
    # maximize number of failing donors, tie-break by lowest median fold change
    all_ct <- all_ct[with(all_ct, order(gene, -n_donors_fail, median_fc)), ]
    worst <- all_ct[!duplicated(all_ct$gene), ]

    df <- data.frame(
        gene = worst$gene,
        annotation_column = annotation_column,
        target_annotation = targetAnnotation,
        closest_annotation = worst$other_annotation,
        median_cpm_target = as.numeric(median_cpm_target[match(worst$gene, names(median_cpm_target))]),
        median_fold_change_across_donors = as.numeric(worst$median_fc),
        n_donors_expr = as.integer(worst$n_donors_expr),
        n_donors_fc_pass = as.integer(worst$n_donors_fc_pass),
        n_donors_pass = as.integer(worst$n_donors_pass),
        n_donors_fail = as.integer(worst$n_donors_fail),
        stringsAsFactors = FALSE
    )

    df <- df[
        with(df, order(-n_donors_pass, -median_fold_change_across_donors)),
    ]

    df
}


#simplify the original dge to metacells summarized by donor and cell type.
#also return the cpm matrix.
prepare_data_for_marker_genes<-function (data_dir, data_name, cellTypeListFile=NULL, donor_col = "donor") {
    #load the DGEList and prepare the data
    d=bican.mccarroll.differentialexpression::prepare_data_for_differential_expression(data_dir, data_name, randVars=c(), fixedVars=c())
    dge=d$dge

    if (!is.null(cellTypeListFile))
        dge=filter_dgelist_by_celltype_list(dge, cellTypeListFile)

    return (dge)
}

simplify_dge_for_marker_genes<-function (dge, donor_col = "donor") {
    logger::log_info("Simplifying DGEList to donor by cell type observations for marker gene selection.")
    cell_type_list=unique(dge$samples$cell_type)

    dge_list=list()

    for (cellType in cell_type_list) {
        logger::log_info("Processing cell type: {cellType}")
        dge_cell <- dge[, dge$samples$cell_type == cellType, keep.lib.sizes = TRUE]

        # Collapse multiple samples per donor
        dge_cell <- collapse_by_donor(dge_cell, donor_col = donor_col, keep_cols = c(donor_col, "cell_type"))
        dge_list[[cellType]] <- dge_cell

    }

    dge <- do.call(cbind, dge_list)
    return (dge)
}


plot_marker_gene <- function(gene, target_cell_type, dge, marker_df) {
    if (!requireNamespace("edgeR", quietly = TRUE)) {
        stop("edgeR is required.", call. = FALSE)
    }
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 is required.", call. = FALSE)
    }
    if (is.null(dge) || is.null(dge$counts) || is.null(dge$samples)) {
        stop("dge must be a DGEList with $counts and $samples.", call. = FALSE)
    }
    if (!all(c("cell_type", "donor") %in% colnames(dge$samples))) {
        stop("dge$samples must contain 'cell_type' and 'donor'.", call. = FALSE)
    }

    gene <- as.character(gene)
    target_cell_type <- as.character(target_cell_type)

    hit <- marker_df[marker_df$gene == gene & marker_df$target_cell_type == target_cell_type, , drop = FALSE]
    if (nrow(hit) == 0) {
        stop("No row in marker_df matches gene='", gene, "' and target_cell_type='", target_cell_type, "'.", call. = FALSE)
    }
    if (nrow(hit) > 1) {
        hit <- hit[1, , drop = FALSE]
    }

    other_cell_type <- as.character(hit$closest_cell_type)

    samp <- dge$samples
    idx_t <- which(as.character(samp$cell_type) == target_cell_type)
    idx_o <- which(as.character(samp$cell_type) == other_cell_type)

    if (length(idx_t) == 0) stop("No samples found for target_cell_type='", target_cell_type, "'.", call. = FALSE)
    if (length(idx_o) == 0) stop("No samples found for worst_other_cell_type='", other_cell_type, "'.", call. = FALSE)

    cpm_mat <- edgeR::cpm(dge, log = FALSE)

    if (!(gene %in% rownames(cpm_mat))) {
        stop("Gene '", gene, "' not found in rownames(edgeR::cpm(dge)).", call. = FALSE)
    }

    donors_t <- as.character(samp$donor[idx_t])
    donors_o <- as.character(samp$donor[idx_o])

    # If input is properly collapsed, donors are unique per cell type; the match(unique(.), .) keeps one col per donor.
    idx_t1 <- idx_t[match(unique(donors_t), donors_t)]
    idx_o1 <- idx_o[match(unique(donors_o), donors_o)]

    t_donors <- as.character(samp$donor[idx_t1])
    o_donors <- as.character(samp$donor[idx_o1])

    t_vals <- cpm_mat[gene, idx_t1]
    o_vals <- cpm_mat[gene, idx_o1]
    names(t_vals) <- t_donors
    names(o_vals) <- o_donors

    common_donors <- intersect(names(t_vals), names(o_vals))
    if (length(common_donors) == 0) {
        stop("No overlapping donors between target and worst other cell type for this gene.", call. = FALSE)
    }

    plot_df <- data.frame(
        donor = common_donors,
        cpm_other = as.numeric(o_vals[common_donors]),
        cpm_target = as.numeric(t_vals[common_donors]),
        stringsAsFactors = FALSE
    )

    title_str <- paste0(gene, ": ", target_cell_type, " vs ", other_cell_type)
    subtitle_str <- paste0(
        "n_pass=", as.integer(hit$n_donors_pass),
        ", median_fc=", format(as.numeric(hit$median_fold_change_across_donors), digits = 4, scientific = FALSE)
    )

    # Shared axis limits (log10 scale, CPM values)
    # Log10 axes cannot include zeros; floor CPM at a small positive value for plotting only.
    vals_pos <- c(plot_df$cpm_other, plot_df$cpm_target)
    vals_pos <- vals_pos[is.finite(vals_pos) & vals_pos > 0]

    eps <- if (length(vals_pos) == 0) 1e-3 else min(vals_pos) / 10

    plot_df$cpm_other_plot  <- pmax(plot_df$cpm_other,  eps)
    plot_df$cpm_target_plot <- pmax(plot_df$cpm_target, eps)

    lims <- range(plot_df$cpm_other_plot, plot_df$cpm_target_plot, na.rm = TRUE)

    # Use the target-cell-type samples to get one lib.size per donor (since the data are one row per donor/cell type).
    lib_df <- data.frame(
        donor = as.character(samp$donor[idx_t]),
        lib_size = as.numeric(samp$lib.size[idx_t]),
        stringsAsFactors = FALSE
    )

    plot_df <- merge(plot_df, lib_df, by = "donor", all.x = TRUE, sort = FALSE)
    plot_df$log10_lib_size <- log10(plot_df$lib_size)

    # MAKE R CMD CHECK Happy
    cpm_other_plot <- cpm_target_plot <- log10_lib_size <- NULL

    ggplot2::ggplot(plot_df, ggplot2::aes(x = cpm_other_plot, y = cpm_target_plot, color = log10_lib_size)) +
        ggplot2::geom_point() +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        ggplot2::scale_x_continuous(trans = "log10", limits = lims) +
        ggplot2::scale_y_continuous(trans = "log10", limits = lims) +
        ggplot2::coord_equal() +
        ggplot2::labs(
            title = title_str,
            subtitle = subtitle_str,
            x = paste0(other_cell_type, " CPM (log10 scale)"),
            y = paste0(target_cell_type, " CPM (log10 scale)"),
            color = "log10(lib.size)"
        ) +
        ggplot2::theme_classic()

}


validate_writable_file <- function(x) {
    if (is.null(x)) {
        return(invisible(TRUE))
    }

    dirn <- dirname(x)

    if (!dir.exists(dirn)) {
        stop("Directory does not exist for output file: ", x, call. = FALSE)
    }

    testfile <- file.path(
        dirn,
        paste0(".__test_write__", Sys.getpid(), "_", as.integer(Sys.time()), ".tmp")
    )

    con <- NULL
    on.exit({
        if (!is.null(con)) close(con)
        if (file.exists(testfile)) file.remove(testfile)
    }, add = TRUE)

    tryCatch({
        con <- file(testfile, open = "w")
        writeLines("test", con)
    }, error = function(e) {
        stop("Directory is not writable for output file: ", x, call. = FALSE)
    })

    invisible(TRUE)
}

write_marker_genes <- function(df, file,
                               high_expression_threshold,
                               fold_change) {
    con <- file(file, open = "wt")
    on.exit(close(con), add = TRUE)

    # Write parameter header
    writeLines(
        paste0(
            "# high_expression_threshold=", high_expression_threshold,
            ", fold_change=", fold_change
        ),
        con
    )

    # Write the table (tab-delimited, with column names)
    utils::write.table(
        df,
        file = con,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE
    )

    invisible(TRUE)
}


# adhocSteve<-function (dge) {
#     dge=prepare_data_for_marker_genes (data_dir, data_name,  cellTypeListFile=NULL, donor_col = "donor")
#     dge_sst=dge[,dge$samples$cell_type=="GABA_SST"]
#
#     cpm_mat <- edgeR::cpm(dge_sst, log = FALSE)
#     apoe <- as.numeric(cpm_mat["APOE", ])
#
#     plot_df <- data.frame(
#         APOE_CPM = apoe,
#         sex = factor(
#             dge_sst$samples$imputed_sex,
#             levels = c(1, 2),
#             labels = c("Male", "Female")
#         ),
#         stringsAsFactors = FALSE
#     )
#
#     plot_df <- plot_df[is.finite(plot_df$APOE_CPM) & !is.na(plot_df$sex), , drop = FALSE]
#
#     # Floor zeros for plotting only (do NOT use this for the test)
#     vals_pos <- plot_df$APOE_CPM[plot_df$APOE_CPM > 0]
#     eps <- if (length(vals_pos) == 0) 1e-3 else min(vals_pos) / 10
#
#     plot_df$APOE_CPM_plot <- pmax(plot_df$APOE_CPM, eps)
#
#     wilcox_res <- stats::wilcox.test(APOE_CPM ~ sex, data = plot_df)
#
#     ggplot2::ggplot(plot_df, ggplot2::aes(x = sex, y = APOE_CPM_plot)) +
#         ggplot2::geom_boxplot(outlier.shape = NA) +
#         ggplot2::geom_jitter(width = 0.15, height = 0, alpha = 0.6) +
#         ggplot2::scale_y_continuous(trans = "log10") +
#         ggplot2::labs(
#             title = "APOE expression in GABA_SST",
#             subtitle = paste0(
#                 "Wilcoxon p=",
#                 format(wilcox_res$p.value, digits = 3, scientific = TRUE),
#                 "; n(Male)=", sum(plot_df$sex == "Male"),
#                 ", n(Female)=", sum(plot_df$sex == "Female")
#             ),
#             x = NULL,
#             y = "APOE CPM (log10 scale)"
#         ) +
#         ggplot2::theme_classic()
#
# }
#
# adhocSteveB<-function () {
#     d=bican.mccarroll.differentialexpression::prepare_data_for_differential_expression(data_dir, data_name, randVars=c(), fixedVars=c())
#     dge=d$dge
#
#     #subset to a single region.
#     z=dge[,dge$samples$region=="DFC"]
#     counts=z$counts
#     samp=z$samples[, c("donor", "cell_type")]
#
#     dge=edgeR::DGEList(counts=counts, samples=samp)
# }
