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
#' @inheritParams find_marker_genes
select_marker_genes_from_files<-function (data_dir, data_name, cellTypeListFile=NULL, markerGeneReportFile,
                                          high_expression_threshold=10, fold_change=10,
                                          method = c("one_vs_all", "target_vs_outgroup")) {

    #validate the outputs are writable
    validate_writable_file(markerGeneReportFile)
    #validate_writable_file(outPDFFile)

    #read in the data, apply cell type filter if provided
    dge=prepare_data_for_marker_genes (data_dir, data_name, cellTypeListFile, donor_col = "donor")

    #perform the test, write the output.
    marker_gene_df<-select_marker_genes (dge, annotation_column = "cell_type",
                                         high_expression_threshold=high_expression_threshold,
                                         fold_change=fold_change, method=method,
                                         markerGeneReportFile=markerGeneReportFile)

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
#' @inheritParams find_marker_genes
select_marker_genes <- function(dgeList,
                                annotation_column = "cell_type",
                                high_expression_threshold = 10,
                                fold_change = 10,
                                method = c("one_vs_all", "target_vs_outgroup"),
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
            fold_change = fold_change,
            method = method
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



#' Identify marker genes for a target annotation using donor-level comparisons
#'
#' This function identifies genes that robustly distinguish a specified
#' `targetAnnotation` from other annotation values in a donor-collapsed
#' `edgeR::DGEList`. Expression is evaluated per donor using CPM, and genes are
#' ranked based on the number of donors that satisfy both a minimum expression
#' threshold in the target group and a minimum fold-change relative to a
#' comparator group.
#'
#' Two comparison strategies are supported via the `method` argument:
#'
#' \itemize{
#'   \item \strong{`"one_vs_all"`}: For each non-target annotation value, the
#'   target group is compared pairwise against that annotation within donors.
#'   For each gene, the "worst" comparison (the one with the greatest number of
#'   failing donors, breaking ties by lowest median fold-change) is selected.
#'   This emphasizes genes that remain well separated even from their closest
#'   competing annotation.
#'
#'   \item \strong{`"target_vs_outgroup"`}: For each donor, raw counts from all
#'   non-target annotations are summed to form a single aggregated outgroup.
#'   A temporary DGEList is constructed with per-donor ingroup and outgroup
#'   counts, CPM is recomputed, and a single comparison (ingroup vs aggregated
#'   outgroup) is performed. This corresponds to a conventional in-group versus
#'   pooled out-group marker analysis.
#' }
#'
#' In both methods, genes are summarized by the median fold-change across
#' donors and the number of donors passing expression and fold-change
#' thresholds. Results are sorted by decreasing number of passing donors and
#' decreasing median fold-change.
#'
#' @param dge A donor-collapsed \code{edgeR::DGEList} with one column per
#'   donor/annotation combination. The \code{samples} slot must contain
#'   columns for \code{donor} and \code{annotation_column}.
#' @param targetAnnotation A single annotation value identifying the target
#'   group.
#' @param annotation_column Character scalar giving the column in
#'   \code{dge$samples} that defines annotation labels. Default is
#'   \code{"cell_type"}.
#' @param high_expression_threshold Minimum CPM required in the target group
#'   (per donor) to count as expressed. Default is 10.
#' @param fold_change Minimum fold-change (target / comparator) required per
#'   donor to count as passing. Default is 10.
#' @param method Character scalar specifying the comparison strategy. Must be
#'   one of \code{"one_vs_all"} or \code{"target_vs_outgroup"}.
#'
#' @return A data.frame with one row per gene containing:
#'   \itemize{
#'     \item \code{gene}: Gene identifier.
#'     \item \code{annotation_column}: Name of the annotation column used.
#'     \item \code{target_annotation}: The target annotation value.
#'     \item \code{closest_annotation}: The most challenging comparator
#'       (for \code{"one_vs_all"}) or \code{"outgroup"} (for
#'       \code{"target_vs_outgroup"}).
#'     \item \code{median_cpm_target}: Median CPM across donors in the target
#'       group.
#'     \item \code{median_fold_change_across_donors}: Median fold-change across
#'       donors.
#'     \item \code{n_donors_expr}: Number of donors with target CPM above the
#'       expression threshold.
#'     \item \code{n_donors_fc_pass}: Number of donors exceeding the
#'       fold-change threshold.
#'     \item \code{n_donors_pass}: Number of donors satisfying both criteria.
#'     \item \code{n_donors_fail}: Number of donors failing one or both
#'       criteria.
#'   }
#'
#' @export
find_marker_genes <- function(dge,
                              targetAnnotation,
                              annotation_column = "cell_type",
                              high_expression_threshold = 10,
                              fold_change = 10,
                              method = c("one_vs_all", "target_vs_outgroup")) {
    method <- match.arg(method)

    .validate_find_marker_genes_inputs(
        dge,
        targetAnnotation = targetAnnotation,
        annotation_column = annotation_column,
        high_expression_threshold = high_expression_threshold,
        fold_change = fold_change
    )

    if (method == "one_vs_all") {
        return(.find_marker_genes_one_vs_all(
            dge = dge,
            targetAnnotation = targetAnnotation,
            annotation_column = annotation_column,
            high_expression_threshold = high_expression_threshold,
            fold_change = fold_change
        ))
    }

    if (method == "target_vs_outgroup") {
        return(.find_marker_genes_target_vs_outgroup(
            dge = dge,
            targetAnnotation = targetAnnotation,
            annotation_column = annotation_column,
            high_expression_threshold = high_expression_threshold,
            fold_change = fold_change
        ))
    }

    stop("Unknown method.", call. = FALSE)
}

#############################
# MARKER GENE FINDER MAIN METHODS
##############################
.find_marker_genes_one_vs_all <- function(dge,
                                          targetAnnotation,
                                          annotation_column,
                                          high_expression_threshold,
                                          fold_change) {
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

        # Kept for compatibility if donor is duplicated within a ct (should not happen in your intended input).
        donors_other <- as.character(samp$donor[idx_other])
        idx_other_by_donor <- .donor_index(donors_other)
        cpm_other <- cpm_mat[, idx_other[idx_other_by_donor], drop = FALSE]

        common_donors <- intersect(colnames(cpm_target), colnames(cpm_other))
        if (length(common_donors) == 0) next

        tmat <- cpm_target[, common_donors, drop = FALSE]
        omat <- cpm_other[, common_donors, drop = FALSE]

        per_ct[[ct]] <- .summarize_marker_comparison(
            tmat = tmat,
            omat = omat,
            high_expression_threshold = high_expression_threshold,
            fold_change = fold_change,
            other_annotation = ct
        )
    }

    all_ct <- do.call(rbind, per_ct)
    if (is.null(all_ct) || nrow(all_ct) == 0) {
        stop("No target/other comparisons had overlapping donors.", call. = FALSE)
    }

    # Worst-case per gene: maximize failing donors, tie-break by lowest median fold change.
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

.find_marker_genes_target_vs_outgroup <- function(dge,
                                                  targetAnnotation,
                                                  annotation_column,
                                                  high_expression_threshold,
                                                  fold_change) {
    dge_tmp <- .build_target_vs_outgroup_dge_tmp(
        dge = dge,
        targetAnnotation = targetAnnotation,
        annotation_column = annotation_column
    )

    mats <- .split_target_vs_outgroup_mats_from_dge_tmp(
        dge_tmp = dge_tmp,
        annotation_column = annotation_column
    )

    cmp <- .summarize_marker_comparison(
        tmat = mats$tmat,
        omat = mats$omat,
        high_expression_threshold = high_expression_threshold,
        fold_change = fold_change,
        other_annotation = "outgroup"
    )

    median_cpm_target <- apply(mats$tmat, 1, stats::median, na.rm = TRUE)

    df <- data.frame(
        gene = cmp$gene,
        annotation_column = annotation_column,
        target_annotation = targetAnnotation,
        closest_annotation = "outgroup",
        median_cpm_target = as.numeric(median_cpm_target[match(cmp$gene, names(median_cpm_target))]),
        median_fold_change_across_donors = as.numeric(cmp$median_fc),
        n_donors_expr = as.integer(cmp$n_donors_expr),
        n_donors_fc_pass = as.integer(cmp$n_donors_fc_pass),
        n_donors_pass = as.integer(cmp$n_donors_pass),
        n_donors_fail = as.integer(cmp$n_donors_fail),
        stringsAsFactors = FALSE
    )

    df <- df[
        with(df, order(-n_donors_pass, -median_fold_change_across_donors)),
    ]

    df
}

###########################
# MARKER GENE FINDER HELPER FUNCTIONS
###########################

# Helper: summarize a single target/comparator comparison given prepared CPM matrices.
# tmat and omat must be genes x donors with identical colnames (donors) and same ordering.
.summarize_marker_comparison <- function(tmat,
                                         omat,
                                         high_expression_threshold,
                                         fold_change,
                                         other_annotation) {
    if (ncol(tmat) != ncol(omat) || nrow(tmat) != nrow(omat)) {
        stop("tmat and omat must have the same dimensions.", call. = FALSE)
    }
    if (!identical(colnames(tmat), colnames(omat))) {
        stop("tmat and omat must have identical column names in the same order.", call. = FALSE)
    }
    if (!identical(rownames(tmat), rownames(omat))) {
        stop("tmat and omat must have identical row names in the same order.", call. = FALSE)
    }

    ratio <- tmat / omat
    ratio[omat == 0 & tmat > 0] <- Inf
    ratio[omat == 0 & tmat == 0] <- NA_real_

    expr_pass <- tmat > high_expression_threshold
    fc_pass <- ratio > fold_change
    pass <- expr_pass & fc_pass
    fail <- !pass

    data.frame(
        gene = rownames(tmat),
        other_annotation = other_annotation,
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

.split_target_vs_outgroup_mats_from_dge_tmp <- function(dge_tmp,
                                                        annotation_column) {
    if (is.null(dge_tmp$samples) || !(annotation_column %in% colnames(dge_tmp$samples))) {
        stop("Temporary DGEList is missing annotation_column in samples.", call. = FALSE)
    }
    if (!("donor" %in% colnames(dge_tmp$samples))) {
        stop("Temporary DGEList is missing 'donor' in samples.", call. = FALSE)
    }

    ann2 <- as.character(dge_tmp$samples[[annotation_column]])
    idx_in <- which(ann2 == "ingroup")
    idx_out <- which(ann2 == "outgroup")

    if (length(idx_in) == 0 || length(idx_out) == 0) {
        stop("Temporary DGEList must contain both ingroup and outgroup samples.", call. = FALSE)
    }
    if (length(idx_in) != length(idx_out)) {
        stop("Temporary DGEList must have the same number of ingroup and outgroup samples.", call. = FALSE)
    }

    cpm_tmp <- edgeR::cpm(dge_tmp, log = FALSE)

    tmat <- cpm_tmp[, idx_in, drop = FALSE]
    omat <- cpm_tmp[, idx_out, drop = FALSE]

    colnames(tmat) <- as.character(dge_tmp$samples$donor[idx_in])
    colnames(omat) <- as.character(dge_tmp$samples$donor[idx_out])

    common_donors <- intersect(colnames(tmat), colnames(omat))
    if (length(common_donors) == 0) {
        stop("No donors had both ingroup and outgroup after construction.", call. = FALSE)
    }

    tmat <- tmat[, common_donors, drop = FALSE]
    omat <- omat[, common_donors, drop = FALSE]

    list(
        tmat = tmat,
        omat = omat
    )
}

.build_target_vs_outgroup_dge_tmp <- function(dge,
                                              targetAnnotation,
                                              annotation_column) {
    samp <- dge$samples
    ann <- as.character(samp[[annotation_column]])
    donors <- as.character(samp$donor)

    counts <- dge$counts
    genes <- rownames(counts)
    if (is.null(genes)) {
        stop("dge$counts must have rownames (gene identifiers).", call. = FALSE)
    }

    elig <- .get_eligible_donors_target_vs_outgroup(
        ann = ann,
        donors = donors,
        targetAnnotation = targetAnnotation
    )

    keep_donors <- elig$keep_donors
    target_col_by_donor <- elig$target_col_by_donor

    counts_target <- counts[, target_col_by_donor, drop = FALSE]
    colnames(counts_target) <- keep_donors

    counts_outgroup <- matrix(0, nrow = nrow(counts), ncol = length(keep_donors))
    rownames(counts_outgroup) <- genes
    colnames(counts_outgroup) <- keep_donors

    for (i in seq_along(keep_donors)) {
        d <- keep_donors[i]
        idx_d <- which(donors == d)
        idx_o <- idx_d[ann[idx_d] != targetAnnotation]
        counts_outgroup[, i] <- rowSums(counts[, idx_o, drop = FALSE])
    }

    tmp_counts <- cbind(counts_target, counts_outgroup)
    tmp_colnames <- c(
        paste(keep_donors, "ingroup", sep = "_"),
        paste(keep_donors, "outgroup", sep = "_")
    )
    colnames(tmp_counts) <- tmp_colnames

    tmp_samples <- data.frame(
        donor = rep(keep_donors, times = 2),
        stringsAsFactors = FALSE
    )
    tmp_samples[[annotation_column]] <- rep(c("ingroup", "outgroup"), each = length(keep_donors))

    edgeR::DGEList(counts = tmp_counts, samples = tmp_samples)
}

.get_eligible_donors_target_vs_outgroup <- function(ann,
                                                    donors,
                                                    targetAnnotation) {
    donors <- as.character(donors)
    ann <- as.character(ann)

    donor_levels <- unique(donors)

    keep_donors <- character(0)
    target_col_by_donor <- integer(0)

    for (d in donor_levels) {
        idx_d <- which(donors == d)
        if (length(idx_d) == 0) next

        idx_t <- idx_d[ann[idx_d] == targetAnnotation]
        idx_o <- idx_d[ann[idx_d] != targetAnnotation]

        if (length(idx_t) != 1) next
        if (length(idx_o) < 1) next

        keep_donors <- c(keep_donors, d)
        target_col_by_donor <- c(target_col_by_donor, idx_t)
    }

    if (length(keep_donors) == 0) {
        stop("No donors had both a target sample and at least one outgroup sample.", call. = FALSE)
    }

    list(
        keep_donors = keep_donors,
        target_col_by_donor = target_col_by_donor
    )
}

# This returns a named integer vector mapping donor -> column index.
# (Kept for compatibility; in donor-collapsed inputs this is typically a no-op.)
.donor_index <- function(donors) {
    donors <- as.character(donors)
    idx <- match(unique(donors), donors)
    names(idx) <- donors[idx]
    idx
}

######################
# OTHER HELPERS
######################

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
        stop("dge$samples must contain annotation_column.", call. = FALSE)
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
    if (!(targetAnnotation %in% ann)) {
        stop("targetAnnotation not present in annotation_column.", call. = FALSE)
    }
}


