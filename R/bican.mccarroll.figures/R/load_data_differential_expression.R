# in_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results/LEVEL_3/sex_age/cell_type"
# cellTypeListFile="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metadata/mash_cell_type_list_simple.txt"
# file_pattern="age"
# fdrThreshold=0.01


#' Gather differential expression results and construct mash-ready matrices
#'
#' @description
#' Read per-condition (or per-cell-type) differential expression (DE) result tables
#' from `in_dir`, optionally restrict to a provided list of cell types, and
#' construct mash-compatible matrices (`Bhat`, `Shat`, `Tstat`, `FDR`, etc.) using
#' `make_mash_inputs()` with `gene_mode = "union"`.
#'
#' The function then filters the matrices to genes that pass an FDR threshold in
#' at least one condition (row-wise across the `FDR` matrix).
#'
#' Gene x conditions that were not tested (missing in the input DE tables) are
#' flagged by the `missing_mask` output matrix.  Those entries have been filled with
#' neutral values (beta = 0, SE = large value) so they do not influence mash results.
#'
#' @param in_dir Character scalar. Directory containing DE result files.
#' @param file_pattern Character scalar. Pattern used to select DE files within
#'   `in_dir` (passed to `list.files(..., pattern = file_pattern)`). Default `"age"`.
#' @param cellTypeListFile Character scalar or `NULL`. Path to a text file with one
#'   cell type per line used to filter the DE inputs. If provided, entries are
#'   matched as prefixes of file basenames after appending `paste0("_", file_pattern)`.
#'   Default `NULL`.
#' @param fdrThreshold Numeric scalar. Genes are retained if they have
#'   `FDR < fdrThreshold` in at least one condition. Default `0.01`.
#'
#' @return
#' A list of mash input matrices after filtering:
#' \itemize{
#'   \item `Bhat` - numeric matrix [genes x conditions] of effect sizes.
#'   \item `Shat` - numeric matrix [genes x conditions] of standard errors.
#'   \item `Tstat` - numeric matrix [genes x conditions] of t-statistics.
#'   \item `P_val` - numeric matrix [genes x conditions] of p-values (if present in inputs).
#'   \item `FDR` - numeric matrix [genes x conditions] of FDR values.
#'   \item `missing_mask` - logical matrix [genes x conditions] indicating filled/missing entries.
#'   \item `idxPassFDRAny` - integer vector of row indices (in the unfiltered union matrices)
#'         that passed the FDR criterion.
#' }
#'
#' @export
gather_de_results<-function (in_dir, file_pattern="age", cellTypeListFile=NULL, fdrThreshold=0.01) {

    d=parse_de_inputs(in_dir, file_pattern, cellTypeListFile)
    #make mash inputs
    inputs_union<-make_mash_inputs(d, coef_col = "logFC", t_col = "t", fdr_col = "adj.P.Val", gene_mode="union")
    #Filter results to genes that pass an FDR threshold in at least one condition
    idxPassFDRAny<-which(apply(inputs_union$FDR, 1, function (x) any(x<fdrThreshold)))
    #how many genes pass an FDR < 0.05 in any condition?
    nGenesFDR=length(idxPassFDRAny)
    logger::log_info("Number of genes with at least one condition that passes [", fdrThreshold, "] FDR [", nGenesFDR, "]")

    #Filter all the gene matrixes to only those that pass the FDR threshold in at least one condition
    is_gene_matrix <- vapply(inputs_union, function(x) is.matrix(x) && nrow(x) == nrow(inputs_union$FDR), logical(1))

    inputs_union[is_gene_matrix] <- lapply(
        inputs_union[is_gene_matrix],
        function(m) m[idxPassFDRAny, , drop = FALSE]
    )

    return (inputs_union)
}

#' Build Bhat/Shat/Tstat/FDR matrices from per-cell-type DE tables
#'
#' @description
#' Construct matrices for `mashr` from a list of differential expression tables
#' (one per cell type or condition). You can choose to include the **union** of all
#' genes seen across tables (filling missing entries with beta = 0 and a very large
#' SE so they have no influence on mash), or restrict to the **intersection**
#' (genes present in every table).
#'
#' For each table, the function reads a coefficient column (for example, `logFC`),
#' a t-statistic column (for example, `t`), and an FDR column (for example, `adj.P.Val`), and
#' computes `Shat = abs(beta / t)` (so `t = beta / se` as in limma).
#'
#' @param lst Named `list` of data frames (one per cell type). Each must have
#'   row names as gene IDs and the columns specified by `coef_col`, `t_col`,
#'   and `fdr_col`.
#' @param coef_col Character scalar. Column name for the effect size (beta),
#'   for example `"logFC"`. Default `"logFC"`.
#' @param t_col Character scalar. Column name for the t-statistic, for example `"t"`.
#'   Default `"t"`.
#' @param p_col Character scalar. Column name for the p-value, for example `"P.Value"`.
#' @param fdr_col Character scalar. Column name for the FDR or q-value,
#'   for example `"adj.P.Val"`. Default `"adj.P.Val"`.
#' @param gene_mode One of `c("union", "intersect")`. Use `"union"` to include
#'   all genes observed in any table (missing entries are replaced with
#'   `beta = 0` and `SE = big_se`). Use `"intersect"` to keep only genes present
#'   in every table. Default `"union"`.
#' @param big_se Numeric scalar. Standard error to use when filling missing or
#'   invalid entries. A very large value ensures such entries contribute ~0
#'   information to mash. Default `1e6`.
#' @param fill_missing_fdr Numeric scalar used to fill missing FDR entries in
#'   the returned `FDR` matrix (for plotting or logic only; not used by mash).
#'   Default `NA_real_`.
#'
#' @return A list with components:
#' \itemize{
#'   \item `Bhat` - matrix [genes x cell types] of betas.
#'   \item `Shat` - matrix [genes x cell types] of standard errors (`abs(beta / t)`).
#'   \item `Tstat` - matrix [genes x cell types] of moderated t-statistics.
#'   \item `FDR`  - matrix [genes x cell types] of FDR values.
#'   \item `missing_mask` - logical matrix [genes x cell types] indicating where
#'         input was missing or invalid and replaced by `beta = 0, SE = big_se`.
#' }
#'
#' @details
#' - With `gene_mode = "union"`, the function guarantees rectangular matrices by
#'   inserting neutral placeholders for missing cells (0, `big_se`). This allows
#'   mash to ingest more genes without biasing fits.
#' - With `gene_mode = "intersect"`, all retained genes have entries in every
#'   column (aside from pathological `t = 0`; those are still treated as missing
#'   and assigned `SE = big_se`).
#'
#' @export
make_mash_inputs <- function(
        lst,
        coef_col = "logFC",
        t_col    = "t",
        p_col = "P.Value",
        fdr_col  = "adj.P.Val",
        gene_mode = c("union", "intersect"),
        big_se   = 1e6,
        fill_missing_fdr = NA_real_
) {
    gene_mode <- match.arg(gene_mode)

    stopifnot(is.list(lst), length(lst) >= 2L)
    if (is.null(names(lst)) || any(names(lst) == "")) {
        names(lst) <- paste0("CT", seq_along(lst))
    }

    # determine gene set
    gene_lists <- lapply(lst, rownames)
    if (gene_mode == "union") {
        genes <- sort(unique(unlist(gene_lists, use.names = FALSE)))
    } else {
        genes <- Reduce(intersect, gene_lists)
        genes <- sort(unique(genes))
    }
    G <- length(genes); C <- length(lst)
    if (G == 0L) stop("No gene names found after applying gene_mode = '", gene_mode, "'.")

    # preallocate outputs
    Bhat <- matrix(NA_real_, G, C, dimnames = list(genes, names(lst)))
    Shat <- matrix(NA_real_, G, C, dimnames = list(genes, names(lst)))
    Tstat<- matrix(NA_real_, G, C, dimnames = list(genes, names(lst)))
    FDR  <- matrix(NA_real_, G, C, dimnames = list(genes, names(lst)))
    P_val<- matrix(NA_real_, G, C, dimnames = list(genes, names(lst)))

    # fill matrices column by column
    for (j in seq_along(lst)) {
        df <- lst[[j]]
        required <- c(coef_col, t_col, fdr_col)
        if (!all(required %in% colnames(df))) {
            stop(paste("Missing required columns in lst[['", names(lst)[j], "']]: ",
                       paste(setdiff(required, colnames(df)), collapse = ", ")))
        }

        g <- intersect(genes, rownames(df))
        if (!length(g)) next

        beta <- as.numeric(df[g, coef_col])
        tt   <- as.numeric(df[g, t_col])
        fdr  <- as.numeric(df[g, fdr_col])
        pval <- as.numeric (df[g, p_col])
        se <- abs(beta / tt)  # SE from beta/t

        Bhat[g, j]  <- beta
        Shat[g, j]  <- se
        P_val[g,j] <- pval
        Tstat[g, j] <- tt
        FDR [g, j]  <- fdr
    }

    # handle missing/invalid
    missing_mask <- !is.finite(Bhat) | !is.finite(Shat) | (Shat <= 0)
    if (any(missing_mask)) {
        Bhat[missing_mask]  <- 0
        Shat[missing_mask]  <- big_se
        Tstat[missing_mask] <- 0
    }

    nas_fdr <- !is.finite(FDR)
    if (any(nas_fdr)) FDR[nas_fdr] <- fill_missing_fdr

    list(Bhat = Bhat, Shat = Shat, Tstat = Tstat, P_val=P_val,
         FDR = FDR, missing_mask = missing_mask)
}



#' Drop common suffix from file names
#' Useful to infer the data set names from differential expression files
#' @param files A list of file names to process
#' @return A vector with the same length as the input files containing the name of the file without
#' the common suffix.
drop_common_suffix <- function(files) {
    # Reverse each filename so suffix becomes a prefix
    rev_files <- sapply(files, function(x) paste0(rev(strsplit(x, NULL)[[1]]), collapse = ""))

    # Find common prefix of reversed strings
    common_rev_prefix <- Reduce(function(a, b) {
        i <- 1
        while (i <= nchar(a) && i <= nchar(b) && substr(a, i, i) == substr(b, i, i)) {
            i <- i + 1
        }
        substr(a, 1, i - 1)
    }, rev_files)

    # Reverse back to get the suffix
    common_suffix <- paste0(rev(strsplit(common_rev_prefix, NULL)[[1]]), collapse = "")

    # Drop suffix from each filename
    sub(paste0(common_suffix, "$"), "", files)
}

#' Remove shared tokens across strings
#'
#' Split each string by a delimiter and drop any tokens that appear in **all**
#' inputs. Useful for shortening condition names with common prefixes/suffixes.
#'
#' @param x Character vector of strings to simplify.
#' @param delim Single-character delimiter used to split and rejoin tokens
#'   (default `"_"`).
#'
#' @return Character vector with common tokens removed. The vector carries an
#'   attribute `"common_tokens"` listing the tokens removed.
#' @noRd
remove_common_tokens <- function(x, delim = "_") {
    stopifnot(is.character(x))
    toks <- strsplit(x, delim, fixed = TRUE)
    common <- Reduce(intersect, lapply(toks, unique))
    if (length(common) == 0L) return(x)
    out <- vapply(
        toks,
        function(v) paste(v[!v %in% common], collapse = delim),
        character(1)
    )
    attr(out, "common_tokens") <- common
    out
}

find_prefix_matches <- function(full_names, partial_names) {
    # indexes in full_names that match any prefix in partial_names
    idx <- which(sapply(full_names, function(x) any(startsWith(x, partial_names))))

    # prefixes in partial_names that had no matches
    unmatched <- partial_names[!sapply(partial_names, function(pn) any(startsWith(full_names, pn)))]

    list(indexes = idx, unmatched = unmatched)
}

# parse in many DE inputs, optionally filter by a list of cell types
# File names are cell type, age, region, DE_results.txt
parse_de_inputs<-function (in_dir, file_pattern="age", cellTypeListFile=NULL) {
    logger::log_info("Reading DE files from {in_dir} matching pattern '{file_pattern}'")
    #get the list of files from a directory that match some pattern.
    f=list.files(in_dir, pattern = paste0(file_pattern), full.names = TRUE)
    d=lapply(f, utils::read.table, sep="\t", header=TRUE,
             stringsAsFactors = FALSE, check.names = FALSE)

    #Set the file names (without directories) for each list element
    names(d)=basename(f)

    # look for the cell type in the names of the files
    if (!is.null(cellTypeListFile)) {
        cellTypeList=utils::read.table(cellTypeListFile, header=FALSE, stringsAsFactors = FALSE)$V1
        #add the variable tested to make this more specific.
        #this differentiates MSN_D1_age from MSN_D1_striosome_age.
        cellTypeList=paste0(cellTypeList, "_", file_pattern)
        r=find_prefix_matches(names(d), cellTypeList)
        if (length(r$unmatched)>0) {
            logger::log_warn("The following cell types in the cell type list had no matches in the data: {paste(r$unmatched, collapse=', ')}")
        }
        d=d[r$indexes]
        logger::log_info("Filtered to {length(d)} cell types based on provided cell type list.")
    }

    #simplifiy names
    names(d)=drop_common_suffix(names(d))
    #drop common substrings from the names
    names(d)=remove_common_tokens(names(d))

    return (d)
}
