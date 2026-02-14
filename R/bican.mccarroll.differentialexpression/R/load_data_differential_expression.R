# in_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results/LEVEL_3/sex_age/cell_type"
# in_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results/LEVEL_3/sex_age/cell_type_region_interaction_absolute_effects"
# cellTypeListFile="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metadata/mash_cell_type_list_simple.txt"
# file_pattern="age"
# fdrThreshold=0.01


#' Gather differential expression results and construct mash-ready matrices
#'
#' @description
#' Read per-condition (or per-cell-type) differential expression (DE) result tables
#' from `in_dir`, optionally restrict to a provided list of cell types, and
#' construct mash-compatible matrices (`Bhat`, `Shat`, `Tstat`, `FDR`, etc.) using
#' `make_joint_inputs()` with `gene_mode = "union"`.
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
#' }
#'
#' @export
gather_de_results<-function (in_dir, file_pattern="age", cellTypeListFile=NULL, fdrThreshold=0.01) {

    d=parse_de_inputs(in_dir, file_pattern, cellTypeListFile)
    #make mash inputs
    inputs_union<-make_joint_inputs(d, coef_col = "logFC", t_col = "t", fdr_col = "adj.P.Val", gene_mode="union")
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


#' Build Bhat/Shat/Tstat/FDR matrices from a long DE results data frame
#'
#' @description
#' Construct rectangular matrices suitable for \code{mashr} or clustering from a
#' single "long" differential expression (DE) result data frame. The input must
#' contain one row per \code{gene} and condition, where the condition is defined
#' by \code{cell_type} and (optionally) \code{interaction}.
#'
#' The input is assumed to be pre-filtered to a single contrast (for example,
#' all rows correspond to \code{contrast == "age"}). This function enforces that
#' exactly one unique \code{contrast} value is present, and then drops contrast
#' from the matrix column names.
#'
#' For each condition, the function reads a coefficient column (for example,
#' \code{logFC}), a t-statistic column (for example, \code{t}), a p-value column,
#' and an FDR column. Standard errors are computed as \code{Shat = abs(beta / t)}
#' (so \code{t = beta / se} as in limma).
#'
#' @param df Data frame in long format. Must contain columns \code{gene},
#'   \code{cell_type}, \code{interaction}, and \code{contrast}, plus the columns
#'   specified by \code{coef_col}, \code{t_col}, \code{p_col}, and \code{fdr_col}.
#'   Each \code{(gene, condition)} combination must appear at most once.
#' @param coef_col Character scalar. Column name for the effect size (beta),
#'   for example \code{"logFC"}. Default \code{"logFC"}.
#' @param t_col Character scalar. Column name for the t-statistic, for example
#'   \code{"t"}. Default \code{"t"}.
#' @param p_col Character scalar. Column name for the p-value, for example
#'   \code{"P.Value"}. Default \code{"P.Value"}.
#' @param fdr_col Character scalar. Column name for the FDR or q-value, for example
#'   \code{"adj.P.Val"}. Default \code{"adj.P.Val"}.
#' @param gene_mode One of \code{c("union", "intersect")}. Use \code{"union"} to
#'   include all genes observed in any condition (missing entries are replaced with
#'   \code{beta = 0} and \code{SE = big_se}). Use \code{"intersect"} to keep only
#'   genes present in every condition. Default \code{"union"}.
#' @param big_se Numeric scalar. Standard error to use when filling missing or
#'   invalid entries. A very large value ensures such entries contribute ~0
#'   information to mash. Default \code{1e6}.
#' @param fill_missing_fdr Numeric scalar used to fill missing FDR entries in
#'   the returned \code{FDR} matrix (for plotting or logic only; not used by mash).
#'   Default \code{NA_real_}.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{Bhat}}{Numeric matrix \code{[genes x conditions]} of betas.}
#'   \item{\code{Shat}}{Numeric matrix \code{[genes x conditions]} of standard errors
#'     (\code{abs(beta / t)}).}
#'   \item{\code{Tstat}}{Numeric matrix \code{[genes x conditions]} of t-statistics.}
#'   \item{\code{P_val}}{Numeric matrix \code{[genes x conditions]} of p-values.}
#'   \item{\code{FDR}}{Numeric matrix \code{[genes x conditions]} of FDR values.}
#'   \item{\code{missing_mask}}{Logical matrix \code{[genes x conditions]} indicating
#'     entries that were missing or invalid and replaced by \code{beta = 0},
#'     \code{SE = big_se}, \code{t = 0}.}
#' }
#'
#' @details
#' Conditions (matrix columns) are defined as:
#' \describe{
#'   \item{No interaction}{If \code{interaction} is entirely \code{NA}, columns are
#'     named by \code{cell_type}.}
#'   \item{With interaction}{If \code{interaction} is entirely non-\code{NA}, columns
#'     are named by \code{paste(cell_type, interaction, sep = "__")}.}
#' }
#' Mixed \code{NA}/non-\code{NA} values in \code{interaction} are treated as an error.
#'
#' With \code{gene_mode = "union"}, rectangular matrices are formed by inserting
#' neutral placeholders for missing entries (\code{beta = 0}, \code{SE = big_se}).
#' With \code{gene_mode = "intersect"}, only genes present in every condition are
#' retained (pathological entries such as \code{t = 0} are still treated as missing
#' and assigned \code{SE = big_se}).
#'
#' @export
make_joint_inputs <- function(
        df,
        coef_col = "logFC",
        t_col    = "t",
        p_col    = "P.Value",
        fdr_col  = "adj.P.Val",
        gene_mode = c("union", "intersect"),
        big_se   = 1e6,
        fill_missing_fdr = NA_real_
) {
    gene_mode <- match.arg(gene_mode)

    stopifnot(is.data.frame(df))

    meta_cols <- c("gene", "cell_type", "interaction", "contrast")
    if (!all(meta_cols %in% colnames(df))) {
        stop(
            "Input df must contain columns: ",
            paste(meta_cols, collapse = ", "),
            call. = FALSE
        )
    }

    required <- c(coef_col, t_col, p_col, fdr_col)
    if (!all(required %in% colnames(df))) {
        stop(
            "Input df is missing required columns: ",
            paste(setdiff(required, colnames(df)), collapse = ", "),
            call. = FALSE
        )
    }

    ## Enforce exactly one unique contrast in the input
    u_contrast <- unique(as.character(df[["contrast"]]))
    u_contrast <- u_contrast[!is.na(u_contrast) & nzchar(u_contrast)]
    if (length(u_contrast) != 1L) {
        stop(
            "Expected exactly one unique contrast in df$contrast; found: ",
            paste(u_contrast, collapse = ", "),
            call. = FALSE
        )
    }

    interaction_chr <- as.character(df[["interaction"]])
    interaction_chr[is.na(interaction_chr) | interaction_chr == ""] <- NA_character_

    ## Enforce your invariant: interaction is either entirely NA or entirely filled.
    all_na_interaction <- all(is.na(interaction_chr))
    any_na_interaction <- any(is.na(interaction_chr))
    if (!all_na_interaction && any_na_interaction) {
        stop("Mixed NA/non-NA values found in 'interaction'. Input must be consistent.", call. = FALSE)
    }

    ## Build condition names WITHOUT contrast (since it is constant)
    if (all_na_interaction) {
        condition <- as.character(df[["cell_type"]])
    } else {
        condition <- paste(df[["cell_type"]], interaction_chr, sep = "__")
    }

    genes_all <- as.character(df[["gene"]])
    cond_levels <- unique(condition)

    ## Determine gene set for matrices
    genes_by_cond <- split(genes_all, condition)
    if (gene_mode == "union") {
        genes <- sort(unique(unlist(genes_by_cond, use.names = FALSE)))
    } else {
        genes <- sort(unique(Reduce(intersect, genes_by_cond)))
    }
    G <- length(genes)
    C <- length(cond_levels)
    if (G == 0L) stop("No genes found after applying gene_mode = '", gene_mode, "'.", call. = FALSE)
    if (C < 2L) stop("Need at least two conditions to build matrices.", call. = FALSE)

    ## Check uniqueness of (gene, condition)
    key <- paste(condition, genes_all, sep = "\r")
    if (any(duplicated(key))) {
        stop("Duplicate (condition, gene) rows found in input df.", call. = FALSE)
    }

    ## Preallocate outputs
    Bhat <- matrix(NA_real_, G, C, dimnames = list(genes, cond_levels))
    Shat <- matrix(NA_real_, G, C, dimnames = list(genes, cond_levels))
    Tstat<- matrix(NA_real_, G, C, dimnames = list(genes, cond_levels))
    FDR  <- matrix(NA_real_, G, C, dimnames = list(genes, cond_levels))
    P_val<- matrix(NA_real_, G, C, dimnames = list(genes, cond_levels))

    ## Fill matrices condition-by-condition
    for (j in seq_along(cond_levels)) {
        cj <- cond_levels[[j]]
        idx <- which(condition == cj)
        if (!length(idx)) next

        g <- genes_all[idx]

        keep <- g %in% genes
        if (!any(keep)) next
        idx <- idx[keep]
        g <- g[keep]

        beta <- as.numeric(df[idx, coef_col])
        tt   <- as.numeric(df[idx, t_col])
        pval <- as.numeric(df[idx, p_col])
        fdr  <- as.numeric(df[idx, fdr_col])

        se <- abs(beta / tt)

        Bhat[g, j]  <- beta
        Shat[g, j]  <- se
        Tstat[g, j] <- tt
        P_val[g, j] <- pval
        FDR[g, j]   <- fdr
    }

    ## Handle missing/invalid
    missing_mask <- !is.finite(Bhat) | !is.finite(Shat) | (Shat <= 0)
    if (any(missing_mask)) {
        Bhat[missing_mask]  <- 0
        Shat[missing_mask]  <- big_se
        Tstat[missing_mask] <- 0
    }

    nas_fdr <- !is.finite(FDR)
    if (any(nas_fdr)) FDR[nas_fdr] <- fill_missing_fdr

    list(
        Bhat = Bhat,
        Shat = Shat,
        Tstat = Tstat,
        P_val = P_val,
        FDR = FDR,
        missing_mask = missing_mask
    )
}






# parse in many DE inputs, optionally filter by a list of cell types
# File names are cell type, age, region, DE_results.txt

#' Parse differential expression result tables from a directory
#'
#' @description
#' Read multiple differential expression (DE) result files from a directory,
#' combine them into a single data frame, and annotate each row with metadata
#' derived from the file name.
#'
#' Filenames must follow one of these basename patterns:
#' \describe{
#'   \item{\code{celltype__contrast_DE_results.txt}}{No interaction term; \code{interaction} is \code{NA}.}
#'   \item{\code{celltype__interaction__contrast_DE_results.txt}}{Includes an interaction term.}
#' }
#'
#' Files are selected using \code{list.files(..., pattern = file_pattern)}, read
#' with \code{utils::read.table()}, and row-bound into one data frame. For each
#' file, three constant columns are prepended to its table:
#' \code{cell_type}, \code{interaction}, and \code{contrast}.
#'
#' Optionally, the combined result can be restricted to a supplied list of
#' cell types.
#'
#' @param in_dir Character scalar. Directory containing DE result files.
#' @param file_pattern Character scalar. Pattern used to select files from
#'   \code{in_dir} via \code{list.files(..., pattern = file_pattern)}.
#' @param cellTypeListFile Character scalar or \code{NULL}. Optional path to a text
#'   file containing one cell type per line. If provided, the combined data frame
#'   is filtered to rows whose \code{cell_type} is in this list. Default \code{NULL}.
#'
#' @return
#' A data frame produced by row-binding all retained DE result tables, with three
#' additional columns prepended:
#' \describe{
#'   \item{\code{cell_type}}{Cell type parsed from the file basename.}
#'   \item{\code{interaction}}{Interaction token parsed from the file basename, or \code{NA}.}
#'   \item{\code{contrast}}{Contrast parsed from the file basename.}
#' }
#'
#' @details
#' Parsing of \code{cell_type}, \code{interaction}, and \code{contrast} is performed
#' by \code{parse_de_result_filenames()}, which expects the common suffix
#' \code{_DE_results.txt}.
#'
#' @export
parse_de_inputs<-function (in_dir, file_pattern="age", cellTypeListFile=NULL) {
    logger::log_info("Reading DE files from {in_dir} matching pattern '{file_pattern}'")
    #get the list of files from a directory that match some pattern.
    f=list.files(in_dir, pattern = paste0(file_pattern), full.names = TRUE)

    d=lapply(f, utils::read.table, sep="\t", header=TRUE,
             stringsAsFactors = FALSE, check.names = FALSE)

    #Set the file names (with5out directories) for each list element
    names(d)=basename(f)

    #extract the cell type and contrast from the file name.
    cell_type_contrast_df<- parse_de_result_filenames(names(d))

    #add the cell type and contrast to the list of data frames
    # Add cell_type and contrast as FIRST two columns
    for (i in seq_along(d)) {
        df <- d[[i]]
        df <- cbind(
            gene = rownames(df),
            cell_type = cell_type_contrast_df$celltype[i],
            interaction = cell_type_contrast_df$interaction[i],
            contrast  = cell_type_contrast_df$contrast[i],
            df
        )
        d[[i]] <- df
    }

    #make one big dataframe
    df<-do.call(rbind, d)
    rownames (df)<-NULL


    # look for the cell type in the names of the files
    if (!is.null(cellTypeListFile)) {
        cellTypeList=utils::read.table(cellTypeListFile, header=FALSE, stringsAsFactors = FALSE)$V1
        #add the variable tested to make this more specific.
        #this differentiates MSN_D1_age from MSN_D1_striosome_age.
        df=df[df$cell_type %in% cellTypeList,]
        num_unique_cell_types=length(unique(df$cell_type))
        logger::log_info("Filtered to {num_unique_cell_types} cell types based on provided cell type list.")
    }

    return (df)
}

#' Parse DE result filenames into structured metadata
#'
#' @description
#' Parse differential expression (DE) result filenames and extract structured
#' metadata fields based on a strict naming convention.
#'
#' Filenames must follow one of the two supported basename formats:
#' \describe{
#'   \item{\code{celltype__contrast_DE_results.txt}}{
#'     No interaction term; \code{interaction} is returned as \code{NA}.
#'   }
#'   \item{\code{celltype__interaction__contrast_DE_results.txt}}{
#'     Includes an interaction term between cell type and contrast.
#'   }
#' }
#'
#' The common suffix \code{_DE_results.txt} is removed before parsing.
#' Components are separated by double underscores (\code{"__"}).
#'
#' @param files Character vector of file paths or basenames to parse.
#'
#' @return
#' A data frame with one row per input file and the following columns:
#' \describe{
#'   \item{\code{file}}{Original input file path or name.}
#'   \item{\code{basename}}{Basename of the file.}
#'   \item{\code{celltype}}{Cell type extracted from the filename.}
#'   \item{\code{interaction}}{Interaction term if present, otherwise \code{NA}.}
#'   \item{\code{contrast}}{Contrast extracted from the filename.}
#' }
#'
#' @details
#' An error is thrown if any filename does not match one of the two supported
#' formats (i.e., does not split into exactly two or three fields after removing
#' the suffix).
#'
#' @export
parse_de_result_filenames <- function(files) {

    bn <- basename(files)

    suffix <- "_DE_results.txt"
    core <- sub(paste0(suffix, "$"), "", bn)

    parts <- strsplit(core, "__", fixed = TRUE)

    n_parts <- vapply(parts, length, integer(1))
    ok <- n_parts %in% c(2L, 3L)
    if (!all(ok)) {
        bad <- bn[!ok]
        stop(
            "Unexpected filename format for: ",
            paste(bad, collapse = ", "),
            call. = FALSE
        )
    }

    celltype <- vapply(parts, function(x) x[[1L]], character(1))

    interaction <- vapply(
        parts,
        function(x) if (length(x) == 3L) x[[2L]] else NA_character_,
        character(1)
    )

    contrast <- vapply(
        parts,
        function(x) if (length(x) == 3L) x[[3L]] else x[[2L]],
        character(1)
    )

    data.frame(
        file = files,
        basename = bn,
        celltype = celltype,
        interaction = interaction,
        contrast = contrast,
        stringsAsFactors = FALSE
    )
}

# find_prefix_matches <- function(full_names, partial_names) {
#     # indexes in full_names that match any prefix in partial_names
#     idx <- which(sapply(full_names, function(x) any(startsWith(x, partial_names))))
#
#     # prefixes in partial_names that had no matches
#     unmatched <- partial_names[!sapply(partial_names, function(pn) any(startsWith(full_names, pn)))]
#
#     list(indexes = idx, unmatched = unmatched)
# }

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
# remove_common_tokens <- function(x, delim = "_") {
#     stopifnot(is.character(x))
#     toks <- strsplit(x, delim, fixed = TRUE)
#     common <- Reduce(intersect, lapply(toks, unique))
#     if (length(common) == 0L) return(x)
#     out <- vapply(
#         toks,
#         function(v) paste(v[!v %in% common], collapse = delim),
#         character(1)
#     )
#     attr(out, "common_tokens") <- common
#     out
# }

#' Drop common suffix from file names
#' Useful to infer the data set names from differential expression files
#' @param files A list of file names to process
#' @return A vector with the same length as the input files containing the name of the file without
#' the common suffix.
#' @noRd
# drop_common_suffix <- function(files) {
#     # Reverse each filename so suffix becomes a prefix
#     rev_files <- sapply(files, function(x) paste0(rev(strsplit(x, NULL)[[1]]), collapse = ""))
#
#     # Find common prefix of reversed strings
#     common_rev_prefix <- Reduce(function(a, b) {
#         i <- 1
#         while (i <= nchar(a) && i <= nchar(b) && substr(a, i, i) == substr(b, i, i)) {
#             i <- i + 1
#         }
#         substr(a, 1, i - 1)
#     }, rev_files)
#
#     # Reverse back to get the suffix
#     common_suffix <- paste0(rev(strsplit(common_rev_prefix, NULL)[[1]]), collapse = "")
#
#     # Drop suffix from each filename
#     sub(paste0(common_suffix, "$"), "", files)
# }
