# inDir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/metacells_eqtl/LEVEL_3"
# outDir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/metacells_eqtl/LEVEL_3_EUR"
# donorListFile="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/eqtl/2025-09-24_0_excluded_donors_all_villages/covariates/eur_donors.txt"

# 8
# filter_covariates_to_subpopulation(inDir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/metacells_eqtl/LEVEL_3", outDir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/metacells_eqtl/LEVEL_3_EUR", donorListFile="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/eqtl/2025-09-24_0_excluded_donors_all_villages/covariates/eur_donors.txt")

#' Filter metacell matrices to a donor subpopulation
#'
#' Read one or more metacell matrices from `inDir` and write filtered versions to
#' `outDir`, keeping only the `gene_symbol` column and donor columns listed in
#' `donorListFile` that are present in each file.
#'
#' Input files are matched by the pattern `metacells\\.txt\\.gz$`. Output is
#' written compressed (gz) to `outDir` using the input basenames.
#'
#' @param inDir Character scalar. Input directory containing `*metacells.txt.gz`
#'   files.
#' @param outDir Character scalar. Output directory to write filtered files. Will
#'   be created if it does not exist.
#' @param donorListFile Character scalar. Path to a tab-delimited file whose
#'   first column contains donor IDs matching donor column names in the metacell
#'   matrices. The file may have no header.
#'
#' @return Invisibly returns a character vector of output file paths written.
#' @export
#'
filter_metacells_to_subpopulation <- function(inDir, outDir, donorListFile) {
    .filter_files_to_subpopulation(
        inDir = inDir,
        outDir = outDir,
        donorListFile = donorListFile,
        file_pattern = "metacells\\.txt\\.gz$",
        id_col = "gene_symbol",
        compress_output = TRUE,
        log_every = 10L
    )
}

#' Filter covariate tables to a donor subpopulation
#'
#' Read one or more covariate tables from `inDir` and write filtered versions to
#' `outDir`, keeping only the `id` column and donor columns listed in
#' `donorListFile` that are present in each file.
#'
#' Input files are matched by the pattern `covariates\\.txt$`. Output is written
#' uncompressed (plain text) to `outDir` using the input basenames.
#'
#' @param inDir Character scalar. Input directory containing `*covariates.txt`
#'   files.
#' @param outDir Character scalar. Output directory to write filtered files. Will
#'   be created if it does not exist.
#' @param donorListFile Character scalar. Path to a tab-delimited file whose
#'   first column contains donor IDs matching donor column names in the covariate
#'   tables. The file may have no header.
#'
#' @return Invisibly returns a character vector of output file paths written.
#' @export
#'
filter_covariates_to_subpopulation <- function(inDir, outDir, donorListFile) {
    .filter_files_to_subpopulation(
        inDir = inDir,
        outDir = outDir,
        donorListFile = donorListFile,
        file_pattern = "covariates\\.txt$",
        id_col = "id",
        compress_output = FALSE,
        log_every = 10L
    )
}


# filter_metacells_to_subpopulation <- function(inDir, outDir, donorListFile) {
#
#     stopifnot(
#         is.character(inDir), length(inDir) == 1L, nzchar(inDir),
#         is.character(outDir), length(outDir) == 1L, nzchar(outDir),
#         is.character(donorListFile), length(donorListFile) == 1L, nzchar(donorListFile)
#     )
#
#     if (!dir.exists(inDir)) {
#         stop("inDir does not exist: ", inDir)
#     }
#     if (!file.exists(donorListFile)) {
#         stop("donorListFile does not exist: ", donorListFile)
#     }
#     if (!dir.exists(outDir)) {
#         dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
#     }
#
#     donors <- utils::read.table(
#         donorListFile,
#         header = FALSE,
#         stringsAsFactors = FALSE,
#         sep = "\t",
#         comment.char = "",
#         quote = ""
#     )[[1L]]
#
#     donors <- unique(as.character(donors))
#     donors <- donors[nzchar(donors)]
#     if (length(donors) == 0L) {
#         stop("No donor IDs found in donorListFile: ", donorListFile)
#     }
#
#     fileList <- list.files(
#         inDir,
#         pattern = "metacells\\.txt\\.gz$",
#         full.names = TRUE
#     )
#     if (length(fileList) == 0L) {
#         stop("No input files matching '*metacells.txt.gz' in: ", inDir)
#     }
#
#     outFiles <- character(0)
#
#     for (index in seq_along(fileList)) {
#         f <- fileList[index]
#
#         if (index %% 10 == 0L)
#             logger::log_info("Processing file ", index, " of ", length(fileList), "...")
#
#         a <- utils::read.table(
#             f,
#             header = TRUE,
#             stringsAsFactors = FALSE,
#             sep = "\t",
#             comment.char = "",
#             quote = ""
#         )
#
#         if (!("gene_symbol" %in% colnames(a))) {
#             stop("Missing required column 'gene_symbol' in file: ", f)
#         }
#
#         donorCols <- intersect(donors, colnames(a))
#         if (length(donorCols) == 0L) {
#             warning("No donor columns matched in file (writing gene_symbol only): ", f)
#         }
#
#         keepCols <- c("gene_symbol", donorCols)
#         aa <- a[, keepCols, drop = FALSE]
#
#         #write file, supporting .gz.
#         outFile <- file.path(outDir, basename(f))
#         con <- gzfile(outFile, open = "wt")
#         tryCatch(
#             utils::write.table(
#                 aa,
#                 file = con,
#                 quote = FALSE,
#                 sep = "\t",
#                 row.names = FALSE,
#                 col.names = TRUE
#             ),
#             finally = {
#                 close(con)
#             }
#         )
#     }
#
#     invisible(outFiles)
# }
#
# filter_covariates_to_subpopulation <- function(inDir, outDir, donorListFile) {
#
#     stopifnot(
#         is.character(inDir), length(inDir) == 1L, nzchar(inDir),
#         is.character(outDir), length(outDir) == 1L, nzchar(outDir),
#         is.character(donorListFile), length(donorListFile) == 1L, nzchar(donorListFile)
#     )
#
#     if (!dir.exists(inDir)) {
#         stop("inDir does not exist: ", inDir)
#     }
#     if (!file.exists(donorListFile)) {
#         stop("donorListFile does not exist: ", donorListFile)
#     }
#     if (!dir.exists(outDir)) {
#         dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
#     }
#
#     donors <- utils::read.table(
#         donorListFile,
#         header = FALSE,
#         stringsAsFactors = FALSE,
#         sep = "\t",
#         comment.char = "",
#         quote = ""
#     )[[1L]]
#
#     donors <- unique(as.character(donors))
#     donors <- donors[nzchar(donors)]
#     if (length(donors) == 0L) {
#         stop("No donor IDs found in donorListFile: ", donorListFile)
#     }
#
#     fileList <- list.files(
#         inDir,
#         pattern = "covariates.txt$",
#         full.names = TRUE
#     )
#     if (length(fileList) == 0L) {
#         stop("No input files matching '*covariates.txt' in: ", inDir)
#     }
#
#     outFiles <- character(0)
#
#     for (index in seq_along(fileList)) {
#         f <- fileList[index]
#
#         if (index %% 10 == 0L)
#             logger::log_info("Processing file ", index, " of ", length(fileList), "...")
#
#         a <- utils::read.table(
#             f,
#             header = TRUE,
#             stringsAsFactors = FALSE,
#             sep = "\t",
#             comment.char = "",
#             quote = ""
#         )
#
#         if (!("id" %in% colnames(a))) {
#             stop("Missing required column 'gene_symbol' in file: ", f)
#         }
#
#         donorCols <- intersect(donors, colnames(a))
#         if (length(donorCols) == 0L) {
#             warning("No donor columns matched in file (writing gene_symbol only): ", f)
#         }
#
#         keepCols <- c("id", donorCols)
#         aa <- a[, keepCols, drop = FALSE]
#
#         #write file, supporting .gz.
#         outFile <- file.path(outDir, basename(f))
#
#         utils::write.table(
#             aa,
#             file = outFile,
#             quote = FALSE,
#             sep = "\t",
#             row.names = FALSE,
#             col.names = TRUE
#         )
#
#     }
#
#     invisible(outFiles)
# }


# Helper ---------------------------------------------------------------

.filter_files_to_subpopulation <- function(inDir,
                                           outDir,
                                           donorListFile,
                                           file_pattern,
                                           id_col,
                                           input_sep = "\t",
                                           output_sep = "\t",
                                           compress_output = TRUE,
                                           log_every = 10L) {
    stopifnot(
        is.character(inDir), length(inDir) == 1L, nzchar(inDir),
        is.character(outDir), length(outDir) == 1L, nzchar(outDir),
        is.character(donorListFile), length(donorListFile) == 1L, nzchar(donorListFile),
        is.character(file_pattern), length(file_pattern) == 1L, nzchar(file_pattern),
        is.character(id_col), length(id_col) == 1L, nzchar(id_col),
        is.character(input_sep), length(input_sep) == 1L,
        is.character(output_sep), length(output_sep) == 1L,
        is.logical(compress_output), length(compress_output) == 1L,
        is.numeric(log_every), length(log_every) == 1L, !is.na(log_every)
    )

    if (!dir.exists(inDir)) {
        stop("inDir does not exist: ", inDir)
    }
    if (!file.exists(donorListFile)) {
        stop("donorListFile does not exist: ", donorListFile)
    }
    if (!dir.exists(outDir)) {
        dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
    }

    donors <- utils::read.table(
        donorListFile,
        header = FALSE,
        stringsAsFactors = FALSE,
        sep = "\t",
        comment.char = "",
        quote = ""
    )[[1L]]

    donors <- unique(as.character(donors))
    donors <- donors[nzchar(donors)]
    if (length(donors) == 0L) {
        stop("No donor IDs found in donorListFile: ", donorListFile)
    }

    fileList <- list.files(
        inDir,
        pattern = file_pattern,
        full.names = TRUE
    )
    if (length(fileList) == 0L) {
        stop("No input files matching pattern '", file_pattern, "' in: ", inDir)
    }

    outFiles <- character(0)

    for (index in seq_along(fileList)) {
        f <- fileList[index]

        if (log_every > 0L && index %% log_every == 0L) {
            logger::log_info("Processing file ", index, " of ", length(fileList), "...")
        }

        a <- utils::read.table(
            f,
            header = TRUE,
            stringsAsFactors = FALSE,
            sep = input_sep,
            comment.char = "",
            quote = ""
        )

        if (!(id_col %in% colnames(a))) {
            stop("Missing required column '", id_col, "' in file: ", f)
        }

        donorCols <- intersect(donors, colnames(a))
        if (length(donorCols) == 0L) {
            warning("No donor columns matched in file (writing ", id_col, " only): ", f)
        }

        keepCols <- c(id_col, donorCols)
        aa <- a[, keepCols, drop = FALSE]

        outFile <- file.path(outDir, basename(f))

        if (compress_output) {
            con <- gzfile(outFile, open = "wt")
            tryCatch(
                utils::write.table(
                    aa,
                    file = con,
                    quote = FALSE,
                    sep = output_sep,
                    row.names = FALSE,
                    col.names = TRUE
                ),
                finally = {
                    close(con)
                }
            )
        } else {
            utils::write.table(
                aa,
                file = outFile,
                quote = FALSE,
                sep = output_sep,
                row.names = FALSE,
                col.names = TRUE
            )
        }

        outFiles <- c(outFiles, outFile)
    }

    invisible(outFiles)
}


