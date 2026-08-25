# Convert Kana et al. 2026 (bioRxiv) Supplementary Table 4 (global CPS) into
# the canonical "voom-like" DE result format used elsewhere in this repo, so
# it can be loaded interchangeably with native limma/voom/dream outputs via
# bican.mccarroll.differentialexpression::parse_de_inputs().
#
# Supplementary Table 4 reports a per-(cell_type, gene) linear model with
# terms (Intercept), Age_at_Death_binned_codes, Sex_codes, method10xV3.1_HT,
# APOE4_Status_codes, and CPS (Continuous Pseudo-progression Score), each with
# logFC_/se_/p_ columns (and fdr_p_ for every non-intercept term). Only the
# CPS term is extracted here -- it is the CaH analogue of the Gabitto 2024
# MTG "ad_cps" contrast already converted by PMID_39402379.R in this
# directory, so the two are given the same contrast token and differ only in
# their region token (CaH vs MTG).
#
# The source file's column *names* are misleading for two columns: despite
# being named "gene_id" and "var_name"/"var_index", column 1 (row name) is
# "<cell_type>_<gene>", "gene_id" is a 1-based row index *within* each cell
# type (not a gene identifier), and "var_name"/"var_index" are identical gene
# symbol columns. The gene symbol used throughout this script is "var_name".

inDegTable <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/external_comparison_kana_2026_biorxiv/original_data/supplemental_table_4_global_cps.csv"
outDir     <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/external_comparison_kana_2026_biorxiv/voom-like"
region     <- "CaH"
contrast   <- "ad_cps"

# ---------------------------------------------------------------------------
# Column selection / reading
# ---------------------------------------------------------------------------

#' Columns read from the source CSV, and the canonical adapter names they map
#' to. Only the CPS term columns plus the gene/cell-type keys are needed; the
#' source file is 612MB with 28 columns, so selecting these six via
#' \code{data.table::fread(select = ...)} avoids materializing the rest.
kana_source_columns <- c(
    gene = "var_name",
    logFC = "logFC_CPS",
    se = "se_CPS",
    p_value = "p_CPS",
    fdr_p_value = "fdr_p_CPS",
    cell_type = "cell_type"
)

#' Read and rename the CPS adapter columns from a Kana 2026 supplementary
#' table.
#'
#' @param in_file Character scalar. Path to the source CSV.
#' @return A plain data.frame (not a data.table) with columns \code{gene},
#'   \code{logFC}, \code{se}, \code{p_value}, \code{fdr_p_value},
#'   \code{cell_type}, in source row order. Converting away from data.table
#'   deliberately avoids a data.table NSE footgun: subsetting a data.table
#'   with an \code{i} expression like \code{full$cell_type == cell_type}
#'   resolves the bare \code{cell_type} symbol against the table's own
#'   \code{cell_type} column first, silently shadowing an identically named
#'   R variable/argument and turning the filter into a no-op (every row
#'   matches). Plain data.frame subsetting uses ordinary lexical scoping and
#'   has no such hazard.
read_kana_2026_table <- function(in_file) {
    missing_cols <- setdiff(kana_source_columns, colnames(data.table::fread(in_file, nrows = 0)))
    if (length(missing_cols) > 0L) {
        stop(
            "Source file '", in_file, "' is missing expected column(s): ",
            paste(missing_cols, collapse = ", "),
            call. = FALSE
        )
    }

    raw <- data.table::fread(in_file, select = unname(kana_source_columns))
    data.table::setnames(raw, unname(kana_source_columns), names(kana_source_columns))
    data.table::setDF(raw)
    raw
}

# ---------------------------------------------------------------------------
# Per-cell-type conversion
# ---------------------------------------------------------------------------

#' Construct the canonical output filename for one cell type.
#'
#' @param cell_type Character scalar. Source cell type name, kept verbatim.
#' @param region Character scalar. Brain region token.
#' @param contrast Character scalar. Canonical contrast name.
#' @return Character scalar basename, e.g. "Astro_1__CaH__ad_cps_DE_results.txt".
kana_output_filename <- function(cell_type, region, contrast) {
    paste0(cell_type, "__", region, "__", contrast, "_DE_results.txt")
}

#' Validate the adapter rows for one cell type against the invariants
#' verified against the full source file (non-empty/unique gene names,
#' finite logFC/se/p, se > 0, p in [0, 1]). The source data has none of
#' these problems, so a violation here indicates the input file changed and
#' should stop the conversion rather than silently drop rows.
#'
#' @param z Data.table with columns \code{gene}, \code{logFC}, \code{se},
#'   \code{p_value}, \code{fdr_p_value} for a single cell type.
#' @param cell_type Character scalar, used only for error messages.
#' @return Invisibly \code{TRUE} if all checks pass; otherwise throws.
validate_kana_2026_cell_type <- function(z, cell_type) {
    if (any(is.na(z$gene) | trimws(z$gene) == "")) {
        stop("Cell type '", cell_type, "': blank or missing gene name(s) found.", call. = FALSE)
    }
    if (anyDuplicated(z$gene) != 0L) {
        stop("Cell type '", cell_type, "': duplicate gene name(s) found.", call. = FALSE)
    }
    if (any(!is.finite(z$logFC)) || any(!is.finite(z$se)) || any(!is.finite(z$p_value))) {
        stop("Cell type '", cell_type, "': non-finite logFC/se/p_value found.", call. = FALSE)
    }
    if (any(z$se <= 0)) {
        stop("Cell type '", cell_type, "': non-positive standard error(s) found.", call. = FALSE)
    }
    if (any(z$p_value < 0 | z$p_value > 1)) {
        stop("Cell type '", cell_type, "': p-value(s) outside [0, 1] found.", call. = FALSE)
    }
    invisible(TRUE)
}

#' Build the canonical seven-column voom-like table for one cell type, and
#' recompute the BH adjustment fresh from p_value (rather than copying the
#' source fdr_p_value) so the output is internally consistent, matching the
#' convention used in PMID_39402379.R.
#'
#' @param z Data.table with columns \code{gene}, \code{logFC}, \code{se},
#'   \code{p_value}, \code{fdr_p_value} for a single cell type, in source row
#'   order.
#' @return A list with \code{de_table} (a data.frame with columns
#'   \code{logFC}, \code{AveExpr}, \code{t}, \code{P.Value}, \code{adj.P.Val},
#'   \code{B}, \code{z.std}, row names set to \code{z$gene}) and
#'   \code{max_fdr_diff} (max absolute difference between the recomputed and
#'   source FDR, for the manifest).
build_kana_2026_de_table <- function(z) {
    adj_p_val <- stats::p.adjust(z$p_value, method = "BH")

    de_table <- data.frame(
        logFC = z$logFC,
        AveExpr = NA_real_,
        t = z$logFC / z$se,
        P.Value = z$p_value,
        adj.P.Val = adj_p_val,
        B = NA_real_,
        z.std = NA_real_,
        stringsAsFactors = FALSE
    )
    rownames(de_table) <- z$gene

    list(
        de_table = de_table,
        max_fdr_diff = max(abs(adj_p_val - z$fdr_p_value))
    )
}

#' Build an empty (zero-count) manifest row template for one cell type.
#'
#' @param cell_type Character scalar. Source cell type name.
#' @param region Character scalar. Brain region token.
#' @param contrast Character scalar. Canonical contrast name.
#' @param output_file Character scalar. Expected/actual output path.
#' @return A one-row manifest data frame.
empty_kana_manifest_row <- function(cell_type, region, contrast, output_file) {
    data.frame(
        cell_type = cell_type,
        region = region,
        contrast = contrast,
        output_file = output_file,
        n_source_rows = NA_integer_,
        n_output_rows = 0L,
        max_fdr_diff = NA_real_,
        status = "skipped",
        message = "",
        stringsAsFactors = FALSE
    )
}

#' Convert exactly one cell type's rows into one canonical DE output file.
#'
#' @param full A data.table as returned by \code{read_kana_2026_table()},
#'   containing every cell type.
#' @param cell_type Character scalar. Cell type to extract and convert.
#' @param region Character scalar. Brain region token.
#' @param contrast Character scalar. Canonical contrast name.
#' @param out_dir Character scalar. Directory into which the DE output file
#'   is written.
#' @return A one-row manifest data frame; see \code{convert_kana_2026()} for
#'   the manifest column list.
convert_kana_2026_cell_type <- function(full, cell_type, region, contrast, out_dir) {
    output_file <- file.path(out_dir, kana_output_filename(cell_type, region, contrast))
    manifest_row <- empty_kana_manifest_row(cell_type, region, contrast, output_file)

    z <- full[full$cell_type == cell_type, ]
    manifest_row$n_source_rows <- nrow(z)

    validate_kana_2026_cell_type(z, cell_type)

    built <- build_kana_2026_de_table(z)
    manifest_row$max_fdr_diff <- built$max_fdr_diff

    if (built$max_fdr_diff > 1e-8) {
        msg <- paste0(
            "Recomputed BH adjustment disagrees with source fdr_p_CPS by up to ",
            built$max_fdr_diff, "; source FDR may not be per-cell-type. Skipping."
        )
        logger::log_warn("Cell type '{cell_type}': {msg}")
        manifest_row$message <- msg
        return(manifest_row)
    }

    utils::write.table(
        built$de_table,
        output_file,
        row.names = TRUE,
        col.names = TRUE,
        quote = FALSE,
        sep = "\t"
    )

    manifest_row$n_output_rows <- nrow(built$de_table)
    manifest_row$status <- "written"

    manifest_row
}

# ---------------------------------------------------------------------------
# Top-level conversion
# ---------------------------------------------------------------------------

#' Convert every cell type in a Kana 2026 Supplementary Table 4-style CSV
#' into canonical voom-like DE output files, one per cell type, for a single
#' contrast term (CPS by default).
#'
#' All output paths are computed and checked for collisions, and round-tripped
#' through \code{parse_de_result_filenames()}, before any file is written
#' ("preflight"), so a naming problem is detected before a partial result set
#' is produced.
#'
#' @param in_file Character scalar. Path to the source CSV
#'   (e.g. supplemental_table_4_global_cps.csv).
#' @param region Character scalar. Brain region token for the output
#'   filenames (e.g. "CaH").
#' @param contrast Character scalar. Canonical contrast token for the output
#'   filenames (e.g. "ad_cps").
#' @param out_dir Character scalar. Directory into which DE output files and
#'   the conversion manifest are written. Created recursively if needed.
#' @return Invisibly, a data frame manifest with one row per cell type and
#'   columns \code{cell_type}, \code{region}, \code{contrast},
#'   \code{output_file}, \code{n_source_rows}, \code{n_output_rows},
#'   \code{max_fdr_diff}, \code{status}, and \code{message}. The same
#'   manifest is written to \code{<out_dir>/conversion_manifest.txt}.
#' @export
convert_kana_2026 <- function(in_file, region, contrast, out_dir) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    logger::log_info("Reading '{in_file}'.")
    full <- read_kana_2026_table(in_file)

    cell_types <- sort(unique(full$cell_type))
    logger::log_info("Found {length(cell_types)} cell type(s) in '{in_file}'.")

    # Preflight: compute every expected output path before converting or
    # writing anything.
    output_files <- vapply(
        cell_types,
        function(ct) file.path(out_dir, kana_output_filename(ct, region, contrast)),
        character(1)
    )

    dup_files <- unique(output_files[duplicated(output_files)])
    if (length(dup_files) > 0L) {
        stop(
            "Output filename collision(s) detected across ", length(dup_files),
            " output path(s). Aborting before writing any DE files.",
            call. = FALSE
        )
    }

    # Round-trip every expected output path through the existing DE filename
    # parser, to confirm the generated filenames stay compatible with
    # downstream infrastructure before any file is written.
    parsed_back <- bican.mccarroll.differentialexpression::parse_de_result_filenames(output_files)
    for (i in seq_along(cell_types)) {
        if (!identical(parsed_back$celltype[i], cell_types[i]) ||
            !identical(parsed_back$interaction[i], region) ||
            !identical(parsed_back$contrast[i], contrast)) {
            stop(
                "Generated output filename for cell type '", cell_types[i],
                "' does not round-trip through parse_de_result_filenames(): '",
                output_files[i], "'",
                call. = FALSE
            )
        }
    }

    logger::log_info("Preflight complete: {length(cell_types)} cell type(s), no output collisions.")

    manifest_rows <- lapply(
        cell_types,
        function(ct) convert_kana_2026_cell_type(full, ct, region, contrast, out_dir)
    )

    manifest_columns <- c(
        "cell_type", "region", "contrast", "output_file",
        "n_source_rows", "n_output_rows", "max_fdr_diff", "status", "message"
    )

    manifest <- do.call(rbind, manifest_rows)
    manifest <- manifest[, manifest_columns]

    n_written <- sum(manifest$status == "written")
    n_skipped <- sum(manifest$status == "skipped")
    logger::log_info(
        "Kana 2026 conversion complete: {n_written} file(s) written, {n_skipped} cell type(s) skipped."
    )

    manifest_file <- file.path(out_dir, "conversion_manifest.txt")
    utils::write.table(
        manifest,
        manifest_file,
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE,
        sep = "\t"
    )

    invisible(manifest)
}

# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

# manifest <- convert_kana_2026(inDegTable, region, contrast, outDir)
# table(manifest$status)
# stopifnot(all(manifest$status == "written"))
