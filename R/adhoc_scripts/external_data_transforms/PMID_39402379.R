# Convert Gabitto 2024 (PMID 39402379) Supplementary Table 5 into the
# canonical "voom-like" DE result format used elsewhere in this repo, so it
# can be loaded interchangeably with native limma/voom/dream outputs via
# bican.mccarroll.differentialexpression::parse_de_inputs().
#
# See Gabitto_2024_Supplementary_Table_5_conversion_spec.md in this directory
# for the full specification this script implements.
#
# Supplementary Table 5 is a single ~3.5GB workbook with 555 DE-result
# sheets. Gabitto-specific sheet-name and column interpretation is confined to
# this script; the files it writes use the same three-field naming convention
# and seven-column schema as other DE results in this project.

inWorkbook <- "/downloads/Gabitto_2024/original/gabitto_2024_supplementary_table_5.xlsx"
outDir     <- "/downloads/Gabitto_2024/voom-like"

# ---------------------------------------------------------------------------
# Sheet-name parsing
# ---------------------------------------------------------------------------

#' Normalize a SEA-AD supertype string for use in an output filename.
#'
#' Excel sheet names cannot contain "/", so layer supertypes such as
#' "L2/3 IT" are stored in the workbook as "L2 3 IT". This restores the
#' slash before applying the documented normalization rules (remove "/",
#' replace spaces with "_"), so that e.g. "L2 3 IT_1" -> "L23_IT_1" and
#' "L5 6 NP_2" -> "L56_NP_2", matching the conversion spec's examples.
#'
#' @param x Character vector of supertype strings.
#' @return Character vector of normalized, filesystem-safe supertype strings.
normalize_gabitto_celltype <- function(x) {
    x <- gsub("(L[0-9]) ([0-9])", "\\1/\\2", x)
    x <- gsub("/", "", x, fixed = TRUE)
    x <- gsub(" ", "_", x, fixed = TRUE)
    x
}

#' Parse a single Gabitto Supplementary Table 5 sheet name into canonical
#' metadata.
#'
#' Sheet-name parsing is strict: every sheet must match one of the four known
#' patterns (MTG/early/late CPS, or MTG versus_all) and must decompose
#' unambiguously into a subclass and a supertype (for CPS sheets) or a bare
#' supertype (for versus_all sheets). An unrecognized sheet name is a hard
#' error.
#'
#' For CPS sheets, the sheet-name middle is always "<subclass>_<supertype>".
#' No fixed subclass vocabulary is needed to split it: SEA-AD subclass names
#' use spaces or hyphens for internal multi-word structure (e.g. "L2 3 IT",
#' "Micro-PVM") and never contain an underscore, so the first underscore in
#' the middle is always exactly the subclass/supertype boundary, even when
#' the supertype does not repeat the subclass token (e.g.
#' "Micro-PVM_Lymphocyte", "VLMC_Pericyte_1").
#'
#' @param sheet_name Character scalar. A single Excel sheet name.
#' @return A one-row data frame with columns \code{supertype}, \code{cell_type}
#'   (normalized supertype), \code{region}, and \code{contrast}.
parse_gabitto_sheet_name <- function(sheet_name) {
    region <- "MTG"

    # (sheet-name prefix, sheet-name suffix, canonical contrast) for the
    # three Continuous Pseudo-progression Score analyses.
    cps_patterns <- list(
        list(prefix = "MTG ",   suffix = "_across_Continuous_Pseudo-progression_Score_DE.csv CPS", contrast = "ad_cps"),
        list(prefix = "early ", suffix = "_across_Continuous_Pseudo-progression_Score_DE.csv CPS", contrast = "early_ad_cps"),
        list(prefix = "late ",  suffix = "_across_Continuous_Pseudo-progression_Score_DE.csv CPS", contrast = "late_ad_cps")
    )

    # (sheet-name prefix, sheet-name suffix, canonical contrast) for the
    # versus-all-other-supertypes analysis.
    vs_all_pattern <- list(
        prefix = "MTG ", suffix = "_versus_all_DE.csv vs_all", contrast = "versus_all"
    )

    for (pat in cps_patterns) {
        if (startsWith(sheet_name, pat$prefix) && endsWith(sheet_name, pat$suffix)) {
            mid <- substr(
                sheet_name,
                nchar(pat$prefix) + 1L,
                nchar(sheet_name) - nchar(pat$suffix)
            )

            sep_pos <- regexpr("_", mid, fixed = TRUE)[[1L]]
            if (sep_pos < 0L) {
                stop(
                    "Could not find a subclass/supertype separator in Gabitto sheet name: '",
                    sheet_name, "' (middle = '", mid, "')",
                    call. = FALSE
                )
            }

            supertype <- substr(mid, sep_pos + 1L, nchar(mid))

            return(data.frame(
                supertype = supertype,
                cell_type = normalize_gabitto_celltype(supertype),
                region = region,
                contrast = pat$contrast,
                stringsAsFactors = FALSE
            ))
        }
    }

    pat <- vs_all_pattern
    if (startsWith(sheet_name, pat$prefix) && endsWith(sheet_name, pat$suffix)) {
        supertype <- substr(
            sheet_name,
            nchar(pat$prefix) + 1L,
            nchar(sheet_name) - nchar(pat$suffix)
        )

        return(data.frame(
            supertype = supertype,
            cell_type = normalize_gabitto_celltype(supertype),
            region = region,
            contrast = pat$contrast,
            stringsAsFactors = FALSE
        ))
    }

    stop("Unrecognized Gabitto Supplementary Table 5 sheet name: '", sheet_name, "'", call. = FALSE)
}

#' Construct the canonical output filename for a parsed sheet.
#'
#' @param cell_type Character scalar. Normalized cell type.
#' @param region Character scalar. Region (always "MTG" for this table).
#' @param contrast Character scalar. Canonical contrast name.
#' @return Character scalar basename, e.g. "L23_IT_1__MTG__ad_cps_DE_results.txt".
gabitto_output_filename <- function(cell_type, region, contrast) {
    paste0(cell_type, "__", region, "__", contrast, "_DE_results.txt")
}

# ---------------------------------------------------------------------------
# Single-sheet conversion
# ---------------------------------------------------------------------------

#' Build an empty (zero-count) manifest row template for a parsed sheet.
#'
#' @param sheet_name Character scalar. Source sheet name.
#' @param parsed One-row data frame as returned by \code{parse_gabitto_sheet_name()}.
#' @param output_file Character scalar. Expected/actual output path.
#' @return A one-row manifest data frame with all count fields set to zero.
empty_gabitto_manifest_row <- function(sheet_name, parsed, output_file) {
    data.frame(
        sheet_name = sheet_name,
        cell_type = parsed$cell_type,
        region = parsed$region,
        contrast = parsed$contrast,
        output_file = output_file,
        n_source_rows = NA_integer_,
        n_output_rows = 0L,
        n_blank_gene_rows_removed = 0L,
        n_exact_duplicates_removed = 0L,
        n_ambiguous_gene_rows_removed = 0L,
        n_invalid_numeric_rows_removed = 0L,
        n_invalid_se_rows_removed = 0L,
        n_invalid_p_rows_removed = 0L,
        status = "skipped",
        message = "",
        stringsAsFactors = FALSE
    )
}

#' Select and rename the adapter columns (gene, beta, se, p_value) for one
#' Gabitto sheet, given its contrast.
#'
#' @param raw A data frame/tibble as read by \code{readxl::read_excel()}.
#' @param contrast Character scalar. Canonical contrast name.
#' @return A list with \code{adapter} (data frame with columns \code{gene},
#'   \code{beta}, \code{se}, \code{p_value}) or, if required source columns
#'   are missing, \code{adapter = NULL} and \code{missing_cols} listing what
#'   was missing.
select_gabitto_adapter_columns <- function(raw, contrast) {
    gene_col <- colnames(raw)[[1L]]

    if (contrast == "versus_all") {
        beta_col <- "logFC_comparison1"
        se_col <- "se_comparison1"
        p_col <- "p_comparison1"
    } else {
        beta_col <- "logFC_Continuous_Pseudo-progression_Score"
        se_col <- "se_Continuous_Pseudo-progression_Score"
        p_col <- "p_Continuous_Pseudo-progression_Score"
    }

    required <- c(beta_col, se_col, p_col)
    missing_cols <- setdiff(required, colnames(raw))

    if (length(missing_cols) > 0L) {
        return(list(adapter = NULL, missing_cols = missing_cols))
    }

    adapter <- data.frame(
        gene = as.character(raw[[gene_col]]),
        beta = raw[[beta_col]],
        se = raw[[se_col]],
        p_value = raw[[p_col]],
        stringsAsFactors = FALSE
    )

    list(adapter = adapter, missing_cols = character(0))
}

#' Remove rows with missing, empty, or whitespace-only gene names.
#'
#' @param adapter Data frame with a \code{gene} column.
#' @param sheet_name Character scalar, used only for the warning message.
#' @return A list with the filtered \code{adapter} and \code{n_removed}.
remove_blank_gene_rows <- function(adapter, sheet_name) {
    is_blank <- is.na(adapter$gene) | trimws(adapter$gene) == ""
    n_removed <- sum(is_blank)

    if (n_removed > 0L) {
        logger::log_warn(
            "Sheet '{sheet_name}': removing {n_removed} row(s) with blank gene names."
        )
    }

    list(adapter = adapter[!is_blank, , drop = FALSE], n_removed = n_removed)
}

#' Coerce beta/se/p_value to numeric and remove rows that fail to convert
#' cleanly or that have a non-finite beta.
#'
#' @param adapter Data frame with columns \code{gene}, \code{beta}, \code{se},
#'   \code{p_value}.
#' @param sheet_name Character scalar, used only for the warning message.
#' @return A list with the filtered \code{adapter} and \code{n_removed}.
remove_invalid_numeric_rows <- function(adapter, sheet_name) {
    to_numeric <- function(x) {
        if (is.numeric(x)) return(x)
        suppressWarnings(as.numeric(x))
    }

    beta <- to_numeric(adapter$beta)
    se <- to_numeric(adapter$se)
    p_value <- to_numeric(adapter$p_value)

    is_invalid <- is.na(beta) | is.na(se) | is.na(p_value) | !is.finite(beta)
    n_removed <- sum(is_invalid)

    adapter$beta <- beta
    adapter$se <- se
    adapter$p_value <- p_value

    if (n_removed > 0L) {
        logger::log_warn(
            "Sheet '{sheet_name}': removing {n_removed} row(s) with non-numeric or non-finite beta/se/p_value."
        )
    }

    list(adapter = adapter[!is_invalid, , drop = FALSE], n_removed = n_removed)
}

#' Remove rows with invalid standard errors (NA, NaN, infinite, zero, or
#' negative).
#'
#' @param adapter Data frame with a numeric \code{se} column.
#' @param sheet_name Character scalar, used only for the warning message.
#' @return A list with the filtered \code{adapter} and \code{n_removed}.
remove_invalid_se_rows <- function(adapter, sheet_name) {
    is_invalid <- is.na(adapter$se) | !is.finite(adapter$se) | adapter$se <= 0
    n_removed <- sum(is_invalid)

    if (n_removed > 0L) {
        logger::log_warn(
            "Sheet '{sheet_name}': removing {n_removed} row(s) with invalid standard errors (NA/NaN/Inf/<=0)."
        )
    }

    list(adapter = adapter[!is_invalid, , drop = FALSE], n_removed = n_removed)
}

#' Remove rows with invalid p-values (NA, NaN, infinite, <0, or >1).
#'
#' @param adapter Data frame with a numeric \code{p_value} column.
#' @param sheet_name Character scalar, used only for the warning message.
#' @return A list with the filtered \code{adapter} and \code{n_removed}.
remove_invalid_p_rows <- function(adapter, sheet_name) {
    is_invalid <- is.na(adapter$p_value) | !is.finite(adapter$p_value) |
        adapter$p_value < 0 | adapter$p_value > 1
    n_removed <- sum(is_invalid)

    if (n_removed > 0L) {
        logger::log_warn(
            "Sheet '{sheet_name}': removing {n_removed} row(s) with invalid p-values (NA/NaN/Inf/<0/>1)."
        )
    }

    list(adapter = adapter[!is_invalid, , drop = FALSE], n_removed = n_removed)
}

#' Resolve duplicate gene rows using only the adapter fields.
#'
#' For each gene appearing more than once, compares \code{beta}, \code{se},
#' and \code{p_value} across its rows using an \code{all.equal()}-style
#' floating-point tolerant comparison. If all rows for a gene are numerically
#' equivalent, only the first occurrence (in original sheet order) is kept.
#' If any rows for a gene disagree, all rows for that gene are removed.
#'
#' @param adapter Data frame with columns \code{gene}, \code{beta}, \code{se},
#'   \code{p_value}, in original sheet order.
#' @param sheet_name Character scalar, used only for the warning message.
#' @return A list with the filtered \code{adapter} (original row order
#'   preserved), \code{n_exact_duplicates_removed}, and
#'   \code{n_ambiguous_gene_rows_removed}.
resolve_duplicate_genes <- function(adapter, sheet_name) {
    n_exact_duplicates_removed <- 0L
    n_ambiguous_gene_rows_removed <- 0L

    dup_genes <- unique(adapter$gene[duplicated(adapter$gene)])

    if (length(dup_genes) == 0L) {
        return(list(
            adapter = adapter,
            n_exact_duplicates_removed = n_exact_duplicates_removed,
            n_ambiguous_gene_rows_removed = n_ambiguous_gene_rows_removed
        ))
    }

    keep <- rep(TRUE, nrow(adapter))
    ambiguous_genes <- character(0)

    for (g in dup_genes) {
        idx <- which(adapter$gene == g)

        equivalent <- vapply(idx[-1L], function(i) {
            isTRUE(all.equal(adapter$beta[idx[1L]], adapter$beta[i], check.attributes = FALSE)) &&
                isTRUE(all.equal(adapter$se[idx[1L]], adapter$se[i], check.attributes = FALSE)) &&
                isTRUE(all.equal(adapter$p_value[idx[1L]], adapter$p_value[i], check.attributes = FALSE))
        }, logical(1))

        if (all(equivalent)) {
            keep[idx[-1L]] <- FALSE
            n_exact_duplicates_removed <- n_exact_duplicates_removed + length(idx) - 1L
        } else {
            keep[idx] <- FALSE
            n_ambiguous_gene_rows_removed <- n_ambiguous_gene_rows_removed + length(idx)
            ambiguous_genes <- c(ambiguous_genes, g)
        }
    }

    if (length(ambiguous_genes) > 0L) {
        logger::log_warn(
            "Sheet '{sheet_name}': removing {n_ambiguous_gene_rows_removed} row(s) for ",
            "{length(ambiguous_genes)} gene(s) with conflicting duplicate DE results: ",
            "{paste(ambiguous_genes, collapse = ', ')}"
        )
    }

    list(
        adapter = adapter[keep, , drop = FALSE],
        n_exact_duplicates_removed = n_exact_duplicates_removed,
        n_ambiguous_gene_rows_removed = n_ambiguous_gene_rows_removed
    )
}

#' Build the canonical seven-column voom-like table from a cleaned adapter
#' table.
#'
#' @param adapter Data frame with columns \code{gene}, \code{beta}, \code{se},
#'   \code{p_value}, one row per retained gene, in original sheet order.
#' @return A data frame with exactly the columns \code{logFC}, \code{AveExpr},
#'   \code{t}, \code{P.Value}, \code{adj.P.Val}, \code{B}, \code{z.std}, in
#'   that order, with row names set to \code{adapter$gene}.
build_canonical_de_table <- function(adapter) {
    df <- data.frame(
        logFC = adapter$beta,
        AveExpr = NA_real_,
        t = adapter$beta / adapter$se,
        P.Value = adapter$p_value,
        adj.P.Val = stats::p.adjust(adapter$p_value, method = "BH"),
        B = NA_real_,
        z.std = NA_real_,
        stringsAsFactors = FALSE
    )
    rownames(df) <- adapter$gene
    df
}

#' Convert exactly one Gabitto Supplementary Table 5 sheet into at most one
#' canonical DE output file.
#'
#' @param in_file Character scalar. Path to the Gabitto Supplementary Table 5
#'   workbook.
#' @param sheet_name Character scalar. Name of the sheet to convert.
#' @param out_dir Character scalar. Directory into which the DE output file
#'   is written.
#' @return A one-row manifest data frame; see \code{convert_gabitto_2024()}
#'   for the manifest column list.
#' @export
convert_gabitto_2024_sheet <- function(in_file, sheet_name, out_dir) {
    parsed <- parse_gabitto_sheet_name(sheet_name)
    output_file <- file.path(
        out_dir,
        gabitto_output_filename(parsed$cell_type, parsed$region, parsed$contrast)
    )

    manifest_row <- empty_gabitto_manifest_row(sheet_name, parsed, output_file)

    raw <- tryCatch(
        # The first column has a blank header in every sheet (it is the gene
        # column, handled positionally below), which makes readxl emit a
        # "New names" message on every read; that message is expected and
        # not useful once per sheet, so it is suppressed here.
        suppressMessages(readxl::read_excel(in_file, sheet = sheet_name)),
        error = function(e) {
            stop(
                "Failed to read Gabitto sheet '", sheet_name, "' from '", in_file, "': ",
                conditionMessage(e),
                call. = FALSE
            )
        }
    )

    manifest_row$n_source_rows <- nrow(raw)

    selected <- select_gabitto_adapter_columns(raw, parsed$contrast)
    if (is.null(selected$adapter)) {
        msg <- paste0(
            "Missing required source column(s): ",
            paste(selected$missing_cols, collapse = ", ")
        )
        logger::log_warn("Sheet '{sheet_name}': {msg}. Skipping.")
        manifest_row$message <- msg
        return(manifest_row)
    }

    adapter <- selected$adapter

    step1 <- remove_blank_gene_rows(adapter, sheet_name)
    adapter <- step1$adapter
    manifest_row$n_blank_gene_rows_removed <- step1$n_removed

    step2 <- remove_invalid_numeric_rows(adapter, sheet_name)
    adapter <- step2$adapter
    manifest_row$n_invalid_numeric_rows_removed <- step2$n_removed

    step3 <- remove_invalid_se_rows(adapter, sheet_name)
    adapter <- step3$adapter
    manifest_row$n_invalid_se_rows_removed <- step3$n_removed

    step4 <- remove_invalid_p_rows(adapter, sheet_name)
    adapter <- step4$adapter
    manifest_row$n_invalid_p_rows_removed <- step4$n_removed

    step5 <- resolve_duplicate_genes(adapter, sheet_name)
    adapter <- step5$adapter
    manifest_row$n_exact_duplicates_removed <- step5$n_exact_duplicates_removed
    manifest_row$n_ambiguous_gene_rows_removed <- step5$n_ambiguous_gene_rows_removed

    if (nrow(adapter) == 0L) {
        msg <- "No usable gene rows remained after filtering."
        logger::log_warn("Sheet '{sheet_name}': {msg}")
        manifest_row$message <- msg
        return(manifest_row)
    }

    de_table <- build_canonical_de_table(adapter)

    utils::write.table(
        de_table,
        output_file,
        row.names = TRUE,
        col.names = TRUE,
        quote = FALSE,
        sep = "\t"
    )

    manifest_row$n_output_rows <- nrow(de_table)
    manifest_row$status <- "written"

    manifest_row
}

# ---------------------------------------------------------------------------
# Top-level conversion
# ---------------------------------------------------------------------------

#' Convert every DE-result sheet in a Gabitto 2024 Supplementary Table 5
#' workbook into canonical voom-like DE output files.
#'
#' All sheet names are parsed and all expected output paths are computed
#' before any file is written ("preflight"), so that an unparsable sheet name
#' or an output filename collision is detected and stops the conversion
#' before a partial result set is produced.
#'
#' @param in_file Character scalar. Path to the Gabitto Supplementary Table 5
#'   workbook.
#' @param out_dir Character scalar. Directory into which DE output files and
#'   the conversion manifest are written. Created recursively if needed.
#' @return Invisibly, a data frame manifest with one row per source sheet and
#'   columns \code{sheet_name}, \code{cell_type}, \code{region},
#'   \code{contrast}, \code{output_file}, \code{n_source_rows},
#'   \code{n_output_rows}, \code{n_blank_gene_rows_removed},
#'   \code{n_exact_duplicates_removed}, \code{n_ambiguous_gene_rows_removed},
#'   \code{n_invalid_numeric_rows_removed}, \code{n_invalid_se_rows_removed},
#'   \code{n_invalid_p_rows_removed}, \code{status}, and \code{message}. The
#'   same manifest is written to \code{<out_dir>/conversion_manifest.txt}.
#' @export
convert_gabitto_2024 <- function(in_file, out_dir) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    sheet_names <- tryCatch(
        readxl::excel_sheets(in_file),
        error = function(e) {
            stop(
                "Failed to enumerate sheets in workbook '", in_file, "': ",
                conditionMessage(e),
                call. = FALSE
            )
        }
    )

    logger::log_info("Found {length(sheet_names)} sheet(s) in '{in_file}'.")

    # Preflight: parse every sheet name and compute every expected output
    # path before converting or writing anything.
    parsed_list <- lapply(sheet_names, parse_gabitto_sheet_name)

    output_files <- vapply(parsed_list, function(parsed) {
        file.path(out_dir, gabitto_output_filename(parsed$cell_type, parsed$region, parsed$contrast))
    }, character(1))

    dup_files <- unique(output_files[duplicated(output_files)])
    if (length(dup_files) > 0L) {
        for (dup in dup_files) {
            colliding_sheets <- sheet_names[output_files == dup]
            logger::log_error(
                "Output filename collision for '{dup}': sheets {paste(colliding_sheets, collapse = ', ')}"
            )
        }
        stop(
            "Output filename collision(s) detected across ", length(dup_files),
            " output path(s); see log for details. Aborting before writing any DE files.",
            call. = FALSE
        )
    }

    # Round-trip every expected output path through the existing DE filename
    # parser, to confirm the generated filenames stay compatible with
    # downstream infrastructure before any file is written.
    parsed_back <- bican.mccarroll.differentialexpression::parse_de_result_filenames(output_files)
    for (i in seq_along(parsed_list)) {
        expected <- parsed_list[[i]]
        if (!identical(parsed_back$celltype[i], expected$cell_type) ||
            !identical(parsed_back$interaction[i], expected$region) ||
            !identical(parsed_back$contrast[i], expected$contrast)) {
            stop(
                "Generated output filename for sheet '", sheet_names[i],
                "' does not round-trip through parse_de_result_filenames(): '",
                output_files[i], "'",
                call. = FALSE
            )
        }
    }

    logger::log_info("Preflight complete: {length(sheet_names)} sheet name(s) parsed, no output collisions.")

    manifest_rows <- vector("list", length(sheet_names))
    for (i in seq_along(sheet_names)) {
        manifest_rows[[i]] <- convert_gabitto_2024_sheet(in_file, sheet_names[[i]], out_dir)
    }

    manifest_columns <- c(
        "sheet_name", "cell_type", "region", "contrast", "output_file",
        "n_source_rows", "n_output_rows",
        "n_blank_gene_rows_removed", "n_exact_duplicates_removed",
        "n_ambiguous_gene_rows_removed", "n_invalid_numeric_rows_removed",
        "n_invalid_se_rows_removed", "n_invalid_p_rows_removed",
        "status", "message"
    )

    manifest <- do.call(rbind, manifest_rows)
    manifest <- manifest[, manifest_columns]

    n_written <- sum(manifest$status == "written")
    n_skipped <- sum(manifest$status == "skipped")
    logger::log_info(
        "Gabitto 2024 conversion complete: {n_written} file(s) written, {n_skipped} sheet(s) skipped."
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
# Validation (QC helper, separately callable; not required by the spec but
# useful for confirming the conversion round-trips correctly)
# ---------------------------------------------------------------------------

#' Re-check every "written" row of a conversion manifest against its source
#' sheet and its output file on disk.
#'
#' @param in_file Character scalar. Path to the Gabitto Supplementary Table 5
#'   workbook.
#' @param manifest Data frame as returned by \code{convert_gabitto_2024()}.
#' @return A data frame with one row per "written" manifest entry, containing
#'   \code{sheet_name}, \code{output_file}, \code{ok} (logical), and
#'   \code{message} describing any failure.
validate_gabitto_2024 <- function(in_file, manifest) {
    canonical_de_columns <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "z.std")

    written <- manifest[manifest$status == "written", , drop = FALSE]

    results <- vector("list", nrow(written))

    for (i in seq_len(nrow(written))) {
        row <- written[i, ]
        problems <- character(0)

        parsed_name <- tryCatch(
            bican.mccarroll.differentialexpression::parse_de_result_filenames(row$output_file),
            error = function(e) NULL
        )
        if (is.null(parsed_name)) {
            problems <- c(problems, "output filename failed to parse")
        } else {
            if (!identical(parsed_name$celltype, row$cell_type)) problems <- c(problems, "celltype mismatch")
            if (!identical(parsed_name$interaction, row$region)) problems <- c(problems, "region/interaction mismatch")
            if (!identical(parsed_name$contrast, row$contrast)) problems <- c(problems, "contrast mismatch")
        }

        out_df <- utils::read.table(
            row$output_file, sep = "\t", header = TRUE,
            row.names = 1, check.names = FALSE, stringsAsFactors = FALSE
        )

        if (!identical(colnames(out_df), canonical_de_columns)) {
            problems <- c(problems, "output columns do not match canonical schema")
        }
        if (any(is.na(rownames(out_df))) || any(trimws(rownames(out_df)) == "")) {
            problems <- c(problems, "blank gene row names present")
        }
        if (length(rownames(out_df)) != length(unique(rownames(out_df)))) {
            problems <- c(problems, "duplicate gene row names present")
        }
        if (nrow(out_df) != row$n_output_rows) {
            problems <- c(problems, "row count does not match manifest")
        }

        expected_adj <- stats::p.adjust(out_df$P.Value, method = "BH")
        if (!isTRUE(all.equal(expected_adj, out_df$adj.P.Val, check.attributes = FALSE))) {
            problems <- c(problems, "adj.P.Val does not match fresh BH adjustment")
        }

        reconstructed_se <- abs(out_df$logFC / out_df$t)
        if (any(!is.finite(reconstructed_se))) {
            problems <- c(problems, "logFC / t did not reconstruct a finite standard error for all rows")
        }

        results[[i]] <- data.frame(
            sheet_name = row$sheet_name,
            output_file = row$output_file,
            ok = length(problems) == 0L,
            message = paste(problems, collapse = "; "),
            stringsAsFactors = FALSE
        )
    }

    do.call(rbind, results)
}

# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

# manifest <- convert_gabitto_2024(inWorkbook, outDir)
# validation <- validate_gabitto_2024(inWorkbook, manifest)
# table(manifest$status)
# table(manifest$contrast)
# stopifnot(all(validation$ok))
