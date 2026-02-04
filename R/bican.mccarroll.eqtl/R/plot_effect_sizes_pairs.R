# Pairwise cell-type effect-size scatter plots within LEVEL_3.
# This file is intended to live in the same package as compare_filtering.R,
# so read_index_file() and read_all_pairs_file() are assumed to be available.

# -----------------------------
# Main driver (top-level)
# -----------------------------

#' Plot pairwise eQTL effect sizes across cell types within a region
#'
#' Generates scatter plots comparing eQTL effect sizes (slopes) for all unordered
#' pairs of cell types listed in \code{cell_types_file}, restricted to a single
#' \code{region}. All data are read from a single results directory.
#'
#' For a single pair of cell types (A, B), the procedure is:
#' \enumerate{
#'   \item Read each cell type's index eQTL table (\code{*.cis_qtl_ann.txt.gz})
#'         and filter to significant eGenes with \code{qval < fdr_threshold}.
#'   \item Form the union of significant genes across A and B, and select a
#'         single reference SNP per gene using the index tables:
#'         choose the most significant (smallest) \code{qval}; if tied, choose
#'         the larger \code{abs(slope)}.
#'   \item For each chosen (gene, SNP) pair, read each cell type's all-pairs table
#'         (\code{*.cis_qtl_pairs.txt.gz}) and extract the slope for that exact
#'         (gene, SNP), yielding comparable effect sizes measured at the same SNP
#'         in both cell types (no allele flipping is required).
#'   \item Plot \code{slope_A} vs \code{slope_B} with a \code{y = x} reference line
#'         and annotate the Pearson correlation.
#' }
#'
#' Plots are written to a single multi-page PDF with \code{nrow} x \code{ncol}
#' plots per page. Pair ordering follows the order of \code{cell_types_file}.
#'
#' @param data_dir Character scalar. Path to the results root directory.
#' @param region Character scalar. Region identifier used in directory and file
#'   naming (e.g., \code{"DFC"}, \code{"CaH"}).
#' @param cell_types_file Character scalar. Path to a text file containing one
#'   cell type prefix per line (e.g., \code{"GABA_CGE_DFC"}). Region is supplied
#'   separately via \code{region}.
#' @param outPDF Character scalar. Output PDF filepath.
#' @param fdr_threshold Numeric scalar. FDR/q-value threshold used to define
#'   significant eGenes in the index tables. Default is \code{0.05}.
#' @param nrow Integer scalar. Number of plot rows per PDF page. Default is \code{2}.
#' @param ncol Integer scalar. Number of plot columns per PDF page. Default is \code{2}.
#'
#' @export
plot_effect_sizes <- function(data_dir,
                              region,
                              cell_types_file,
                              outPDF,
                              fdr_threshold = 0.05,
                              nrow = 2,
                              ncol = 2) {

    cell_types <- read_cell_types_file(cell_types_file)
    if (length(cell_types) < 2) {
        stop("Need at least two cell types.")
    }

    paths <- build_paths_map(data_dir = data_dir, region = region, cell_types = cell_types)
    validate_eqtl_paths(paths)

    grDevices::pdf(outPDF, width = 11, height = 11)
    on.exit(grDevices::dev.off(), add = TRUE)

    plots_per_page <- nrow * ncol
    plot_buf <- list()

    for (pair in cell_type_pairs(cell_types)) {
        ctA <- pair[[1]]
        ctB <- pair[[2]]

        plt <- plot_one_pair(
            paths = paths,
            cell_type_A = ctA,
            cell_type_B = ctB,
            region = region,
            fdr_threshold = fdr_threshold
        )

        plot_buf <- append(plot_buf, list(plt))
        plot_buf <- plot_page(
            plots = plot_buf,
            nrow = nrow,
            ncol = ncol,
            flush_when_full = TRUE,
            plots_per_page = plots_per_page
        )
    }

    plot_page(
        plots = plot_buf,
        nrow = nrow,
        ncol = ncol,
        flush_when_full = FALSE,
        plots_per_page = plots_per_page
    )

    invisible(TRUE)
}

# -----------------------------
# One pair: load -> select refs -> join slopes -> plot
# -----------------------------
plot_one_pair <- function(paths,
                          cell_type_A,
                          cell_type_B,
                          region,
                          fdr_threshold) {

    pA <- paths[[cell_type_A]]
    pB <- paths[[cell_type_B]]

    idx_A <- read_index_file(pA$index_file, colsToKeep = c("gene_name", "variant_id", "slope", "qval"))
    idx_B <- read_index_file(pB$index_file, colsToKeep = c("gene_name", "variant_id", "slope", "qval"))
    ap_A <- read_all_pairs_file(pA$all_pairs_file, colsToKeep = c("phenotype_id", "variant_id", "slope"))
    ap_B <- read_all_pairs_file(pB$all_pairs_file, colsToKeep = c("phenotype_id", "variant_id", "slope"))

    sig_A <- filter_significant_index(idx_A, fdr_threshold = fdr_threshold)
    sig_B <- filter_significant_index(idx_B, fdr_threshold = fdr_threshold)

    ref <- select_reference_pairs(sig_A = sig_A, sig_B = sig_B)

    eff <- build_pair_effect_table(
        ref_dt = ref,
        all_pairs_A = ap_A,
        all_pairs_B = ap_B
    )

    plot_pair_effects(
        effect_dt = eff,
        cell_type_A = cell_type_A,
        cell_type_B = cell_type_B,
        region = region
    )
}

# -----------------------------
# Page writer: prints plots in an nrow x ncol grid; returns remaining buffer
# -----------------------------
plot_page <- function(plots,
                      nrow,
                      ncol,
                      flush_when_full,
                      plots_per_page) {

    if (length(plots) == 0) {
        return(plots)
    }

    if (flush_when_full && length(plots) < plots_per_page) {
        return(plots)
    }

    to_print <- if (length(plots) >= plots_per_page) {
        plots[seq_len(plots_per_page)]
    } else {
        plots
    }

    .print_plot_grid(to_print, nrow = nrow, ncol = ncol)

    if (length(plots) >= plots_per_page) {
        return(plots[-seq_len(plots_per_page)])
    }

    list()
}

# -----------------------------
# Helpers: filesystem / inputs
# -----------------------------
read_cell_types_file <- function(cell_types_file) {
    x <- readLines(cell_types_file, warn = FALSE)
    x <- trimws(x)
    x <- x[nzchar(x)]
    x
}

build_paths_map <- function(data_dir, region, cell_types) {
    paths <- lapply(cell_types, function(ct) build_eqtl_paths(data_dir = data_dir, cell_type_prefix = ct, region = region))
    names(paths) <- cell_types
    paths
}

build_eqtl_paths <- function(data_dir, cell_type_prefix, region) {
    file_separator <- "__"
    dir_name <- paste0(cell_type_prefix, file_separator, region)
    base <- paste0(data_dir, "/", dir_name, "/", dir_name)

    list(
        cell_type_prefix = cell_type_prefix,
        region = region,
        dir_name = dir_name,
        index_file = paste0(base, ".cis_qtl_ann.txt.gz"),
        all_pairs_file = paste0(base, ".cis_qtl_pairs.txt.gz")
    )
}

cell_type_pairs <- function(cell_types) {
    pairs <- list()
    k <- 0L
    n <- length(cell_types)

    if (n < 2) {
        return(pairs)
    }

    for (i in seq_len(n - 1)) {
        for (j in (i + 1):n) {
            k <- k + 1L
            pairs[[k]] <- list(cell_types[[i]], cell_types[[j]])
        }
    }

    pairs
}

# -----------------------------
# Helpers: data filtering / selection
# -----------------------------
filter_significant_index <- function(index_dt, fdr_threshold) {
    # MAKE R CMD CHECK happy (data.table NSE)
    qval <- NULL

    dt <- data.table::as.data.table(index_dt)

    required <- c("gene_name", "variant_id", "slope", "qval")
    missing_cols <- setdiff(required, names(dt))
    if (length(missing_cols) > 0) {
        stop(paste0("Index table missing required column(s): ", paste(missing_cols, collapse = ", ")))
    }

    dt[qval < fdr_threshold]
}

select_reference_pairs <- function(sig_A, sig_B) {

    # MAKE R CMD CHECK happy (data.table NSE)
    gene_name <- variant_id <- qval <- slope <- abs_slope <- source_cell <- .SD <- NULL

    sig_A_dt <- data.table::as.data.table(sig_A)
    sig_B_dt <- data.table::as.data.table(sig_B)

    keep_cols <- c("gene_name", "variant_id", "qval", "slope")
    sig_A_dt <- sig_A_dt[, intersect(keep_cols, names(sig_A_dt)), with = FALSE]
    sig_B_dt <- sig_B_dt[, intersect(keep_cols, names(sig_B_dt)), with = FALSE]

    sig_A_dt[, source_cell := "A"]
    sig_B_dt[, source_cell := "B"]

    both <- data.table::rbindlist(list(sig_A_dt, sig_B_dt), use.names = TRUE, fill = TRUE)
    if (nrow(both) == 0) {
        return(data.table::data.table(gene_name = character(0), variant_id = character(0)))
    }

    both[, abs_slope := abs(slope)]
    ref <- both[order(qval, -abs_slope), .SD[1], by = gene_name]
    ref[, abs_slope := NULL]

    ref[, c("gene_name", "variant_id", "source_cell", "qval", "slope"), with = FALSE]
}

build_pair_effect_table <- function(ref_dt, all_pairs_A, all_pairs_B) {

    ref <- data.table::as.data.table(ref_dt)[, c("gene_name", "variant_id"), with = FALSE]
    if (nrow(ref) == 0) {
        return(data.table::data.table(
            gene_name = character(0),
            variant_id = character(0),
            slope_A = numeric(0),
            slope_B = numeric(0)
        ))
    }

    apA <- data.table::as.data.table(all_pairs_A)
    apB <- data.table::as.data.table(all_pairs_B)

    need_ap <- c("phenotype_id", "variant_id", "slope")
    apA <- apA[, intersect(need_ap, names(apA)), with = FALSE]
    apB <- apB[, intersect(need_ap, names(apB)), with = FALSE]

    if ("phenotype_id" %in% names(apA)) {
        data.table::setnames(apA, "phenotype_id", "gene_name")
    }
    if ("phenotype_id" %in% names(apB)) {
        data.table::setnames(apB, "phenotype_id", "gene_name")
    }

    data.table::setnames(apA, "slope", "slope_A")
    data.table::setnames(apB, "slope", "slope_B")

    # Left joins: preserve row count of ref, fill missing slopes as NA
    m1 <- merge(ref, apA, by = c("gene_name", "variant_id"), all.x = TRUE, sort = FALSE)
    m2 <- merge(m1, apB, by = c("gene_name", "variant_id"), all.x = TRUE, sort = FALSE)

    m2
}

compute_pair_metrics <- function(effect_dt) {

    if (nrow(effect_dt) == 0) {
        return(list(
            n_total = 0L,
            n_complete = 0L,
            pearson_r = NA_real_,
            n_tested_sign = 0L,
            n_agree_sign = 0L,
            frac_agree_sign = NA_real_,
            sign_pct_txt = "NA"
        ))
    }

    n_total <- nrow(effect_dt)

    ok <- stats::complete.cases(effect_dt[, c("slope_A", "slope_B"), with = FALSE])
    n_complete <- sum(ok)

    pearson_r <- if (n_complete > 0) {
        suppressWarnings(stats::cor(effect_dt$slope_A[ok], effect_dt$slope_B[ok], method = "pearson"))
    } else {
        NA_real_
    }

    s1 <- sign(effect_dt$slope_A[ok])
    s2 <- sign(effect_dt$slope_B[ok])
    ok_sign <- (s1 != 0) & (s2 != 0)

    n_tested <- sum(ok_sign)
    n_agree <- if (n_tested > 0) sum(s1[ok_sign] == s2[ok_sign]) else 0L
    frac_agree <- if (n_tested > 0) n_agree / n_tested else NA_real_

    sign_pct_txt <- if (!is.na(frac_agree)) sprintf("%.1f%%", 100 * frac_agree) else "NA"

    list(
        n_total = as.integer(n_total),
        n_complete = as.integer(n_complete),
        pearson_r = pearson_r,
        n_tested_sign = as.integer(n_tested),
        n_agree_sign = as.integer(n_agree),
        frac_agree_sign = frac_agree,
        sign_pct_txt = sign_pct_txt
    )
}


# -----------------------------
# Helpers: plotting
# -----------------------------
plot_pair_effects <- function(effect_dt, cell_type_A, cell_type_B, region) {

    slope_A <- slope_B <- NULL

    base_empty_plot <- function(subtitle) {
        ggplot2::ggplot() +
            ggplot2::theme_bw() +
            ggplot2::labs(
                title = paste0(cell_type_A, " vs ", cell_type_B, " (", region, ")"),
                subtitle = subtitle,
                x = paste0(cell_type_A, " slope"),
                y = paste0(cell_type_B, " slope")
            )
    }

    if (nrow(effect_dt) == 0) {
        return(base_empty_plot("No genes after filtering/joining"))
    }

    metrics <- compute_pair_metrics(effect_dt)

    plot_dt <- effect_dt[
        stats::complete.cases(effect_dt[, c("slope_A", "slope_B"), with = FALSE])
    ]

    if (nrow(plot_dt) == 0) {
        return(base_empty_plot("No complete cases (all slopes missing on at least one axis)"))
    }

    r_lab <- if (is.finite(metrics$pearson_r)) sprintf("Pearson r = %.3f", metrics$pearson_r) else "Pearson r = NA"
    sign_lab <- paste0(
        "Sign agree = ", metrics$sign_pct_txt,
        " (", metrics$n_agree_sign, "/", metrics$n_tested_sign, ")"
    )

    annot_lab <- paste(
        r_lab,
        sign_lab,
        paste0("Total genes = ", metrics$n_total),
        paste0("Complete cases = ", metrics$n_complete),
        sep = "\n"
    )

    xr <- range(plot_dt$slope_A, plot_dt$slope_B, finite = TRUE)
    x_min <- xr[1]
    x_max <- xr[2]

    ggplot2::ggplot(plot_dt, ggplot2::aes(x = slope_A, y = slope_B)) +
        ggplot2::geom_point(size = 0.8, alpha = 0.7) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linewidth = 0.6) +
        ggplot2::annotate(
            geom = "text",
            x = x_min,
            y = x_max,
            hjust = 0,
            vjust = 1,
            label = annot_lab,
            size = 3.5
        ) +
        ggplot2::theme_bw() +
        ggplot2::labs(
            title = paste0(cell_type_A, " vs ", cell_type_B, " (", region, ")"),
            x = paste0(cell_type_A, " slope"),
            y = paste0(cell_type_B, " slope")
        )
}

.print_plot_grid <- function(plotlist, nrow, ncol) {
    if (length(plotlist) == 0) {
        return(invisible(TRUE))
    }

    print(cowplot::plot_grid(plotlist = plotlist, nrow = nrow, ncol = ncol))
    invisible(TRUE)
}


# -----------------------------
# Validate expected filesystem paths
# -----------------------------
validate_eqtl_paths <- function(paths) {

    # paths: named list as returned by build_paths_map()

    missing <- list()

    for (ct in names(paths)) {
        p <- paths[[ct]]

        if (!dir.exists(paste0(dirname(p$index_file)))) {
            missing[[ct]] <- c(missing[[ct]], "directory")
        }

        if (!file.exists(p$index_file)) {
            missing[[ct]] <- c(missing[[ct]], "index_file")
        }

        if (!file.exists(p$all_pairs_file)) {
            missing[[ct]] <- c(missing[[ct]], "all_pairs_file")
        }
    }

    if (length(missing) > 0) {
        msg <- c(
            "One or more required paths are missing:",
            unlist(
                lapply(
                    names(missing),
                    function(ct) {
                        paste0(
                            "  ", ct, ": ",
                            paste(missing[[ct]], collapse = ", ")
                        )
                    }
                ),
                use.names = FALSE
            )
        )
        stop(paste(msg, collapse = "\n"))
    }

    invisible(TRUE)
}

