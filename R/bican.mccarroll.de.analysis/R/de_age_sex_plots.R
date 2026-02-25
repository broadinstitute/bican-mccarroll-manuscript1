#########################
# Basic DE plots  - volcano plots, scatter plots comparing effect sizes,
# and correlation heatmaps across cell type x region groups.
######################

#' Read cell types from a text file
#'
#' Reads one cell type label per line (or whitespace-separated tokens).
#'
#' @param ct_file Path to a text file of cell types.
#' @return A character vector of cell type names.
#' @export
read_cell_types <- function(ct_file) {
    base::scan(ct_file, what = character(), quiet = TRUE)
}

#' Read gene-to-chromosome mapping
#'
#' The input file must contain at least columns "gene" and "chr".
#' Chromosomes are filtered to 1:22, X, Y, and M by default.
#'
#' @param gene_to_chr_file Path to a tabular file readable by data.table::fread.
#' @return A data.table with columns including gene and chr.
#' @export
#' @import data.table
read_gene_to_chr <- function(gene_to_chr_file) {
  chr <- gene <- NULL

  dt <- data.table::fread(gene_to_chr_file)

  keep <- c(as.character(1:22), "X", "Y", "M")
  dt <- dt[dt[["chr"]] %in% keep,]

  dt[]
}

#' Read and format differential expression results
#'
#' This reads DE results via bican.mccarroll.differentialexpression::parse_de_inputs
#' and standardizes column names, merges chromosome annotations, and computes
#' log_fc standard errors from the t-statistic.
#'
#' @param de_dir Directory containing DE results.
#' @param test Test name passed to parse_de_inputs (e.g., "age", "female_vs_male").
#' @param ct_file Path to cell types file.
#' @param gene_to_chr A data.table with at least "gene" and "chr".
#' @return A data.table with standardized DE results.
#' @export
read_de_results <- function(de_dir, test, ct_file, gene_to_chr) {
    df <- bican.mccarroll.differentialexpression::parse_de_inputs(de_dir, test, ct_file)
    prep_de(df, gene_to_chr)
}

#' Prepare DE results for plotting and downstream analysis
#' (originally format_de_results)
#' @param df A data.frame-like object returned by parse_de_inputs.
#' @param gene_to_chr A data.table with at least "gene" and "chr".
#' @return A data.table with standardized columns and annotations.
#' @export
prep_de <- function(df, gene_to_chr) {
    gene <- chr <- log_fc <- t <- log_fc_se<- NULL

    dt <- data.table::as.data.table(df)

    data.table::setnames(
        dt,
        c("gene", "cell_type", "region", "test", "log_fc", "ave_expr",
          "t", "p_value", "adj_p_val", "b", "z_std")
    )

    dt <- merge(gene_to_chr, dt, by = "gene", all.y = TRUE)
    dt[, log_fc_se := log_fc / t]

    return (dt)
}

#' Volcano plot for DE results
#'
#' @param de_dt A prepared DE data.table from prep_de/read_de_results.
#' @param cell_type_use Cell type label to plot.
#' @param region_use Region label to plot, or NA for region-combined results.
#' @param fdr_cutoff Adjusted p-value threshold.
#' @param abs_log_fc_cutoff Absolute log2 fold-change threshold.
#' @return Invisibly returns NULL.
#' @export
plot_de_volcano <- function(de_dt,
                            cell_type_use,
                            region_use,
                            fdr_cutoff = 0.05,
                            abs_log_fc_cutoff = base::log2(1.05)) {

    cell_type <- region <- log_fc <- adj_p_val <- NULL

    if (is.na(region_use)) {
        dt <- de_dt[cell_type == cell_type_use & is.na(region)]
    } else {
        dt <- de_dt[cell_type == cell_type_use & region == region_use]
    }

    rng <- dt[, base::max(base::abs(log_fc))]
    rng <- c(-rng, rng)

    graphics::par(pty = "s", xpd = FALSE)

    dt[, graphics::plot(
        log_fc,
        -base::log10(adj_p_val),
        pch = 20,
        xlim = rng,
        col = "lightgrey",
        xlab = "Effect size, log2",
        ylab = "Adjusted p-value, -log10",
        main = ""
    )]

    dt[adj_p_val < fdr_cutoff, graphics::points(
        log_fc,
        -base::log10(adj_p_val),
        pch = 20,
        col = "cornflowerblue"
    )]

    graphics::title(main = cell_type_use, adj = 0)
    graphics::abline(h = -base::log10(fdr_cutoff), lty = 2)
    graphics::abline(v = c(-abs_log_fc_cutoff, abs_log_fc_cutoff), lty = 2)

    # Make R CMD CHECK Happy
    .N <- NULL

    up <- dt[log_fc > abs_log_fc_cutoff & adj_p_val < fdr_cutoff, .N]
    down <- dt[log_fc < -abs_log_fc_cutoff & adj_p_val < fdr_cutoff, .N]

    p <- stats::binom.test(up, up + down, 0.5)$p.value

    graphics::par(xpd = NA)
    if (p < 0.05) {
        p_txt <- base::formatC(p, format = "e", digits = 1)
        graphics::legend(
            "topright",
            inset = c(0, -0.16),
            c(
                base::paste("proportion up =", base::round(up / (up + down), 2)),
                base::paste("p-value =", p_txt)
            ),
            bty = "n"
        )
    }

    invisible(NULL)
}

#' Scatter plot comparing DE effect sizes
#'
#' @param de_dt A prepared DE data.table from prep_de/read_de_results.
#' @param cell_type_a First cell type.
#' @param cell_type_b Second cell type.
#' @param region_a Region for A, or NA for region-combined results.
#' @param region_b Region for B, or NA for region-combined results.
#' @param fdr_cutoff Adjusted p-value threshold.
#' @param add_fit Whether to add a robust linear fit when rho^2 is high.
#' @return Invisibly returns NULL.
#' @export
plot_de_scatter <- function(de_dt,
                            cell_type_a,
                            cell_type_b,
                            region_a = NA,
                            region_b = NA,
                            fdr_cutoff = 0.05,
                            add_fit = TRUE) {

    cell_type <- region <- chr <- gene <- adj_p_val <- log_fc <- NULL

    if (is.na(region_a)) {
        x <- de_dt[cell_type == cell_type_a & is.na(region)]
    } else {
        x <- de_dt[cell_type == cell_type_a & region == region_a]
    }

    if (is.na(region_b)) {
        y <- de_dt[cell_type == cell_type_b & is.na(region)]
    } else {
        y <- de_dt[cell_type == cell_type_b & region == region_b]
    }

    name_a <- base::paste0(base::toupper(base::substr(cell_type_a, 1, 1)),
                           base::substr(cell_type_a, 2, base::nchar(cell_type_a)))
    name_b <- base::paste0(base::toupper(base::substr(cell_type_b, 1, 1)),
                           base::substr(cell_type_b, 2, base::nchar(cell_type_b)))

    m <- merge(x, y, by = c("chr", "gene"))

    # Make R CMD CHECK Happy
    adj_p_val.x <- adj_p_val.y <- log_fc.x <- log_fc.y <- NULL

    rng <- m[adj_p_val.x < fdr_cutoff | adj_p_val.y < fdr_cutoff,
             base::max(base::abs(c(log_fc.x, log_fc.y)))]
    rng <- c(-rng, rng)

    graphics::par(pty = "s", mar = c(6, 6, 4, 2))

    m[, graphics::plot(
        log_fc.x,
        log_fc.y,
        pch = 20,
        col = "lightgrey",
        xlim = rng,
        ylim = rng,
        xlab = c("Effect size, log2", base::paste(name_a, region_a, sep = ", ")),
        ylab = c(base::paste(name_b, region_b, sep = ", "), "Effect size, log2")
    )]

    m[(adj_p_val.x < fdr_cutoff | adj_p_val.y < fdr_cutoff),
      graphics::points(log_fc.x, log_fc.y, pch = 20, col = "cornflowerblue")]

    graphics::abline(h = 0, v = 0, lty = 2)
    graphics::abline(0, 1, lty = 2)

    ct <- m[adj_p_val.x < fdr_cutoff | adj_p_val.y < fdr_cutoff,
            stats::cor.test(log_fc.x, log_fc.y, method = "spearman")]

    rho_sqrd <- round(ct$estimate^2, 2)

    graphics::legend(
      "topleft",
      legend = bquote(rho^2 == .(rho_sqrd)),
      bty = "n"
    )

    if (rho_sqrd > 0.5 && isTRUE(add_fit)) {
        fit <- stats::lm(log_fc.y ~ log_fc.x,
                         data = m[adj_p_val.x < fdr_cutoff & adj_p_val.y < fdr_cutoff])

        coef_tab <- lmtest::coeftest(fit, vcov = sandwich::vcovHC(fit, type = "HC1"))

        b <- coef_tab["log_fc.x", "Estimate"]
        b_se <- coef_tab["log_fc.x", "Std. Error"]
        b_ci <- b + c(-1, 1) * 1.96 * b_se

        fit_color <- "tomato"

        if (b_ci[1] >= 1 & b_ci[1] <= 1) {
            fit_color <- "lightgrey"
        }

        graphics::abline(fit, lty = 2, col = fit_color)
        graphics::legend(
          "bottomright",
          legend = bquote(beta == .(base::round(b_ci[1], 2)) * " - " * .(base::round(b_ci[2], 2))),
          bty = "n"
        )
    }

    invisible(NULL)
}

#' Compute correlation matrix across cell_type x region groups
#'
#' Correlations are computed on log_fc among genes passing an FDR filter in either
#' of the compared groups. The returned value is signed rho^2.
#'
#' @param de_dt A prepared DE data.table from prep_de/read_de_results.
#' @param cell_types_use Character vector of cell types to include.
#' @param regions_use Character vector of regions to include.
#' @param non_neuron_types Character vector of non-neuronal cell types.
#' @param fdr_cutoff Adjusted p-value threshold.
#' @return A square numeric matrix with dimnames "cell_type__region".
#' @export
compute_de_cor_mat <- function(de_dt,
                               cell_types_use,
                               regions_use,
                               non_neuron_types,
                               fdr_cutoff = 0.05) {

    cell_type <- region <- cr <- gene <- adj_p_val <- log_fc <- NULL

    dt <- data.table::copy(de_dt)

    dt <- dt[cell_type %in% cell_types_use]
    dt <- dt[region %in% regions_use]

    # Exclude ic for neurons (preserved from original)
    dt <- dt[region != "ic" | cell_type %in% non_neuron_types]

    dt[, cell_type := base::factor(cell_type, levels = cell_types_use)]
    dt[, region := base::factor(region, levels = regions_use)]
    dt[, cr := base::paste(cell_type, region, sep = "__")]

    data.table::setorderv(dt, c("cell_type", "region"))

    keys <- unique(dt$cr)
    n <- base::length(keys)

    out_mat <- base::matrix(NA_real_,
                            nrow = n,
                            ncol = n,
                            dimnames = list(keys, keys))

    for (i in keys) {
        base::message(i)
        for (j in keys) {
            if (identical(i, j)) {
                out_mat[i, j] <- 1
                next
            }

            a <- dt[cr == i]
            b <- dt[cr == j]
            m <- merge(a, b, by = "gene")

            # Make R CMD CHECK Happy
            adj_p_val.x <- adj_p_val.y <- log_fc.x <- log_fc.y <- NULL

            ctest <- m[adj_p_val.x < fdr_cutoff | adj_p_val.y < fdr_cutoff,
                       stats::cor.test(log_fc.x, log_fc.y, method = "spearman")]

            r <- base::sign(ctest$estimate) * ctest$estimate^2
            out_mat[i, j] <- r
        }
    }

    out_mat
}

#' Plot a correlation heatmap with pheatmap
#'
#' @param cor_mat A numeric matrix produced by compute_de_cor_mat.
#' @param clustering_method Clustering method passed to pheatmap.
#' @param breaks Numeric breakpoints for the color scale.
#' @param palette_colors Vector of colors used for the palette.
#' @return The value returned by pheatmap::pheatmap.
#' @export
plot_de_cor_heatmap <- function(cor_mat,
                                clustering_method = "complete",
                                breaks = base::seq(-1, 1, length.out = 101),
                                palette_colors = c("steelblue", "white", "darkorange")) {

    pal_fun <- grDevices::colorRampPalette(palette_colors)

    pheatmap::pheatmap(
        cor_mat,
        breaks = breaks,
        color = pal_fun(base::length(breaks) - 1),
        clustering_method = clustering_method
    )
}

###############################################################################
# META CELL PARSING AND KMEANS CLUSTERING
###############################################################################

#' Read cell metadata
#'
#' @param cell_metadata_file Path to annotated cell metadata file.
#' @return A data.table of cell metadata.
#' @export
read_cell_metadata <- function(cell_metadata_file) {
  data.table::fread(cell_metadata_file)
}

#' Extract donor ages from cell metadata
#'
#' @param cell_metadata A data.table containing donor_external_id and age.
#' @return A named numeric vector of donor ages.
#' @export
extract_donor_ages <- function(cell_metadata) {
  donor_external_id <- age <- NULL

  donor_ages_dt <- unique(cell_metadata[, list(donor_external_id, age)])

  donor_ages <- as.numeric(donor_ages_dt$age)
  names(donor_ages) <- donor_ages_dt$donor_external_id

  sort(donor_ages)
}

#' Read and aggregate metacells
#' (originally load_metacells)
#'
#' @param path Path to metacell count matrix file.
#' @param cell_types_use Character vector of cell types to retain.
#' @param regions_use Character vector of regions to retain.
#' @return A list with elements:
#'   - metacells: dgCMatrix of aggregated TPM
#'   - col_metadata: data.table of donor/cell_type/region
#' @export
read_metacells <- function(path,
                           cell_types_use,
                           regions_use = c("CaH", "Pu", "NAC", "ic", "DFC")) {

  gene <- donor <- village <- cell_type <- region <- single_cell_assay <- group <- NULL

  dt <- data.table::fread(path)
  names(dt)[1] <- "gene"

  mat <- as.matrix(dt, rownames = "gene")

  stopifnot(requireNamespace("Matrix", quietly = TRUE))
  mat <- methods::as(mat, "dgCMatrix")

  col_metadata <- data.table::as.data.table(
    do.call(base::rbind, base::strsplit(colnames(mat), "__"))
  )
  data.table::setnames(
    col_metadata,
    c("donor", "village", "cell_type", "region", "single_cell_assay")
  )

  keep <- col_metadata[, cell_type %in% cell_types_use & region %in% regions_use]
  col_metadata <- col_metadata[keep]
  mat <- mat[, keep]

  col_metadata[, group := base::paste(donor, cell_type, region, sep = "__")]
  f <- base::factor(col_metadata$group)

  dm <- Matrix::sparseMatrix(
    i = seq_along(f),
    j = as.integer(f),
    x = 1,
    dims = c(length(f), nlevels(f)),
    dimnames = list(NULL, levels(f))
  )

  mat <- mat %*% dm

  cs <- Matrix::colSums(mat)
  d <- Matrix::Diagonal(x = 1e6 / cs)
  mat_tpm <- mat %*% d
  colnames(mat_tpm) <- colnames(mat)

  col_metadata <- data.table::as.data.table(
    do.call(base::rbind, base::strsplit(colnames(mat_tpm), "__"))
  )
  data.table::setnames(col_metadata, c("donor", "cell_type", "region"))

  list(metacells = mat_tpm, col_metadata = col_metadata)
}


#' Fast row-wise summary statistics for dense matrices
#'
#' Computes median, MAD, selected quantiles, and non-missing counts
#' for each row of a numeric matrix.
#'
#' @param mat Numeric matrix (genes x samples).
#' @return A list with elements:
#'   - median
#'   - mad
#'   - q_10
#'   - q_25
#'   - q_75
#'   - q_90
#'   - n (number of non-NA values per row)
#'
#' @keywords internal
row_stats_block_fast <- function(mat) {

  med <- matrixStats::rowMedians(mat, na.rm = TRUE)

  mad <- matrixStats::rowMads(
    mat,
    na.rm = TRUE,
    constant = 1.4826
  )

  q <- matrixStats::rowQuantiles(
    mat,
    probs = c(0.10, 0.25, 0.75, 0.90),
    na.rm = TRUE,
    type = 7
  )

  n <- base::rowSums(!is.na(mat))

  list(
    median = med,
    mad = mad,
    q_10 = q[, 1],
    q_25 = q[, 2],
    q_75 = q[, 3],
    q_90 = q[, 4],
    n = n
  )
}

#' Summarize metacells by gene and cell_type-region
#'
#' @param metacells dgCMatrix of TPM values.
#' @param col_metadata data.table returned by read_metacells.
#' @param donor_ages Named numeric vector of donor ages.
#' @return A data.table summary.
#' @export
summarize_metacells <- function(metacells,
                                col_metadata,
                                donor_ages) {

  donor <- age <- age_bin <- cell_type <- region <- cr <- genes <- NULL

  col_metadata <- data.table::copy(col_metadata)

  col_metadata[, age := unname(donor_ages[donor])]

  col_metadata[, age_bin := base::cut(
    age,
    breaks = c(-Inf, 39, 49, 59, 69, 79, 89, Inf),
    labels = c("30", "40", "50", "60", "70", "80", "90"),
    right = TRUE
  )]

  col_metadata[, cr := base::paste(cell_type, region, sep = "__")]

  groups <- split(seq_len(nrow(col_metadata)), col_metadata$cr)

  genes <- rownames(metacells)
  age_bins <- c("30", "40", "50", "60", "70", "80", "90")

  res_list <- vector("list", length(groups))
  names(res_list) <- names(groups)

  i <- 0L

  for (gn in names(groups)) {

    base::message(gn)

    i <- i + 1L
    cols <- groups[[gn]]

    sub_sp <- metacells[, cols, drop = FALSE]
    sub_dense <- base::as.matrix(sub_sp)

    st <- row_stats_block_fast(sub_dense)

    cr_parts <- data.table::tstrsplit(gn, "__", fixed = TRUE)
    ct <- cr_parts[[1]][1]
    rg <- cr_parts[[2]][1]

    dt_out <- data.table::data.table(
      gene = genes,
      cell_type = ct,
      region = rg,
      median = as.numeric(st$median),
      mad = as.numeric(st$mad),
      q_10 = as.numeric(st$q_10),
      q_25 = as.numeric(st$q_25),
      q_75 = as.numeric(st$q_75),
      q_90 = as.numeric(st$q_90),
      n_donors = as.integer(st$n)
    )

    bins_for_cols <- col_metadata$age_bin[cols]

    for (b in age_bins) {

      idx <- base::which(bins_for_cols == b)
      colname <- base::paste0("median_", b)

      if (length(idx) == 0L) {
        dt_out[[colname]] <- NA_real_
      } else if (length(idx) == 1L) {
        dt_out[[colname]] <- as.numeric(sub_dense[, idx])
      } else {
        dt_out[[colname]] <- matrixStats::rowMedians(
          sub_dense[, idx, drop = FALSE]
        )
      }
    }

    res_list[[i]] <- dt_out
  }

  #Make R CMD CHECK Happy
  gene <- cell_type <- region <- NULL

  out <- data.table::rbindlist(res_list, use.names = TRUE, fill = TRUE)
  data.table::setorder(out, gene, cell_type, region)

  out
}

#' Split metacells by cell_type and region
#'
#' @param metacells dgCMatrix of TPM values.
#' @param col_metadata data.table returned by read_metacells.
#' @param donor_ages Named numeric vector of donor ages.
#' @return A named list of dgCMatrix objects, sorted by donor age.
#' @export
split_metacells_by_cell_type_region <- function(metacells,
                                                col_metadata,
                                                donor_ages) {

  # Make R CMD CHECK Happy
  donor <- cell_type <- region <- cr <- age<- NULL

  col_metadata <- data.table::copy(col_metadata)

  col_metadata[, age := unname(donor_ages[donor])]
  col_metadata[, cr := base::paste(cell_type, region, sep = "__")]

  groups <- split(seq_len(nrow(col_metadata)), col_metadata$cr)

  out <- vector("list", length(groups))
  names(out) <- names(groups)

  i <- 0L

  for (gn in names(groups)) {

    i <- i + 1L
    cols <- groups[[gn]]

    sub_sp <- metacells[, cols, drop = FALSE]

    donors <- col_metadata$donor[cols]
    ages <- unname(donor_ages[donors])

    ord <- base::order(ages, donors)

    sub_sp <- sub_sp[, ord, drop = FALSE]
    donors <- donors[ord]

    colnames(sub_sp) <- donors

    out[[i]] <- sub_sp
  }

  out
}

#' Prepare DE matrices for clustering and heatmaps
#' (originally prepare_de_matrices)
#'
#' @param de_dt Prepared DE data.table.
#' @param metacell_summary Output from summarize_metacells.
#' @param cell_types_use Character vector of cell types.
#' @param fdr_cutoff Adjusted p-value threshold.
#' @param abs_lfc_cutoff Absolute log2 fold-change threshold.
#' @param min_tpm Minimum median TPM in either CaH or DFC.
#' @param regions_use Regions to consider (currently CaH and DFC logic retained).
#' @return A list with lfc_mat, fdr_mat, lfc_mat_z.
#' @export
prep_de_matrices <- function(de_dt,
                             metacell_summary,
                             cell_types_use,
                             fdr_cutoff = 0.01,
                             abs_lfc_cutoff = log2(1.05),
                             min_tpm = 10,
                             regions_use = c("CaH", "DFC")) {

  gene <- cell_type <- region <- median <- log_fc <- adj_p_val <- NULL

  dt <- data.table::copy(de_dt)
  dt <- dt[cell_type %in% cell_types_use]

  dt <- merge(
    dt,
    metacell_summary[region == "CaH", list(gene, cell_type, median_tpm_ca = median)],
    by = c("gene", "cell_type"),
    all.x = TRUE
  )

  dt <- merge(
    dt,
    metacell_summary[region == "DFC", list(gene, cell_type, median_tpm_dfc = median)],
    by = c("gene", "cell_type"),
    all.x = TRUE
  )

  # Make R CMD CHECK Happy
  median_tpm_ca <- median_tpm_dfc <- NULL

  dt <- dt[median_tpm_ca > min_tpm | median_tpm_dfc > min_tpm]

  lfc_mat <- as.matrix(
    data.table::dcast(dt, gene ~ cell_type, value.var = "log_fc"),
    rownames = "gene"
  )[, cell_types_use, drop = FALSE]

  fdr_mat <- as.matrix(
    data.table::dcast(dt, gene ~ cell_type, value.var = "adj_p_val"),
    rownames = "gene"
  )[, cell_types_use, drop = FALSE]

  n_sig <- base::rowSums(
    (fdr_mat < fdr_cutoff) & (abs(lfc_mat) > abs_lfc_cutoff),
    na.rm = TRUE
  )

  lfc_mat <- lfc_mat[n_sig > 0, , drop = FALSE]
  fdr_mat <- fdr_mat[n_sig > 0, , drop = FALSE]

  lfc_mat[is.na(lfc_mat)] <- 0
  fdr_mat[is.na(fdr_mat)] <- 1

  lfc_mat_z <- base::scale(lfc_mat, scale = TRUE, center = TRUE)

  list(lfc_mat = lfc_mat, fdr_mat = fdr_mat, lfc_mat_z = lfc_mat_z)
}


#' Prepare region-specific LFC matrix
#' (originally prepare_region_lfc_matrix)
#'
#' @param de_dt Prepared DE data.table.
#' @param genes_use Genes to retain.
#' @param cell_types_use Character vector of cell types.
#' @param regions_use Character vector of regions.
#' @return A numeric matrix of log fold-changes.
#' @export
prep_region_lfc_matrix <- function(de_dt,
                                   genes_use,
                                   cell_types_use,
                                   regions_use) {

  gene <- cell_type <- region <- log_fc <- NULL

  dt <- data.table::copy(de_dt)

  dt <- dt[
    cell_type %in% cell_types_use &
      region %in% regions_use
  ]

  lfc_dt <- data.table::dcast(
    dt,
    gene ~ cell_type + region,
    value.var = "log_fc",
    sep = "__"
  )

  ## IMPORTANT: avoid coercion to character by using gene as rownames
  lfc_mat <- as.matrix(lfc_dt, rownames = "gene")

  ## optional: stop early if genes are missing
  missing_genes <- setdiff(genes_use, rownames(lfc_mat))
  if (length(missing_genes) > 0L) {
    stop(
      "prep_region_lfc_matrix: ", length(missing_genes),
      " genes in genes_use are missing from the cast matrix. Example: ",
      missing_genes[[1]]
    )
  }

  lfc_mat <- lfc_mat[genes_use, , drop = FALSE]

  tmp <- data.table::as.data.table(
    do.call(base::rbind, base::strsplit(colnames(lfc_mat), "__"))
  )

  # Make R CMD CHECK Happy
  V1 <- V2 <- NULL

  tmp[, V1 := base::factor(V1, levels = cell_types_use)]
  tmp[, V2 := base::factor(V2, levels = regions_use)]

  col_order <- base::order(tmp$V1, tmp$V2)
  lfc_mat <- lfc_mat[, col_order, drop = FALSE]

  lfc_mat
}
#' Plot average silhouette width across k for k-means clustering
#' (originally plot_k_means_silhouette)
#'
#' Computes k-means clustering across a range of cluster numbers and
#' visualizes the average silhouette width for each k.
#'
#' @param mat Numeric matrix used for clustering.
#'   Rows represent the units being clustered (typically genes),
#'   and columns represent features (e.g., scaled log fold-changes).
#'   The matrix should be numeric and is typically centered and scaled
#'   prior to calling this function.
#'
#' @param ks Integer vector of cluster numbers (k) to evaluate. Default is 10:30.
#'
#' @return Invisibly returns NULL. Produces a base R plot as a side effect.
#'
#' @details
#' - Euclidean distance is computed on rows of `mat`.
#' - K-means clustering is performed for k = 10:30.
#' - For each k, clustering is run with nstart = 200 and iter.max = 20.
#' - Average silhouette width is computed using cluster::silhouette.
#' - A line plot of k versus average silhouette width is produced.
#'
#' @export
plot_kmeans_silhouette <- function(mat, ks = 10:30) {

  d <- stats::dist(mat)

  avg_sil <- base::sapply(ks, function(k) {

    set.seed(42)
    km <- stats::kmeans(mat, centers = k, nstart = 200, iter.max = 20)

    sil <- cluster::silhouette(km$cluster, d)

    base::mean(sil[, "sil_width"])
  })

  graphics::plot(
    ks,
    avg_sil,
    type = "b",
    xlab = "Number of clusters",
    ylab = "Average silhouette width"
  )

  invisible(NULL)
}

#' Plot k-means clustered heatmap of log fold-changes
#' (originally plot_k_means_heatmap)
#'
#' Performs k-means clustering on a gene-by-feature matrix and visualizes
#' the corresponding log fold-change matrix as a heatmap, ordered by cluster.
#'
#' @param k_means_mat Numeric matrix used for clustering.
#'   Rows represent genes (or other units being clustered). Rownames must be
#'   non-NULL and must correspond to gene identifiers. Typically this matrix is
#'   a scaled (e.g., z-scored) version of `lfc_mat`, since k-means is
#'   scale-sensitive.
#'
#' @param lfc_mat Numeric matrix of log fold-changes used for visualization.
#'   Rows represent genes and must include all genes present in `k_means_mat`
#'   (at least those that remain after optional cluster dropping/releveling).
#'   Columns typically represent cell types or cell_type-region combinations.
#'
#' @param scaling_factor Numeric multiplier applied to `lfc_mat` prior to
#'   plotting. This is used purely for visualization scaling and does not
#'   affect clustering.
#'
#' @param k Integer number of clusters (centers) for k-means.
#'
#' @param cluster_level_order Optional integer vector defining the desired
#'   ordering of k-means cluster labels. If NULL, cluster labels are used as
#'   returned by k-means (no reordering and no cluster dropping). If non-NULL,
#'   clusters not present in this vector will be dropped unless
#'   `allow_drop_clusters = FALSE`.
#'
#' @param allow_drop_clusters Logical; if FALSE, stop when `cluster_level_order`
#'   would drop any cluster labels present in the k-means solution.
#'
#' @return A named integer vector of cluster assignments for genes, where names
#'   correspond to gene identifiers (after any optional dropping implied by
#'   `cluster_level_order`).
#'
#' @details
#' - Euclidean k-means clustering is performed with `stats::kmeans()` using
#'   `centers = k`, `nstart = 200`, and `iter.max = 20`.
#' - If `cluster_level_order` is provided, cluster labels are reordered
#'   according to that order; optionally, clusters not included may be dropped.
#' - The heatmap is generated using `pheatmap::pheatmap()` without column or row
#'   clustering. Columns are shown in the existing order of `lfc_mat`.
#'
#' @export
plot_kmeans_heatmap <- function(k_means_mat,
                                lfc_mat,
                                scaling_factor,
                                k = 19,
                                cluster_level_order = c(2, 10, 6, 3, 5, 14, 9, 13, 4, 1, 15, 19, 18, 7, 17, 8, 11, 12),
                                allow_drop_clusters = TRUE) {

  ## -----------------------
  ## Assertions / validation
  ## -----------------------

  if (!is.matrix(k_means_mat) || !is.numeric(k_means_mat)) {
    stop("k_means_mat must be a numeric matrix.")
  }

  if (is.null(rownames(k_means_mat))) {
    stop("k_means_mat must have non-NULL rownames (gene identifiers).")
  }

  if (!is.matrix(lfc_mat) || !is.numeric(lfc_mat)) {
    stop("lfc_mat must be a numeric matrix.")
  }

  if (is.null(rownames(lfc_mat))) {
    stop("lfc_mat must have non-NULL rownames (gene identifiers).")
  }

  if (!is.numeric(scaling_factor) || length(scaling_factor) != 1L || is.na(scaling_factor)) {
    stop("scaling_factor must be a single non-NA numeric value.")
  }

  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 2) {
    stop("k must be a single numeric/integer value >= 2.")
  }
  k <- as.integer(k)

  if (!is.null(cluster_level_order)) {

    if (!is.numeric(cluster_level_order) || length(cluster_level_order) < 1L) {
      stop("cluster_level_order must be NULL or a non-empty numeric/integer vector.")
    }

    if (anyDuplicated(cluster_level_order)) {
      stop("cluster_level_order must not contain duplicate entries.")
    }
  }

  ## -----------------------
  ## K-means
  ## -----------------------

  set.seed(42)
  km <- stats::kmeans(
    k_means_mat,
    centers = k,
    nstart = 200,
    iter.max = 20
  )

  gn <- rownames(k_means_mat)

  missing_in_lfc <- setdiff(gn, rownames(lfc_mat))
  if (length(missing_in_lfc) > 0L) {
    stop(
      "lfc_mat is missing ", length(missing_in_lfc),
      " genes present in k_means_mat. Example: ",
      missing_in_lfc[[1]]
    )
  }

  present_clusters <- sort(unique(km$cluster))

  ## -----------------------
  ## Cluster labeling / ordering
  ## -----------------------

  if (is.null(cluster_level_order)) {

    ## No releveling; keep k-means labels as-is
    k_use <- as.integer(km$cluster)
    names(k_use) <- gn

  } else {

    dropped_clusters <- setdiff(present_clusters, cluster_level_order)

    if (length(dropped_clusters) > 0L && !isTRUE(allow_drop_clusters)) {
      stop(
        "cluster_level_order would drop these cluster labels present in the k-means solution: ",
        paste(dropped_clusters, collapse = ", "),
        ". Set allow_drop_clusters=TRUE to allow this."
      )
    }

    ## Relevel and optionally drop clusters not in the ordering vector
    k_use <- base::factor(km$cluster, levels = cluster_level_order)
    k_use <- as.numeric(k_use)
    names(k_use) <- gn

    keep <- !is.na(k_use)
    k_use <- k_use[keep]

    if (length(k_use) == 0L) {
      stop("After applying cluster_level_order, no genes remain (all clusters were dropped).")
    }
  }

  gene_order <- names(k_use)[base::order(k_use)]
  boundaries <- base::cumsum(base::table(k_use))

  ## -----------------------
  ## Plot
  ## -----------------------

  pheatmap::pheatmap(
    base::t(lfc_mat[gene_order, , drop = FALSE] * scaling_factor),
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    breaks = base::seq(-1, 1, length.out = 101),
    color = grDevices::colorRampPalette(c("steelblue", "white", "darkorange"))(100),
    show_colnames = FALSE,
    gaps_col = boundaries
  )

  k_use
}

#' Write lightweight DE outputs (summary + top up/down gene tables)
#'
#' Writes three files:
#' 1) A plain-text summary of inputs and DE results
#' 2) A TSV of top upregulated genes (by t-statistic)
#' 3) A TSV of top downregulated genes (by t-statistic)
#'
#' The DE table is subset to `cell_type_use`. Expression summaries are taken
#' from `metacell_summary` for the specified `cell_type_use` and `region_use`.
#'
#' @param de_dt data.table containing DE results with at least columns:
#'   gene, cell_type, log_fc, log_fc_se, t, adj_p_val.
#' @param metacell_summary data.table containing per-gene expression summaries
#'   with at least columns: gene, cell_type, region, n_donors, and
#'   median_30 ... median_90.
#' @param donor_ages Named numeric vector of donor ages (used for summary text).
#' @param cell_type_use Character scalar; cell type to extract.
#' @param region_use Character scalar; region to extract expression medians from.
#' @param out_name Basename for output files (no extension).
#' @param out_dir Output directory path.
#' @param n_top Maximum number of top up and top down genes to write.
#' @param fdr_thresh Adjusted p-value threshold.
#'
#' @return Invisibly returns a list with file paths and the two output tables.
#' @export
write_de_lite <- function(de_dt,
                          metacell_summary,
                          donor_ages,
                          cell_type_use,
                          region_use,
                          out_name,
                          out_dir,
                          n_top = 200,
                          fdr_thresh = 0.05) {

  gene <- cell_type <- region <- adj_p_val <- log_fc <- log_fc_se <- t <- median <- n_donors <- NULL
  median_30 <- median_40 <- median_50 <- median_60 <- median_70 <- median_80 <- median_90 <- NULL

  if (!dir.exists(out_dir)) {
    stop("out_dir does not exist: ", out_dir)
  }

  if (!is.numeric(n_top) || length(n_top) != 1L || is.na(n_top) || n_top < 1) {
    stop("n_top must be a positive scalar.")
  }
  n_top <- as.integer(n_top)

  if (!is.numeric(fdr_thresh) || length(fdr_thresh) != 1L || is.na(fdr_thresh) || fdr_thresh <= 0 || fdr_thresh > 1) {
    stop("fdr_thresh must be a scalar in (0, 1].")
  }

  de_dt <- data.table::copy(de_dt)
  metacell_summary <- data.table::copy(metacell_summary)

  ## subset DE to requested cell type
  de_dt <- de_dt[cell_type == cell_type_use]

  if (nrow(de_dt) == 0L) {
    stop("No DE rows found for cell_type_use='", cell_type_use, "'.")
  }

  ## medians for requested cell type + region
  medians <- metacell_summary[
    cell_type == cell_type_use & region == region_use,
    list(gene, median_30, median_40, median_50, median_60, median_70, median_80, median_90)
  ]

  if (nrow(medians) == 0L) {
    stop(
      "No metacell_summary rows found for cell_type_use='",
      cell_type_use,
      "' and region_use='",
      region_use,
      "'."
    )
  }

  ## rename medians columns:
  ## median_30 -> median_30s, etc. (preserving the collaborator’s intent)
  idx_dec <- which(!(names(medians) %in% c("gene")))
  data.table::setnames(medians, idx_dec, paste0(names(medians)[idx_dec], "s"))

  ## build up table (up)
  up_genes <- de_dt[
    adj_p_val < fdr_thresh & log_fc > 0,
    list(gene, log_fc, log_fc_se, t_stat = t, fdr = adj_p_val)
  ]

  up_genes <- merge(up_genes, medians, by = "gene", all.x = TRUE)
  data.table::setorderv(up_genes, "t_stat", -1)
  if (nrow(up_genes) > 0L) {
    up_genes <- up_genes[seq_len(min(n_top, nrow(up_genes)))]
  }

  ## build down table (down)
  down_genes <- de_dt[
    adj_p_val < fdr_thresh & log_fc < 0,
    list(gene, log_fc, log_fc_se, t_stat = t, fdr = adj_p_val)
  ]

  down_genes <- merge(down_genes, medians, by = "gene", all.x = TRUE)
  data.table::setorderv(down_genes, "t_stat", 1)
  if (nrow(down_genes) > 0L) {
    down_genes <- down_genes[seq_len(min(n_top, nrow(down_genes)))]
  }

  ## capture summary messages as text
  summary_lines <- character(0)

  add_line <- function(...) {
    summary_lines <<- c(summary_lines, paste0(...))
  }

  add_line("Input summary")
  add_line("cell type: ", cell_type_use)
  add_line("brain region: ", region_use)

  n_donors_val <- metacell_summary[
    cell_type == cell_type_use & region == region_use,
    n_donors
  ][1]

  add_line("N donors: ", n_donors_val)

  donor_ages_num <- as.numeric(donor_ages)
  donor_ages_num <- donor_ages_num[!is.na(donor_ages_num)]

  if (length(donor_ages_num) == 0L) {
    add_line("donor age range: NA - NA")
    add_line("donors median age: NA")
    add_line("donors per decade bin: NA")
  } else {
    add_line("donor age range: ", min(donor_ages_num), " - ", max(donor_ages_num))
    add_line("donors median age: ", stats::median(donor_ages_num))
    str <- convert_ages_to_decade_string(donor_ages_num)
    add_line("donors per decade bin: ", str)
  }

  # Make R CMD CHECK Happy
  .N <- NULL

  add_line("")
  add_line("Results summary")
  add_line("genes tested: ", de_dt[, .N])

  n_sig <- de_dt[adj_p_val < fdr_thresh, .N]
  n_up <- de_dt[adj_p_val < fdr_thresh & log_fc > 0, .N]
  n_dn <- de_dt[adj_p_val < fdr_thresh & log_fc < 0, .N]

  add_line("FDR < ", fdr_thresh, " total: ", n_sig)
  add_line("FDR < ", fdr_thresh, " up: ", n_up)
  add_line("FDR < ", fdr_thresh, " down: ", n_dn)

  pct_down <- (n_dn / n_sig) * 100
  add_line(
    "percent down FDR < ",
    fdr_thresh,
    ": ",
    ifelse(is.na(pct_down), "NA", round(pct_down, 3))
  )

  ## out paths
  summary_path <- base::file.path(out_dir, paste0(out_name, "_summary.txt"))
  up_path <- base::file.path(out_dir, paste0(out_name, "_up_genes.txt"))
  down_path <- base::file.path(out_dir, paste0(out_name, "_down_genes.txt"))

  summary_text <- paste(
    paste(summary_lines, collapse = "\n"),
    "",
    "Additional files",
    base::basename(up_path),
    base::basename(down_path),
    sep = "\n"
  )

  base::writeLines(summary_text, con = summary_path)

  data.table::fwrite(up_genes, file = up_path, sep = "\t", quote = FALSE, na = "NA")
  data.table::fwrite(down_genes, file = down_path, sep = "\t", quote = FALSE, na = "NA")

  invisible(list(
    summary_path = summary_path,
    up_path = up_path,
    down_path = down_path,
    up_genes = up_genes,
    down_genes = down_genes,
    summary_lines = summary_lines
  ))
}


#' Convert ages to a compact decade-bin summary string
#' (originally convert_ages_to_decade_string)
#'
#' @param donor_ages Numeric vector of donor ages.
#' @return Character scalar like "30s N = 12, 40s N = 20, ..."
#' @keywords internal
convert_ages_to_decade_string <- function(donor_ages) {

  donor_ages_chr <- as.character(donor_ages)
  dd <- substr(donor_ages_chr, 1, 1)

  ## preserve original behavior: treat 20s as 30s
  dd[dd == "2"] <- "3"

  dd <- as.numeric(dd)
  dd <- dd * 10

  tbl <- table(dd)
  paste0(names(tbl), "s N = ", tbl, collapse = ", ")
}


#' Run GSEA across cell types and GMT files
#'
#'
#' @param de_dt data.table with at least columns gene, cell_type, t.
#' @param fgsea_cell_types Character vector of cell types to run.
#' @param gmt_files Character vector of GMT file paths.
#' @param seed Random seed for reproducibility.
#' @return data.table of fgsea results with added columns cell_type, gmt, and id.
#' @export
run_gsea <- function(de_dt, fgsea_cell_types, gmt_files, seed=42) {
  set.seed(seed)
  gene <- cell_type <- t <- NULL

  gsea_results <- list()
  for (i in fgsea_cell_types) {

    message(i)

    ranks <- de_dt[cell_type == i, t] # using t-statistic to rank genes
    names(ranks) <- de_dt[cell_type == i, gene]
    ranks <- sort(ranks)

    for (j in gmt_files) {

      message(j)

      pathways <- fgsea::gmtPathways(j) # EDIT: would be faster to load all of these once; but fgsea is rate limiting

      new_results <- fgsea::fgsea(pathways = pathways,
                           stats = ranks,
                           minSize = 15,
                           maxSize = 500)

      # colnames(new_results) = tolower(colnames(new_results)) # leaving these alone for now; would be nice if they were snake_case like everything else

      new_results$cell_type <- i
      new_results$gmt <- basename(j)

      data.table::setorderv(new_results, "pval")

      gsea_results[[length(gsea_results) + 1]] <- new_results

    }

  }

  gsea_results <- data.table::rbindlist(gsea_results)

  # Make R CMD CHECK Happy
  id <- NULL
  gsea_results[, id := 1:nrow(gsea_results)]
}

#' Write lightweight GSEA outputs per cell type
#' (originally write_gsea_lite)
#'
#' Writes two TSV files per cell type:
#' - "<cell_type>_pos_gsea.txt": positive NES, padj < threshold
#' - "<cell_type>_neg_gsea.txt": negative NES, padj < threshold
#'
#' Currently filters to GMTs whose basename contains "c5.go".
#'
#' @param gsea_results data.table returned by run_gsea.
#' @param out_dir Output directory.
#' @param padj_thresh Adjusted p-value threshold.
#' @param gmt_pattern Fixed-pattern string applied to the `gmt` column.
#' @export
write_gsea_lite <- function(gsea_results,
                            out_dir,
                            padj_thresh = 0.05,
                            gmt_pattern = "c5.go") {

  cell_type <- gmt <- NES <- padj <- pathway <- size <- NULL

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed.")
  }

  if (!data.table::is.data.table(gsea_results)) {
    gsea_results <- data.table::as.data.table(gsea_results)
  }

  req_cols <- c("cell_type", "gmt", "NES", "padj", "pathway", "size")
  missing_cols <- setdiff(req_cols, names(gsea_results))
  if (length(missing_cols) > 0L) {
    stop("gsea_results is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!base::dir.exists(out_dir)) {
    stop("out_dir does not exist: ", out_dir)
  }

  if (!is.numeric(padj_thresh) || length(padj_thresh) != 1L || is.na(padj_thresh) || padj_thresh <= 0 || padj_thresh > 1) {
    stop("padj_thresh must be a scalar in (0, 1].")
  }

  for (ct in base::sort(unique(gsea_results$cell_type))) {

    base::message(ct)

    pos <- gsea_results[
      cell_type == ct &
        base::grepl(gmt_pattern, gmt, fixed = TRUE) &
        NES > 0 &
        padj < padj_thresh
    ]

    pos <- pos[, list(pathway, p_adj = padj, nes = NES, size, gmt)]

    utils::write.table(
      pos,
      file = base::paste0(out_dir, "/", ct, "_pos_gsea.txt"),
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE,
      sep = "\t"
    )

    neg <- gsea_results[
      cell_type == ct &
        base::grepl(gmt_pattern, gmt, fixed = TRUE) &
        NES < 0 &
        padj < padj_thresh
    ]

    neg <- neg[, list(pathway, p_adj = padj, nes = NES, size, gmt)]

    utils::write.table(
      neg,
      file = base::paste0(out_dir, "/", ct, "_neg_gsea.txt"),
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE,
      sep = "\t"
    )
  }

  invisible(NULL)
}

#' Donor-level GEX heatmap ordered by donor age
#' (originally plot_donor_gex_age_heatmap)
#'
#' Plots a heatmap of donor-level expression (TPM) for a gene set, with donors
#' implicitly ordered by column order in `metacells` and decade boundaries
#' indicated by gaps.
#'
#' The expression is scaled per gene using the 10th and 90th percentiles across donors:
#' values are shifted by p10, divided by (p90 - p10), and clipped to [0, 1].
#'
#' @param metacells Numeric matrix (or Matrix) with rows = genes and columns = donors.
#'   Column names must be donor IDs that match names(donor_ages).
#' @param gs Character vector of gene symbols to plot.
#' @param donor_ages Named numeric vector of donor ages (names are donor IDs).
#' @param gs_gaps Optional integer vector giving gap positions between gene groups.
#' @param cluster_gs Logical; if TRUE, cluster genes using Spearman correlation.
#' @param transpose Logical; if TRUE, plot donors as rows and genes as columns.
#'
#' @return Invisibly returns NULL.
#' @export
plot_donor_gex_age_heatmap <- function(metacells,
                                       gs,
                                       donor_ages,
                                       gs_gaps = NULL,
                                       cluster_gs = FALSE,
                                       transpose = FALSE) {

  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' is required but not installed.")
  }
  if (!requireNamespace("viridisLite", quietly = TRUE)) {
    stop("Package 'viridisLite' is required but not installed.")
  }

  if (is.null(rownames(metacells))) {
    stop("metacells must have non-NULL rownames (gene identifiers).")
  }
  if (is.null(colnames(metacells))) {
    stop("metacells must have non-NULL colnames (donor identifiers).")
  }
  if (is.null(names(donor_ages))) {
    stop("donor_ages must be a named numeric vector (names are donor IDs).")
  }

  missing_genes <- gs[!(gs %in% rownames(metacells))]
  if (length(missing_genes) > 0L) {
    message("WARNING: missing ", paste(missing_genes, collapse = ", "))
  }

  gs <- gs[gs %in% rownames(metacells)]
  if (length(gs) == 0L) {
    stop("No genes from `gs` were found in rownames(metacells).")
  }

  metacells <- metacells[gs, , drop = FALSE]

  donors <- colnames(metacells)
  if (!all(donors %in% names(donor_ages))) {
    bad <- donors[!(donors %in% names(donor_ages))]
    stop("Some metacells colnames are missing from donor_ages. Example: ", bad[[1]])
  }

  ## Scale each gene: (x - p10) / (p90 - p10), clamp to [0,1]
  ## Use apply to preserve original behavior exactly.
  p_10 <- apply(metacells, 1, stats::quantile, probs = 0.10, na.rm = TRUE)
  p_90 <- apply(metacells, 1, stats::quantile, probs = 0.90, na.rm = TRUE)

  range <- p_90 - p_10

  ## Avoid divide-by-zero: if range==0, set to 1 so the row becomes all 0 after subtraction.
  range[range == 0] <- 1

  metacells <- sweep(metacells, 1, p_10, "-")
  metacells <- sweep(metacells, 1, range, "/")

  metacells[metacells < 0] <- 0
  metacells[metacells > 1] <- 1

  ## Decade transition gaps
  donor_decades <- donor_ages[donors] %/% 10
  idx <- which(diff(donor_decades) != 0)

  if (isTRUE(cluster_gs)) {
    ## cor() needs dense; keep conversion local to this branch
    dist_mat <- stats::as.dist(1 - stats::cor(t(as.matrix(metacells)), method = "spearman"))
    cluster_arg <- stats::hclust(dist_mat, method = "average")
  } else {
    cluster_arg <- FALSE
  }

  if (isTRUE(transpose)) {
    pheatmap::pheatmap(
      t(metacells),
      cluster_rows = FALSE,
      cluster_cols = cluster_arg,
      color = viridisLite::viridis(100),
      show_rownames = FALSE,
      gaps_col = gs_gaps,
      gaps_row = idx
    )
  } else {
    pheatmap::pheatmap(
      metacells,
      cluster_rows = cluster_arg,
      cluster_cols = FALSE,
      color = viridisLite::viridis(100),
      show_colnames = FALSE,
      gaps_col = idx,
      gaps_row = gs_gaps
    )
  }

  invisible(NULL)
}


#' Donor-level GEX scatterplot vs age (ggplot2)
#' (originally plot_donor_gex_age_scatterplot)
#'
#' @param exp_vector Named numeric vector of expression values (names = donor IDs).
#' @param donor_ages Named numeric vector of donor ages (names = donor IDs).
#' @param main Character plot title.
#' @param show_spearman Logical; if TRUE, compute Spearman correlation and display it above the panel.
#' @param size If show_spearman is TRUE, dictates the size of the correlation text.
#' @param rho_threshold Numeric threshold for |rho| above which points are black; otherwise light grey.
#'   Default is 0.2.
#' @param y_axis_floor Numeric minimum upper Y-axis limit. If NULL, no Y limits are imposed.
#'   Default is 10.
#' @return A ggplot object.
#' @export
plot_donor_gex_age_scatterplot <- function(exp_vector,
                                           donor_ages,
                                           main = "",
                                           show_spearman = FALSE,
                                           size = 6,
                                           rho_threshold = 0.2,
                                           y_axis_floor = 10) {

  age <- expression <- NULL

  if (is.null(names(exp_vector))) {
    stop("exp_vector must be a named numeric vector (names are donor IDs).")
  }
  if (is.null(names(donor_ages))) {
    stop("donor_ages must be a named numeric vector (names are donor IDs).")
  }

  donors <- names(exp_vector)

  if (!all(donors %in% names(donor_ages))) {
    bad <- donors[!(donors %in% names(donor_ages))]
    stop("Some exp_vector names are missing from donor_ages. Example: ", bad[[1]])
  }

  df <- data.frame(
    age = as.numeric(donor_ages[donors]),
    expression = as.numeric(exp_vector),
    stringsAsFactors = FALSE
  )

  point_color <- "black"
  rho <- NA_real_
  subtitle_txt <- NULL

  if (isTRUE(show_spearman)) {

    ok <- stats::complete.cases(df$age, df$expression)

    if (sum(ok) >= 3L) {

      ct <- stats::cor.test(
        df$age[ok],
        df$expression[ok],
        method = "spearman",
        exact = FALSE
      )

      rho <- unname(ct$estimate)
      subtitle_txt <- paste0("rho = ", formatC(rho, format = "f", digits = 2))

      if (is.finite(rho) && abs(rho) < rho_threshold) {
        point_color <- "lightgrey"
      }

    } else {
      warning("Not enough non-missing points to compute Spearman correlation.", call. = FALSE)
    }
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = age, y = expression)) +
    ggplot2::geom_point(size = 1.0, color = point_color) +
    ggplot2::labs(
      x = "Age",
      y = "Expression, TPM",
      title = main,
      subtitle = subtitle_txt
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.subtitle = ggplot2::element_text(size = size, hjust = 0)
    )

  if (!is.null(y_axis_floor)) {
    y_max <- max(y_axis_floor, max(df$expression, na.rm = TRUE))
    p <- p + ggplot2::scale_y_continuous(limits = c(0, y_max))
  }

  p
}
