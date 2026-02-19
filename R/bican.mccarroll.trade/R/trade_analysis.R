#' Run TRADE on differential expression results
#'
#' This function loads differential expression results from a directory
#' produced by the differential expression pipeline, formats them for TRADEtools,
#' optionally filters by region, and computes TRADE statistics for
#' each combination of test, cell type, and region.
#'
#' Cell types may optionally be filtered upstream by providing a
#' \code{cellTypeListFile}, which is passed directly to
#' \code{bican.mccarroll.differentialexpression::parse_de_inputs}.
#'
#' Regions can also be filtered before running TRADE. If
#' \code{regions_use = NULL}, all regions present in the DE results
#' are retained.
#'
#' @param data_path Character scalar. Directory containing differential
#'   expression results.
#' @param contrast Character scalar. Pattern identifying the contrast
#'   (e.g. \code{"age"}, \code{"female_vs_male"}).
#' @param gene_to_chr_path Character scalar. Path to a two-column table
#'   mapping genes to chromosomes.
#' @param cellTypeListFile Optional character scalar. Path to a file
#'   containing one cell type per line. If \code{NULL}, all cell types
#'   are used.
#' @param regions_use Optional character vector of region names to retain.
#'   If \code{NULL}, no region filtering is applied.
#' @param include_na_region Logical. If \code{TRUE}, rows with
#'   \code{region == NA} are retained when region filtering is applied.
#'
#' @return A \code{data.table} with one row per test × cell type × region,
#'   containing TRADE summary statistics.
#' @importFrom data.table as.data.table setnames data.table
#' @export
trade_run_from_de_dir <- function(data_path,
                                  contrast,
                                  gene_to_chr_path,
                                  cellTypeListFile = NULL,
                                  regions_use = c("CaH", "Pu", "NAC", "ic", "DFC"),
                                  include_na_region = TRUE) {

    gene_to_chr <- read_gene_to_chr(gene_to_chr_path)

    de_raw <- bican.mccarroll.differentialexpression::parse_de_inputs(
        in_dir = data_path,
        file_pattern = contrast,
        cellTypeListFile = cellTypeListFile
    )

    de_dt <- format_de_results(de_raw, gene_to_chr)
    de_dt <- filter_de_by_region(de_dt, regions_use = regions_use, include_na_region = include_na_region)

    run_trade(de_dt)
}

#' Plot TRADE results as a heatmap across regions
#'
#' This function visualizes TRADE statistics across regions for each
#' cell type. The input should typically contain multiple regions;
#' heatmaps are most informative when region-resolved TRADE results
#' are provided.
#'
#' Cell types may optionally be filtered for plotting only. If
#' \code{cell_types_use = NULL}, all cell types in the input are shown.
#'
#' @param trade_results A \code{data.table} produced by
#'   \code{trade_run_from_de_dir()}.
#' @param cell_types_use Optional character vector specifying the order
#'   and subset of cell types to display.
#' @param region_order Optional character vector specifying column order
#'   for regions. If \code{NULL}, regions are displayed in their current order.
#' @param non_neuron_types Character vector identifying non-neuronal
#'   cell types. Used to reproduce the original behavior of masking
#'   the \code{"ic"} region for neuronal types.
#' @param value_var Character scalar naming the TRADE statistic to plot.
#' @param na_region_label Character scalar used to label rows where
#'   \code{region == NA}.
#'
#' @return Invisibly returns the \code{pheatmap} object.
#' @importFrom data.table copy dcast setorderv
#' @export
trade_heatmap <- function(trade_results,
                          cell_types_use = NULL,
                          region_order = NULL,
                          non_neuron_types = c("astrocyte", "OPC", "oligodendrocyte", "microglia"),
                          value_var = "trade_twi",
                          na_region_label = "all") {

    dt <- data.table::copy(trade_results)

    dt[, region_plot := region]
    dt[is.na(region_plot), region_plot := na_region_label]

    wide <- data.table::dcast(dt, cell_type ~ region_plot, value.var = value_var)

    if (!is.null(cell_types_use)) {
        wide <- wide[cell_type %in% cell_types_use]
        wide[, cell_type := factor(cell_type, levels = cell_types_use)]
        data.table::setorderv(wide, "cell_type")
    }

    # If region_order not provided, use all region columns (excluding "cell_type") in current order
    if (is.null(region_order)) {
        region_cols <- setdiff(colnames(wide), "cell_type")
        region_order <- region_cols
    }

    # Preserve the original behavior: neuron-only gets NA for "ic" (if present)
    if ("ic" %in% region_order) {
        wide[!(cell_type %in% non_neuron_types), ic := NA]
    }

    mat <- as.matrix(wide, rownames = "cell_type")
    mat <- mat[, region_order, drop = FALSE]

    pheatmap::pheatmap(
        mat,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        na_col = "white",
        color = viridisLite::mako(100)
    )
}

#' Plot TRADE results as a barplot across cell types
#'
#' This function visualizes a single TRADE statistic across cell types.
#' It is intended for results where regions have been collapsed
#' (i.e., one row per cell type).
#'
#' If multiple regions are present in \code{trade_results}, the user
#' should filter to a single region or to region-combined results
#' before calling this function.
#'
#' @param trade_results A \code{data.table} produced by
#'   \code{trade_run_from_de_dir()}.
#' @param cell_types_use Optional character vector specifying the order
#'   and subset of cell types to display.
#' @param value_var Character scalar naming the TRADE statistic to plot.
#' @param na_region_label Character scalar used to label rows where
#'   \code{region == NA}.
#'
#' @return Invisibly returns \code{NULL}. A base graphics plot is drawn.
#' @importFrom data.table copy
#' @export
trade_barplot <- function(trade_results,
                          cell_types_use = NULL,
                          value_var = "trade_twi",
                          na_region_label = "all") {

    dt <- data.table::copy(trade_results)

    dt[, region_plot := region]
    dt[is.na(region_plot), region_plot := na_region_label]

    # Barplot assumes one row per cell_type. If multiple regions exist, user should filter first.
    dup_ct <- dt[, .N, by = .(cell_type)][N > 1]
    if (nrow(dup_ct) > 0) {
        stop("trade_barplot() requires one row per cell_type (e.g. region combined only). Filter trade_results first.")
    }

    tmp <- dt[, .(cell_type, value = get(value_var))]
    tmp <- as.matrix(tmp, rownames = "cell_type")

    if (!is.null(cell_types_use)) {
        tmp <- tmp[cell_types_use[length(cell_types_use):1], , drop = FALSE]
    }

    graphics::par(mar = c(7, 10, 2, 2), mgp = c(5, 1, 0))
    graphics::barplot(
        tmp,
        col = "black",
        las = 2,
        horiz = TRUE,
        xlab = paste0("TRADE ", value_var)
    )
}

read_gene_to_chr <- function(gene_to_chr_path) {

    gene_to_chr <- data.table::fread(gene_to_chr_path)
    gene_to_chr <- gene_to_chr[chr %in% c(1:22, "X", "Y", "M")]
    gene_to_chr
}

format_de_results <- function(
        df,
        gene_to_chr,
        expected_cols = c(
            "gene", "cell_type", "region", "test",
            "log_fc", "ave_expr", "t",
            "p_value", "adj_p_val", "b", "z_std"
        )) {

    dt <- data.table::as.data.table(df)

    if (length(expected_cols) != ncol(dt)) {
        stop(
            "Number of columns in df (", ncol(dt),
            ") does not match length of expected_cols (",
            length(expected_cols), ")."
        )
    }

    data.table::setnames(dt, expected_cols)

    dt <- merge(gene_to_chr, dt, by = "gene", all.y = TRUE)

    dt[, log_fc_se := log_fc / t]

    dt
}


filter_de_by_region <- function(de_dt, regions_use = c("CaH", "Pu", "NAC", "ic", "DFC"), include_na_region = TRUE) {

    if (is.null(regions_use)) {
        return(de_dt)
    }

    if (include_na_region) {
        de_dt[is.na(region) | region %in% regions_use]
    } else {
        de_dt[region %in% regions_use]
    }
}

run_trade <- function(de_dt) {

    combinations <- unique(de_dt[, .(test, cell_type, region)])

    res_list <- vector("list", nrow(combinations))

    for (i in seq_len(nrow(combinations))) {

        test_use <- combinations[i, test]
        cell_type_use <- combinations[i, cell_type]
        region_use <- combinations[i, region]

        if (is.na(region_use)) {
            tmp <- de_dt[test == test_use & cell_type == cell_type_use & is.na(region)]
        } else {
            tmp <- de_dt[test == test_use & cell_type == cell_type_use & region == region_use]
        }

        if (nrow(tmp) == 0) {
            res_list[[i]] <- NULL
            next
        }

        tmp_a <- tmp[chr %in% 1:22]
        out_a <- TRADEtools::TRADE(
            mode = "univariate",
            results1 = tmp_a,
            results2 = NULL,
            annot_table = NULL,
            log2FoldChange = "log_fc",
            lfcSE = "log_fc_se",
            pvalue = "p_value",
            model_significant = TRUE,
            genes_exclude = NULL,
            estimate_sampling_covariance = FALSE,
            covariance_matrix_set = "combined",
            component_varexplained_threshold = 0,
            weight_nocorr = 1,
            n_sample = NULL,
            verbose = FALSE
        )

        tmp_x <- tmp[chr == "X"]
        out_x <- TRADEtools::TRADE(
            mode = "univariate",
            results1 = tmp_x,
            log2FoldChange = "log_fc",
            lfcSE = "log_fc_se",
            pvalue = "p_value"
        )

        res_list[[i]] <- data.table::data.table(
            test = test_use,
            cell_type = cell_type_use,
            region = region_use,
            trade_twi = out_a$distribution_summary$transcriptome_wide_impact,
            trade_degs = out_a$distribution_summary$Me,
            trade_twi_x = out_x$distribution_summary$transcriptome_wide_impact,
            trade_degs_x = out_x$distribution_summary$Me
        )
    }

    data.table::rbindlist(res_list, use.names = TRUE, fill = TRUE)
}

###################################
# Main function to run all TRADE analyses and generate plots, preserving the original behavior of the script
###################################
run_all_trade <- function() {

    # --------------------------------------------------------------------------
    # Hard-coded paths (edit for your environment)
    # --------------------------------------------------------------------------

    # DE results
    de_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results/LEVEL_3/sex_age/cell_type"
    de_region_subset_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results/LEVEL_3/sex_age/cell_type_subset_region"
    de_region_interaction_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results/LEVEL_3/sex_age/cell_type_region_interaction_absolute_effects"

    # Gene to chromosome table
    gene_to_chr_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/metadata/gene_to_chromosome.txt"

    # Cell types file (optional upstream filter for parse_de_inputs)
    ct_file <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/metadata/cell_types_for_de_filtering_plot.txt"

    # Plot ordering defaults (matches the original script)
    region_order <- c("CaH", "Pu", "NAC", "ic", "DFC")
    cell_types_use <- scan(ct_file, what = character(), quiet = TRUE)

    # Optional region filter for TRADE runs:
    # - Set to NULL to keep all regions present in the DE results
    # this preseves the original logic.
    regions_use <- NULL
    include_na_region <- TRUE

    # --------------------------------------------------------------------------
    # Run TRADE
    # --------------------------------------------------------------------------

    trade_age <- trade_run_from_de_dir(
        data_path = de_dir,
        contrast = "age",
        gene_to_chr_path = gene_to_chr_path,
        cellTypeListFile = ct_file,
        regions_use = regions_use,
        include_na_region = include_na_region
    )

    trade_sex <- trade_run_from_de_dir(
        data_path = de_dir,
        contrast = "female_vs_male",
        gene_to_chr_path = gene_to_chr_path,
        cellTypeListFile = ct_file,
        regions_use = regions_use,
        include_na_region = include_na_region
    )

    trade_rs_age <- trade_run_from_de_dir(
        data_path = de_region_subset_dir,
        contrast = "age",
        gene_to_chr_path = gene_to_chr_path,
        cellTypeListFile = ct_file,
        regions_use = regions_use,
        include_na_region = include_na_region
    )

    trade_rs_sex <- trade_run_from_de_dir(
        data_path = de_region_subset_dir,
        contrast = "female_vs_male",
        gene_to_chr_path = gene_to_chr_path,
        cellTypeListFile = ct_file,
        regions_use = regions_use,
        include_na_region = include_na_region
    )

    trade_ri_age <- trade_run_from_de_dir(
        data_path = de_region_interaction_dir,
        contrast = "age",
        gene_to_chr_path = gene_to_chr_path,
        cellTypeListFile = ct_file,
        regions_use = regions_use,
        include_na_region = include_na_region
    )

    trade_ri_sex <- trade_run_from_de_dir(
        data_path = de_region_interaction_dir,
        contrast = "female_vs_male",
        gene_to_chr_path = gene_to_chr_path,
        cellTypeListFile = ct_file,
        regions_use = regions_use,
        include_na_region = include_na_region
    )

    trade_list <- list(
        age = trade_age,
        sex = trade_sex,
        region_subset_age = trade_rs_age,
        region_subset_sex = trade_rs_sex,
        region_interaction_age = trade_ri_age,
        region_interaction_sex = trade_ri_sex
    )

    # --------------------------------------------------------------------------
    # Plots
    # --------------------------------------------------------------------------

    # Region-resolved runs only
    trade_heatmap(trade_list$region_subset_age, cell_types_use = cell_types_use, region_order = region_order)
    trade_heatmap(trade_list$region_subset_sex, cell_types_use = cell_types_use, region_order = region_order)

    trade_heatmap(trade_list$region_interaction_age, cell_types_use = cell_types_use, region_order = region_order)
    trade_heatmap(trade_list$region_interaction_sex, cell_types_use = cell_types_use, region_order = region_order)

    # Regions-combined runs
    trade_barplot(trade_list$age, cell_types_use = cell_types_use)
    trade_barplot(trade_list$sex, cell_types_use = cell_types_use)

    invisible(trade_list)
}
