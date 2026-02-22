#' Load and format differential expression results for TRADE
#'
#' This function loads differential expression results from a directory,
#' formats them for downstream TRADE analysis, and optionally filters by region.
#' The returned table is suitable for input to \code{run_trade()}.
#'
#' Cell types may optionally be filtered upstream by providing a
#' \code{cellTypeListFile}, which is passed directly to
#' \code{bican.mccarroll.differentialexpression::parse_de_inputs}.
#'
#' Region filtering is applied only when \code{regions_use} is non-\code{NULL}.
#' If all values of the \code{region} column are \code{NA} and \code{regions_use}
#' is non-\code{NULL}, this function errors, which usually indicates that a
#' regions-combined dataset was provided when region-resolved data was expected.
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
#'
#' @return A \code{data.table} of formatted DE results suitable for TRADE.
#'
#' @importFrom data.table as.data.table setnames
#' @export
load_trade_data <- function(data_path,
                            contrast,
                            gene_to_chr_path,
                            cellTypeListFile = NULL,
                            regions_use = NULL) {

    gene_to_chr <- read_gene_to_chr(gene_to_chr_path)

    de_raw <- bican.mccarroll.differentialexpression::parse_de_inputs(
        in_dir = data_path,
        file_pattern = contrast,
        cellTypeListFile = cellTypeListFile
    )

    de_dt <- format_de_results(de_raw, gene_to_chr)

    # Unconditional call; filter_de_by_region() decides what to do.
    de_dt <- filter_de_by_region(de_dt, regions_use = regions_use)

    de_dt
}

#' Run TRADE on formatted differential expression results
#'
#' This function computes TRADE summary statistics for each combination of
#' test, cell type, and region present in the formatted DE results.
#'
#' The input is expected to be produced by \code{load_trade_data()} or an
#' equivalent table with the required columns. This function runs on the input
#' rows exactly as provided and does not perform any chromosome subsetting.
#'
#' @param de_dt A \code{data.table} of formatted DE results.
#'
#' @return A \code{data.table} with one row per test × cell type × region,
#'   containing TRADE summary statistics.
#'
#' @importFrom data.table data.table rbindlist
#' @export
run_trade <- function(de_dt) {

    #Make R CMD CHECK happy
    test <- cell_type <- region <- NULL

    combinations <- unique(de_dt[, list(test, cell_type, region)])
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

        out <- TRADEtools::TRADE(
            mode = "univariate",
            results1 = tmp,
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

        res_list[[i]] <- data.table::data.table(
            test = test_use,
            cell_type = cell_type_use,
            region = region_use,
            trade_twi = out$distribution_summary$transcriptome_wide_impact,
            trade_degs = out$distribution_summary$Me
        )
    }

    data.table::rbindlist(res_list, use.names = TRUE, fill = TRUE)
}


#' Plot TRADE results as a barplot across cell types
#'
#' @param trade_results A data.table produced by run_trade().
#' @param cell_types_use Optional character vector specifying the order and subset of cell types to display.
#' @param value_var Character scalar naming the TRADE statistic to plot.
#'
#' @return A ggplot object.
#' @export
trade_barplot <- function(trade_results,
                          cell_types_use = NULL,
                          value_var = "trade_twi") {

    # Make R CMD CHECK happy
    cell_type <- value <- N <- .N <- NULL

    dt <- data.table::as.data.table(trade_results)

    dup_ct <- dt[, .N, by = list(cell_type)][N > 1]
    if (nrow(dup_ct) > 0) {
        stop("trade_barplot(): requires one row per cell_type. Filter trade_results to a single region (or region-combined) before plotting.")
    }

    dt <- dt[, list(cell_type, value = get(value_var))]

    if (!is.null(cell_types_use)) {
        dt <- dt[cell_type %in% cell_types_use]
        dt[, cell_type := factor(cell_type, levels = rev(cell_types_use))]
    } else {
        dt <- dt[order(value)]
        dt[, cell_type := factor(cell_type, levels = cell_type)]
    }

    ggplot2::ggplot(dt, ggplot2::aes(x = value, y = cell_type)) +
        ggplot2::geom_col() +
        ggplot2::labs(x = paste0("TRADE ", value_var), y = NULL) +
        ggplot2::theme_classic()
}

#' Plot TRADE results as a heatmap across regions (ggplot2 version)
#'
#' Row ordering matches the legacy pheatmap implementation: the first row
#' of the wide table appears at the top of the heatmap.
#'
#' @param trade_results A data.table produced by run_trade().
#' @param cell_types_use Optional character vector specifying order/subset.
#' @param region_order Optional character vector specifying column order.
#' @param value_var Character scalar naming the TRADE statistic to plot.
#' @param na_region_label Label for region == NA.
#'
#' @return A ggplot object.
#'
#' @importFrom data.table copy dcast melt setorderv
#' @export
trade_heatmap<- function(trade_results,
                                  cell_types_use = NULL,
                                  region_order = NULL,
                                  value_var = "trade_twi",
                                  na_region_label = "all") {

    #Make R CMD CHECK happy
    cell_type <- region_plot <- value <- region <- NULL

    dt <- data.table::copy(trade_results)

    # Relabel NA regions
    dt[, region_plot := region]
    dt[is.na(region_plot), region_plot := na_region_label]

    # Construct wide table to define canonical ordering
    wide <- data.table::dcast(dt, cell_type ~ region_plot, value.var = value_var)

    if (!is.null(cell_types_use)) {
        wide <- wide[cell_type %in% cell_types_use]
        wide[, cell_type := factor(cell_type, levels = cell_types_use)]
        data.table::setorderv(wide, "cell_type")
    }

    if (is.null(region_order)) {
        region_cols <- setdiff(colnames(wide), "cell_type")
        region_order <- region_cols
    }

    # Melt to long format
    long <- data.table::melt(
        wide,
        id.vars = "cell_type",
        variable.name = "region_plot",
        value.name = "value",
        variable.factor = FALSE
    )

    # Match pheatmap row display (first row at top)
    ct_order <- as.character(wide$cell_type)
    long[, cell_type := factor(cell_type, levels = rev(ct_order))]

    # Column order
    long[, region_plot := factor(region_plot, levels = region_order)]

    ggplot2::ggplot(long, ggplot2::aes(x = region_plot, y = cell_type, fill = value)) +
        ggplot2::geom_tile(color = "grey70") +
        ggplot2::scale_fill_viridis_c(
            option = "mako",
            na.value = "white",
            name = NULL,
            guide = ggplot2::guide_colorbar(
                barheight = grid::unit(0.35, "npc"),
                barwidth  = grid::unit(0.025, "npc")
            )
        ) +
        ggplot2::scale_y_discrete(position = "right") +
        ggplot2::labs(x = NULL, y = NULL) +
        ggplot2::theme_classic() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            axis.ticks.y = ggplot2::element_blank(),
            axis.line = ggplot2::element_blank(),
            legend.position = "right",
            legend.justification = "top",
            legend.box.just = "top"
        )
}


read_gene_to_chr <- function(gene_to_chr_path) {
    #Make R CMD CHECK happy
    chr <- NULL

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

    #Make R CMD CHECK happy
    log_fc <- log_fc_se <- t<- NULL
    dt[, log_fc_se := log_fc / t]

    dt
}


filter_de_by_region <- function(de_dt, regions_use = NULL) {

    if (is.null(regions_use)) {
        return(de_dt)
    }

    if (!("region" %in% colnames(de_dt))) {
        stop("filter_de_by_region(): expected column 'region' was not found in de_dt.")
    }

    if (all(is.na(de_dt$region))) {
        stop(
            "Requested regions were provided (regions_use = {",
            paste(regions_use, collapse = ", "),
            "}) but all values of de_dt$region are NA. ",
            "This usually indicates that you pointed at a regions-combined DE dataset."
        )
    }

    present_regions <- unique(de_dt$region[!is.na(de_dt$region)])
    missing_regions <- setdiff(regions_use, present_regions)

    if (length(missing_regions) > 0) {
        stop(
            "Requested regions were not found in de_dt$region. Missing: {",
            paste(missing_regions, collapse = ", "),
            "}. Present: {",
            paste(present_regions, collapse = ", "),
            "}."
        )
    }

    #MAKE R CMD CHECK Happy
    region <- NULL

    de_dt[region %in% regions_use]
}

###################################
# Main function to run all TRADE analyses and generate plots, preserving the original behavior of the script
###################################
#' Run all hard-coded TRADE analyses and generate plots
#'
#' This convenience function reproduces the original workflow using hard-coded
#' paths and contrasts, but runs each dataset as a serial series of steps:
#' load -> filter -> split autosome/X -> run -> plot.
#'
#' @return Invisibly returns \code{NULL}. Plots are produced as side effects.
#' @export
run_all_trade <- function() {

    # --------------------------------------------------------------------------
    # Hard-coded paths (edit for your environment)
    # --------------------------------------------------------------------------

    de_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results/LEVEL_3/sex_age/cell_type"
    de_region_subset_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results/LEVEL_3/sex_age/cell_type_subset_region"
    de_region_interaction_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results/LEVEL_3/sex_age/cell_type_region_interaction_absolute_effects"

    gene_to_chr_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/metadata/gene_to_chromosome.txt"
    ct_file <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/metadata/cell_types_for_de_filtering_plot.txt"

    region_order <- c("CaH", "Pu", "NAC", "ic", "DFC")
    cell_types_use <- scan(ct_file, what = character(), quiet = TRUE)

    # Helper: split and run TRADE for autosomes and X
    run_trade_by_chr <- function(de_dt) {

        #Make R CMD CHECK happy
        chr <- NULL
        de_auto <- de_dt[chr %in% 1:22]

        de_x <- de_dt[chr == "X"]

        trade_auto <- run_trade(de_auto)
        trade_x <- run_trade(de_x)

        list(trade_auto = trade_auto, trade_x = trade_x)
    }

    # Helper: filter out non-neuronal cell types from the "ic" region
    # before running TRADE.
    filter_ic_to_non_neurons <- function(de_dt) {

        #Make R CMD CHECK happy
        region <- cell_type <- NULL

        non_neuron_types <- c("astrocyte", "OPC", "oligodendrocyte", "microglia")

        de_dt[!(region == "ic" & !(cell_type %in% non_neuron_types))]
    }

    # --------------------------------------------------------------------------
    # Dataset 1: regions combined, age
    # --------------------------------------------------------------------------

    de_dt <- load_trade_data(
        data_path = de_dir,
        contrast = "age",
        gene_to_chr_path = gene_to_chr_path,
        cellTypeListFile = ct_file,
        regions_use = NULL
    )
    chr_res <- run_trade_by_chr(de_dt)

    # Regions-combined runs: barplots
    trade_barplot(chr_res$trade_auto, cell_types_use = cell_types_use, value_var = "trade_twi")
    trade_barplot(chr_res$trade_x, cell_types_use = cell_types_use, value_var = "trade_twi")

    # --------------------------------------------------------------------------
    # Dataset 2: regions combined, sex
    # --------------------------------------------------------------------------

    de_dt <- load_trade_data(
        data_path = de_dir,
        contrast = "female_vs_male",
        gene_to_chr_path = gene_to_chr_path,
        cellTypeListFile = ct_file,
        regions_use = NULL
    )

    chr_res <- run_trade_by_chr(de_dt)

    trade_barplot(chr_res$trade_auto, cell_types_use = cell_types_use, value_var = "trade_twi")
    trade_barplot(chr_res$trade_x, cell_types_use = cell_types_use, value_var = "trade_twi")

    # --------------------------------------------------------------------------
    # Dataset 3: region subset, age
    # --------------------------------------------------------------------------

    de_dt <- load_trade_data(
        data_path = de_region_subset_dir,
        contrast = "age",
        gene_to_chr_path = gene_to_chr_path,
        cellTypeListFile = ct_file,
        regions_use = NULL
    )

    # Optional user-injected filtering goes here, EG: filtering out neurons from the "ic" region
    de_dt <- filter_ic_to_non_neurons(de_dt)
    chr_res <- run_trade_by_chr(de_dt)

    # Region-resolved runs: heatmaps
    trade_heatmap(chr_res$trade_auto, cell_types_use = cell_types_use, region_order = region_order, value_var = "trade_twi")
    trade_heatmap(chr_res$trade_x, cell_types_use = cell_types_use, region_order = region_order, value_var = "trade_twi")

    # --------------------------------------------------------------------------
    # Dataset 4: region subset, sex
    # --------------------------------------------------------------------------

    de_dt <- load_trade_data(
        data_path = de_region_subset_dir,
        contrast = "female_vs_male",
        gene_to_chr_path = gene_to_chr_path,
        cellTypeListFile = ct_file,
        regions_use = NULL
    )

    # Optional user-injected filtering goes here, EG: filtering out neurons from the "ic" region
    de_dt <- filter_ic_to_non_neurons(de_dt)
    chr_res <- run_trade_by_chr(de_dt)

    trade_heatmap(chr_res$trade_auto, cell_types_use = cell_types_use, region_order = region_order, value_var = "trade_twi")
    trade_heatmap(chr_res$trade_x, cell_types_use = cell_types_use, region_order = region_order, value_var = "trade_twi")

    # --------------------------------------------------------------------------
    # Dataset 5: region interaction, age
    # --------------------------------------------------------------------------

    de_dt <- load_trade_data(
        data_path = de_region_interaction_dir,
        contrast = "age",
        gene_to_chr_path = gene_to_chr_path,
        cellTypeListFile = ct_file,
        regions_use = NULL
    )

    # Optional user-injected filtering goes here, EG: filtering out neurons from the "ic" region
    de_dt <- filter_ic_to_non_neurons(de_dt)
    chr_res <- run_trade_by_chr(de_dt)

    trade_heatmap(chr_res$trade_auto, cell_types_use = cell_types_use, region_order = region_order, value_var = "trade_twi")
    trade_heatmap(chr_res$trade_x, cell_types_use = cell_types_use, region_order = region_order, value_var = "trade_twi")

    # --------------------------------------------------------------------------
    # Dataset 6: region interaction, sex
    # --------------------------------------------------------------------------

    de_dt <- load_trade_data(
        data_path = de_region_interaction_dir,
        contrast = "female_vs_male",
        gene_to_chr_path = gene_to_chr_path,
        cellTypeListFile = ct_file,
        regions_use = NULL
    )

    # Optional user-injected filtering goes here, EG: filtering out neurons from the "ic" region
    de_dt <- filter_ic_to_non_neurons(de_dt)
    chr_res <- run_trade_by_chr(de_dt)

    trade_heatmap(chr_res$trade_auto, cell_types_use = cell_types_use, region_order = region_order, value_var = "trade_twi")
    trade_heatmap(chr_res$trade_x, cell_types_use = cell_types_use, region_order = region_order, value_var = "trade_twi")

    invisible(NULL)
}
