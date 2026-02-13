library (ggplot2)
library (svglite)



#########################################
# EQTL GENE DISCOCVERY AND FILTERING PLOTS
#########################################

#' Clustered trajectories of eGene yield across eQTL filtering levels
#'
#' Generate a summary figure showing how eGene discovery
#' changes as increasingly stringent filtering is applied across multiple eQTL
#' runs. The function loads or computes a table of yields across filtering
#' levels, subsets to a configured set of cell types, clusters cell types by
#' their yield trajectories, then draws a faceted panel plot suitable for the
#' manuscript.
#'
#' To avoid repeatedly recomputing the underlying summary data, this function
#' maintains a tab-delimited cache in `data_cache_dir`. When the cache file is
#' present, it is read and used directly. When missing, the function calls
#' `bican.mccarroll.eqtl::compare_all_eQTL_runs()` to regenerate the data and
#' writes the cache before plotting.
#'
#' When `outDir` is not `NULL`, the figure is saved as
#' `eqtl_filtering_plot.svg` in `outDir`. When `outDir` is `NULL`, the plot is
#' constructed but not saved.
#'
#' @param data_dir Directory containing per-level eQTL results used by
#'   `bican.mccarroll.eqtl::compare_all_eQTL_runs()`.
#' @param outDir Output directory for the SVG. If `NULL`, no file is written.
#' @param filter_levels Integer vector of filtering levels to compare.
#' @param fdr_threshold FDR cutoff used when computing eGene yields.
#' @param num_clusters Number of clusters to use when grouping cell types by
#'   trajectory.
#' @param cellTypeListFile Path to a one-column text file listing cell types to
#'   include in the plot. Cell types not present in the data are reported via a
#'   warning.
#' @param data_cache_dir Directory used to store and read cached intermediate
#'   data for this figure.
#'
#' @return Returns `NULL` invisibly. This function is called for its side
#'   effects (plot generation and optional SVG export).
#'
#' @seealso
#'   \code{\link[bican.mccarroll.eqtl]{compare_all_eQTL_runs}}
#'   \code{\link[bican.mccarroll.eqtl]{eqtl_cluster_filtering_trajectories}}
#' @export
#' @family sample filtering figures
plot_eqtl_filtering_trajectories<-function (data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results",
                              outDir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/figure_repository",
                              filter_levels=c(0,1,2,3,4), fdr_threshold=0.05, num_clusters=4,
                              cellTypeListFile="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/metadata/cell_types_for_de_filtering_plot.txt",
                              data_cache_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/figure_repository/data_cache") {


    #if the expected cache file is present, use it instead of rebuilding the data.
    cache_file <- file.path(data_cache_dir, "eqtl_filtering_plot_cache.txt")
    if (file.exists(cache_file)) {
        logger::log_info("Using cached data from {cache_file}")
        df=read.table(cache_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    } else { # process the data as usual, write to the cache
        #clustering_min_genes and num_clusters don't do anything when outDir is NULL
        logger::log_info("No cached data from {cache_file} regenerating data from sources.  This can be extremely slow.")
        eqtl_cache_dir=paste(data_cache_dir, "eqtl_filtering_plot", sep="/")
        df= bican.mccarroll.eqtl::compare_all_eQTL_runs(data_dir=data_dir, outDir=NULL, filter_levels=filter_levels, fdr_threshold=fdr_threshold, cache_dir=eqtl_cache_dir)
        write.table(df, file=cache_file, sep="\t", row.names=FALSE, quote=FALSE)
    }


    cell_types_to_use <- read.table(cellTypeListFile, header=FALSE, stringsAsFactors=FALSE)[,1]
    df <- df[df$cell_type %in% cell_types_to_use, ]
    missing=setdiff(cell_types_to_use, df$cell_type)
    if (length(missing)>0) {
        missing_list=paste(missing, collapse=', ')
        logger::log_warn("The following cell types were requested but not found in the data: {missing_list}")
    }

    #change the baseline name and comparison name to match the DE analysis
    df$baseline_level <- as.numeric(gsub("LEVEL_", "", df$baseline_name))
    df$comparison_level <- as.numeric(gsub("LEVEL_", "", df$comparison_name))


    # Steve has asked to include level 0 - this is a temporary hack until that's approved as the
    # best way to visualize data.
    df2=.add_baseline_comparison_level_eqtl (df_filtered)

    res <- bican.mccarroll.eqtl::eqtl_cluster_filtering_trajectories(df2, value_col = "yield",
                                                                     comparison_col="comparison_level",
                                                                     K = num_clusters, title="Change in number of eGenes discovered compared to baseline")

    #make compatible with plotting code.
    clusters=res$clusters
    colnames(clusters)[1]="cell_type"

    df_plot=data.frame(
        cell_type=paste(df2$cell_type, df2$region, sep="__"),
        base_level=as.numeric(df2$baseline_level),
        comparison_level=as.numeric(df2$comparison_level),
        frac_genes_discovered=df2$yield
    )

    if (length(setdiff(df_plot$cell_type, clusters$cell_type)))
        stop("Data merge issue.")

    if (length(setdiff(clusters$cell_type, df_plot$cell_type)))
        stop("Data merge issue.")

    p<- make_filtering_cluster_panels(df_plot, clusters, legend_scale=0.8)

    if (!is.null(outDir)) {
        output_svg <- file.path(outDir, "eqtl_filtering_plot.svg")
        ggplot2::ggsave(filename = output_svg, plot = p, device = "svg", width = 8, height = 8)
    }


}

#' Example eQTL changes between two analysis levels for selected cell types
#'
#' Create a manuscript-formatted figure illustrating how eQTL results differ
#' between a baseline and comparison analysis level for a small set of cell
#' types in one region. For each requested cell type, the function draws an
#' effect-size comparison plot and a corresponding significance comparison plot,
#' stacks these into a labeled row, and then assembles all rows into a single
#' combined panel figure.
#'
#' The required per-variant summary data are collected via
#' `get_all_eqtl_filtering_example_data()`, which may use `data_cache_dir` to
#' speed up repeated calls. Plotting is performed by
#' `bican.mccarroll.eqtl::plot_pair_effects()` and
#' `bican.mccarroll.eqtl::plot_pair_pvals()`, followed by layout with `cowplot`.
#'
#' The figure is saved as `eqtl_filtering_cell_type_examples.svg` in `outDir`.
#'
#' @param cell_type_list Character vector of cell types to include.
#' @param region Region label to subset the example data.
#' @param baseline_data_dir Directory containing the baseline-level eQTL results.
#' @param comparison_data_dir Directory containing the comparison-level eQTL results.
#' @param baseline_name Label used in plot annotations for the baseline level.
#' @param comparison_name Label used in plot annotations for the comparison level.
#' @param outDir Output directory for the SVG file.
#' @param data_cache_dir Directory used by helper routines to cache intermediate
#'   data for the example panels.
#'
#' @return Returns `NULL` invisibly. This function is called for its side effect
#'   of writing an SVG file.
#'
#' @seealso
#'   \code{\link[bican.mccarroll.eqtl]{plot_pair_effects}}
#'   \code{\link[bican.mccarroll.eqtl]{plot_pair_pvals}}
#'
#' @export
#' @family sample filtering figures
plot_eqtl_filtering_examples<-function (cell_type_list=c("astrocyte", "microglia", "MSN_D1", "OPC"),
                                        region="CaH",
                                        baseline_data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/LEVEL_3",
                                        comparison_data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/LEVEL_4",
                                        baseline_name="LEVEL 3", comparison_name="LEVEL 4",
                                        outDir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/figure_repository",
                                        data_cache_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/figure_repository/data_cache") {

    #this is slightly different than other internal functions with the same names to tune
    #these particular plots.
    make_celltype_row <- function(p_left, p_right, label, strip_size = 14, gutter = 8) {
        p_left  <- p_left  + ggplot2::labs(title = NULL) +
            ggplot2::theme(plot.margin = ggplot2::margin(0, gutter, 0, 0))  # add right margin

        p_right <- p_right + ggplot2::labs(title = NULL) +
            ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, gutter))  # optional: add left margin

        panels <- cowplot::plot_grid(p_left, p_right, ncol = 2, align = "hv", axis = "tblr")

        strip <- cowplot::ggdraw() +
            cowplot::draw_label(label, fontface = "bold", size = strip_size, x = 0.5)

        cowplot::plot_grid(strip, panels, ncol = 1, rel_heights = c(0.1, 1))
    }

    df=get_all_eqtl_filtering_example_data(cell_type_list=cell_type_list, region=region, baseline_data_dir=baseline_data_dir,
                                         comparison_data_dir=comparison_data_dir, baseline_name=baseline_name,
                                         comparison_name=comparison_name, data_cache_dir=data_cache_dir)

    plot_list <- list()
    for (ct in cell_type_list) {

        d=df[df$cell_type == ct,]

        p1<-bican.mccarroll.eqtl::plot_pair_effects(
            effect_dt = d,
            cell_type_A = baseline_name,
            cell_type_B = comparison_name,
            region = region,
            annot_size=2.5
        )

        p2<-bican.mccarroll.eqtl::plot_pair_pvals(
            effect_dt = d,
            cell_type_A = baseline_name,
            cell_type_B = comparison_name,
            region = region,
            annot_size=2.5
        )

        #bind plots together
        ct_lab <- gsub("_", " ", ct, fixed = TRUE)
        plot_list[[ct]] <- make_celltype_row(p1, p2, ct_lab)
    }

    #make one big plot and save as SVG
    # If you have N plots, give each row less vertical real estate by increasing output height,
    # or explicitly control rel_heights:
    combined_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = 2)

    ggplot2::ggsave(
        filename = file.path(outDir, "eqtl_filtering_cell_type_examples.svg"),
        plot = combined_plot, device = "svg", width = 12, height = 6
    )



}

get_all_eqtl_filtering_example_data<-function (cell_type_list=c("astrocyte", "microglia", "MSN_D1", "OPC"),
                                        region="CaH",
                                        baseline_data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/LEVEL_3",
                                        comparison_data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/LEVEL_4",
                                        baseline_name="LEVEL 3", comparison_name="LEVEL 4",
                                        data_cache_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/figure_repository/data_cache") {


    df_list <- list()
    for (cell_type in cell_type_list) {
        logger::log_info("Getting data for cell type {cell_type} for eQTL filtering example plot")
        df_list[[cell_type]] <- get_eqtl_filtering_example_data(cell_type=cell_type, region=region, baseline_data_dir=baseline_data_dir,
                                                               comparison_data_dir=comparison_data_dir, baseline_name=baseline_name,
                                                                    comparison_name=comparison_name, data_cache_dir=data_cache_dir)
    }
    df_all <- do.call(rbind, df_list)

    #will need to be a data table later.
    data.table::setDT(df_all)
    return(df_all)
}

get_eqtl_filtering_example_data<-function (cell_type="astrocyte", region="CaH",
                                       baseline_data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/LEVEL_3",
                                       comparison_data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/LEVEL_4",
                                       baseline_name="LEVEL 3", comparison_name="LEVEL 4",
                                       file_separator="__", fdr_threshold=0.05,
                                       data_cache_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/figure_repository/data_cache") {


    cache_name<-paste0(paste("eqtl_filtering_example_cache_", cell_type, region, baseline_name, comparison_name, sep="_"),".txt")
    cache_name<-gsub(" ", "_", cache_name, fixed=TRUE)
    cache_file <- file.path(data_cache_dir, cache_name)

    if (file.exists(cache_file)) {
        logger::log_info("Using cached data from {cache_file}")
        df=read.table(cache_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
        return (df)
    }
    logger::log_info("No cached data file found {cache_file}")

    index_file=paste0(baseline_data_dir,"/", cell_type, file_separator, region, "/", cell_type, file_separator, region, ".cis_qtl_ann.txt.gz")
    index_file_comparison=paste0(comparison_data_dir,"/", cell_type, file_separator, region, "/", cell_type, file_separator, region, ".cis_qtl_ann.txt.gz")
    all_pairs_file=paste0(baseline_data_dir,"/", cell_type, file_separator, region, "/", cell_type, file_separator, region, ".cis_qtl_pairs.txt.gz")
    all_pairs_file_comparison=paste0(comparison_data_dir,"/", cell_type, file_separator, region, "/", cell_type, file_separator, region, ".cis_qtl_pairs.txt.gz")

    idx_A <- bican.mccarroll.eqtl::read_index_file(index_file, colsToKeep = c("gene_name", "variant_id", "slope", "qval", "pval_nominal"))
    idx_B <- bican.mccarroll.eqtl::read_index_file(index_file_comparison, colsToKeep = c("gene_name", "variant_id", "slope", "qval", "pval_nominal"))
    ap_A <- bican.mccarroll.eqtl::read_all_pairs_file(all_pairs_file, colsToKeep = c("phenotype_id", "variant_id", "slope", "pval_nominal"))
    ap_B <- bican.mccarroll.eqtl::read_all_pairs_file(all_pairs_file_comparison, colsToKeep = c("phenotype_id", "variant_id", "slope", "pval_nominal"))

    sig_A <- bican.mccarroll.eqtl::filter_significant_index(idx_A, fdr_threshold = fdr_threshold)
    sig_B <- bican.mccarroll.eqtl::filter_significant_index(idx_B, fdr_threshold = fdr_threshold)

    ref <- bican.mccarroll.eqtl::select_reference_pairs(sig_A = sig_A, sig_B = sig_B)

    eff <- bican.mccarroll.eqtl::build_pair_effect_table(
        ref_dt = ref,
        all_pairs_A = ap_A,
        all_pairs_B = ap_B,
        value_col = "slope"
    )

    pvals<-bican.mccarroll.eqtl::build_pair_effect_table(
        ref_dt = ref,
        all_pairs_A = ap_A,
        all_pairs_B = ap_B,
        value_col = "pval_nominal"
    )

    df=.merge_eqtl_effects_and_pvals(
        eff = eff,
        pvals = pvals,
        cell_type = cell_type,
        region = region,
        baseline_name = baseline_name,
        comparison_name = comparison_name
    )

    write.table(df, file=cache_file, sep="\t", row.names=FALSE, quote=FALSE)
    return (df)


}

.merge_eqtl_effects_and_pvals <- function(eff,
                                         pvals,
                                         cell_type,
                                         region,
                                         baseline_name,
                                         comparison_name) {

    key_cols <- c("gene_name", "variant_id")

    if (!all(key_cols %in% names(eff))) {
        stop("eff must contain: gene_name, variant_id")
    }
    if (!all(key_cols %in% names(pvals))) {
        stop("pvals must contain: gene_name, variant_id")
    }

    eff_dt <- data.table::as.data.table(eff)
    pval_dt <- data.table::as.data.table(pvals)

    merged <- merge(
        eff_dt,
        pval_dt,
        by = key_cols,
        all.x = TRUE,
        sort = FALSE
    )

    # Add fixed metadata columns
    merged[, cell_type := cell_type]
    merged[, region := region]
    merged[, baseline_name := baseline_name]
    merged[, comparison_name := comparison_name]

    # Put metadata first
    meta_cols <- c("cell_type", "region",
                   "baseline_name", "comparison_name",
                   key_cols)
    other_cols <- setdiff(names(merged), meta_cols)

    merged <- merged[, c(meta_cols, other_cols), with = FALSE]

    merged
}


#########################################
# DE GENE DISCOVERY AND FILTERING PLOTS
#########################################

# This produces the change in the number of DE genes discovered at each level of filtering

#' Clustered trajectories of DE gene discovery across filtering levels
#'
#' Generate a summary figure showing how the number of
#' significantly age-associated genes changes as increasingly stringent filtering
#' is applied across differential expression runs. The function loads or computes
#' a table of gene discovery counts across filtering levels, subsets to a
#' configured list of cell types, filters to cell types with sufficient baseline
#' discovery for clustering, clusters trajectories, and then draws a faceted
#' panel plot.
#'
#' To avoid repeatedly recomputing the underlying summary data, this function
#' maintains a tab-delimited cache in `data_cache_dir`. When the cache file is
#' present, it is read and used directly. When missing, the function calls
#' `bican.mccarroll.differentialexpression::compare_all_age_de_runs()` to
#' regenerate the data and writes the cache before plotting.
#'
#' When `outDir` is not `NULL`, the figure is saved as `de_filtering_plot.svg` in
#' `outDir`. When `outDir` is `NULL`, the plot is constructed but not saved.
#'
#' @param data_dir Directory containing per-level DE results used by
#'   `bican.mccarroll.differentialexpression::compare_all_age_de_runs()`.
#' @param data_name Included for interface symmetry with other figure functions;
#'   not used by the current implementation.
#' @param data_cache_dir Directory used to store and read cached intermediate
#'   data for this figure.
#' @param cellTypeListFile Path to a one-column text file listing cell types to
#'   include in the plot. Cell types not present in the data are reported via a
#'   warning.
#' @param outDir Output directory for the SVG. If `NULL`, no file is written.
#' @param clustering_min_genes Minimum baseline discovery required for a cell
#'   type to be included in clustering and plotting.
#' @param num_clusters Number of clusters to use when grouping cell types by
#'   trajectory.
#'
#' @return Returns `NULL` invisibly. This function is called for its side
#'   effects (plot generation and optional SVG export).
#'
#' @seealso
#'   \code{\link[bican.mccarroll.differentialexpression]{compare_all_age_de_runs}}
#'   \code{\link[bican.mccarroll.differentialexpression]{cluster_filtering_trajectories}}
#' @export
#' @family sample filtering figures
plot_de_filtering_trajectories<-function (
        data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results",
        data_name="donor_rxn_DGEList",
        data_cache_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/figure_repository/data_cache",
        cellTypeListFile="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/metadata/cell_types_for_de_filtering_plot.txt",
        outDir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/figure_repository",
        clustering_min_genes=150, num_clusters=3) {

    #if the expected cache file is present, use it instead of rebuilding the data.
    cache_file <- file.path(data_cache_dir, "de_filtering_plot_cache.txt")
    if (file.exists(cache_file)) {
        logger::log_info("Using cached data from {cache_file}")
        df=read.table(cache_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    } else { # process the data as usual, write to the cache
        #clustering_min_genes and num_clusters don't do anything when outDir is NULL
        logger::log_info("No cached data from {cache_file} regenerating data from sources")
        df=bican.mccarroll.differentialexpression::compare_all_age_de_runs (data_dir, outDir=NULL, filter_levels=c(0,1,2,3,4), fdr_cutoff=0.05, clustering_min_genes=clustering_min_genes, num_clusters=num_clusters)
        write.table(df, file=cache_file, sep="\t", row.names=FALSE, quote=FALSE)
    }

    #filter by the cell type list before clustering
    cell_types_to_use <- read.table(cellTypeListFile, header=FALSE, stringsAsFactors=FALSE)[,1]
    df <- df[df$cell_type %in% cell_types_to_use, ]
    missing=setdiff(cell_types_to_use, df$cell_type)
    if (length(missing)>0) {
        missing_list=paste(missing, collapse=', ')
        logger::log_warn("The following cell types were requested but not found in the data: {missing_list}")
    }

    df_filtered=df[df$num_genes_significant_old>=clustering_min_genes, ]

    # Steve has asked to include level 0 - this is a temporary hack until that's approved as the
    # best way to visualize data.
    df2=.add_baseline_comparison_level_de (df_filtered)
    res <- bican.mccarroll.differentialexpression::cluster_filtering_trajectories(df2, K = num_clusters)
    p<- make_filtering_cluster_panels(df2, res$clusters, legend_scale=0.8)

    if (!is.null(outDir)) {
        output_svg <- file.path(outDir, "de_filtering_plot.svg")
        ggplot2::ggsave(filename = output_svg, plot = p, device = "svg", width = 8, height = 8)
    }

    #I need a plot of the initial number of discoveries.
    #z=df_filtered[df_filtered$base_level==0 & df_filtered$comparison_level==1, c("cell_type", "num_genes_significant_old")]
    #p_num_de_genes<-plot_num_de_genes(z)

    #plot

}

#' Example DE changes between two analysis levels for selected cell types
#'
#' Create a figure illustrating how differential
#' expression results differ between a baseline and comparison analysis level
#' for a small set of cell types. For each requested cell type, the function
#' computes a paired comparison using
#' `bican.mccarroll.differentialexpression::compare_age_de_run()`, extracts the
#' effect-size scatter plot and the FDR scatter plot.
#'
#' The figure is saved as `de_filtering_cell_type_examples.svg` in `outDir`.
#'
#' @param cell_type_list Character vector of cell types to include.
#' @param baseline_data_dir Directory containing the baseline-level DE results.
#' @param comparison_data_dir Directory containing the comparison-level DE results.
#' @param baseline_name Label used in plot annotations for the baseline level.
#' @param comparison_name Label used in plot annotations for the comparison level.
#' @param outDir Output directory for the SVG file.
#'
#' @return Returns `NULL` invisibly. This function is called for its side effect
#'   of writing an SVG file.
#'
#' @seealso
#'   \code{\link[bican.mccarroll.differentialexpression]{compare_age_de_run}}
#' @family sample filtering figures
#' @export
plot_de_filtering_examples<-function (cell_type_list=c("astrocyte", "microglia", "MSN_D1", "glutamatergic_IT"),
                                      baseline_data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results/LEVEL_3/sex_age/cell_type",
                                      comparison_data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results/LEVEL_4/sex_age/cell_type",
                                      baseline_name="LEVEL 3", comparison_name="LEVEL 4",
                                      outDir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/figure_repository") {

    make_celltype_row <- function(p_left, p_right, label, strip_size = 14) {
        p_left  <- p_left  + ggplot2::labs(title = NULL) +
            ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
        p_right <- p_right + ggplot2::labs(title = NULL) +
            ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))

        panels <- cowplot::plot_grid(p_left, p_right, ncol = 2, align = "hv", axis = "tblr")

        strip <- cowplot::ggdraw() +
            cowplot::draw_label(label, fontface = "bold", size = strip_size, x = 0.5)

        # Small strip on top, panels below; strip gets very little height
        cowplot::plot_grid(strip, panels, ncol = 1, rel_heights = c(0.1, 1))
    }

    plot_list <- list()
    for (cell_type in cell_type_list) {
        z <- bican.mccarroll.differentialexpression::compare_age_de_run(
            cell_type = cell_type,
            old_data_dir = baseline_data_dir,
            new_data_dir = comparison_data_dir,
            baseline_name = baseline_name,
            comparison_name = comparison_name,
            fdr_cutoff = 0.05
        )

        ct_lab <- gsub("_", " ", cell_type, fixed = TRUE)
        plot_list[[cell_type]] <- make_celltype_row(z$scatter_effect, z$scatter_fdr, ct_lab)
    }

    # If you have N plots, give each row less vertical real estate by increasing output height,
    # or explicitly control rel_heights:
    combined_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = 2)

    ggplot2::ggsave(
        filename = file.path(outDir, "de_filtering_cell_type_examples.svg"),
        plot = combined_plot, device = "svg", width = 12, height = 6  # try smaller height first
    )

}

plot_num_de_genes<-function (df) {
    cell_type <- num_genes_significant_old <- NULL

    df$cell_type <- factor(
        z$cell_type,
        levels = z$cell_type[order(z$num_genes_significant_old)]
    )

    ggplot(df, aes(x = cell_type, y = num_genes_significant_old)) +
        geom_col() +
        geom_text(
            aes(label = num_genes_significant_old),
            hjust = -0.1,
            size = 3
        ) +
        scale_y_log10() +
        coord_flip() +
        labs(
            x = "cell type",
            y = "Number of age DE genes"
        ) +
        theme_bw()

}

#this plot requires:
# df columns: [cell_type, base_level, comparison_level, frac_genes_discovered]
# clusters_df columns: [cell_type, cluster]
make_filtering_cluster_panels <- function(
        df,
        clusters_df,
        ncol = 2,
        legend_position = c(0.02, 0.02),   # default: lower left
        legend_justification = c(0, 0),
        panel_titles = NULL,
        legend_scale = 1                   # NEW: scale factor for inset legend
) {
    cell_type <- comparison_level <- frac_genes_discovered <- base_level <- cluster <- NULL

    df2 <- merge(df, clusters_df, by = "cell_type", all.x = TRUE, sort = FALSE)
    df2 <- df2[df2$base_level == 0 & !is.na(df2$cluster), , drop = FALSE]
    df2$cluster <- as.factor(df2$cluster)

    y_lim <- range(df2$frac_genes_discovered, na.rm = TRUE)
    ks <- sort(unique(df2$cluster))

    # Titles
    if (is.null(panel_titles)) {
        title_map <- setNames(paste0("Cluster ", ks), ks)
    } else {
        if (!is.null(names(panel_titles)) && all(as.character(ks) %in% names(panel_titles))) {
            title_map <- panel_titles[as.character(ks)]
        } else {
            if (length(panel_titles) != length(ks)) {
                stop("panel_titles must be NULL, named for all clusters, or same length as clusters.")
            }
            title_map <- setNames(panel_titles, ks)
        }
    }

    make_one <- function(k) {
        sub <- df2[df2$cluster == k, , drop = FALSE]

        max_x <- max(sub$comparison_level, na.rm = TRUE)
        end_df <- sub[sub$comparison_level == max_x,
                      c("cell_type", "frac_genes_discovered")]
        end_df <- end_df[!is.na(end_df$frac_genes_discovered), , drop = FALSE]
        end_df <- end_df[order(-end_df$frac_genes_discovered, end_df$cell_type), , drop = FALSE]

        sub$cell_type <- factor(sub$cell_type, levels = end_df$cell_type)

        ggplot2::ggplot(
            sub,
            ggplot2::aes(
                x = comparison_level,
                y = frac_genes_discovered,
                group = cell_type,
                color = cell_type
            )
        ) +
            ggplot2::geom_hline(yintercept = 1, linewidth = 0.3) +
            ggplot2::geom_line(alpha = 0.7, linewidth = 0.6) +
            ggplot2::geom_point(alpha = 0.7, size = 1.2) +
            ggplot2::scale_x_continuous(breaks = sort(unique(df2$comparison_level))) +
            ggplot2::coord_cartesian(ylim = y_lim) +
            ggplot2::labs(
                title = title_map[[as.character(k)]],
                x = "Filtering level",
                y = "Fraction of genes discovered"
            ) +
            ggplot2::theme_bw() +
            ggplot2::theme(
                plot.title = ggplot2::element_text(hjust = 0.5),

                legend.title = ggplot2::element_blank(),
                legend.position = legend_position,
                legend.justification = legend_justification,

                legend.text = ggplot2::element_text(size = 9 * legend_scale),
                legend.key.height = ggplot2::unit(0.35 * legend_scale, "cm"),
                legend.key.width  = ggplot2::unit(0.35 * legend_scale, "cm"),
                legend.spacing.y  = ggplot2::unit(0.2 * legend_scale, "cm"),

                legend.background = ggplot2::element_rect(
                    fill = "white",
                    color = "grey70"
                )
            )
    }

    plots <- lapply(ks, make_one)
    cowplot::plot_grid(plotlist = plots, ncol = ncol, align = "hv")
}


.add_baseline_comparison_level_de <- function(df) {
    ## Expect exactly one non-zero comparison level per (cell_type, base_level)
    ## Baseline rows are copied from comparison_level == 1

    df_lvl1 <- df[df$comparison_level == 1, , drop = FALSE]

    baseline <- df_lvl1

    baseline$comparison_level <- 0

    baseline$logFC_correlation <- 1
    baseline$logFC_sign_agreement <- 1
    baseline$FDR_correlation <- 1
    baseline$frac_genes_discovered <- 1

    baseline$n_dropped_significant <- 0
    baseline$n_dropped_non_significant <- 0

    baseline$num_genes_significant_new <- baseline$num_genes_significant_old

    rbind(baseline, df)
}

.add_baseline_comparison_level_eqtl <- function(df) {
    ## Expect exactly one non-zero comparison level per (cell_type, base_level)
    ## Baseline rows are copied from comparison_level == 1

    df_lvl1 <- df[df$comparison_level == 1, , drop = FALSE]

    baseline <- df_lvl1

    baseline$comparison_level <- 0

    baseline$egene_jaccard_index <- 1
    baseline$abs_slope_cor_val <- 1
    baseline$qvalue_cor_val <- 1
    baseline$yield <- 1

    baseline$num_genes_significant_new <- baseline$num_genes_significant_old

    rbind(baseline, df)
}
