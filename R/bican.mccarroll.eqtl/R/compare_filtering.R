# library (data.table)
# library (ggplot2)
# library (ggvenn)
# library (cowplot)
# library (digest)
# library (logger)

# cell_type="MSN";region="CaH"
# index_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/LEVEL_0/MSN__CaH/MSN__CaH.cis_qtl_ann.txt.gz"
# index_file_comparison="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/LEVEL_2/MSN__CaH/MSN__CaH.cis_qtl_ann.txt.gz"
# all_pairs_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/LEVEL_0/MSN__CaH/MSN__CaH.cis_qtl_pairs.txt.gz"
# all_pairs_file_comparison="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/LEVEL_2/MSN__CaH/MSN__CaH.cis_qtl_pairs.txt.gz"

# baseline_name="LEVEL_0"
# comparison_name="LEVEL_2"


# For the full pipeline run.
#cache_dir="/downloads/eqtl_cache"
#fdr_threshold=0.05;filter_levels=c(0,1,2,3)
#outDir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/compare_filtering"

#Looking at a pair of eQTL results at different filter levels, comparisons:
#1. Number of eQTLs found in the baseline and comparison levels at FDR<0.05.
    # How many eQTLs are found in both levels (FDR<0.05)
    # Recall = both levels FDR < 0.05 / baseline level FDR <0.05
    # Yield - (filtered level FDR <0.05) / (baseline level FDR <0.05) [are there more eQLTs found after filtering?]
#2. Effect size correlation between levels for SNP/genes that are found at both levels (FDR<0.05)
#3. Sign test agreement between levels for SNP/genes that are found in the baseline level (baseline FDR <0.05)
#4. Empiric pvalue comparison between levels for SNP/genes that are found at the baseline level (baseline FDR <0.05)
#5. Lambda GC (or equivalent)?




#compare all levels to each other.

#compare_all_eQTL_runs(data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results", outDir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/compare_filtering", filter_levels=c(0,1,2,3), fdr_threshold=0.05, cache_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/compare_filtering/cache")

#' Compare eQTL results across multiple filtering levels
#'
#' Iterates over a set of filtering levels (e.g., LEVEL_0, LEVEL_1, ...) and
#' performs pairwise comparisons between the first level (baseline) and each
#' subsequent level. For each comparison, summary statistics are collected and
#' written to a tab-delimited file in \code{outDir}.
#'
#' @param data_dir Character scalar. Base directory containing \code{LEVEL_<n>}
#'   subdirectories.
#' @param outDir Character scalar. Output directory for per-comparison PDFs and
#'   summary statistics files.
#' @param filter_levels Integer vector. Filtering levels to compare. The first
#'   element is used as the baseline level; each subsequent element is compared
#'   against the baseline.
#' @param fdr_threshold Numeric scalar. FDR threshold used by downstream
#'   comparison routines (default 0.05).
#' @param cache_dir Character scalar. Cache directory used by downstream routines
#'   to avoid re-reading large files.
#'
#' @return Invisibly returns a data.frame of per-dataset summary statistics
#'   (one row per cell_type/region per comparison). Also writes a combined
#'   summary file \code{compare_eQTL_all_levels_stats.txt} to \code{outDir}.
#'
#' @importFrom logger log_info
#' @importFrom utils write.table
#' @export
compare_all_eQTL_runs<-function (data_dir, outDir, filter_levels=c(0,1,2,3), fdr_threshold=0.05, cache_dir) {
    base_level=filter_levels[1]
    results=list()
    for (i in 1:(length(filter_levels)-1)) {
        comparison_level=filter_levels[i+1]
        baseline_data_dir=paste(data_dir,"/LEVEL_", base_level, sep="")
        comparison_data_dir=paste(data_dir,"/LEVEL_", comparison_level, sep="")
        logger::log_info(paste0("Comparing eQTL results between LEVEL ", base_level, " and LEVEL ", comparison_level, "\n"))
        outPDF=paste(outDir, "/compare_eQTL_LEVEL_", base_level, "_vs_LEVEL_", comparison_level, ".pdf", sep="")
        outFile=paste(outDir, "/compare_eQTL_LEVEL_", base_level, "_vs_LEVEL_", comparison_level, ".txt", sep="")
        z=compare_eqtl_runs_ctr(baseline_data_dir, comparison_data_dir, fdr_threshold=fdr_threshold, outPDF, cache_dir=cache_dir)
        df=z$df
        results[[i]]=df
    }
    df=do.call(rbind, results)
    #Save file summary statistics across all runs for trend plotting.
    write.table(df, file=paste(outDir,"/compare_eQTL_all_levels_stats.txt", sep=""), sep="\t", quote=F, row.names=F)

}


#' Compare eQTL runs for one baseline/comparison directory pair
#'
#' For a given baseline and comparison directory, iterates over all
#' \code{cell_type__region} subdirectories in \code{baseline_data_dir},
#' runs per-dataset comparisons via \code{compare_eqtl_runs()}, and optionally
#' writes plots to a PDF and summary statistics to a tab-delimited file.
#'
#' @param baseline_data_dir Character scalar. Directory for the baseline eQTL
#'   results (e.g., \code{.../LEVEL_0}).
#' @param comparison_data_dir Character scalar. Directory for the comparison eQTL
#'   results (e.g., \code{.../LEVEL_2}).
#' @param fdr_threshold Numeric scalar. FDR threshold passed to
#'   \code{compare_eqtl_runs()} (default 0.05).
#' @param outPDF Character scalar or NULL. If non-NULL, a PDF is created and all
#'   per-dataset plots are written to it.
#' @param outFile Character scalar or NULL. If non-NULL, per-dataset summary
#'   statistics are written to this file as tab-delimited text.
#' @param cache_dir Character scalar or NULL. Cache directory passed to
#'   \code{compare_eqtl_runs()}.
#'
#' @return A data.frame of per-dataset summary statistics (one row per
#'   cell_type/region).
#'
#' @importFrom logger log_info
#' @importFrom utils write.table
#' @importFrom grDevices pdf dev.off
#' @export
compare_eqtl_runs_ctr<-function (baseline_data_dir, comparison_data_dir, fdr_threshold=0.05, outPDF=NULL, outFile=NULL, cache_dir=NULL) {

    file_separator="__"

    data_sets=list.files(baseline_data_dir, include.dirs=T)
    #split into cell and region to drive individual comparisons
    z=strsplit (data_sets, file_separator)
    cell_types=sapply (z, function (x) x[1])
    regions=sapply (z, function (x) x[2])
    df=data.frame(
        cell_type=cell_types,
        region=regions,
        stringsAsFactors = FALSE
    )
    #so that the comparisons are always in the same order
    df=df[order(df$cell_type, df$region),]

    baseline_name= basename (baseline_data_dir)
    comparison_name= basename (comparison_data_dir)
    results=list()

    if (!is.null(outPDF)) {
        grDevices::pdf(outPDF, width = 11, height = 11)
    }

    for (i in 1:nrow(df)) {
        cell_type=df$cell_type[i]
        region=df$region[i]
        logger::log_info(paste0("Comparing eQTL results for ", cell_type, " ", region, "\n"))
        index_file=paste0(baseline_data_dir,"/", cell_type, file_separator, region, "/", cell_type, file_separator, region, ".cis_qtl_ann.txt.gz")
        index_file_comparison=paste0(comparison_data_dir,"/", cell_type, file_separator, region, "/", cell_type, file_separator, region, ".cis_qtl_ann.txt.gz")
        all_pairs_file=paste0(baseline_data_dir,"/", cell_type, file_separator, region, "/", cell_type, file_separator, region, ".cis_qtl_pairs.txt.gz")
        all_pairs_file_comparison=paste0(comparison_data_dir,"/", cell_type, file_separator, region, "/", cell_type, file_separator, region, ".cis_qtl_pairs.txt.gz")

        z=compare_eqtl_runs(
            cell_type=cell_type,
            region=region,
            baseline_name=baseline_name,
            comparison_name=comparison_name,
            index_file=index_file,
            index_file_comparison=index_file_comparison,
            all_pairs_file=all_pairs_file,
            all_pairs_file_comparison=all_pairs_file_comparison,
            fdr_threshold=fdr_threshold,
            cache_dir=cache_dir
        )

        results[[i]]=z$stats_df
        if (!is.null(outPDF)) {
            print(z$final_p1)
            print(z$final_p2)
        }

    }

    if (!is.null(outPDF)) {
        grDevices::dev.off()
    }

    #merge summary stats write to file and return
    df=do.call(rbind, results)

    if (!is.null(outFile)) {
        write.table(df, file=outFile, sep="\t", quote=F, row.names=F)
    }

    return (df)
}


#' Compare eQTL runs
#'
#' Compares two eQTL runs (baseline and comparison) for a given \code{cell_type}
#' and \code{region}. The function loads the index eQTL results for each run and
#' augments them with cross-experiment slopes derived from the corresponding
#' all-pairs (all SNP-gene tests) files. The augmented slopes are used to compute
#' sign-test concordance between runs.
#'
#' The function returns two multi-panel plots (intended as two PDF pages):
#' (1) discovery and effect-size / q-value comparisons, and (2) sign-test and
#' p-value inflation diagnostics. It also returns a single-row statistics
#' data.frame summarizing overlap, correlation, yield, and sign-test concordance.
#'
#' Cross-slope augmentation is cached to avoid repeatedly reading large all-pairs
#' files. The cache location defaults to a persistent user cache directory and
#' can be overridden either by setting the \code{EQTL_CACHE_DIR} environment
#' variable or by providing \code{cache_dir}.
#'
#' @param cell_type Character scalar. Cell type label used for plot titles and
#'   for locating input files in higher-level wrappers.
#' @param region Character scalar. Region label used for plot titles and for
#'   locating input files in higher-level wrappers.
#' @param baseline_name Character scalar. Label for the baseline run (used in
#'   plot labels and titles).
#' @param comparison_name Character scalar. Label for the comparison run (used
#'   in plot labels and titles).
#' @param index_file Character scalar. Path to the baseline index eQTL file
#'   (one top SNP per gene).
#' @param index_file_comparison Character scalar. Path to the comparison index
#'   eQTL file (one top SNP per gene).
#' @param all_pairs_file Character scalar. Path to the baseline all-pairs file
#'   (all SNP-gene tests).
#' @param all_pairs_file_comparison Character scalar. Path to the comparison
#'   all-pairs file (all SNP-gene tests).
#' @param fdr_threshold Numeric scalar. FDR threshold used for defining eGenes
#'   and for sign-test evaluation (default 0.05).
#' @param cache_dir Character scalar or NULL. Optional override for the cache
#'   directory used to store augmented index files with cross slopes. If NULL,
#'   the default cache directory is used (see Details).
#'
#' @return A named list with:
#'   \describe{
#'     \item{final_p1}{A cowplot object containing discovery and effect/q-value
#'       comparison plots.}
#'     \item{final_p2}{A cowplot object containing sign-test and p-value inflation
#'       diagnostic plots.}
#'     \item{stats_df}{A single-row data.frame of summary statistics for this
#'       \code{cell_type}/\code{region} comparison.}
#'   }
#'
#' @importFrom data.table fread setDF as.data.table is.data.table
#' @importFrom ggplot2 ggplot aes geom_point geom_abline labs ggtitle theme theme_bw element_text scale_color_manual
#' @importFrom cowplot plot_grid ggdraw draw_label
#' @importFrom logger log_info
#' @export
compare_eqtl_runs<-function (cell_type="MSN", region="CaH", baseline_name, comparison_name,
                             index_file, index_file_comparison,
                             all_pairs_file, all_pairs_file_comparison, fdr_threshold=0.05,
                             cache_dir=NULL) {


    cache_dir_final<-get_eqtl_cache_dir(cache_dir = cache_dir)

    aug <- get_or_build_augmented_indices_for_sign(
        index_file = index_file,
        index_file_comparison = index_file_comparison,
        all_pairs_file = all_pairs_file,
        all_pairs_file_comparison = all_pairs_file_comparison,
        cache_dir = cache_dir_final
    )

    index_dt <- aug$index_dt
    index_comparison_dt <- aug$index_comparison_dt

    results_index= compare_index_eqtls(index_dt, index_comparison_dt,
                                       baseline_name=baseline_name,
                                       comparison_name=comparison_name,
                                       fdr_threshold=fdr_threshold)

    # Create first page plot of eQTL comparisons.
    p1=cowplot::plot_grid(
        results_index$venn_plot_genes_tested,
        results_index$venn_plot_eGenes,
        results_index$effect_size_scatter_plot,
        results_index$qval_scatter_plot,
        rel_heights = c(0.4, 0.6),
        ncol=2
    )

    final_p1 <- add_supertitle_cowplot(p = p1,
        title = paste("eQTL Comparison:", cell_type, region),
        title_size = 16,
        rel_height_title = 0.08
    )

    # Create second page of eQTL sign test comparisons + qq plots
    sign_test_result_baseline <- plot_sign_test(
        index_dt = index_dt,
        baseline_name = baseline_name,
        comparison_name = comparison_name,
        fdr_threshold = fdr_threshold
    )

    sign_test_result_comparison <- plot_sign_test(
        index_dt = index_comparison_dt,
        baseline_name = comparison_name,
        comparison_name = baseline_name,
        fdr_threshold = fdr_threshold
    )

    p2=cowplot::plot_grid(
        sign_test_result_baseline$plot,
        sign_test_result_comparison$plot,
        results_index$pval_inflation_plots_baseline,
        results_index$pval_inflation_plots_comparison,
        rel_heights = c(0.4, 0.6),
        ncol=2
    )

    final_p2 <- add_supertitle_cowplot(p = p2,
        title = paste("eQTL Sign Test and P-value Inflation:", cell_type, region),
        title_size = 16,
        rel_height_title = 0.08
    )

    stats_df=results_index$stats_df
    #add the sign test results to the stats df
    stats_df$sign_test_baseline_n_tested=sign_test_result_baseline$stats$n_tested
    stats_df$sign_test_baseline_n_concordant=sign_test_result_baseline$stats$n_agree
    stats_df$sign_test_baseline_concordance_rate=sign_test_result_baseline$stats$frac_agree
    stats_df$sign_test_comparison_n_tested=sign_test_result_comparison$stats$n_tested
    stats_df$sign_test_comparison_n_concordant=sign_test_result_comparison$stats$n_agree
    stats_df$sign_test_comparison_concordance_rate=sign_test_result_comparison$stats$frac_agree

    #add the cell type and region to the front of the stats df
    stats_df=data.frame(
        cell_type=cell_type,
        region=region,
        stats_df,
        stringsAsFactors = FALSE
    )

    return (list(
        final_p1=final_p1,
        final_p2=final_p2,
        stats_df=stats_df
    ))


}

compare_index_eqtls<-function (index_dt, index_comparison_dt, baseline_name, comparison_name, fdr_threshold=0.05) {
    #how many total genes tested? (barplot 1)
    z=plot_gene_venn(index_dt, index_comparison_dt, text_size=6, title_size=16,
                     baseline_name=baseline_name, comparison_name=comparison_name,
                     title="Number of genes tested")

    venn_plot_genes_tested=z$plot
    stats_genes_tested=z$stats

    #how many eQTLs are found in both levels (FDR<0.05)
    eQTLs_base=index_dt[index_dt$qval < fdr_threshold, ]
    eQTLs_comp=index_comparison_dt[index_comparison_dt$qval < fdr_threshold, ]

    z=plot_gene_venn(eQTLs_base, eQTLs_comp, text_size=6, title_size=16,
                     baseline_name=baseline_name, comparison_name=comparison_name,
                     title="Number of eQTLs discovered")

    venn_plot_eGenes=z$plot
    stats_eGenes=z$stats

    cowplot::plot_grid(venn_plot_genes_tested, venn_plot_eGenes, ncol=2)

    #what are the effect sizes of the eGenes?
    #Use the absolute effect sizes.
    z=plot_egene_effect_scatter(index_dt, index_comparison_dt,
                              baseline_name=baseline_name, comparison_name=comparison_name,
                              fdr_threshold = 0.05, abs_slopes = TRUE,
                              title = "Effect sizes [union sig eGenes]",
                              point_size = 1.2)
    effect_size_scatter_plot=z$plot
    effect_size_stats=z$stats

    z=plot_egene_qval_compare(index_dt, index_comparison_dt,
                            baseline_name=baseline_name, comparison_name=comparison_name,
                            fdr_threshold = 0.05, title = "eGene q-value comparison",
                            point_size = 1.2, title_size = 16)
    qval_scatter_plot=z$plot
    qval_stats=z$stats

    z=test_eqtl_pval_inflation(index_dt,
                             main_prefix = paste(baseline_name),
                             pval_col = "pval_beta",
                             qq_point_size = 1,
                             qq_point_alpha = 0.6,
                             hist_bins = 40)
    pval_inflation_plots_baseline=z$combined_plot
    pval_inflation_stats_baseline=z$stats

    z=test_eqtl_pval_inflation(index_comparison_dt,
                               main_prefix = paste(comparison_name),
                               pval_col = "pval_beta",
                               qq_point_size = 1,
                               qq_point_alpha = 0.6,
                               hist_bins = 40)
    pval_inflation_plots_comparison=z$combined_plot
    pval_inflation_stats_comparison=z$stats

    #create a single stats df
    stats_df=make_compare_index_eqtls_stats_df(stats_genes_tested,
                                               stats_eGenes,
                                               effect_size_stats,
                                               qval_stats)
    return (list(
        venn_plot_genes_tested=venn_plot_genes_tested,
        venn_plot_eGenes=venn_plot_eGenes,
        effect_size_scatter_plot=effect_size_scatter_plot,
        qval_scatter_plot=qval_scatter_plot,
        pval_inflation_plots_baseline=pval_inflation_plots_baseline,
        pval_inflation_plots_comparison=pval_inflation_plots_comparison,
        stats_df=stats_df
    ))
}


##########################
# PLOTS
##########################


add_supertitle_cowplot <- function(p,
                                   title,
                                   title_size = 16,
                                   rel_height_title = 0.08) {
    title_grob <- cowplot::ggdraw() +
        cowplot::draw_label(
            title,
            fontface = "bold",
            size = title_size,
            x = 0.5,
            hjust = 0.5
        )

    cowplot::plot_grid(
        title_grob,
        p,
        ncol = 1,
        rel_heights = c(rel_height_title, 1)
    )
}


#' Venn diagram comparing overlap between two eQTL result sets
#'
#' Creates a two-set Venn diagram showing overlap of items defined by
#' \code{gene_name} (or any pre-filtered subset of rows, e.g. FDR < 0.05).
#' Also computes union and intersection counts for tracking across runs.
#'
#' @param index_dt data.frame for the baseline run.
#'   Must include a \code{gene_name} column.
#' @param index_comparison_dt data.frame for the comparison run.
#'   Must include a \code{gene_name} column.
#' @param baseline_name character scalar. Label for the baseline dataset.
#' @param comparison_name character scalar. Label for the comparison dataset.
#' @param title character scalar. Plot title.
#' @param text_size numeric scalar. Text size for counts and percentages in the Venn.
#' @param title_size numeric scalar. Title font size.
#'
#' @return A list with:
#' \describe{
#' \item{\code{plot}}{A ggplot object (the Venn diagram).}
#' \item{\code{stats}}{A one-row data.frame with union and intersection counts.}
#' }
#'
#' @importFrom ggvenn ggvenn
#' @importFrom ggplot2 ggtitle theme element_text
plot_gene_venn <- function(index_dt,
                           index_comparison_dt,
                           baseline_name = "baseline",
                           comparison_name = "comparison",
                           title = "Number of genes tested",
                           text_size = 6,
                           title_size = 16) {
    gene_name <- NULL

    genes_baseline <- unique(index_dt[["gene_name"]])
    genes_comp <- unique(index_comparison_dt[["gene_name"]])

    n_intersect <- length(intersect(genes_baseline, genes_comp))
    n_union <- length(union(genes_baseline, genes_comp))

    stats <- data.frame(
        baseline_name = baseline_name,
        comparison_name = comparison_name,
        n_union = n_union,
        n_intersect = n_intersect,
        stringsAsFactors = FALSE
    )

    sets <- list(
        baseline = genes_baseline,
        comparison = genes_comp
    )
    names(sets) <- c(baseline_name, comparison_name)

    p <- ggvenn::ggvenn(
        sets,
        fill_color = c("#0072B2", "#009E73"),
        text_size = text_size
    ) +
        ggplot2::ggtitle(title) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(
                size = title_size,
                hjust = 0.5
            )
        )

    list(plot = p, stats = stats)
}



# Helper: prepare per-gene merged data + union categories (shared by multiple plots)
prepare_egene_union_df <- function(index_dt,
                                   index_comparison_dt,
                                   baseline_name = "baseline",
                                   comparison_name = "comparison",
                                   fdr_threshold = 0.05) {
    gene_name <- NULL

    d1 <- index_dt[, c("gene_name", "slope", "qval")]
    d2 <- index_comparison_dt[, c("gene_name", "slope", "qval")]
    colnames(d1) <- c("gene_name", "slope_base", "qval_base")
    colnames(d2) <- c("gene_name", "slope_comp", "qval_comp")

    if (any(duplicated(d1$gene_name))) {
        stop("index_dt has duplicated gene_name values; expected one row per gene.")
    }
    if (any(duplicated(d2$gene_name))) {
        stop("index_comparison_dt has duplicated gene_name values; expected one row per gene.")
    }

    m <- merge(d1, d2, by = "gene_name", all = TRUE)

    base_sig <- !is.na(m$qval_base) & (m$qval_base < fdr_threshold)
    comp_sig <- !is.na(m$qval_comp) & (m$qval_comp < fdr_threshold)

    in_union_sig <- base_sig | comp_sig
    df <- m[in_union_sig, , drop = FALSE]

    df$category <- NA_character_
    df$category[base_sig[in_union_sig] & comp_sig[in_union_sig]] <- "significant_in_both"
    df$category[base_sig[in_union_sig] & !comp_sig[in_union_sig]] <- "baseline_only"
    df$category[!base_sig[in_union_sig] & comp_sig[in_union_sig]] <- "comparison_only"

    stats <- data.frame(
        baseline_name = baseline_name,
        comparison_name = comparison_name,
        fdr_threshold = fdr_threshold,
        n_union_egenes = nrow(df),
        n_intersect_egenes = sum(df$category == "significant_in_both", na.rm = TRUE),
        n_baseline_only_egenes = sum(df$category == "baseline_only", na.rm = TRUE),
        n_comparison_only_egenes = sum(df$category == "comparison_only", na.rm = TRUE),
        stringsAsFactors = FALSE
    )

    list(df = df, stats = stats)
}

#' Effect size scatter plot for union of discovered eGenes
#'
#' Scatter plot of per-gene slopes for eGenes significant in either run
#' at \code{qval < fdr_threshold}, colored by Venn category:
#' significant in both runs, baseline only, or comparison only.
#'
#' Inputs should be unfiltered per-gene result tables (one row per gene).
#'
#' @param index_dt data.frame for the baseline run. Must include columns:
#'   \code{gene_name}, \code{slope}, \code{qval}.
#' @param index_comparison_dt data.frame for the comparison run. Must include columns:
#'   \code{gene_name}, \code{slope}, \code{qval}.
#' @param baseline_name character scalar. Label for the baseline run.
#' @param comparison_name character scalar. Label for the comparison run.
#' @param fdr_threshold numeric scalar. Threshold applied to \code{qval} to define the union set.
#' @param abs_slopes logical scalar. If TRUE, plot \code{abs(slope)} for both runs.
#' @param title character scalar. Plot title.
#' @param point_size numeric scalar. Point size.
#' @param title_size numeric scalar. Title font size.
#'
#' @return A list with:
#' \describe{
#' \item{\code{plot}}{A ggplot object.}
#' \item{\code{stats}}{A one-row data.frame with union/intersection counts and \code{n_plotted}.}
#' \item{\code{df}}{The underlying merged data.frame used for plotting/statistics.}
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point labs ggtitle theme element_text scale_color_manual
plot_egene_effect_scatter <- function(index_dt,
                                      index_comparison_dt,
                                      baseline_name = "baseline",
                                      comparison_name = "comparison",
                                      fdr_threshold = 0.05,
                                      abs_slopes = FALSE,
                                      title = "Effect sizes for union of discovered eGenes",
                                      point_size = 1.2,
                                      title_size = 16) {
    slope_base <- slope_comp <- category <- NULL

    prep <- prepare_egene_union_df(
        index_dt = index_dt,
        index_comparison_dt = index_comparison_dt,
        baseline_name = baseline_name,
        comparison_name = comparison_name,
        fdr_threshold = fdr_threshold
    )

    df <- prep$df

    ok <- is.finite(df$slope_base) & is.finite(df$slope_comp)
    df_plot <- df[ok, , drop = FALSE]

    if (abs_slopes) {
        df_plot$slope_base <- abs(df_plot$slope_base)
        df_plot$slope_comp <- abs(df_plot$slope_comp)
        xlab_txt <- paste0(baseline_name, " abs(slope)")
        ylab_txt <- paste0(comparison_name, " abs(slope)")
    } else {
        xlab_txt <- paste0(baseline_name, " slope")
        ylab_txt <- paste0(comparison_name, " slope")
    }

    p <- ggplot2::ggplot(
        df_plot,
        ggplot2::aes(x = slope_base, y = slope_comp, color = category)
    ) +
        ggplot2::geom_point(size = point_size, alpha = 0.8) +
        ggplot2::scale_color_manual(
            values = c(
                significant_in_both = "#000000",
                baseline_only = "#0072B2",
                comparison_only = "#D55E00"
            ),
            labels = c(
                significant_in_both = "sig both",
                baseline_only = paste0(baseline_name),
                comparison_only = paste0(comparison_name)
            )
        ) +
        ggplot2::labs(x = xlab_txt, y = ylab_txt, color = NULL) +
        ggplot2::ggtitle(title) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(size = title_size, hjust = 0.5),
            legend.position = "top",
            legend.key.height = grid::unit(0.6, "lines"),
            legend.key.width = grid::unit(1.2, "lines")
        )

    stats <- prep$stats
    stats$n_plotted <- nrow(df_plot)

    list(plot = p, stats = stats, df = df)
}

#' Compare eGene q-values between baseline and comparison runs
#'
#' Scatter plot of \code{-log10(qval)} in baseline vs comparison for the union
#' of significant eGenes (qval < fdr_threshold in either run), colored by Venn category.
#' Includes a black dashed y=x reference line and annotates correlation and yield.
#'
#' Inputs should be unfiltered per-gene result tables (one row per gene).
#'
#' @param index_dt data.frame for the baseline run. Must include columns:
#'   \code{gene_name}, \code{slope}, \code{qval}.
#' @param index_comparison_dt data.frame for the comparison run. Must include columns:
#'   \code{gene_name}, \code{slope}, \code{qval}.
#' @param baseline_name character scalar. Label for the baseline run.
#' @param comparison_name character scalar. Label for the comparison run.
#' @param fdr_threshold numeric scalar. Threshold applied to \code{qval} to define the union set.
#' @param title character scalar. Plot title.
#' @param point_size numeric scalar. Point size.
#' @param title_size numeric scalar. Title font size.
#' @param cor_method character scalar. Correlation method passed to \code{cor()}.
#'
#' @return A list with:
#' \describe{
#' \item{\code{plot}}{A ggplot object.}
#' \item{\code{stats}}{A one-row data.frame with union/intersection counts, correlation, and yield.}
#' \item{\code{df}}{The underlying merged data.frame used for plotting/statistics.}
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_abline annotate labs ggtitle theme element_text scale_color_manual
plot_egene_qval_compare <- function(index_dt,
                                    index_comparison_dt,
                                    baseline_name = "baseline",
                                    comparison_name = "comparison",
                                    fdr_threshold = 0.05,
                                    title = "eGene q-value comparison",
                                    point_size = 1.2,
                                    title_size = 16,
                                    cor_method = "pearson") {
    qval_base <- qval_comp <- base_log10q <- comp_log10q <- category <- NULL

    prep <- prepare_egene_union_df(
        index_dt = index_dt,
        index_comparison_dt = index_comparison_dt,
        baseline_name = baseline_name,
        comparison_name = comparison_name,
        fdr_threshold = fdr_threshold
    )

    df <- prep$df

    df$base_log10q <- -log10(df$qval_base)
    df$comp_log10q <- -log10(df$qval_comp)

    ok <- is.finite(df$base_log10q) & is.finite(df$comp_log10q)
    df_plot <- df[ok, , drop = FALSE]

    cor_val <- if (nrow(df_plot) >= 2) {
        stats::cor(df_plot$base_log10q, df_plot$comp_log10q, method = cor_method)
    } else {
        NA_real_
    }

    total_eGenes_base <- sum(!is.na(df$qval_base) & (df$qval_base < fdr_threshold))
    total_eGenes_comp <- sum(!is.na(df$qval_comp) & (df$qval_comp < fdr_threshold))
    yield <- if (total_eGenes_base > 0) (total_eGenes_comp / total_eGenes_base) else NA_real_

    x_min <- min(df_plot$base_log10q, na.rm = TRUE)
    x_max <- max(df_plot$base_log10q, na.rm = TRUE)
    y_min <- min(df_plot$comp_log10q, na.rm = TRUE)
    y_max <- max(df_plot$comp_log10q, na.rm = TRUE)

    x_anno <- x_min + 0.02 * (x_max - x_min)
    y_anno1 <- y_max - 0.05 * (y_max - y_min)
    y_anno2 <- y_anno1 - 0.06 * (y_max - y_min)

    cor_label <- if (is.finite(cor_val)) {
        paste0("cor (", cor_method, ") = ", sprintf("%.3f", cor_val))
    } else {
        paste0("cor (", cor_method, ") = NA")
    }

    yield_label <- if (is.finite(yield)) {
        paste0("yield = ", sprintf("%.3f", yield),
               " (", total_eGenes_comp, " / ", total_eGenes_base, ")")
    } else {
        paste0("yield = NA (", total_eGenes_comp, " / ", total_eGenes_base, ")")
    }

    # Make R CMD CHECK happy
    base_log10q <- comp_log10q <- NULL

    p <- ggplot2::ggplot(
        df_plot,
        ggplot2::aes(x = base_log10q, y = comp_log10q, color = category)
    ) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.7) +
        ggplot2::geom_point(size = point_size, alpha = 0.8) +
        ggplot2::scale_color_manual(
            values = c(
                significant_in_both = "#000000",
                baseline_only = "#0072B2",
                comparison_only = "#D55E00"
            ),
            labels = c(
                significant_in_both = "sig both",
                baseline_only = paste0(baseline_name),
                comparison_only = paste0(comparison_name)
            )
        ) +
        ggplot2::annotate("text", x = x_anno, y = y_anno1, hjust = 0, label = cor_label) +
        ggplot2::annotate("text", x = x_anno, y = y_anno2, hjust = 0, label = yield_label) +
        ggplot2::labs(
            x = paste0(baseline_name, " -log10(qval)"),
            y = paste0(comparison_name, " -log10(qval)"),
            color = NULL
        ) +
        ggplot2::ggtitle(title) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(size = title_size, hjust = 0.5),
            legend.position = "top",
            legend.key.height = grid::unit(0.6, "lines"),
            legend.key.width = grid::unit(1.2, "lines")
        )

    stats <- prep$stats
    stats$n_complete_pairs <- nrow(df_plot)
    stats$cor_method <- cor_method
    stats$cor_val <- cor_val
    stats$total_eGenes_base <- total_eGenes_base
    stats$total_eGenes_comp <- total_eGenes_comp
    stats$yield <- yield

    list(plot = p, stats = stats, df = df)
}

#' Diagnose p-value inflation/deflation for eQTL/eGene results
#'
#' Computes a genomic-control style inflation factor (lambda) from p-values
#' and produces a QQ plot and histogram. Optionally reports a tail-uniformity
#' summary for p-values >= tail_min_p.
#'
#' Interpretation:
#' - lambda_gc summarizes global signal inflation and is only meaningful on the
#'   full, untruncated p-value distribution.
#' - Tail diagnostics use a normalized mean ratio where 1.0 indicates the
#'   expected mean under Uniform(tail_min_p, 1).
#'
#' The histogram includes a red dashed horizontal segment over [tail_min_p, 1]
#' at the expected per-bin count under Uniform(tail_min_p, 1).
#'
#' @param eqtl_dataframe data.frame containing p-values to diagnose.
#' @param main_prefix character scalar. Used in plot titles.
#' @param pval_col character scalar. Column name holding p-values.
#' @param qq_point_size numeric scalar. Point size for QQ plot.
#' @param qq_point_alpha numeric scalar. Alpha for QQ plot points.
#' @param hist_bins integer scalar. Number of histogram bins.
#' @param tail_min_p numeric scalar in [0,1] or NULL. If not NULL, compute a
#'   tail-uniformity summary on p-values satisfying p >= tail_min_p, and draw
#'   the tail reference segment on the histogram.
#'
#' @return A list with:
#' \describe{
#' \item{\code{combined_plot}}{Combined QQ + histogram plot (cowplot object).}
#' \item{\code{stats}}{List containing lambda and tail summary.}
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_abline annotate labs theme_bw
#' @importFrom ggplot2 geom_histogram geom_segment
#' @importFrom stats qchisq
#' @importFrom cowplot plot_grid
test_eqtl_pval_inflation <- function(eqtl_dataframe,
                                     main_prefix = "eQTL Analysis",
                                     pval_col = "pval_beta",
                                     qq_point_size = 1,
                                     qq_point_alpha = 0.6,
                                     hist_bins = 40,
                                     tail_min_p = 0.5) {
    expected <- observed <- pvals <- NULL

    stopifnot(is.data.frame(eqtl_dataframe))
    stopifnot(pval_col %in% names(eqtl_dataframe))

    pvals_raw <- eqtl_dataframe[[pval_col]]
    ok <- is.finite(pvals_raw) & !is.na(pvals_raw) &
        (pvals_raw >= 0) & (pvals_raw <= 1)
    pvals <- pvals_raw[ok]
    n <- length(pvals)

    if (n < 2) {
        stop("Not enough valid p-values to run diagnostics (need at least 2).")
    }

    ## ---- Global lambda (genomic control style) ----
    chi2 <- stats::qchisq(1 - pvals, df = 1)
    lambda_gc <- stats::median(chi2, na.rm = TRUE) /
        stats::qchisq(0.5, df = 1)

    ## ---- QQ plot data ----
    p_sorted <- sort(pvals)
    expected <- -log10((seq_len(n) - 0.5) / n)
    observed <- -log10(p_sorted)
    dfqq <- data.frame(expected = expected, observed = observed)

    ## ---- Annotation placement ----
    x_range <- range(dfqq$expected, finite = TRUE)
    y_range <- range(dfqq$observed, finite = TRUE)
    x_annot <- x_range[1] + 0.05 * diff(x_range)
    y_annot <- y_range[2] - 0.05 * diff(y_range)

    ## ---- Annotation text ----
    label_lines <- c(sprintf("lambda = %.3f", lambda_gc))

    ## ---- Tail-uniformity summary ----
    tail_stats <- NULL
    if (!is.null(tail_min_p)) {
        stopifnot(is.numeric(tail_min_p), length(tail_min_p) == 1,
                  tail_min_p >= 0, tail_min_p <= 1)

        p_tail <- pvals[pvals >= tail_min_p]
        if (length(p_tail) >= 2) {

            #Compare the mean of p-values in the tail to the expected mean under uniform
            #Values < 1 indicate some skew towards smaller p-values in the tail
            tail_mean <- mean(p_tail)
            tail_mean_expected <- (tail_min_p + 1) / 2
            tail_mean_ratio <- tail_mean / tail_mean_expected

            label_lines <- c(
                label_lines,
                "",
                sprintf("tail p>=%.2f:", tail_min_p),
                sprintf("  n = %d", length(p_tail)),
                sprintf("  mean ratio = %.3f", tail_mean_ratio)
            )

            tail_stats <- list(
                tail_min_p = tail_min_p,
                n_used = length(p_tail),
                mean_p = tail_mean,
                mean_p_expected = tail_mean_expected,
                mean_ratio = tail_mean_ratio
            )
        } else {
            label_lines <- c(
                label_lines,
                "",
                sprintf("tail p>=%.2f: n<2", tail_min_p)
            )
        }
    }

    annot_label <- paste(label_lines, collapse = "\n")

    ## ---- QQ plot ----
    qqplot <- ggplot2::ggplot(dfqq, ggplot2::aes(x = expected, y = observed)) +
        ggplot2::geom_point(alpha = qq_point_alpha, size = qq_point_size) +
        ggplot2::geom_abline(intercept = 0, slope = 1,
                             linetype = "dashed", linewidth = 0.7) +
        ggplot2::annotate(
            "text",
            x = x_annot, y = y_annot,
            label = annot_label,
            hjust = 0, vjust = 1,
            size = 3
        ) +
        ggplot2::labs(
            title = paste(main_prefix, "QQ plot of p-values"),
            x = "Expected -log10(p)",
            y = "Observed -log10(p)"
        ) +
        ggplot2::theme_bw()

    ## ---- Histogram + tail uniform reference segment ----
    histplot <- ggplot2::ggplot(
        data.frame(pvals = pvals),
        ggplot2::aes(x = pvals)
    ) +
        ggplot2::geom_histogram(
            bins = hist_bins,
            color = "black",
            fill = "grey80"
        )

    expected_tail_count_per_bin <- NA_real_
    n_tail <- NA_integer_
    n_bins_tail <- NA_integer_

    if (!is.null(tail_min_p)) {
        n_tail <- sum(pvals >= tail_min_p)
        n_bins_tail <- max(1L, as.integer(round(hist_bins * (1 - tail_min_p))))
        expected_tail_count_per_bin <- n_tail / n_bins_tail

        histplot <- histplot +
            ggplot2::geom_segment(
                x = tail_min_p,
                xend = 1,
                y = expected_tail_count_per_bin,
                yend = expected_tail_count_per_bin,
                linetype = "dashed",
                linewidth = 0.9,
                color = "red"
            )
    }

    histplot <- histplot +
        ggplot2::labs(
            title = paste(main_prefix, "P-value histogram"),
            x = "P-value",
            y = "Count"
        ) +
        ggplot2::theme_bw()

    ## ---- Combine panels ----
    combined_plot <- cowplot::plot_grid(
        qqplot,
        histplot,
        ncol = 1,
        rel_heights = c(2, 1),
        align = "v"
    )

    list(
        combined_plot = combined_plot,
        stats = list(
            n_used = n,
            lambda_gc = lambda_gc,
            tail = tail_stats,
            hist_tail_reference = list(
                n_tail = n_tail,
                n_bins_tail = n_bins_tail,
                expected_tail_count_per_bin = expected_tail_count_per_bin
            )
        )
    )
}


plot_sign_test <- function(index_dt,
                                      baseline_name,
                                      comparison_name,
                                      slope_col = "slope",
                                      slope_cross_col = "slope_cross",
                                      qval_col = "qval",
                                      fdr_threshold = 0.05,
                                      point_size = 0.6,
                                      point_alpha = 0.6,
                                      title_size = 12,
                                      axis_title_size = 11,
                                      axis_text_size = 10) {
    stopifnot(is.data.frame(index_dt) || data.table::is.data.table(index_dt))
    stopifnot(all(c(slope_col, slope_cross_col, qval_col) %in% names(index_dt)))

    df <- index_dt

    ## Same filter used for sign test + plotting
    idx_fdr <- !is.na(df[[qval_col]]) & (df[[qval_col]] < fdr_threshold)
    df_fdr <- df[idx_fdr, , drop = FALSE]

    ok <- !is.na(df_fdr[[slope_col]]) & !is.na(df_fdr[[slope_cross_col]])
    n_compared <- sum(ok)

    ## Sign agreement on the same plotted set (exclude zeros to match earlier logic)
    s1 <- sign(df_fdr[[slope_col]][ok])
    s2 <- sign(df_fdr[[slope_cross_col]][ok])
    ok_sign <- (s1 != 0) & (s2 != 0)

    n_tested <- sum(ok_sign)
    n_agree <- if (n_tested > 0) sum(s1[ok_sign] == s2[ok_sign]) else 0L
    frac_agree <- if (n_tested > 0) n_agree / n_tested else NA_real_

    sign_pct_txt <- if (!is.na(frac_agree)) sprintf("%.1f%%", 100 * frac_agree) else "NA"

    title_txt <- paste0(
        baseline_name, " vs ", comparison_name,
        "\nfeatures compared (q<", fdr_threshold, "): ", n_compared,
        " | sign agreement: ", sign_pct_txt
    )

    xlab_txt <- paste0("slope (", baseline_name, ")")
    ylab_txt <- paste0("slope (", comparison_name, ")")

    p <- ggplot2::ggplot(
        df_fdr[ok, , drop = FALSE],
        ggplot2::aes(x = .data[[slope_col]], y = .data[[slope_cross_col]])
    ) +
        ggplot2::geom_point(size = point_size, alpha = point_alpha) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
        ggplot2::labs(title = title_txt, x = xlab_txt, y = ylab_txt) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(size = title_size, hjust = 0.5),
            axis.title = ggplot2::element_text(size = axis_title_size),
            axis.text = ggplot2::element_text(size = axis_text_size)
        )

    list(
        plot = p,
        stats = data.frame(
            baseline_name = baseline_name,
            comparison_name = comparison_name,
            fdr_threshold = fdr_threshold,
            n_compared = as.integer(n_compared),
            n_tested = as.integer(n_tested),
            n_agree = as.integer(n_agree),
            frac_agree = as.numeric(frac_agree)
        )
    )
}



#################################
# SIGN TEST
# Because the files are so large, cache the augmented index files with cross slopes
#################################

.add_slope_cross_join <- function(index_dt,
                                  all_pairs_dt,
                                  out_col,
                                  all_pairs_slope_col) {
    need_left <- c("gene_name", "variant_id")
    need_right <- c("phenotype_id", "variant_id", all_pairs_slope_col)

    if (!(all(need_left %in% names(index_dt)) &&
          all(need_right %in% names(all_pairs_dt)))) {
        logger::log_info("Sign-test: missing join columns; leaving slope_cross as NA.")
        return(invisible(NULL))
    }

    index_dt[
        all_pairs_dt,
        (out_col) := get(paste0("i.", all_pairs_slope_col)),
        on = c("gene_name" = "phenotype_id", "variant_id" = "variant_id")
    ]

    invisible(NULL)
}

add_cross_slopes_for_sign_test <- function(index_dt,
                                           index_comparison_dt,
                                           all_pairs,
                                           all_pairs_comparison,
                                           out_col = "slope_cross",
                                           all_pairs_slope_col = "slope") {
    stopifnot(is.data.frame(index_dt) || data.table::is.data.table(index_dt))
    stopifnot(is.data.frame(index_comparison_dt) || data.table::is.data.table(index_comparison_dt))
    stopifnot(is.data.frame(all_pairs) || data.table::is.data.table(all_pairs))
    stopifnot(is.data.frame(all_pairs_comparison) || data.table::is.data.table(all_pairs_comparison))

    ## Track original classes
    idx_is_dt <- data.table::is.data.table(index_dt)
    idxc_is_dt <- data.table::is.data.table(index_comparison_dt)

    ## Work in data.table
    index_dt <- data.table::as.data.table(index_dt)
    index_comparison_dt <- data.table::as.data.table(index_comparison_dt)
    all_pairs <- data.table::as.data.table(all_pairs)
    all_pairs_comparison <- data.table::as.data.table(all_pairs_comparison)

    ## Default: NA (covers missing matches or missing columns)
    index_dt[[out_col]] <- NA_real_
    index_comparison_dt[[out_col]] <- NA_real_

    ## Baseline index_dt <- comparison all_pairs_comparison
    .add_slope_cross_join(
        index_dt = index_dt,
        all_pairs_dt = all_pairs_comparison,
        out_col = out_col,
        all_pairs_slope_col = all_pairs_slope_col
    )

    ## Comparison index_comparison_dt <- baseline all_pairs
    .add_slope_cross_join(
        index_dt = index_comparison_dt,
        all_pairs_dt = all_pairs,
        out_col = out_col,
        all_pairs_slope_col = all_pairs_slope_col
    )

    ## Restore original class
    if (!idx_is_dt) {
        data.table::setDF(index_dt)
    }
    if (!idxc_is_dt) {
        data.table::setDF(index_comparison_dt)
    }

    list(
        index_dt = index_dt,
        index_comparison_dt = index_comparison_dt
    )
}


compute_sign_concordance_baseline_sig <- function(index_dt,
                                                  fdr_threshold = 0.05,
                                                  qval_col = "qval",
                                                  slope_col = "slope",
                                                  slope_cross_col = "slope_cross") {
    stopifnot(is.data.frame(index_dt))
    stopifnot(all(c(qval_col, slope_col, slope_cross_col) %in% names(index_dt)))

    idx <- !is.na(index_dt[[qval_col]]) & (index_dt[[qval_col]] < fdr_threshold) &
        !is.na(index_dt[[slope_col]]) & !is.na(index_dt[[slope_cross_col]])

    if (sum(idx) == 0) {
        return(list(n_tested = 0L, n_agree = 0L, frac_agree = NA_real_))
    }

    s1 <- sign(index_dt[[slope_col]][idx])
    s2 <- sign(index_dt[[slope_cross_col]][idx])

    ok <- (s1 != 0) & (s2 != 0)
    n_tested <- sum(ok)
    if (n_tested == 0) {
        return(list(n_tested = 0L, n_agree = 0L, frac_agree = NA_real_))
    }

    n_agree <- sum(s1[ok] == s2[ok])
    list(
        n_tested = as.integer(n_tested),
        n_agree = as.integer(n_agree),
        frac_agree = as.numeric(n_agree / n_tested)
    )
}

#' Build cache file paths for sign-test augmented eQTL index files
#'
#' Constructs deterministic cache file paths for baseline and comparison
#' eQTL index tables augmented with cross-experiment slopes for sign tests.
#' The cache key is derived from file metadata (normalized path, file size,
#' and modification time) for the input files, avoiding any need to read
#' large file contents.
#'
#' @param index_file Character scalar. Path to the baseline eQTL annotation file.
#' @param index_file_comparison Character scalar. Path to the comparison eQTL
#'   annotation file.
#' @param all_pairs_file Character scalar. Path to the baseline all-pairs file.
#' @param all_pairs_file_comparison Character scalar. Path to the comparison
#'   all-pairs file.
#' @param cache_dir Character scalar. Directory in which cache files should be
#'   created.
#' @param prefix Character scalar. Filename prefix used for cache files
#'   (default \code{"sign_cache"}).
#' @param path_suffix_depth Integer scalar. Number of path components from the
#'   end of each input file path to include in the cache key (default 4).
#' @return A list with two elements:
#'   \describe{
#'     \item{index_path}{Cache path for the baseline augmented index file.}
#'     \item{comparison_path}{Cache path for the comparison augmented index file.}
#'   }
#'
#' @importFrom digest digest
# build_sign_cache_paths <- function(index_file,
#                                    index_file_comparison,
#                                    all_pairs_file,
#                                    all_pairs_file_comparison,
#                                    cache_dir,
#                                    prefix = "sign_cache") {
#     files <- c(index_file, index_file_comparison, all_pairs_file, all_pairs_file_comparison)
#     stopifnot(all(file.exists(files)))
#
#     fi <- file.info(files)
#     stopifnot(!any(is.na(fi$size)), !any(is.na(fi$mtime)))
#
#     norm_paths <- vapply(
#         files,
#         function(x) normalizePath(x, winslash = "/", mustWork = TRUE),
#         character(1)
#     )
#
#     meta <- paste(
#         norm_paths,
#         fi$size,
#         as.numeric(fi$mtime),
#         sep = "|",
#         collapse = "||"
#     )
#
#     key <- digest::digest(meta, algo = "md5")
#
#     list(
#         index_path = file.path(cache_dir, paste0(prefix, "_index_", key, ".tsv.gz")),
#         comparison_path = file.path(cache_dir, paste0(prefix, "_comparison_", key, ".tsv.gz"))
#     )
# }
build_sign_cache_paths <- function(index_file,
                                   index_file_comparison,
                                   all_pairs_file,
                                   all_pairs_file_comparison,
                                   cache_dir,
                                   prefix = "sign_cache",
                                   path_suffix_depth = 4L) {
    .build_path_suffix <- function(p, n = 4L) {
        parts <- strsplit(p, "/")[[1]]
        parts <- parts[nzchar(parts)]
        utils::tail(parts, n)
    }

    files <- c(index_file, index_file_comparison, all_pairs_file, all_pairs_file_comparison)
    stopifnot(all(file.exists(files)))

    fi <- file.info(files)
    stopifnot(!any(is.na(fi$size)), !any(is.na(fi$mtime)))

    suffixes <- vapply(
        files,
        function(x) paste(.build_path_suffix(x, n = path_suffix_depth), collapse = "/"),
        character(1)
    )

    meta <- paste(
        suffixes,
        fi$size,
        as.numeric(fi$mtime),
        sep = "|",
        collapse = "||"
    )

    key <- digest::digest(meta, algo = "md5")

    list(
        index_path = file.path(cache_dir, paste0(prefix, "_index_", key, ".tsv.gz")),
        comparison_path = file.path(cache_dir, paste0(prefix, "_comparison_", key, ".tsv.gz"))
    )
}

get_or_build_augmented_indices_for_sign <- function(index_file,
                                                    index_file_comparison,
                                                    all_pairs_file,
                                                    all_pairs_file_comparison,
                                                    cache_dir,
                                                    force = FALSE) {
    paths <- build_sign_cache_paths(
        index_file = index_file,
        index_file_comparison = index_file_comparison,
        all_pairs_file = all_pairs_file,
        all_pairs_file_comparison = all_pairs_file_comparison,
        cache_dir = cache_dir
    )

    ## ---- Cache hit (unless forced) ----
    if (!force &&
        file.exists(paths$index_path) &&
        file.exists(paths$comparison_path)) {

        log_info("Sign-test cache hit. Using cached augmented index data.")

        index_dt <- read_index_file(paths$index_path)
        index_comparison_dt <- read_index_file(paths$comparison_path)

        stopifnot("slope_cross" %in% names(index_dt))
        stopifnot("slope_cross" %in% names(index_comparison_dt))

        return(list(
            index_dt = index_dt,
            index_comparison_dt = index_comparison_dt
        ))
    }

    ## ---- Cache miss or forced rebuild ----
    if (force) {
        log_info("Sign-test cache force rebuild.")
    } else {
        log_info("Sign-test cache miss. Reading source files.")
    }

    index_dt <- read_index_file(index_file)
    index_comparison_dt <- read_index_file(index_file_comparison)

    log_info("Reading SNP/Gene pairs baseline")
    all_pairs <- read_all_pairs_file(all_pairs_file)
    log_info("Reading SNP/Gene pairs comparison")
    all_pairs_comparison <- read_all_pairs_file(all_pairs_file_comparison)

    log_info("Adding cross-experiment slopes for sign test.")

    augmented <- add_cross_slopes_for_sign_test(
        index_dt = index_dt,
        index_comparison_dt = index_comparison_dt,
        all_pairs = all_pairs,
        all_pairs_comparison = all_pairs_comparison,
        all_pairs_slope_col = "slope",
        out_col = "slope_cross"
    )

    log_info("Writing sign-test cache.")

    utils::write.table(augmented$index_dt, file = paths$index_path,
                       sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

    utils::write.table(augmented$index_comparison_dt, file = paths$comparison_path,
                       sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

    augmented
}

get_eqtl_cache_dir <- function(cache_dir = NULL) {
    ## 1) Explicit argument override
    if (!is.null(cache_dir)) {
        base <- cache_dir
    } else {
        ## 2) Environment variable override
        env_override <- Sys.getenv("EQTL_CACHE_DIR", unset = NA_character_)
        if (!is.na(env_override) && nzchar(env_override)) {
            base <- env_override
        } else {
            ## 3) Default user cache directory
            base <- tools::R_user_dir("eqtl_analysis", which = "cache")
        }
    }

    cache_dir_final <- file.path(base)

    if (!dir.exists(cache_dir_final)) {
        dir.create(cache_dir_final, recursive = TRUE, showWarnings = FALSE)
    }

    cache_dir_final
}




#################################
# Aggregate summary
#################################
make_compare_index_eqtls_stats_df <- function(stats_genes_tested,
                                              stats_eGenes,
                                              effect_size_stats,
                                              qval_stats) {
    stopifnot(is.data.frame(stats_genes_tested), nrow(stats_genes_tested) == 1)
    stopifnot(is.data.frame(stats_eGenes), nrow(stats_eGenes) == 1)
    stopifnot(is.data.frame(effect_size_stats), nrow(effect_size_stats) == 1)
    stopifnot(is.data.frame(qval_stats), nrow(qval_stats) == 1)

    data.frame(
        baseline_name = as.character(qval_stats$baseline_name[1]),
        comparison_name = as.character(qval_stats$comparison_name[1]),
        n_union_genes = as.numeric(stats_genes_tested$n_union[1]),
        n_intersect_genes = as.numeric(stats_genes_tested$n_intersect[1]),
        n_union_egenes = as.numeric(stats_eGenes$n_union[1]),
        n_intersect_egenes = as.numeric(stats_eGenes$n_intersect[1]),
        n_baseline_only_egenes = as.numeric(effect_size_stats$n_baseline_only_egenes[1]),
        n_comparison_only_egenes = as.numeric(effect_size_stats$n_comparison_only_egenes[1]),
        cor_val = as.numeric(qval_stats$cor_val[1]),
        yield = as.numeric(qval_stats$yield[1])
    )
}


#################################
# I/O functions
#################################


read_index_file<-function (index_file, colsToKeep=c("gene_name", "variant_id", "slope", "slope_cross", "pval_nominal", "pval_beta", "qval")) {
    index_dt <- data.table::fread(index_file)

    ## Keep only requested columns that actually exist in the file
    cols_present <- intersect(colsToKeep, names(index_dt))
    index_dt <- index_dt[, cols_present, with = FALSE]

    data.table::setDF(index_dt)
    index_dt
}

read_all_pairs_file<-function (all_pairs_file, colsToKeep=c("phenotype_id", "variant_id", "slope", "pval_nominal")) {
    # cmd <- sprintf(
    #     "gzip -cd %s | grep -v -E '^[[:space:]]*#'",
    #     shQuote(all_pairs_file)
    # )
    # all_pairs <- data.table::fread(cmd = cmd)

    #simplified!
    all_pairs <- data.table::fread(all_pairs_file)

    ## Keep only requested columns that actually exist in the file
    cols_present <- intersect(colsToKeep, names(all_pairs))
    all_pairs <- all_pairs[, cols_present, with = FALSE]

    data.table::setDF(all_pairs)
    return (all_pairs)
}
