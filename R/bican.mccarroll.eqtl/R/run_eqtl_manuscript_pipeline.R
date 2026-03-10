#' Run the manuscript eQTL analysis pipeline
#'
#' Runs the manuscript-specific end-to-end eQTL pipeline used to generate the
#' intermediate matrices and figures for the paper. The pipeline performs
#' preprocessing of tensorQTL outputs, clustering of eQTL effect-size profiles,
#' and generation of summary visualizations used in the manuscript.
#'
#' The function assumes that intermediate and output file names follow the
#' manuscript pipeline naming convention. Only the top-level input locations,
#' core analysis settings, and gene-SNP plotting cases are exposed as arguments;
#' all other paths are inferred from these inputs.
#'
#' When \code{force = FALSE}, a step is skipped if its expected output file
#' already exists. When \code{force = TRUE}, all steps are rerun and existing
#' outputs may be overwritten by the underlying functions.
#'
#' The K-means clustering step is implemented in Python and is executed via the
#' external command \code{test-eqtl-pipeline}. See the Python component of the
#' repository for installation instructions.
#'
#' @details
#' The pipeline proceeds through the following stages:
#'
#' \strong{Step 1: Identify the union of eGene-variant pairs.}
#' For each cell type-region combination, the pipeline reads the tensorQTL
#' cis-eQTL index file, filters to significant eGenes using the specified
#' q-value threshold, and collects the corresponding
#' \code{(phenotype_id, variant_id)} pairs. These pairs are combined across all
#' groups and deduplicated so that each unique gene-variant pair appears once,
#' retaining the smallest observed q-value.
#'
#' \strong{Steps 2-4: Construct cross-cell-type eQTL statistic matrices.}
#' Using the union of eGene-variant pairs defined in Step 1, the pipeline builds
#' parallel matrices of effect sizes (slopes), nominal p-values, and nominal
#' significance thresholds across all cell type-region combinations. Each row
#' corresponds to a gene-variant pair and each column corresponds to a specific
#' cell type and brain region.
#'
#' \strong{Step 5: Select index variants and construct the clustering matrix.}
#' The full slope matrix is reduced to one representative index variant per
#' gene by selecting the variant with the largest absolute effect across cell
#' types. Effect directions are reoriented so that the strongest effect for each
#' gene is positive, and missing values are replaced with zero to produce a
#' complete matrix suitable for clustering.
#'
#' \strong{Steps 6-7: Quantify and visualize similarity of eQTL effects.}
#' Pairwise Spearman correlations of eQTL effect sizes are computed between
#' cell type-region combinations using gene-variant pairs that are nominally
#' significant in at least one of the two groups being compared. The resulting
#' R-squared matrix summarizes cross-group similarity in regulatory effects and
#' is visualized as a clustered heatmap.
#'
#' \strong{Step 8: Plot example eQTL genotype-expression relationships.}
#' For selected gene-variant pairs, donor genotypes are extracted from the VCF
#' and combined with expression measurements across cell types. Expression is
#' plotted by genotype for each cell type-region combination and annotated with
#' a multiple-testing-adjusted association p-value.
#'
#' \strong{Step 9: Cluster genes by eQTL effect-size profile.}
#' The sign-adjusted index-SNP slope matrix from Step 5 is clustered using
#' K-means to group genes with similar cross-cell-type regulatory patterns.
#' Clusters are then manually reordered for visualization and displayed as a
#' heatmap to highlight recurring patterns of regulatory effects across cell
#' types and brain regions.
#'
#' @param eqtl_dir Character scalar giving the directory containing the
#'   tensorQTL eQTL result directories.
#' @param out_dir Character scalar giving the directory where all intermediate
#'   outputs and plots will be written.
#' @param region_cell_type_path Character scalar giving the path to a
#'   tab-delimited file describing the cell type-region combinations to include
#'   in the analysis. The file must contain columns \code{cell_type} and
#'   \code{region}. Each row defines one group, and these values are combined as
#'   \code{<cell_type>__<region>} to locate tensorQTL result directories within
#'   \code{eqtl_dir}.
#' @param vcf_path Character scalar giving the path to the VCF file used to
#'   obtain donor genotypes for gene-SNP plots.
#' @param K Integer scalar giving the number of K-means clusters used to group
#'   genes by eQTL effect-size profile.
#' @param cluster_order Character scalar giving the desired cluster order as a
#'   comma-separated string used for final heatmap visualization.
#' @param qval Numeric scalar giving the q-value threshold used to define the
#'   union set of eGenes.
#' @param force Logical scalar indicating whether to rerun steps even if their
#'   expected output files already exist.
#' @param gene_snp_cases List of named lists defining gene-SNP plots. Each
#'   element must contain \code{gene}, \code{chr}, and \code{pos}.
#'
#' @return A named list containing the inferred output paths for all pipeline
#'   products generated during the run.
#'
#' @export
run_eqtl_manuscript_pipeline <- function(
        eqtl_dir,
        out_dir,
        region_cell_type_path,
        vcf_path,
        K,
        cluster_order,
        qval = 0.01,
        force = FALSE,
        gene_snp_cases = list(
            list(gene = "XRRA1", chr = "chr11", pos = 74935168),
            list(gene = "NPAS3", chr = "chr14", pos = 32935820),
            list(gene = "CEP112", chr = "chr17", pos = 66192315)
        )) {

    .validate_run_eqtl_manuscript_pipeline_inputs(
        eqtl_dir = eqtl_dir,
        out_dir = out_dir,
        region_cell_type_path = region_cell_type_path,
        vcf_path = vcf_path,
        K = K,
        cluster_order = cluster_order,
        qval = qval,
        force = force,
        gene_snp_cases = gene_snp_cases
    )

    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }

    paths <- .get_run_eqtl_manuscript_pipeline_paths(
        out_dir = out_dir,
        qval = qval,
        K = K
    )

    .run_eqtl_manuscript_pipeline_step(
        step_label = "Step 1: get_egene_union_pairs",
        output_path = paths$egene_path,
        force = force,
        fun = function() {
            bican.mccarroll.eqtl::get_egene_union_pairs(
                eqtl_dir = eqtl_dir,
                region_cell_type_path = region_cell_type_path,
                qval_threshold = qval,
                output_path = paths$egene_path
            )
        }
    )

    .run_eqtl_manuscript_pipeline_step(
        step_label = "Step 2: get_slope_matrix",
        output_path = paths$slope_path,
        force = force,
        fun = function() {
            bican.mccarroll.eqtl::get_slope_matrix(
                eqtl_dir = eqtl_dir,
                region_cell_type_path = region_cell_type_path,
                egene_union_pairs_path = paths$egene_path,
                output_path = paths$slope_path
            )
        }
    )

    .run_eqtl_manuscript_pipeline_step(
        step_label = "Step 3: get_pval_nominal_matrix",
        output_path = paths$pval_path,
        force = force,
        fun = function() {
            bican.mccarroll.eqtl::get_pval_nominal_matrix(
                eqtl_dir = eqtl_dir,
                region_cell_type_path = region_cell_type_path,
                egene_union_pairs_path = paths$egene_path,
                output_path = paths$pval_path
            )
        }
    )

    .run_eqtl_manuscript_pipeline_step(
        step_label = "Step 4: get_pval_nominal_threshold_matrix",
        output_path = paths$pval_thresh_path,
        force = force,
        fun = function() {
            bican.mccarroll.eqtl::get_pval_nominal_threshold_matrix(
                eqtl_dir = eqtl_dir,
                region_cell_type_path = region_cell_type_path,
                egene_union_pairs_path = paths$egene_path,
                output_path = paths$pval_thresh_path
            )
        }
    )

    .run_eqtl_manuscript_pipeline_step(
        step_label = "Step 5: get_index_snp_slope_matrix_with_impute",
        output_path = paths$index_snp_path,
        force = force,
        fun = function() {
            bican.mccarroll.eqtl::get_index_snp_slope_matrix_with_impute(
                slope_matrix_path = paths$slope_path,
                output_path = paths$index_snp_path
            )
        }
    )

    .run_eqtl_manuscript_pipeline_step(
        step_label = "Step 6: get_cell_type_pairwise_cor_matrix",
        output_path = paths$r_squared_path,
        force = force,
        fun = function() {
            bican.mccarroll.eqtl::get_cell_type_pairwise_cor_matrix(
                slope_matrix_path = paths$slope_path,
                pval_nominal_matrix_path = paths$pval_path,
                pval_nominal_threshold_matrix_path = paths$pval_thresh_path,
                egene_union_pairs_path = paths$egene_path,
                region_cell_type_path = region_cell_type_path,
                output_path = paths$r_squared_path
            )
        }
    )

    .run_eqtl_manuscript_pipeline_step(
        step_label = "Step 7: plot_cell_type_pairwise_cor",
        output_path = paths$cor_plot_path,
        force = force,
        fun = function() {
            bican.mccarroll.eqtl::plot_cell_type_pairwise_cor(
                r_squared_path = paths$r_squared_path,
                output_path = paths$cor_plot_path
            )
        }
    )

    #this does it's own logging, and is step 8.
    .run_eqtl_manuscript_pipeline_gene_snp_plots(
        out_dir = out_dir,
        vcf_path = vcf_path,
        expression_path = paths$combined_expression_path,
        gene_snp_cases = gene_snp_cases,
        force = force
    )

    # .run_eqtl_manuscript_pipeline_step(
    #     step_label = "Step 8: get_index_snp_start_distance",
    #     output_path = paths$start_distance_path,
    #     force = force,
    #     fun = function() {
    #         bican.mccarroll.eqtl::get_index_snp_start_distance(
    #             eqtl_dir = eqtl_dir,
    #             region_cell_type_path = region_cell_type_path,
    #             index_snp_matrix_path = paths$index_snp_path,
    #             output_path = paths$start_distance_path
    #         )
    #     }
    # )

    # .run_eqtl_manuscript_pipeline_step(
    #     step_label = "Step 9: combine_expression_across_cell_types",
    #     output_path = paths$combined_expression_path,
    #     force = force,
    #     fun = function() {
    #         bican.mccarroll.eqtl::combine_expression_across_cell_types(
    #             eqtl_dir = eqtl_dir,
    #             region_cell_type_path = region_cell_type_path,
    #             output_path = paths$combined_expression_path
    #         )
    #     }
    # )

    .run_eqtl_manuscript_pipeline_kmeans(
        out_dir = out_dir,
        K = K,
        cluster_order = cluster_order,
        force = force
    )

    # .run_eqtl_manuscript_pipeline_step(
    #     step_label = "Step 11: plot_eqtl_distance_to_tss_boxplot",
    #     output_path = paths$boxplot_output,
    #     force = force,
    #     fun = function() {
    #         bican.mccarroll.eqtl::plot_eqtl_distance_to_tss_boxplot(
    #             index_snp_matrix_path = paths$index_snp_path,
    #             cluster_assignments_path = paths$cluster_assignments_path,
    #             start_distance_path = paths$start_distance_path,
    #             output_path = paths$boxplot_output,
    #             orientation = "horizontal"
    #         )
    #     }
    # )

    cat("\n===== All steps completed! =====\n")

    return(paths)
}

#' Run the manuscript pipeline with default manuscript settings
#'
#' Convenience wrapper around `run_eqtl_manuscript_pipeline()` that fills in the
#' manuscript paths and default analysis settings used for the paper. Any of the
#' defaults can be overridden, which is useful for regression tests or writing
#' outputs to an alternate directory.
#'
#' @inheritParams run_eqtl_manuscript_pipeline
#' @export
run_eqtl_manuscript_pipeline_defaults <- function(
        eqtl_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/LEVEL_3",
        out_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/eqtl_analysis_pipeline_run_jim",
        region_cell_type_path = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/manuscript_data/region_cell_type.tsv",
        vcf_path = "/broad/bican_um1_mccarroll/vcfs/2025-05-05/gvs_concat_outputs_2025-05-05T14-10-02.donors_renamed_filtered_norm.vcf.gz",
        K = 13,
        cluster_order = "11,0,5,4,2,12,9,1,6,3,7,10,8",
        qval = 0.01,
        force = TRUE,
        gene_snp_cases = list(
            list(gene = "XRRA1", chr = "chr11", pos = 74935168),
            list(gene = "NPAS3", chr = "chr14", pos = 32935820),
            list(gene = "CEP112", chr = "chr17", pos = 66192315)
        )) {

    run_eqtl_manuscript_pipeline(
        eqtl_dir = eqtl_dir,
        out_dir = out_dir,
        region_cell_type_path = region_cell_type_path,
        vcf_path = vcf_path,
        K = K,
        cluster_order = cluster_order,
        qval = qval,
        force = force,
        gene_snp_cases = gene_snp_cases
    )
}

.run_eqtl_manuscript_pipeline_step <- function(
        step_label,
        output_path,
        force,
        fun) {

    cat("\n=====", step_label, "=====\n")

    if (!isTRUE(force) && !is.null(output_path) && file.exists(output_path)) {
        cat("  SKIPPED: file already exists at", output_path, "\n")
        return(invisible(NULL))
    }

    out <- fun()

    if (!is.null(output_path)) {
        cat("  Written to:", output_path, "\n")
    }

    invisible(out)
}

.run_eqtl_manuscript_pipeline_kmeans <- function(
        out_dir,
        K,
        cluster_order,
        force) {

    cat("\n===== Step 10: K-means clustering (Python) =====\n")

    args <- c(
        "--out-dir", out_dir,
        sprintf("--K=%d", K),
        "--desired-order", cluster_order
    )

    if (isTRUE(force)) {
        args <- c(args, "--force")
    }

    cat("Running:", paste("test-eqtl-pipeline", paste(args, collapse = " ")), "\n")

    exit_code <- system2("test-eqtl-pipeline", args)

    if (!identical(exit_code, 0L)) {
        stop("Python K-means pipeline failed")
    }

    invisible(NULL)
}

.run_eqtl_manuscript_pipeline_gene_snp_plots <- function(
        out_dir,
        vcf_path,
        expression_path,
        gene_snp_cases,
        force) {

    cat("\n===== Step 12: plot_gene_snp =====\n")

    for (case in gene_snp_cases) {
        out_file <- file.path(
            out_dir,
            paste0(case$gene, "_", case$chr, "_", case$pos, ".svg")
        )

        if (!isTRUE(force) && file.exists(out_file)) {
            cat("  SKIPPED:", case$gene, "file already exists at", out_file, "\n")
            next
        }

        bican.mccarroll.eqtl::plot_gene_snp(
            gene = case$gene,
            chr = case$chr,
            pos = case$pos,
            vcf_path = vcf_path,
            expression_path = expression_path,
            output_path = out_file
        )

        cat("  ", case$gene, "saved to:", out_file, "\n")
    }

    invisible(NULL)
}

.get_run_eqtl_manuscript_pipeline_paths <- function(
        out_dir,
        qval,
        K) {

    list(
        egene_path = file.path(
            out_dir,
            paste0("egene_union_pairs_qval_", qval, ".tsv")
        ),
        slope_path = file.path(
            out_dir,
            paste0("slope_matrix_qval_", qval, ".tsv")
        ),
        pval_path = file.path(
            out_dir,
            paste0("pval_nominal_matrix_qval_", qval, ".tsv")
        ),
        pval_thresh_path = file.path(
            out_dir,
            paste0("pval_nominal_threshold_matrix_qval_", qval, ".tsv")
        ),
        index_snp_path = file.path(
            out_dir,
            paste0("index_snp_slope_matrix_with_zero_impute_qval_", qval, ".tsv")
        ),
        r_squared_path = file.path(
            out_dir,
            paste0("cell_type_pairwise_r_squared_qval_", qval, ".tsv")
        ),
        cor_plot_path = file.path(
            out_dir,
            paste0("cell_type_cor_plot_qval_", qval, ".svg")
        ),
        start_distance_path = file.path(
            out_dir,
            paste0("index_snp_start_distance_qval_", qval, ".tsv")
        ),
        combined_expression_path = file.path(
            out_dir,
            "combined_tpm_expression_across_cell_types.tsv"
        ),
        cluster_assignments_path = file.path(
            out_dir,
            paste0("cluster_assignments_qval_", qval, "_k", K, ".tsv")
        ),
        boxplot_output = file.path(
            out_dir,
            "eqtl_distance_to_tss_boxplot.svg"
        )
    )
}

.validate_run_eqtl_manuscript_pipeline_inputs <- function(
        eqtl_dir,
        out_dir,
        region_cell_type_path,
        vcf_path,
        K,
        cluster_order,
        qval,
        force,
        gene_snp_cases) {

    if (!is.character(eqtl_dir) || length(eqtl_dir) != 1L || is.na(eqtl_dir)) {
        stop("'eqtl_dir' must be a single character string")
    }

    if (!dir.exists(eqtl_dir)) {
        stop("'eqtl_dir' does not exist: ", eqtl_dir)
    }

    if (!is.character(out_dir) || length(out_dir) != 1L || is.na(out_dir)) {
        stop("'out_dir' must be a single character string")
    }

    if (!is.character(region_cell_type_path) ||
        length(region_cell_type_path) != 1L ||
        is.na(region_cell_type_path)) {
        stop("'region_cell_type_path' must be a single character string")
    }

    if (!file.exists(region_cell_type_path)) {
        stop("'region_cell_type_path' does not exist: ", region_cell_type_path)
    }

    if (!is.character(vcf_path) || length(vcf_path) != 1L || is.na(vcf_path)) {
        stop("'vcf_path' must be a single character string")
    }

    if (!file.exists(vcf_path)) {
        stop("'vcf_path' does not exist: ", vcf_path)
    }

    if (!is.numeric(K) || length(K) != 1L || is.na(K) || K <= 0) {
        stop("'K' must be a single positive number")
    }

    if (K != as.integer(K)) {
        stop("'K' must be an integer value")
    }

    if (!is.character(cluster_order) ||
        length(cluster_order) != 1L ||
        is.na(cluster_order)) {
        stop("'cluster_order' must be a single character string")
    }

    if (!is.numeric(qval) || length(qval) != 1L || is.na(qval)) {
        stop("'qval' must be a single numeric value")
    }

    if (qval < 0 || qval > 1) {
        stop("'qval' must be between 0 and 1")
    }

    if (!is.logical(force) || length(force) != 1L || is.na(force)) {
        stop("'force' must be TRUE or FALSE")
    }

    if (!is.list(gene_snp_cases) || length(gene_snp_cases) == 0L) {
        stop("'gene_snp_cases' must be a non-empty list")
    }

    for (i in seq_along(gene_snp_cases)) {
        case <- gene_snp_cases[[i]]

        if (!is.list(case)) {
            stop("'gene_snp_cases[[", i, "]]' must be a list")
        }

        required_names <- c("gene", "chr", "pos")

        if (!all(required_names %in% names(case))) {
            stop(
                "'gene_snp_cases[[", i,
                "]]' must contain named elements: gene, chr, pos"
            )
        }

        if (!is.character(case$gene) || length(case$gene) != 1L || is.na(case$gene)) {
            stop("'gene_snp_cases[[", i, "]]$gene' must be a single character string")
        }

        if (!is.character(case$chr) || length(case$chr) != 1L || is.na(case$chr)) {
            stop("'gene_snp_cases[[", i, "]]$chr' must be a single character string")
        }

        if (!is.numeric(case$pos) || length(case$pos) != 1L || is.na(case$pos)) {
            stop("'gene_snp_cases[[", i, "]]$pos' must be a single numeric value")
        }

        if (case$pos != as.integer(case$pos) || case$pos < 1) {
            stop("'gene_snp_cases[[", i, "]]$pos' must be a positive integer value")
        }
    }

    invisible(NULL)
}
