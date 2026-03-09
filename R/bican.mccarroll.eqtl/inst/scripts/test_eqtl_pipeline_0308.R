## ============================================================
## Test pipeline 0308
##
## Produces:
##   (1) K-means clustering heatmap
##   (2) Pairwise correlation plot
##   (3) Distance-to-TSS boxplot
##   (4) Gene-snp plots (skip for local testing)
##
## Runs R Part 1 → Python (K-means) → R Part 2 end-to-end.
## ============================================================

## -----------------------
## Config
## -----------------------

force <- TRUE

## -----------------------
## Paths
## -----------------------

base_dir   <- "/Users/tracyyuan/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls"
eqtl_dir   <- file.path(base_dir, "results/LEVEL_3")
out_dir    <- file.path(base_dir, "manuscript_test_0308")
qval       <- 0.01

region_cell_type_path    <- file.path(base_dir, "manuscript_data/region_cell_type.tsv")
combined_expression_path <- file.path(out_dir, "combined_tpm_expression_across_cell_types.tsv")

## Path to standalone scripts (same directory as this file)
script_dir <- "/Users/tracyyuan/Desktop/McCarroll_Lab/BICAN/bican-mccarroll-manuscript1/R/bican.mccarroll.eqtl/inst/scripts"
python_script <- "/Users/tracyyuan/Desktop/McCarroll_Lab/BICAN/bican-mccarroll-manuscript1/python/bican_mccarroll_eqtl/scripts/test_eqtl_pipeline_0308.py"

## VCF (for gene-snp plots)
vcf_path <- "/Users/tracyyuan/broad/bican_um1_mccarroll/vcfs/2025-05-05/gvs_concat_outputs_2025-05-05T14-10-02.donors_renamed_filtered_norm.vcf.gz"

## -----------------------
## Helpers
## -----------------------

.run_step <- function(step_label, output_path, fun) {
    cat("\n=====", step_label, "=====\n")
    if (!isTRUE(force) && !is.null(output_path) && file.exists(output_path)) {
        cat("  SKIPPED: file already exists at", output_path, "\n")
        return(invisible(NULL))
    }
    out <- fun()
    if (!is.null(output_path)) cat("  Written to:", output_path, "\n")
    invisible(out)
}

## ===========================================================
## R Part 1: Generate intermediate matrices
## ===========================================================

## -----------------------
## Step 1: get_egene_union_pairs
## -----------------------

egene_path <- file.path(out_dir, paste0("egene_union_pairs_qval_", qval, ".tsv"))

.run_step(
    "Step 1: get_egene_union_pairs",
    egene_path,
    function() {
        bican.mccarroll.eqtl::get_egene_union_pairs(
            eqtl_dir              = eqtl_dir,
            region_cell_type_path = region_cell_type_path,
            qval_threshold        = qval,
            output_path           = egene_path
        )
    }
)

## -----------------------
## Step 2: get_slope_matrix
## -----------------------

slope_path <- file.path(out_dir, paste0("slope_matrix_qval_", qval, ".tsv"))

.run_step(
    "Step 2: get_slope_matrix",
    slope_path,
    function() {
        bican.mccarroll.eqtl::get_slope_matrix(
            eqtl_dir               = eqtl_dir,
            region_cell_type_path  = region_cell_type_path,
            egene_union_pairs_path = egene_path,
            output_path            = slope_path
        )
    }
)

## -----------------------
## Step 3: get_pval_nominal_matrix
## -----------------------

pval_path <- file.path(out_dir, paste0("pval_nominal_matrix_qval_", qval, ".tsv"))

.run_step(
    "Step 3: get_pval_nominal_matrix",
    pval_path,
    function() {
        bican.mccarroll.eqtl::get_pval_nominal_matrix(
            eqtl_dir               = eqtl_dir,
            region_cell_type_path  = region_cell_type_path,
            egene_union_pairs_path = egene_path,
            output_path            = pval_path
        )
    }
)

## -----------------------
## Step 4: get_pval_nominal_threshold_matrix
## -----------------------

pval_thresh_path <- file.path(out_dir, paste0("pval_nominal_threshold_matrix_qval_", qval, ".tsv"))

.run_step(
    "Step 4: get_pval_nominal_threshold_matrix",
    pval_thresh_path,
    function() {
        bican.mccarroll.eqtl::get_pval_nominal_threshold_matrix(
            eqtl_dir               = eqtl_dir,
            region_cell_type_path  = region_cell_type_path,
            egene_union_pairs_path = egene_path,
            output_path            = pval_thresh_path
        )
    }
)

## -----------------------
## Step 5: get_index_snp_slope_matrix_with_impute
## -----------------------

index_snp_path <- file.path(out_dir, paste0("index_snp_slope_matrix_with_zero_impute_qval_", qval, ".tsv"))

.run_step(
    "Step 5: get_index_snp_slope_matrix_with_impute",
    index_snp_path,
    function() {
        bican.mccarroll.eqtl::get_index_snp_slope_matrix_with_impute(
            slope_matrix_path = slope_path,
            output_path       = index_snp_path
        )
    }
)

## -----------------------
## Step 6: get_cell_type_pairwise_cor_matrix
## -----------------------

r_squared_path <- file.path(out_dir, paste0("cell_type_pairwise_r_squared_qval_", qval, ".tsv"))

.run_step(
    "Step 6: get_cell_type_pairwise_cor_matrix",
    r_squared_path,
    function() {
        bican.mccarroll.eqtl::get_cell_type_pairwise_cor_matrix(
            slope_matrix_path                  = slope_path,
            pval_nominal_matrix_path           = pval_path,
            pval_nominal_threshold_matrix_path = pval_thresh_path,
            egene_union_pairs_path             = egene_path,
            region_cell_type_path              = region_cell_type_path,
            output_path                        = r_squared_path
        )
    }
)

## -----------------------
## Step 7: plot_cell_type_pairwise_cor
## -----------------------

cor_plot_path <- file.path(out_dir, paste0("cell_type_cor_plot_qval_", qval, ".svg"))

.run_step(
    "Step 7: plot_cell_type_pairwise_cor",
    cor_plot_path,
    function() {
        bican.mccarroll.eqtl::plot_cell_type_pairwise_cor(
            r_squared_path = r_squared_path,
            output_path    = cor_plot_path
        )
    }
)

## -----------------------
## Step 8: get_index_snp_start_distance (standalone script)
## -----------------------

start_distance_path <- file.path(out_dir, paste0("index_snp_start_distance_qval_", qval, ".tsv"))

cat("\n===== Step 8: get_index_snp_start_distance =====\n")
if (isTRUE(force) || !file.exists(start_distance_path)) {
    system2("Rscript", c(
        file.path(script_dir, "get_index_snp_start_distance.R"),
        eqtl_dir,
        region_cell_type_path,
        index_snp_path,
        start_distance_path
    ))
    cat("  Written to:", start_distance_path, "\n")
} else {
    cat("  SKIPPED: file already exists at", start_distance_path, "\n")
}

## -----------------------
## Step 9: combine_expression_across_cell_types (needed for gene-snp plots)
## -----------------------

.run_step(
    "Step 9: combine_expression_across_cell_types",
    combined_expression_path,
    function() {
        bican.mccarroll.eqtl::combine_expression_across_cell_types(
            eqtl_dir              = eqtl_dir,
            region_cell_type_path = region_cell_type_path,
            output_path           = combined_expression_path
        )
    }
)

## ===========================================================
## Python: K-means clustering
## ===========================================================

cat("\n===== Step 10: K-means clustering (Python) =====\n")
exit_code <- system2("python3", python_script)
if (exit_code != 0) stop("Python K-means pipeline failed")

## ===========================================================
## R Part 2: Plots that depend on cluster assignments
## ===========================================================

## -----------------------
## Step 11: plot_eqtl_distance_to_tss_boxplot
## -----------------------

cluster_assignments_path <- file.path(out_dir, paste0("cluster_assignments_qval_", qval, "_k11.tsv"))
boxplot_output <- file.path(out_dir, "eqtl_distance_to_tss_boxplot.svg")

cat("\n===== Step 11: plot_eqtl_distance_to_tss_boxplot =====\n")
system2("Rscript", c(
    file.path(script_dir, "plot_eqtl_distance_to_tss_boxplot.R"),
    index_snp_path,
    cluster_assignments_path,
    start_distance_path,
    boxplot_output,
    "horizontal"
))
cat("  Written to:", boxplot_output, "\n")

## -----------------------
## Step 12: plot_gene_snp
## -----------------------

cat("\n===== Step 12: plot_gene_snp =====\n")

gene_snp_cases <- list(
    list(gene = "XRRA1",  chr = "chr11", pos = 74935168),
    list(gene = "NPAS3",  chr = "chr14", pos = 32935820),
    list(gene = "CEP112", chr = "chr17", pos = 66192315)
)

for (case in gene_snp_cases) {
    out_file <- file.path(out_dir, paste0(case$gene, "_", case$chr, "_", case$pos, ".svg"))

    if (!isTRUE(force) && file.exists(out_file)) {
        cat("  SKIPPED:", case$gene, "file already exists at", out_file, "\n")
        next
    }

    bican.mccarroll.eqtl::plot_gene_snp(
        gene            = case$gene,
        chr             = case$chr,
        pos             = case$pos,
        vcf_path        = vcf_path,
        expression_path = combined_expression_path,
        output_path     = out_file
    )
    cat("  ", case$gene, "saved to:", out_file, "\n")
}

## -----------------------
## Done
## -----------------------

cat("\n===== All steps completed! =====\n")
