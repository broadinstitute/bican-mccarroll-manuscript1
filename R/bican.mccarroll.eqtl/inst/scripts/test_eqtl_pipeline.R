## ============================================================
## Test script: Run all bican.mccarroll.eqtl R functions in order
## ============================================================
##
## Usage:
##   1. SSH into Broad server
##   2. Start R
##   3. Install/update the package:
##
##      remotes::install_github(
##        "broadinstitute/bican-mccarroll-manuscript1",
##        subdir = "R/bican.mccarroll.eqtl",
##        ref = "ty_eqtl",
##        dependencies = TRUE
##      )
##
##   4. source("test_eqtl_pipeline.R")
##
## Behavior:
##   - Set force = FALSE to skip steps whose output file already exists.
##   - Set force = TRUE to always run and overwrite output files.
## ============================================================

## -----------------------
## Config
## -----------------------

force <- FALSE  # set TRUE to rerun everything and overwrite existing outputs

## -----------------------
## Paths
## -----------------------

#TODO: these paths should all be mutable.
eqtl_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/LEVEL_3"
region_cell_type_path <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/manuscript_data/region_cell_type.tsv"

# Output directory
out_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/manuscript_data_test"

qval <- 0.01

## VCF and expression paths (for plot_gene_snp)
vcf_path <- "/broad/bican_um1_mccarroll/vcfs/2025-05-05/gvs_concat_outputs_2025-05-05T14-10-02.donors_renamed_filtered_norm.vcf.gz"
combined_expression_path <- file.path(out_dir, "combined_tpm_expression_across_cell_types.tsv")

## Coloc paths (for get_sig_coloc and plot_fisher_exact)
ad_coloc_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/LEVEL_3/coloc/AD_2022"
scz_coloc_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/eqtls/results/LEVEL_3/coloc/SCZ_eur"

## -----------------------
## Helpers
## -----------------------

.run_step_table <- function(step_label, output_path, fun) {
    cat("\n=====", step_label, "=====\n")

    if (!isTRUE(force) && !is.null(output_path) && file.exists(output_path)) {
        cat("  SKIPPED: file already exists at", output_path, "\n")
        return(invisible(NULL))
    }

    out <- fun()

    if (!is.null(output_path)) {
        if (is.data.frame(out)) {
            cat("  Rows:", nrow(out), " Cols:", ncol(out), " Written to:", output_path, "\n")
        } else if (is.matrix(out)) {
            cat("  Matrix:", nrow(out), "x", ncol(out), " Written to:", output_path, "\n")
        } else {
            cat("  Written to:", output_path, "\n")
        }
    }

    invisible(out)
}

.run_step_plot <- function(step_label, output_path, fun) {
    cat("\n=====", step_label, "=====\n")

    if (!isTRUE(force) && file.exists(output_path)) {
        cat("  SKIPPED: file already exists at", output_path, "\n")
        return(invisible(NULL))
    }

    fun()
    cat("  Written to:", output_path, "\n")
    invisible(NULL)
}

## -----------------------
## Step 1: get_egene_union_pairs
## -----------------------

egene_path <- file.path(out_dir, paste0("egene_union_pairs_qval_", qval, ".tsv"))

egene_dt <- .run_step_table(
    step_label = "Step 1: get_egene_union_pairs",
    output_path = egene_path,
    fun = function() {
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

slope_dt <- .run_step_table(
    step_label = "Step 2: get_slope_matrix",
    output_path = slope_path,
    fun = function() {
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

pval_dt <- .run_step_table(
    step_label = "Step 3: get_pval_nominal_matrix",
    output_path = pval_path,
    fun = function() {
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

pval_thresh_dt <- .run_step_table(
    step_label = "Step 4: get_pval_nominal_threshold_matrix",
    output_path = pval_thresh_path,
    fun = function() {
        bican.mccarroll.eqtl::get_pval_nominal_threshold_matrix(
            eqtl_dir               = eqtl_dir,
            region_cell_type_path  = region_cell_type_path,
            egene_union_pairs_path = egene_path,
            output_path            = pval_thresh_path
        )
    }
)

## -----------------------
## Step 5: get_index_snp_slope_matrix_with_median_impute
## -----------------------

index_snp_path <- file.path(out_dir, paste0("index_snp_slope_matrix_with_median_impute_qval_", qval, ".tsv"))

index_snp_dt <- .run_step_table(
    step_label = "Step 5: get_index_snp_slope_matrix_with_median_impute",
    output_path = index_snp_path,
    fun = function() {
        bican.mccarroll.eqtl::get_index_snp_slope_matrix_with_median_impute(
            slope_matrix_path = slope_path,
            output_path       = index_snp_path
        )
    }
)

## -----------------------
## Step 6: get_cell_type_pairwise_cor_matrix
## -----------------------

r_squared_path <- file.path(out_dir, paste0("cell_type_pairwise_r_squared_qval_", qval, ".tsv"))

r_squared <- .run_step_table(
    step_label = "Step 6: get_cell_type_pairwise_cor_matrix",
    output_path = r_squared_path,
    fun = function() {
        bican.mccarroll.eqtl::get_cell_type_pairwise_cor_matrix(
            slope_matrix_path                    = slope_path,
            pval_nominal_matrix_path             = pval_path,
            pval_nominal_threshold_matrix_path   = pval_thresh_path,
            egene_union_pairs_path               = egene_path,
            region_cell_type_path                = region_cell_type_path,
            output_path                          = r_squared_path
        )
    }
)

## -----------------------
## Step 7: plot_cell_type_pairwise_cor
## -----------------------

cor_plot_path <- file.path(out_dir, paste0("cell_type_cor_plot_qval_", qval, ".png"))

.run_step_plot(
    step_label = "Step 7: plot_cell_type_pairwise_cor",
    output_path = cor_plot_path,
    fun = function() {
        bican.mccarroll.eqtl::plot_cell_type_pairwise_cor(
            r_squared_path = r_squared_path,
            output_path    = cor_plot_path
        )
    }
)

## -----------------------
## Step 8: get_heatmap_index_snp_median_expression
## -----------------------

median_expr_path <- file.path(out_dir, paste0("heatmap_index_snp_median_expression_qval_", qval, ".tsv"))

median_expr_dt <- .run_step_table(
    step_label = "Step 8: get_heatmap_index_snp_median_expression",
    output_path = median_expr_path,
    fun = function() {
        bican.mccarroll.eqtl::get_heatmap_index_snp_median_expression(
            index_snp_path        = index_snp_path,
            region_cell_type_path = region_cell_type_path,
            expression_path       = combined_expression_path,
            output_path           = median_expr_path
        )
    }
)

## -----------------------
## Step 9: combine_expression_across_cell_types
## -----------------------

combined_dt <- .run_step_table(
    step_label = "Step 9: combine_expression_across_cell_types",
    output_path = combined_expression_path,
    fun = function() {
        bican.mccarroll.eqtl::combine_expression_across_cell_types(
            eqtl_dir              = eqtl_dir,
            region_cell_type_path = region_cell_type_path,
            output_path           = combined_expression_path
        )
    }
)

## -----------------------
## Step 10: get_sig_coloc (AD and SCZ)
## -----------------------

cat("\n===== Step 10: get_sig_coloc =====\n")

ad_coloc_path <- file.path(out_dir, "AD_2022_coloc_genes_pp_h4_0.9.tsv")
scz_coloc_path <- file.path(out_dir, "SCZ_eur_coloc_genes_pp_h4_0.9.tsv")

if (dir.exists(ad_coloc_dir)) {
    .run_step_table(
        step_label = "  AD: get_sig_coloc",
        output_path = ad_coloc_path,
        fun = function() {
            bican.mccarroll.eqtl::get_sig_coloc(
                coloc_dir       = ad_coloc_dir,
                pp_h4_threshold = 0.9,
                output_path     = ad_coloc_path
            )
        }
    )
} else {
    cat("  SKIPPED: AD coloc dir not found at", ad_coloc_dir, "\n")
}

if (dir.exists(scz_coloc_dir)) {
    .run_step_table(
        step_label = "  SCZ: get_sig_coloc",
        output_path = scz_coloc_path,
        fun = function() {
            bican.mccarroll.eqtl::get_sig_coloc(
                coloc_dir       = scz_coloc_dir,
                pp_h4_threshold = 0.9,
                output_path     = scz_coloc_path
            )
        }
    )
} else {
    cat("  SKIPPED: SCZ coloc dir not found at", scz_coloc_dir, "\n")
}

## -----------------------
## Step 11: plot_fisher_exact (AD and SCZ)
## -----------------------
## NOTE: Requires contingency tables from the Python pipeline (build_fisher_contingency_table).
## Run test_eqtl_pipeline.py steps 2-3 first, then come back to this step.

cat("\n===== Step 11: plot_fisher_exact =====\n")

cluster_order <- c(8, 1, 3, 5, 4, 10, 2, 0, 7, 6, 9)

ad_fisher_path <- file.path(out_dir, "AD_2022_fisher_contingency_counts_gene_clusters.tsv")
scz_fisher_path <- file.path(out_dir, "SCZ_eur_fisher_contingency_counts_gene_clusters.tsv")

ad_enrich_plot_path <- file.path(out_dir, "AD_2022_cluster_enrichment.png")
scz_enrich_plot_path <- file.path(out_dir, "SCZ_eur_cluster_enrichment.png")

if (file.exists(ad_fisher_path)) {
    .run_step_plot(
        step_label = "  AD: plot_fisher_exact",
        output_path = ad_enrich_plot_path,
        fun = function() {
            bican.mccarroll.eqtl::plot_fisher_exact(
                fisher_table_path  = ad_fisher_path,
                plot_disease_label = "Alzheimer's disease",
                cluster_order      = cluster_order,
                output_path        = ad_enrich_plot_path
            )
        }
    )
} else {
    cat("  SKIPPED: AD contingency table not found. Run Python pipeline first.\n")
}

if (file.exists(scz_fisher_path)) {
    .run_step_plot(
        step_label = "  SCZ: plot_fisher_exact",
        output_path = scz_enrich_plot_path,
        fun = function() {
            bican.mccarroll.eqtl::plot_fisher_exact(
                fisher_table_path  = scz_fisher_path,
                plot_disease_label = "schizophrenia",
                cluster_order      = cluster_order,
                output_path        = scz_enrich_plot_path
            )
        }
    )
} else {
    cat("  SKIPPED: SCZ contingency table not found. Run Python pipeline first.\n")
}

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
    out_file <- file.path(out_dir, paste0(case$gene, "_", case$chr, "_", case$pos, ".png"))

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
## Done with R steps
## -----------------------

cat("\n===== All R steps completed successfully! =====\n")
cat("\nNext: run the Python pipeline.\n")
cat("Run the following in a terminal:\n\n")
cat("  python3 -m venv ~/test_eqtl_env\n")
cat("  source ~/test_eqtl_env/bin/activate\n")
cat('  pip install "bican_mccarroll_eqtl @ git+https://github.com/broadinstitute/bican-mccarroll-manuscript1.git@ty_eqtl#subdirectory=python/bican_mccarroll_eqtl"\n')
cat("  python3 test_eqtl_pipeline.py\n\n")
