
# bican.mccarroll.eqtl

<!-- badges: start -->
<!-- badges: end -->

R package for eQTL analysis in the first McCarroll lab BICAN manuscript. Provides functions for building eQTL effect-size matrices, K-means clustering of eQTL patterns, colocalization enrichment testing, and gene-by-genotype visualization.

## Installation

``` r
# install.packages("remotes")
remotes::install_github(
  "broadinstitute/bican-mccarroll-manuscript1",
  subdir = "R/bican.mccarroll.eqtl",
  dependencies = TRUE
)
```

## Running the full pipeline

Test scripts at the repository root run all functions automatically:

1. **`test_eqtl_pipeline.R`** — Runs all 12 R steps (data extraction, matrix building, visualization)
2. **`test_eqtl_pipeline.py`** — Runs all 5 Python steps (K-means clustering, Fisher contingency tables, expression heatmap)

Each script runs all its steps sequentially. Steps that depend on outputs from the other pipeline will be skipped automatically if the required files don't exist yet.

The recommended run order is R first, then Python, then R again for the Fisher plots:

``` bash
# On the Broad server:

# 1. Install packages
Rscript -e 'remotes::install_github("broadinstitute/bican-mccarroll-manuscript1", subdir="R/bican.mccarroll.eqtl", ref="ty_eqtl", dependencies=TRUE)'
python3 -m venv ~/test_eqtl_env && source ~/test_eqtl_env/bin/activate
pip install "bican_mccarroll_eqtl @ git+https://github.com/broadinstitute/bican-mccarroll-manuscript1.git@ty_eqtl#subdirectory=python/bican_mccarroll_eqtl"

# 2. Run R pipeline (steps 1-12; Fisher plots skip if Python hasn't run yet)
Rscript test_eqtl_pipeline.R

# 3. Run Python pipeline (K-means, Fisher contingency tables, expression heatmap)
python3 test_eqtl_pipeline.py

# 4. Re-run R pipeline to generate Fisher enrichment plots (now that contingency tables exist)
Rscript test_eqtl_pipeline.R
```

## Functions

### Data processing
- `get_egene_union_pairs` — Get union of significant gene-variant pairs across cell types
- `get_slope_matrix` — Build matrix of eQTL effect sizes (slopes) across cell types
- `get_pval_nominal_matrix` — Build matrix of nominal p-values across cell types
- `get_pval_nominal_threshold_matrix` — Build matrix of nominal p-value thresholds
- `get_index_snp_slope_matrix_with_median_impute` — Select index SNP per gene and impute missing slopes
- `get_cell_type_pairwise_cor_matrix` — Compute pairwise R-squared of eQTL effect sizes
- `get_heatmap_index_snp_median_expression` — Compute median expression per cell type for heatmap genes
- `combine_expression_across_cell_types` — Combine per-sample gene expression TPM into a single matrix
- `get_sig_coloc` — Extract significant colocalization genes from coloc results

### Visualization
- `plot_cell_type_pairwise_cor` — Heatmap of pairwise R-squared of eQTL effect sizes
- `plot_fisher_exact` — Forest plot of Fisher's exact test enrichment by K-means cluster
- `plot_gene_snp` — Violin + beeswarm plot of gene expression by SNP genotype across cell types
- `plot_effect_sizes` — Scatter plots of eQTL effect sizes between cell types

## Companion Python package

K-means clustering and expression heatmap visualization are in the Python package `bican_mccarroll_eqtl` (see `python/bican_mccarroll_eqtl/`). Functions include:

- `run_k_selection` / `plot_k_selection` — Silhouette-based K selection for K-means
- `run_kmeans_heatmap` — K-means clustering and heatmap of eQTL effect sizes
- `build_fisher_contingency_table` — Build contingency tables for Fisher's exact test
- `order_genes_by_correlation` — Order genes within clusters by expression correlation
- `plot_expression_heatmap` — Median expression heatmap grouped by K-means cluster
