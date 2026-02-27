
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

## Functions

### Data processing
- `get_slope_matrix` — Build matrix of eQTL effect sizes (slopes) across cell types
- `get_pval_nominal_matrix` — Build matrix of nominal p-values across cell types
- `get_pval_nominal_threshold_matrix` — Build matrix of nominal p-value thresholds
- `get_egene_union_pairs` — Get union of significant gene-variant pairs across cell types
- `get_index_snp_slope_matrix_with_median_impute` — Select index SNP per gene and impute missing slopes
- `get_heatmap_index_snp_median_expression` — Compute median expression per cell type for heatmap genes
- `combine_expression_across_cell_types` — Combine per-sample gene expression TPM into a single matrix
- `get_sig_coloc` — Extract significant colocalization genes from coloc results

### Visualization
- `plot_gene_snp` — Violin + beeswarm plot of gene expression by SNP genotype across cell types
- `plot_cell_type_pairwise_cor` — Heatmap of pairwise R-squared of eQTL effect sizes
- `plot_fisher_exact` — Forest plot of Fisher's exact test enrichment by K-means cluster
- `plot_effect_sizes` — Scatter plots of eQTL effect sizes between cell types

## Companion Python package

K-means clustering and expression heatmap visualization are in the Python package `bican_mccarroll_eqtl` (see `python/bican_mccarroll_eqtl/`).
