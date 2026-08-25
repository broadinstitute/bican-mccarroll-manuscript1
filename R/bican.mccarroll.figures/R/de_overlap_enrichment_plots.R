# source("R/paths.R")
#
# options(
#     bican.mccarroll.figures.data_root_dir =
#         "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis",
#
#     bican.mccarroll.figures.out_dir =
#         "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/figure_repository",
#
#     bican.mccarroll.figures.cache_dir =
#         "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/figure_repository/data_cache"
# )
#
# plot_de_overlap_enrichment_bican_sea_ad_vs_pmid_39402379()
# plot_de_overlap_enrichment_bican_vs_kana_2026()

#' Plot BICAN age vs. PMID_39402379 AD directional DE overlap enrichment
#'
#' Wires the BICAN LEVEL_6_sea_ad_mtg_mmc age-effect differential expression
#' results and the external PMID_39402379 (Gabitto et al. 2024) \code{ad_cps}
#' contrast into
#' \code{bican.mccarroll.de.analysis::build_de_overlap_gene_table()} and
#' \code{bican.mccarroll.de.analysis::compute_de_overlap_enrichment()},
#' restricted to the SEA-AD supertype cell types for which BICAN DE results
#' were generated (same three input paths as
#' \code{\link{plot_de_cor_heatmap_bican_sea_ad_vs_pmid_39402379}}).
#'
#' For each matched cell type, tests whether age-upregulated genes are
#' enriched among AD-upregulated genes, and likewise for downregulated
#' genes, via one-sided Fisher exact tests on the intersected tested-gene
#' universe. Also reports the Pearson/Spearman correlation of the continuous
#' effect sizes, and genes significant in both datasets with opposite
#' effect directions.
#'
#' The joined per-gene table (threshold-independent) is cached to a gzipped
#' TSV so that changing \code{min_de_genes}, \code{fdr_cutoff_age}, or
#' \code{fdr_cutoff_ad} does not require re-reading the raw DE result files.
#' The derived result tables (\code{fisher}, \code{counts},
#' \code{discordant}) are small and are always recomputed and rewritten,
#' since they depend on the thresholds.
#'
#' @param min_de_genes Minimum number of significant genes (summed across up
#'   and down) required in \emph{both} datasets for a cell type to be
#'   included in the Fisher results and the heatmap. Defaults to 10; raise
#'   this for a more stringent set of tested cell types.
#' @param fdr_cutoff_age Adjusted P value threshold used to call a gene
#'   significant in the BICAN age dataset. Defaults to 0.05.
#' @param fdr_cutoff_ad Adjusted P value threshold used to call a gene
#'   significant in the PMID_39402379 \code{ad_cps} dataset. Defaults to
#'   0.05.
#' @param force_recompute Logical scalar. If \code{TRUE}, ignore the cached
#'   joined per-gene table, re-read the raw DE result files, and overwrite
#'   the cache. Defaults to \code{FALSE}.
#'
#' @return Invisibly returns the list produced by
#'   \code{bican.mccarroll.de.analysis::compute_de_overlap_enrichment()}
#'   (elements \code{fisher}, \code{counts}, \code{discordant}).
#'
#' @export
plot_de_overlap_enrichment_bican_sea_ad_vs_pmid_39402379 <- function(
  min_de_genes = 10,
  fdr_cutoff_age = 0.05,
  fdr_cutoff_ad = 0.05,
  force_recompute = FALSE
) {
  contrast <- "ad_cps"
  ct_file_rel <- "differential_expression/metadata/cell_types_for_sea_ad_mtg_mmc_plots.txt"

  bican_paths <- resolve_de_cor_paths(
    de_region_interaction_dir = "differential_expression/results/LEVEL_6_sea_ad_mtg_mmc/sex_age/cell_type",
    ct_file = ct_file_rel
  )

  # resolve_de_cor_paths() is a generic in_dir/outDir/cache resolver despite
  # its "region_interaction" parameter name; reused here (with a PMID in_dir
  # override) rather than writing a near-duplicate resolver, matching
  # plot_de_cor_heatmap_bican_sea_ad_vs_pmid_39402379().
  pmid_paths <- resolve_de_cor_paths(
    de_region_interaction_dir = "differential_expression/external_comparison_PMID_39402379/voom-like",
    ct_file = ct_file_rel
  )

  cell_types_use <- bican.mccarroll.de.analysis::read_cell_types(bican_paths$ct_file)

  out_prefix <- sprintf(
    "de_overlap_enrichment_bican_sea_ad_vs_pmid_39402379_%s",
    contrast
  )

  genes_cache_file <- file.path(
    bican_paths$data_cache_dir,
    sprintf("%s_genes.tsv.gz", out_prefix)
  )

  if (!force_recompute && file.exists(genes_cache_file)) {
    gene_dt <- data.table::fread(genes_cache_file)
  } else {
    gene_to_chr <- bican.mccarroll.de.analysis::read_gene_to_chr(bican_paths$gene_to_chr_file)

    de_dt <- bican.mccarroll.de.analysis::read_de_results_datasets(
      de_dir1 = bican_paths$de_region_interaction_dir,
      test1 = "__age_DE_results\\.txt$",
      dataset1 = "BICAN",
      de_dir2 = pmid_paths$de_region_interaction_dir,
      test2 = sprintf("__MTG__%s_DE_results\\.txt$", contrast),
      dataset2 = "PMID_39402379",
      ct_file = bican_paths$ct_file,
      gene_to_chr = gene_to_chr
    )

    gene_dt <- bican.mccarroll.de.analysis::build_de_overlap_gene_table(
      de_dt,
      cell_types_use = cell_types_use,
      dataset_age = "BICAN",
      dataset_ad = "PMID_39402379"
    )

    data.table::fwrite(gene_dt, file = genes_cache_file, sep = "\t")
  }

  result <- bican.mccarroll.de.analysis::compute_de_overlap_enrichment(
    gene_dt,
    fdr_cutoff_age = fdr_cutoff_age,
    fdr_cutoff_ad = fdr_cutoff_ad,
    min_de_genes = min_de_genes
  )

  utils::write.table(
    result$fisher,
    file = file.path(bican_paths$data_cache_dir, sprintf("%s_fisher.tsv", out_prefix)),
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
  )
  utils::write.table(
    result$counts,
    file = file.path(bican_paths$data_cache_dir, sprintf("%s_counts.tsv", out_prefix)),
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
  )
  utils::write.table(
    result$discordant,
    file = file.path(bican_paths$data_cache_dir, sprintf("%s_discordant_genes.tsv", out_prefix)),
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
  )

  p <- bican.mccarroll.de.analysis::plot_de_overlap_enrichment_heatmap(
    result$fisher,
    cell_type_order = cell_types_use
  )

  # The reference figure this is modeled on has only 6-7 rows; we have up
  # to 33 SEA-AD supertypes, so plot height must scale with row count or
  # every row gets crushed into an illegibly short tile.
  plot_height <- max(6, 0.3 * length(unique(result$fisher$cell_type)) + 2)

  save_plot_svg(
    p,
    out_file = sprintf("%s.svg", out_prefix),
    out_dir = bican_paths$outDir,
    width = 11,
    height = plot_height
  )

  logger::log_info("DONE plotting DE overlap enrichment heatmap: {out_prefix}.svg")

  invisible(result)
}


#' Plot BICAN age vs. Kana 2026 CaH CPS directional DE overlap enrichment
#'
#' Wires the BICAN CaH-region age-effect differential expression results and
#' the external Kana et al. 2026 (bioRxiv) CaH \code{ad_cps} contrast into
#' \code{bican.mccarroll.de.analysis::build_de_overlap_gene_table()} and
#' \code{bican.mccarroll.de.analysis::compute_de_overlap_enrichment()},
#' restricted to the six striatal cell types for which a Kana-to-BICAN cell
#' type name mapping was defined (see \code{kana_2026_comparison_setup.R} in
#' \code{adhoc_scripts/external_data_transforms/}; same inputs and cell type
#' renaming as
#' \code{\link{plot_de_cor_heatmap_bican_vs_kana_2026}}).
#'
#' For each matched cell type, tests whether age-upregulated genes are
#' enriched among AD-upregulated genes, and likewise for downregulated
#' genes, via one-sided Fisher exact tests on the intersected tested-gene
#' universe. Also reports the Pearson/Spearman correlation of the continuous
#' effect sizes, and genes significant in both datasets with opposite
#' effect directions.
#'
#' The joined per-gene table (threshold-independent) is cached to a gzipped
#' TSV so that changing \code{min_de_genes}, \code{fdr_cutoff_age}, or
#' \code{fdr_cutoff_ad} does not require re-reading the raw DE result files.
#' The derived result tables (\code{fisher}, \code{counts},
#' \code{discordant}) are small and are always recomputed and rewritten,
#' since they depend on the thresholds.
#'
#' @param min_de_genes Minimum number of significant genes (summed across up
#'   and down) required in \emph{both} datasets for a cell type to be
#'   included in the Fisher results and the heatmap. Defaults to 10; both
#'   datasets have thousands of significant genes for all six cell types, so
#'   no cell type is expected to be excluded at this default.
#' @param fdr_cutoff_age Adjusted P value threshold used to call a gene
#'   significant in the BICAN age dataset. Defaults to 0.05.
#' @param fdr_cutoff_ad Adjusted P value threshold used to call a gene
#'   significant in the Kana 2026 \code{ad_cps} dataset. Defaults to 0.05.
#' @param force_recompute Logical scalar. If \code{TRUE}, ignore the cached
#'   joined per-gene table, re-read the raw DE result files, and overwrite
#'   the cache. Defaults to \code{FALSE}.
#'
#' @return Invisibly returns the list produced by
#'   \code{bican.mccarroll.de.analysis::compute_de_overlap_enrichment()}
#'   (elements \code{fisher}, \code{counts}, \code{discordant}).
#'
#' @export
plot_de_overlap_enrichment_bican_vs_kana_2026 <- function(
  min_de_genes = 10,
  fdr_cutoff_age = 0.05,
  fdr_cutoff_ad = 0.05,
  force_recompute = FALSE
) {
  contrast <- "ad_cps"
  ct_file_rel <- "differential_expression/metadata/cell_types_for_kana_2026_bican_plots.txt"

  bican_paths <- resolve_de_cor_paths(
    de_region_interaction_dir = "differential_expression/results/LEVEL_6/sex_age/cell_type_region_interaction_absolute_effects",
    ct_file = ct_file_rel
  )

  # resolve_de_cor_paths() is a generic in_dir/outDir/cache resolver despite
  # its "region_interaction" parameter name; reused here (with a Kana in_dir
  # override) rather than writing a near-duplicate resolver, matching
  # plot_de_cor_heatmap_bican_vs_kana_2026().
  kana_paths <- resolve_de_cor_paths(
    de_region_interaction_dir = "differential_expression/external_comparison_kana_2026_biorxiv/voom-like",
    ct_file = ct_file_rel
  )

  cell_types_use <- bican.mccarroll.de.analysis::read_cell_types(bican_paths$ct_file)

  cell_type_map <- utils::read.table(
    file.path(bican_paths$data_root_dir, "differential_expression/external_comparison_kana_2026_biorxiv/metadata/kana_2026_bican_cell_type_map.tsv"),
    header = TRUE, sep = "\t", stringsAsFactors = FALSE
  )[, c("kana", "bican")]

  out_prefix <- sprintf(
    "de_overlap_enrichment_bican_vs_kana_2026_%s",
    contrast
  )

  genes_cache_file <- file.path(
    bican_paths$data_cache_dir,
    sprintf("%s_genes.tsv.gz", out_prefix)
  )

  if (!force_recompute && file.exists(genes_cache_file)) {
    gene_dt <- data.table::fread(genes_cache_file)
  } else {
    gene_to_chr <- bican.mccarroll.de.analysis::read_gene_to_chr(bican_paths$gene_to_chr_file)

    de_dt <- .read_de_results_datasets_mapped(
      de_dir1 = bican_paths$de_region_interaction_dir,
      test1 = "__CaH__age_DE_results\\.txt$",
      dataset1 = "BICAN",
      de_dir2 = kana_paths$de_region_interaction_dir,
      test2 = sprintf("__CaH__%s_DE_results\\.txt$", contrast),
      dataset2 = "KANA_2026",
      ct_file = bican_paths$ct_file,
      gene_to_chr = gene_to_chr,
      cell_type_map = cell_type_map
    )

    gene_dt <- bican.mccarroll.de.analysis::build_de_overlap_gene_table(
      de_dt,
      cell_types_use = cell_types_use,
      dataset_age = "BICAN",
      dataset_ad = "KANA_2026"
    )

    data.table::fwrite(gene_dt, file = genes_cache_file, sep = "\t")
  }

  result <- bican.mccarroll.de.analysis::compute_de_overlap_enrichment(
    gene_dt,
    fdr_cutoff_age = fdr_cutoff_age,
    fdr_cutoff_ad = fdr_cutoff_ad,
    min_de_genes = min_de_genes
  )

  utils::write.table(
    result$fisher,
    file = file.path(bican_paths$data_cache_dir, sprintf("%s_fisher.tsv", out_prefix)),
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
  )
  utils::write.table(
    result$counts,
    file = file.path(bican_paths$data_cache_dir, sprintf("%s_counts.tsv", out_prefix)),
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
  )
  utils::write.table(
    result$discordant,
    file = file.path(bican_paths$data_cache_dir, sprintf("%s_discordant_genes.tsv", out_prefix)),
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
  )

  p <- bican.mccarroll.de.analysis::plot_de_overlap_enrichment_heatmap(
    result$fisher,
    cell_type_order = cell_types_use
  )

  plot_height <- max(6, 0.3 * length(unique(result$fisher$cell_type)) + 2)

  save_plot_svg(
    p,
    out_file = sprintf("%s.svg", out_prefix),
    out_dir = bican_paths$outDir,
    width = 11,
    height = plot_height
  )

  logger::log_info("DONE plotting DE overlap enrichment heatmap: {out_prefix}.svg")

  invisible(result)
}
