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


# plot_de_sign_test_snap200_age_bican()
# plot_de_sign_test_snap200_sex_bican()
# plot_de_sign_test_pmid_40903571_bican()
# plot_de_sign_test_pmid_39227716_bican()
# plot_de_sign_test_pmid_39402379_bican()
# plot_de_sign_test_kana_2026_bican()


#' Plot BICAN vs. SNAP age-effect DE sign-test comparison
#'
#' Wires the BICAN/SNAP age-effect overlap manifest into
#' \code{bican.mccarroll.differentialexpression::plot_de_primary_secondary_manifest()}.
#' Writes per-cell-type scatter SVGs, a sign-concordance summary SVG, and a
#' summary TSV to a \code{de_sign_test_snap200} subdirectory of the configured
#' figure output directory (see \code{\link{get_out_dir}}).
#'
#' @param outDir Output directory for the generated SVGs/TSV. If \code{NULL},
#'   resolved via configured output directory options.
#' @param min_num_genes Minimum number of overlapping genes required for a
#'   cell type to be included in the sign-concordance summary.
#' @param alpha Adjusted P-value threshold used to select BICAN-significant
#'   genes.
#'
#' @return Invisibly returns the list produced by
#'   \code{bican.mccarroll.differentialexpression::plot_de_primary_secondary_manifest()}.
#'
#' @export
plot_de_sign_test_snap200_age_bican <- function(outDir = NULL, min_num_genes = 20, alpha = 0.05) {
  .plot_de_sign_test_figure(
    manifest_file = "differential_expression/external_comparison_snap200/metadata/bican_dfc_snap_age_de_overlap_manifest.tsv",
    dataset_id = "snap200",
    primary_dataset = "bican",
    secondary_dataset = "snap",
    primary_label = "Burger et al.",
    secondary_label = "Ling et al.",
    effect_name = "age effects",
    min_num_genes = min_num_genes,
    alpha = alpha,
    outDir = outDir
  )
}


#' Plot PMID_40903571 vs. BICAN age-effect DE sign-test comparison
#'
#' Wires the PMID_40903571/BICAN age-effect overlap manifest into
#' \code{bican.mccarroll.differentialexpression::plot_de_primary_secondary_manifest()}.
#' Writes per-cell-type scatter SVGs, a sign-concordance summary SVG, and a
#' summary TSV to a \code{de_sign_test_pmid_40903571} subdirectory of the
#' configured figure output directory (see \code{\link{get_out_dir}}).
#'
#' @param outDir Output directory for the generated SVGs/TSV. If \code{NULL},
#'   resolved via configured output directory options.
#' @param min_num_genes Minimum number of overlapping genes required for a
#'   cell type to be included in the sign-concordance summary.
#' @param alpha Adjusted P-value threshold used to select PMID_40903571-significant
#'   genes.
#'
#' @return Invisibly returns the list produced by
#'   \code{bican.mccarroll.differentialexpression::plot_de_primary_secondary_manifest()}.
#'
#' @export
plot_de_sign_test_pmid_40903571_bican <- function(outDir = NULL, min_num_genes = 20, alpha = 0.05) {
  .plot_de_sign_test_figure(
    manifest_file = "differential_expression/external_comparison_PMID_40903571/metadata/PMID_40903571_bican_dfc_age_de_overlap_manifest.tsv",
    dataset_id = "pmid_40903571",
    primary_dataset = "PMID_40903571",
    secondary_dataset = "bican",
    primary_label = "Jeffries et al.",
    secondary_label = "Burger et al.",
    effect_name = "age effects",
    min_num_genes = min_num_genes,
    alpha = alpha,
    outDir = outDir
  )
}


#' Plot BICAN vs. PMID_39227716 age-effect DE sign-test comparison
#'
#' Wires the BICAN/PMID_39227716 age-effect overlap manifest into
#' \code{bican.mccarroll.differentialexpression::plot_de_primary_secondary_manifest()}.
#' Writes per-cell-type scatter SVGs, a sign-concordance summary SVG, and a
#' summary TSV to a \code{de_sign_test_pmid_39227716} subdirectory of the
#' configured figure output directory (see \code{\link{get_out_dir}}).
#'
#' \code{min_num_genes} defaults to \code{1} to retain the one cell type that
#' would otherwise be filtered out of the summary.
#'
#' @param outDir Output directory for the generated SVGs/TSV. If \code{NULL},
#'   resolved via configured output directory options.
#' @param min_num_genes Minimum number of overlapping genes required for a
#'   cell type to be included in the sign-concordance summary.
#' @param alpha Adjusted P-value threshold used to select BICAN-significant
#'   genes.
#'
#' @return Invisibly returns the list produced by
#'   \code{bican.mccarroll.differentialexpression::plot_de_primary_secondary_manifest()}.
#'
#' @export
plot_de_sign_test_pmid_39227716_bican <- function(outDir = NULL, min_num_genes = 20, alpha = 0.05) {
  .plot_de_sign_test_figure(
    manifest_file = "differential_expression/external_comparison_PMID_39227716/metadata/PMID_39227716_bican_dfc_age_de_overlap_manifest.tsv",
    dataset_id = "pmid_39227716",
    primary_dataset = "bican",
    secondary_dataset = "PMID_39227716",
    primary_label = "Burger et al.",
    secondary_label = "Fröhlich et al.",
    effect_name = "age effects",
    min_num_genes = min_num_genes,
    alpha = alpha,
    outDir = outDir
  )
}


#' Plot BICAN vs. PMID_39402379 age-effect DE sign-test comparison
#'
#' Wires the BICAN/PMID_39402379 (Gabitto et al. 2024) age-effect overlap
#' manifests into
#' \code{bican.mccarroll.differentialexpression::plot_de_primary_secondary_manifest()}.
#' PMID_39402379 provides four contrasts (\code{ad_cps}, \code{early_ad_cps},
#' \code{late_ad_cps}, \code{versus_all}); each is compared separately
#' against BICAN's \code{age} contrast. Writes per-cell-type scatter SVGs, a
#' sign-concordance summary SVG, and a summary TSV to a
#' \code{<contrast>} subdirectory of the
#' \code{de_sign_test_pmid_39402379} subdirectory of the configured figure
#' output directory (see \code{\link{get_out_dir}}).
#'
#' @param contrast Character vector of one or more of \code{"ad_cps"},
#'   \code{"early_ad_cps"}, \code{"late_ad_cps"}, \code{"versus_all"}.
#' @param outDir Output directory for the generated SVGs/TSVs. If \code{NULL},
#'   resolved via configured output directory options.
#' @param min_num_genes Minimum number of overlapping genes required for a
#'   cell type to be included in the sign-concordance summary. Defaults to
#'   \code{1} (rather than the usual \code{20}) since the sea_ad_mtg_mmc
#'   cell-type list is restricted to the 33 BICAN supertypes, several of
#'   which have few overlapping significant genes.
#' @param alpha Adjusted P-value threshold used to select BICAN-significant
#'   genes.
#' @param data_cache_dir Directory used to store the cached per-gene
#'   comparison data (one TSV per contrast). If \code{NULL}, the directory is
#'   resolved from \code{options("bican.mccarroll.figures.cache_dir")}.
#' @param force_recompute Logical scalar. If \code{TRUE}, ignore existing
#'   cache files, recompute the comparison data, and overwrite the cache.
#'   Defaults to \code{FALSE}.
#' @param scale_effects Logical scalar. If \code{TRUE}, each cell type's
#'   BICAN and PMID_39402379 logFC values are independently rescaled to
#'   \code{[-1, 1]} before plotting, making effect sizes visually comparable
#'   across cell types despite the two datasets' very different native logFC
#'   scales (BICAN age effects vs. AD pseudo-progression score effects).
#'   Does not affect the underlying sign-concordance statistics or the cached
#'   per-gene data, which stay on the raw scale. Defaults to \code{TRUE} for
#'   this comparison.
#'
#' @return Invisibly returns a named list (by \code{contrast}) of the lists
#'   produced by
#'   \code{bican.mccarroll.differentialexpression::plot_de_primary_secondary_from_data()}.
#'
#' @export
plot_de_sign_test_pmid_39402379_bican <- function(
  contrast = c("ad_cps", "early_ad_cps", "late_ad_cps", "versus_all"),
  outDir = NULL, min_num_genes = 20, alpha = 0.05,
  data_cache_dir = NULL, force_recompute = FALSE,
  scale_effects = TRUE
) {
  results <- lapply(contrast, function(ct) {
    .plot_de_sign_test_figure(
      manifest_file = sprintf(
        "differential_expression/external_comparison_PMID_39402379/metadata/PMID_39402379_bican_%s_de_overlap_manifest.tsv",
        ct
      ),
      dataset_id = "pmid_39402379",
      sub_dir = ct,
      primary_dataset = "bican",
      secondary_dataset = "PMID_39402379",
      primary_label = "Burger et al.",
      secondary_label = "Gabitto et al.",
      effect_name = sprintf("age vs %s effects", ct),
      min_num_genes = min_num_genes,
      alpha = alpha,
      cache_name = sprintf("de_sign_test_pmid_39402379_%s.tsv", ct),
      data_cache_dir = data_cache_dir,
      force_recompute = force_recompute,
      scale_effects = scale_effects,
      outDir = outDir
    )
  })
  names(results) <- contrast

  invisible(results)
}


#' Plot BICAN vs. Kana 2026 CaH CPS-effect DE sign-test comparison
#'
#' Wires the BICAN/Kana 2026 (bioRxiv) CaH age-vs-CPS overlap manifest into
#' \code{bican.mccarroll.differentialexpression::plot_de_primary_secondary_manifest()}.
#' The manifest is restricted to the six striatal cell types for which a
#' Kana-to-BICAN cell type name mapping was defined (see
#' \code{kana_2026_comparison_setup.R} in
#' \code{adhoc_scripts/external_data_transforms/}); no renaming is needed
#' here since the manifest already pairs each BICAN file with its
#' differently-named Kana counterpart under one \code{cell_type} label.
#' Writes per-cell-type scatter SVGs, a sign-concordance summary SVG, and a
#' summary TSV to a \code{de_sign_test_kana_2026} subdirectory of the
#' configured figure output directory (see \code{\link{get_out_dir}}).
#'
#' @param outDir Output directory for the generated SVGs/TSV. If \code{NULL},
#'   resolved via configured output directory options.
#' @param min_num_genes Minimum number of overlapping genes required for a
#'   cell type to be included in the sign-concordance summary.
#' @param alpha Adjusted P-value threshold used to select BICAN-significant
#'   genes.
#' @param scale_effects Logical scalar. If \code{TRUE}, each cell type's
#'   BICAN and Kana 2026 logFC values are independently rescaled to
#'   \code{[-1, 1]} before plotting, matching the treatment of the
#'   PMID_39402379 comparison since both compare a BICAN age effect against
#'   an AD pseudo-progression score effect on an unrelated native scale.
#'   Defaults to \code{TRUE}.
#'
#' @return Invisibly returns the list produced by
#'   \code{bican.mccarroll.differentialexpression::plot_de_primary_secondary_manifest()}.
#'
#' @export
plot_de_sign_test_kana_2026_bican <- function(outDir = NULL, min_num_genes = 20, alpha = 0.05, scale_effects = TRUE) {
  .plot_de_sign_test_figure(
    manifest_file = "differential_expression/external_comparison_kana_2026_biorxiv/metadata/kana_2026_bican_ad_cps_de_overlap_manifest.tsv",
    dataset_id = "kana_2026",
    primary_dataset = "bican",
    secondary_dataset = "KANA_2026",
    primary_label = "Burger et al.",
    secondary_label = "Kana et al. 2026",
    effect_name = "age vs ad_cps effects",
    min_num_genes = min_num_genes,
    alpha = alpha,
    scale_effects = scale_effects,
    outDir = outDir
  )
}


#' Plot BICAN vs. SNAP sex-effect DE sign-test comparison
#'
#' Wires the BICAN/SNAP sex-effect overlap manifest into
#' \code{bican.mccarroll.differentialexpression::plot_de_primary_secondary_manifest()},
#' restricted to autosomal genes, sex-chromosome genes, or both. Writes
#' per-cell-type scatter SVGs, a sign-concordance summary SVG, and a summary
#' TSV to the \code{de_sign_test_snap200} subdirectory of the configured
#' figure output directory (see \code{\link{get_out_dir}}) — the same
#' subdirectory used by \code{\link{plot_de_sign_test_snap200_age_bican}},
#' since both compare against the same external dataset.
#'
#' By default all three gene sets (\code{"both"}, \code{"autosome"},
#' \code{"xy"}) are plotted. Pass a single value to generate just one.
#'
#' @param gene_set Character vector of one or more of \code{"both"},
#'   \code{"autosome"}, \code{"xy"}.
#' @param outDir Output directory for the generated SVGs/TSV. If \code{NULL},
#'   resolved via configured output directory options.
#' @param alpha Adjusted P-value threshold used to select BICAN-significant
#'   genes.
#'
#' @return Invisibly returns a named list (by \code{gene_set}) of the lists
#'   produced by
#'   \code{bican.mccarroll.differentialexpression::plot_de_primary_secondary_manifest()}.
#'
#' @export
plot_de_sign_test_snap200_sex_bican <- function(gene_set = c("both", "autosome", "xy"), outDir = NULL, alpha = 0.05) {
  results <- lapply(gene_set, function(gs) {
    .plot_de_sign_test_figure(
      manifest_file = "differential_expression/external_comparison_snap200/metadata/bican_dfc_snap_sex_de_overlap_manifest.tsv",
      dataset_id = "snap200",
      primary_dataset = "bican",
      secondary_dataset = "snap",
      primary_label = "Burger et al.",
      secondary_label = "Ling et al.",
      effect_name = "sex effects",
      gene_set = gs,
      contig_yaml_file = "metadata/GRCh38_ensembl_v43.contig_groups.yaml",
      reduced_gtf_file = "metadata/GRCh38_ensembl_v43.reduced.gtf.gz",
      alpha = alpha,
      outDir = outDir
    )
  })
  names(results) <- gene_set

  invisible(results)
}


.plot_de_sign_test_figure <- function(manifest_file,
                                      dataset_id,
                                      primary_dataset,
                                      secondary_dataset,
                                      primary_label = NULL,
                                      secondary_label = NULL,
                                      effect_name = "age effects",
                                      gene_set = "both",
                                      contig_yaml_file = NULL,
                                      reduced_gtf_file = NULL,
                                      min_num_genes = 20,
                                      alpha = 0.05,
                                      sub_dir = NULL,
                                      cache_name = NULL,
                                      data_cache_dir = NULL,
                                      force_recompute = FALSE,
                                      scale_effects = FALSE,
                                      show_full_concordance_line = TRUE,
                                      order_alphabetically = TRUE,
                                      outDir = NULL) {
  paths <- resolve_de_sign_test_paths(
    manifest_file = manifest_file,
    contig_yaml_file = contig_yaml_file,
    reduced_gtf_file = reduced_gtf_file,
    outDir = outDir
  )

  dataset_out_dir <- file.path(
    paths$outDir,
    paste0("de_sign_test_", dataset_id)
  )
  if (!is.null(sub_dir)) {
    dataset_out_dir <- file.path(dataset_out_dir, sub_dir)
  }
  .ensure_dir(dataset_out_dir)

  if (is.null(cache_name)) {
    de_data <- bican.mccarroll.differentialexpression::compute_de_primary_secondary_data(
      manifest_file = paths$manifest_file,
      primary_dataset = primary_dataset,
      secondary_dataset = secondary_dataset,
      gene_set = gene_set,
      contig_yaml_file = paths$contig_yaml_file,
      reduced_gtf_file = paths$reduced_gtf_file,
      alpha = alpha
    )
  } else {
    cache_dir <- .resolve_cache_dir(data_cache_dir)
    if (is.null(data_cache_dir)) {
      cache_dir <- file.path(cache_dir, "differential_expression")
    }
    .ensure_dir(cache_dir)
    cache_file <- file.path(cache_dir, cache_name)

    if (!force_recompute && file.exists(cache_file)) {
      de_data <- data.table::fread(cache_file)
    } else {
      de_data <- bican.mccarroll.differentialexpression::compute_de_primary_secondary_data(
        manifest_file = paths$manifest_file,
        primary_dataset = primary_dataset,
        secondary_dataset = secondary_dataset,
        gene_set = gene_set,
        contig_yaml_file = paths$contig_yaml_file,
        reduced_gtf_file = paths$reduced_gtf_file,
        alpha = alpha
      )

      utils::write.table(
        de_data,
        file = cache_file,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE
      )
    }
  }

  result <- bican.mccarroll.differentialexpression::plot_de_primary_secondary_from_data(
    de_data,
    out_dir = dataset_out_dir,
    primary_dataset = primary_dataset,
    secondary_dataset = secondary_dataset,
    primary_label = primary_label,
    secondary_label = secondary_label,
    effect_name = effect_name,
    gene_set = gene_set,
    min_num_genes = min_num_genes,
    scale_effects = scale_effects,
    show_full_concordance_line = show_full_concordance_line,
    order_alphabetically = order_alphabetically
  )

  logger::log_info(
    "DONE plotting DE sign test comparison for {dataset_id} ",
    "({effect_name}, gene_set={gene_set})"
  )

  invisible(result)
}


resolve_de_sign_test_paths <- function(manifest_file,
                                       contig_yaml_file = NULL,
                                       reduced_gtf_file = NULL,
                                       outDir = NULL) {
  root <- .resolve_data_root_dir(NULL)

  pick_in <- function(x) {
    if (is.null(x)) {
      return(NULL)
    }
    .resolve_under_root(root, x)
  }

  out <- .resolve_out_dir(outDir)
  .ensure_dir(out)

  list(
    data_root_dir = root,
    manifest_file = .resolve_under_root(root, manifest_file),
    contig_yaml_file = pick_in(contig_yaml_file),
    reduced_gtf_file = pick_in(reduced_gtf_file),
    outDir = out
  )
}
