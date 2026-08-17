# options(
#     bican.mccarroll.figures.data_root_dir =
#         "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis",
#
#     bican.mccarroll.figures.out_dir =
#         "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/figure_repository"
# )
#
# plot_de_sign_test_snap200_age_bican()
# plot_de_sign_test_pmid_40903571_bican()
# plot_de_sign_test_pmid_39227716_bican()
# plot_de_sign_test_pmid_39402379_bican()
# plot_de_sign_test_snap200_sex_bican()

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
    primary_label = "BICAN",
    secondary_label = "SNAP",
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
    primary_label = "PMID 40903571",
    secondary_label = "BICAN",
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
plot_de_sign_test_pmid_39227716_bican <- function(outDir = NULL, min_num_genes = 1, alpha = 0.05) {
  .plot_de_sign_test_figure(
    manifest_file = "differential_expression/external_comparison_PMID_39227716/metadata/PMID_39227716_bican_dfc_age_de_overlap_manifest.tsv",
    dataset_id = "pmid_39227716",
    primary_dataset = "bican",
    secondary_dataset = "PMID_39227716",
    primary_label = "BICAN",
    secondary_label = "PMID 39227716",
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
#'   cell type to be included in the sign-concordance summary.
#' @param alpha Adjusted P-value threshold used to select BICAN-significant
#'   genes.
#'
#' @return Invisibly returns a named list (by \code{contrast}) of the lists
#'   produced by
#'   \code{bican.mccarroll.differentialexpression::plot_de_primary_secondary_manifest()}.
#'
#' @export
plot_de_sign_test_pmid_39402379_bican <- function(
  contrast = c("ad_cps", "early_ad_cps", "late_ad_cps", "versus_all"),
  outDir = NULL, min_num_genes = 20, alpha = 0.05
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
      primary_label = "BICAN",
      secondary_label = "PMID 39402379",
      effect_name = sprintf("age vs %s effects", ct),
      min_num_genes = min_num_genes,
      alpha = alpha,
      outDir = outDir
    )
  })
  names(results) <- contrast

  invisible(results)
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
      primary_label = "BICAN",
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

  result <- bican.mccarroll.differentialexpression::plot_de_primary_secondary_manifest(
    manifest_file = paths$manifest_file,
    out_dir = dataset_out_dir,
    primary_dataset = primary_dataset,
    secondary_dataset = secondary_dataset,
    primary_label = primary_label,
    secondary_label = secondary_label,
    effect_name = effect_name,
    gene_set = gene_set,
    contig_yaml_file = paths$contig_yaml_file,
    reduced_gtf_file = paths$reduced_gtf_file,
    min_num_genes = min_num_genes,
    alpha = alpha
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
