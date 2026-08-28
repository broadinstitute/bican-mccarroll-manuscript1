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

# de_region_interaction_dir <- gene_to_chr_file <- ct_file <- outDir <- data_cache_dir <- NULL

#' Plot differential expression correlation heatmaps for age (main and supplement)
#'
#' Generates two manuscript-ready correlation heatmaps of differential
#' expression effect sizes for the age contrast: a main heatmap for a subset
#' of regions and a supplementary heatmap for all regions.
#' Correlation matrices are cached to plain-text TSV files (one per
#' plot) to speed up subsequent runs and to provide reviewer-inspectable
#' intermediate data.
#'
#' Matrix row and column names are cleaned prior to plotting using
#' \code{clean_cor_mat_names()} (double underscores then single underscores
#' replaced with spaces).
#'
#' All file path arguments may be \code{NULL}, in which case they are resolved
#' using \code{resolve_de_cor_paths()} and the configured data root, cache,
#' and output directory options.
#'
#' @param de_region_interaction_dir Directory containing region-interaction
#'   differential expression results. If \code{NULL}, resolved under the data
#'   root.
#' @param gene_to_chr_file Path to the gene-to-chromosome mapping file. If
#'   \code{NULL}, resolved under the data root.
#' @param ct_file Path to the cell type list file. If \code{NULL}, resolved
#'   under the data root.
#' @param outDir Output directory for generated SVG plots. If \code{NULL},
#'   resolved via configured output directory options.
#' @param data_cache_dir Directory used to store cached correlation matrices
#'   as TSV files. If \code{NULL}, resolved via configured cache directory
#'   options.
#'
#' @export
plot_de_cor_heatmaps_age <- function(de_region_interaction_dir = NULL,
                                     gene_to_chr_file = NULL,
                                     ct_file = NULL,
                                     outDir = NULL,
                                     data_cache_dir = NULL) {
  paths <- resolve_de_cor_paths(
    de_region_interaction_dir = de_region_interaction_dir,
    gene_to_chr_file = gene_to_chr_file,
    ct_file = ct_file,
    outDir = outDir,
    data_cache_dir = data_cache_dir
  )

  ## -----------------------
  ## Hard-coded manuscript parameters
  ## -----------------------

  test <- "age"

  non_neuron_types <- c("astrocyte", "OPC", "oligodendrocyte", "microglia")

  region_order <- c("CaH", "Pu", "NAC", "ic", "DFC")
  regions_main <- c("CaH", "DFC")
  regions_supp <- region_order

  fdr_cutoff <- 0.05

  breaks <- seq(-1, 1, length.out = 101)
  palette_colors <- c("steelblue", "white", "darkorange")
  clustering_method <- "complete"

  ## -----------------------
  ## Small input (load once)
  ## -----------------------

  cell_types_use <- bican.mccarroll.de.analysis::read_cell_types(paths$ct_file)

  ## -----------------------
  ## Plot 1: MAIN (cache = matrix only)
  ## -----------------------


  cache_file <- file.path(paths$data_cache_dir, "de_cor_mat_age_main_CaH_DFC.tsv")

  cor_mat_main <- get_or_build_de_cor_mat_cache(
    cache_file = cache_file,
    de_region_interaction_dir = paths$de_region_interaction_dir,
    test = test,
    ct_file = paths$ct_file,
    gene_to_chr_file = paths$gene_to_chr_file,
    regions_use = regions_main,
    non_neuron_types = non_neuron_types,
    fdr_cutoff = fdr_cutoff
  )

  # Clean up the names.
  cor_mat_main <- clean_cor_mat_names(cor_mat_main)

  out_file <- file.path(paths$outDir, "de_cor_heatmap_age_main_CaH_DFC.svg")
  svglite_manuscript(out_file, width = 6, height = 6)

  # legend_title = "Spearman rho^2\nof age DE logFC"
  legend_title <- NULL

  # bican.mccarroll.de.analysis::plot_de_cor_heatmap(
  #     cor_mat_main,
  #     clustering_method = clustering_method,
  #     breaks = breaks,
  #     palette_colors = palette_colors)

  ht <- bican.mccarroll.de.analysis::plot_de_cor_heatmap_complex(
    cor_mat_main,
    clustering_method = clustering_method,
    breaks = breaks,
    palette_colors = palette_colors,
    legend_title = legend_title,
    show_dendrograms = FALSE
  )

  ComplexHeatmap::draw(ht)
  grDevices::dev.off()

  ## -----------------------
  ## Plot 2: SUPP (cache = matrix only)
  ## -----------------------

  cache_file <- file.path(paths$data_cache_dir, "de_cor_mat_age_supp_all_regions.tsv")

  cor_mat_supp <- get_or_build_de_cor_mat_cache(
    cache_file = cache_file,
    de_region_interaction_dir = paths$de_region_interaction_dir,
    test = test,
    ct_file = paths$ct_file,
    gene_to_chr_file = paths$gene_to_chr_file,
    regions_use = regions_supp,
    non_neuron_types = non_neuron_types,
    fdr_cutoff = fdr_cutoff
  )

  # Clean up the names.
  cor_mat_supp <- clean_cor_mat_names(cor_mat_supp)

  out_file <- file.path(paths$outDir, "de_cor_heatmap_age_supp_all_regions.svg")
  svglite_manuscript(out_file, width = 10, height = 10)

  # bican.mccarroll.de.analysis::plot_de_cor_heatmap(
  #     cor_mat_supp,
  #     clustering_method = clustering_method,
  #     breaks = breaks,
  #     palette_colors = palette_colors)

  ht2 <- bican.mccarroll.de.analysis::plot_de_cor_heatmap_complex(
    cor_mat_supp,
    clustering_method = clustering_method,
    breaks = breaks,
    palette_colors = palette_colors,
    legend_title = legend_title,
    show_dendrograms = FALSE
  )

  ComplexHeatmap::draw(ht2)
  grDevices::dev.off()

  invisible(NULL)
}

#' Plot combined BICAN sea_ad_mtg_mmc vs. PMID_39402379 DE correlation heatmap
#'
#' Builds a correlation heatmap whose matrix combines two datasets: BICAN
#' LEVEL_6_sea_ad_mtg_mmc age-effect differential expression results and the
#' external PMID_39402379 (Gabitto et al. 2024) \code{ad_cps} contrast
#' results, both restricted to the SEA-AD supertype cell types for which
#' BICAN DE results were generated. Each cell type appears twice (once per
#' dataset), so hierarchical clustering can freely reveal whether a BICAN
#' cell type's age-effect profile clusters tightly with its PMID counterpart,
#' with some other PMID cell type, or with neither (no forced grouping/split
#' is applied to the clustering). Row/column labels show the contrast that
#' was tested (\code{"age"} for BICAN, \code{"AD"} for PMID_39402379) followed
#' by the cell type name, with a color rug (blue = age, green = AD) and cell
#' type names aligned vertically -- see
#' \code{bican.mccarroll.de.analysis::plot_de_cor_heatmap_complex()}.
#'
#' The correlation matrix is cached to a plain-text TSV file to speed up
#' subsequent runs and to provide reviewer-inspectable intermediate data.
#' A cell type/dataset combination's correlation with a given partner is
#' only computed if their shared/either-significant gene set reaches
#' \code{min_num_genes} (see
#' \code{bican.mccarroll.de.analysis::compute_de_cor_mat_datasets()}); a
#' combination that never reaches that with any partner is dropped from the
#' matrix entirely, so a cache built at one \code{min_num_genes} value must
#' be regenerated (\code{force_recompute = TRUE}) if that value changes.
#'
#' @param min_num_genes Minimum number of genes significant in either member
#'   of a pair required to compute that pair's correlation. Defaults to
#'   \code{20}; raise (e.g. to \code{100}) for a stricter matrix.
#' @param outDir Output directory for generated SVG plots. If \code{NULL},
#'   resolved via configured output directory options.
#' @param data_cache_dir Directory used to store cached correlation matrices
#'   as TSV files. If \code{NULL}, resolved via configured cache directory
#'   options.
#' @param force_recompute Logical scalar. If \code{TRUE}, ignore an existing
#'   cache file, recompute the matrix, and overwrite the cache. Defaults to
#'   \code{FALSE}.
#'
#' @return Invisibly returns the correlation matrix.
#'
#' @export
plot_de_cor_heatmap_bican_sea_ad_vs_pmid_39402379 <- function(
  min_num_genes = 20,
  outDir = NULL,
  data_cache_dir = NULL,
  force_recompute = FALSE
) {
  contrast <- "ad_cps"

  ct_file_rel <- "differential_expression/metadata/cell_types_for_sea_ad_mtg_mmc_plots.txt"

  bican_paths <- resolve_de_cor_paths(
    de_region_interaction_dir = "differential_expression/results/LEVEL_6_sea_ad_mtg_mmc/sex_age/cell_type",
    ct_file = ct_file_rel,
    outDir = outDir,
    data_cache_dir = data_cache_dir
  )

  # resolve_de_cor_paths() is a generic in_dir/outDir/cache resolver despite
  # its "region_interaction" parameter name; reused here (with a PMID in_dir
  # override) rather than writing a near-duplicate resolver.
  pmid_paths <- resolve_de_cor_paths(
    de_region_interaction_dir = "differential_expression/external_comparison_PMID_39402379/voom-like",
    ct_file = ct_file_rel,
    outDir = outDir,
    data_cache_dir = data_cache_dir
  )

  cell_types_use <- bican.mccarroll.de.analysis::read_cell_types(bican_paths$ct_file)

  cache_file <- file.path(
    bican_paths$data_cache_dir,
    sprintf("de_cor_mat_bican_sea_ad_vs_pmid_39402379_%s.tsv", contrast)
  )

  if (!force_recompute && file.exists(cache_file)) {
    cor_dt <- data.table::fread(cache_file)
    rn <- cor_dt[[1]]
    cor_dt[[1]] <- NULL
    cor_mat <- as.matrix(cor_dt)
    rownames(cor_mat) <- rn
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

    cor_mat <- bican.mccarroll.de.analysis::compute_de_cor_mat_datasets(
      de_dt,
      cell_types_use = cell_types_use,
      datasets_use = c("BICAN", "PMID_39402379"),
      min_num_genes = min_num_genes
    )

    cor_dt <- as.data.frame(cor_mat, check.names = FALSE)
    cor_dt <- cbind(cell_type = rownames(cor_mat), cor_dt)

    utils::write.table(
      cor_dt,
      file = cache_file,
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE
    )
  }

  # Row/column labels are split into a contrast label ("age"/"AD") and a
  # cleaned cell type name, rather than cleaned as a single combined string
  # (see .split_cor_mat_labels()); this must run on the raw cell_type__dataset
  # dimnames, before clean_cor_mat_names() collapses them into one string.
  contrast_map <- data.frame(
    dataset = c("BICAN", "PMID_39402379"),
    contrast = c("age", "AD"),
    stringsAsFactors = FALSE
  )
  cor_mat_labels <- .split_cor_mat_labels(cor_mat, contrast_map)

  cor_mat <- clean_cor_mat_names(cor_mat)

  out_file <- file.path(
    bican_paths$outDir,
    sprintf("de_cor_heatmap_bican_sea_ad_vs_pmid_39402379_%s.svg", contrast)
  )
  svglite_manuscript(out_file, width = 12, height = 12)

  ht <- bican.mccarroll.de.analysis::plot_de_cor_heatmap_complex(
    cor_mat,
    clustering_method = "complete",
    breaks = seq(-1, 1, length.out = 101),
    palette_colors = c("steelblue", "white", "darkorange"),
    legend_title = NULL,
    show_dendrograms = TRUE,
    group_labels = cor_mat_labels$group_labels,
    name_labels = cor_mat_labels$name_labels,
    group_colors = .de_cor_contrast_colors
  )

  ComplexHeatmap::draw(ht)
  grDevices::dev.off()

  logger::log_info("DONE plotting DE correlation heatmap: {out_file}")

  invisible(cor_mat)
}


#' Plot combined BICAN CaH vs. Kana 2026 DE correlation heatmap
#'
#' Builds a correlation heatmap whose matrix combines two datasets: BICAN
#' CaH-region age-effect differential expression results and the external
#' Kana et al. 2026 (bioRxiv) CaH \code{ad_cps} contrast results, restricted
#' to the six striatal cell types for which a Kana-to-BICAN cell type name
#' mapping was defined (see \code{kana_2026_comparison_setup.R} in
#' \code{adhoc_scripts/external_data_transforms/}). Each cell type appears
#' twice (once per dataset), so hierarchical clustering can freely reveal
#' whether a BICAN cell type's age-effect profile clusters tightly with its
#' Kana counterpart, with some other Kana cell type, or with neither (no
#' forced grouping/split is applied to the clustering). Row/column labels
#' show the contrast that was tested (\code{"age"} for BICAN, \code{"AD"} for
#' KANA_2026) followed by the cell type name, with a color rug (blue = age,
#' green = AD) and cell type names aligned vertically -- see
#' \code{bican.mccarroll.de.analysis::plot_de_cor_heatmap_complex()}.
#'
#' Unlike \code{plot_de_cor_heatmap_bican_sea_ad_vs_pmid_39402379()}, the two
#' datasets do not share a cell type vocabulary, so
#' \code{.read_de_results_datasets_mapped()} is used in place of
#' \code{bican.mccarroll.de.analysis::read_de_results_datasets()}, renaming
#' Kana's native cell type labels onto their BICAN counterparts (via the
#' two-column map read from \code{kana_2026_bican_cell_type_map.tsv}) before
#' filtering/matrix construction.
#'
#' The correlation matrix is cached to a plain-text TSV file to speed up
#' subsequent runs and to provide reviewer-inspectable intermediate data.
#'
#' @param min_num_genes Minimum number of genes significant in either member
#'   of a pair required to compute that pair's correlation. Defaults to
#'   \code{20}; both datasets have thousands of significant genes for all six
#'   cell types, so no cell type is expected to be dropped at this default.
#' @param outDir Output directory for generated SVG plots. If \code{NULL},
#'   resolved via configured output directory options.
#' @param data_cache_dir Directory used to store cached correlation matrices
#'   as TSV files. If \code{NULL}, resolved via configured cache directory
#'   options.
#' @param force_recompute Logical scalar. If \code{TRUE}, ignore an existing
#'   cache file, recompute the matrix, and overwrite the cache. Defaults to
#'   \code{FALSE}.
#'
#' @return Invisibly returns the correlation matrix.
#'
#' @export
plot_de_cor_heatmap_bican_vs_kana_2026 <- function(
  min_num_genes = 20,
  outDir = NULL,
  data_cache_dir = NULL,
  force_recompute = FALSE
) {
  contrast <- "ad_cps"

  ct_file_rel <- "differential_expression/metadata/cell_types_for_kana_2026_bican_plots.txt"

  bican_paths <- resolve_de_cor_paths(
    de_region_interaction_dir = "differential_expression/results/LEVEL_6/sex_age/cell_type_region_interaction_absolute_effects",
    ct_file = ct_file_rel,
    outDir = outDir,
    data_cache_dir = data_cache_dir
  )

  # resolve_de_cor_paths() is a generic in_dir/outDir/cache resolver despite
  # its "region_interaction" parameter name; reused here (with a Kana in_dir
  # override) rather than writing a near-duplicate resolver.
  kana_paths <- resolve_de_cor_paths(
    de_region_interaction_dir = "differential_expression/external_comparison_kana_2026_biorxiv/voom-like",
    ct_file = ct_file_rel,
    outDir = outDir,
    data_cache_dir = data_cache_dir
  )

  cell_types_use <- bican.mccarroll.de.analysis::read_cell_types(bican_paths$ct_file)

  cell_type_map <- utils::read.table(
    file.path(bican_paths$data_root_dir, "differential_expression/external_comparison_kana_2026_biorxiv/metadata/kana_2026_bican_cell_type_map.tsv"),
    header = TRUE, sep = "\t", stringsAsFactors = FALSE
  )[, c("kana", "bican")]

  cache_file <- file.path(
    bican_paths$data_cache_dir,
    sprintf("de_cor_mat_bican_vs_kana_2026_%s.tsv", contrast)
  )

  if (!force_recompute && file.exists(cache_file)) {
    cor_dt <- data.table::fread(cache_file)
    rn <- cor_dt[[1]]
    cor_dt[[1]] <- NULL
    cor_mat <- as.matrix(cor_dt)
    rownames(cor_mat) <- rn
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

    cor_mat <- bican.mccarroll.de.analysis::compute_de_cor_mat_datasets(
      de_dt,
      cell_types_use = cell_types_use,
      datasets_use = c("BICAN", "KANA_2026"),
      min_num_genes = min_num_genes
    )

    cor_dt <- as.data.frame(cor_mat, check.names = FALSE)
    cor_dt <- cbind(cell_type = rownames(cor_mat), cor_dt)

    utils::write.table(
      cor_dt,
      file = cache_file,
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE
    )
  }

  # Row/column labels are split into a contrast label ("age"/"AD") and a
  # cleaned cell type name, rather than cleaned as a single combined string
  # (see .split_cor_mat_labels()); this must run on the raw cell_type__dataset
  # dimnames, before clean_cor_mat_names() collapses them into one string.
  contrast_map <- data.frame(
    dataset = c("BICAN", "KANA_2026"),
    contrast = c("age", "AD"),
    stringsAsFactors = FALSE
  )
  cor_mat_labels <- .split_cor_mat_labels(cor_mat, contrast_map)

  cor_mat <- clean_cor_mat_names(cor_mat)

  out_file <- file.path(
    bican_paths$outDir,
    sprintf("de_cor_heatmap_bican_vs_kana_2026_%s.svg", contrast)
  )
  svglite_manuscript(out_file, width = 8, height = 8)

  ht <- bican.mccarroll.de.analysis::plot_de_cor_heatmap_complex(
    cor_mat,
    clustering_method = "complete",
    breaks = seq(-1, 1, length.out = 101),
    palette_colors = c("steelblue", "white", "darkorange"),
    legend_title = NULL,
    show_dendrograms = TRUE,
    group_labels = cor_mat_labels$group_labels,
    name_labels = cor_mat_labels$name_labels,
    group_colors = .de_cor_contrast_colors
  )

  ComplexHeatmap::draw(ht)
  grDevices::dev.off()

  logger::log_info("DONE plotting DE correlation heatmap: {out_file}")

  invisible(cor_mat)
}


# Rug colors shared by the cross-dataset DE correlation heatmaps: blue for
# the age contrast, green for the AD contrast. Chosen to be distinguishable
# from the steelblue/white/darkorange correlation fill scale used by those
# same plots.
.de_cor_contrast_colors <- c(age = "#1F4E9C", AD = "#2E8B57")

#' Split cell_type__dataset matrix dimnames into contrast + cell type labels
#'
#' Cross-dataset correlation matrices (see
#' \code{bican.mccarroll.de.analysis::compute_de_cor_mat_datasets()}) have
#' dimnames of the form \code{"<cell_type>__<dataset>"}. This splits each
#' dimname on the *last* \code{"__"} (cell type names may themselves contain
#' underscores, and some dataset labels, e.g. \code{PMID_39402379}, contain a
#' single underscore too), maps the dataset half to a display contrast label
#' (e.g. \code{"BICAN" -> "age"}) via \code{contrast_map}, and cleans the cell
#' type half with the same \code{__}/\code{_} -> space rule used by
#' \code{clean_cor_mat_names()}.
#'
#' Must be called on the raw (uncleaned) matrix dimnames, before
#' \code{clean_cor_mat_names()} runs -- that function replaces \code{__} and
#' \code{_} with the same space, which destroys the cell_type/dataset
#' boundary irrecoverably.
#'
#' @param cor_mat A matrix whose row and column names are identical and of
#'   the form \code{"<cell_type>__<dataset>"}.
#' @param contrast_map A two-column data frame: column 1 is the dataset label
#'   (as it appears in the dimnames), column 2 is the display contrast label.
#' @return A list with \code{group_labels} (contrast) and \code{name_labels}
#'   (cleaned cell type), both in \code{dimnames(cor_mat)} order.
.split_cor_mat_labels <- function(cor_mat, contrast_map) {
  labels <- rownames(cor_mat)
  if (!identical(labels, colnames(cor_mat))) {
    stop(".split_cor_mat_labels() expects a matrix with identical row/column names.")
  }

  # Cell type names use single underscores only; the "__" inserted by
  # paste(cell_type, dataset, sep = "__") is the sole double-underscore in
  # the string, so a greedy match up to the *last* "__" correctly splits on
  # it even when the dataset half itself contains single underscores (e.g.
  # "PMID_39402379").
  cell_type <- sub("^(.+)__.*$", "\\1", labels)
  dataset <- sub("^.+__", "", labels)

  match_idx <- match(dataset, contrast_map[[1]])
  if (anyNA(match_idx)) {
    stop(
      "Unmapped dataset(s) in contrast_map: ",
      paste(unique(dataset[is.na(match_idx)]), collapse = ", ")
    )
  }
  group_labels <- contrast_map[[2]][match_idx]

  clean_vec <- function(x) {
    x <- gsub("__", " ", x, fixed = TRUE)
    gsub("_", " ", x, fixed = TRUE)
  }

  list(
    group_labels = group_labels,
    name_labels = clean_vec(cell_type)
  )
}


clean_cor_mat_names <- function(cor_mat) {
  if (!is.matrix(cor_mat)) {
    stop("clean_cor_mat_names expects a matrix.")
  }

  rn <- rownames(cor_mat)
  cn <- colnames(cor_mat)

  clean_vec <- function(x) {
    if (is.null(x)) {
      return(x)
    }

    # First replace double underscores
    x <- gsub("__", " ", x, fixed = TRUE)

    # Then replace single underscores
    x <- gsub("_", " ", x, fixed = TRUE)

    x
  }

  rownames(cor_mat) <- clean_vec(rn)
  colnames(cor_mat) <- clean_vec(cn)

  cor_mat
}


#' Read and combine two DE datasets whose cell types use different
#' nomenclatures
#'
#' Like \code{bican.mccarroll.de.analysis::read_de_results_datasets()}, but
#' for a dataset 2 that names its cell types differently than \code{ct_file}
#' (e.g. Kana et al. 2026's own striatal cell type names, vs. BICAN's). Reads
#' dataset 1 and filters it to \code{ct_file} exactly as
#' \code{read_de_results_datasets()} would. Reads dataset 2 unfiltered (its
#' native names aren't in \code{ct_file} yet), renames its \code{cell_type}
#' column via \code{cell_type_map} (a two-column data frame: column 1 is
#' dataset 2's native cell type name, column 2 is the corresponding
#' \code{ct_file}/dataset-1 name), drops any dataset 2 row whose native cell
#' type has no entry in \code{cell_type_map}, and only then filters to
#' \code{ct_file} - matching what \code{read_de_results_datasets()} does for
#' datasets that already share one vocabulary.
#'
#' @param de_dir1,test1,dataset1 Passed to
#'   \code{bican.mccarroll.de.analysis::read_de_results()} for dataset 1.
#' @param de_dir2,test2,dataset2 Passed to
#'   \code{bican.mccarroll.de.analysis::read_de_results()} for dataset 2.
#' @param ct_file Path to the cell types file, in dataset 1's/the target
#'   nomenclature.
#' @param gene_to_chr A data.table with at least "gene" and "chr", shared by
#'   both datasets.
#' @param cell_type_map A two-column data frame mapping dataset 2's native
#'   cell type names (column 1) to their \code{ct_file} counterparts
#'   (column 2).
#' @return A combined data.table (rbind of both datasets'
#'   \code{read_de_results()} output, dataset 2 renamed) with an added
#'   \code{dataset} column, exactly as
#'   \code{bican.mccarroll.de.analysis::read_de_results_datasets()} returns.
.read_de_results_datasets_mapped <- function(de_dir1, test1, dataset1,
                                             de_dir2, test2, dataset2,
                                             ct_file, gene_to_chr,
                                             cell_type_map) {
  dataset <- cell_type <- NULL

  de1 <- bican.mccarroll.de.analysis::read_de_results(de_dir1, test1, ct_file, gene_to_chr)
  de1[, dataset := dataset1]

  de2 <- bican.mccarroll.de.analysis::read_de_results(de_dir2, test2, ct_file = NULL, gene_to_chr)

  native_names <- cell_type_map[[1]]
  target_names <- cell_type_map[[2]]

  match_idx <- match(de2$cell_type, native_names)
  keep <- !is.na(match_idx)
  de2 <- de2[keep]
  de2[, cell_type := target_names[match_idx[keep]]]

  cell_types_use <- bican.mccarroll.de.analysis::read_cell_types(ct_file)
  de2 <- de2[cell_type %in% cell_types_use]
  de2[, dataset := dataset2]

  rbind(de1, de2)
}

get_or_build_de_cor_mat_cache <- function(cache_file,
                                          de_region_interaction_dir,
                                          test,
                                          ct_file,
                                          gene_to_chr_file,
                                          regions_use,
                                          non_neuron_types,
                                          fdr_cutoff) {
  ## Cache is the matrix only.
  if (file.exists(cache_file)) {
    cor_dt <- data.table::fread(cache_file)

    rn <- cor_dt[[1]]
    cor_dt[[1]] <- NULL

    cor_mat <- as.matrix(cor_dt)
    rownames(cor_mat) <- rn

    return(cor_mat)
  }

  ## No cache: read what we need and compute the matrix.
  gene_to_chr <- bican.mccarroll.de.analysis::read_gene_to_chr(gene_to_chr_file)

  de_ri <- bican.mccarroll.de.analysis::read_de_results(
    de_region_interaction_dir,
    test,
    ct_file,
    gene_to_chr
  )

  cell_types_use <- bican.mccarroll.de.analysis::read_cell_types(ct_file)

  cor_mat <- bican.mccarroll.de.analysis::compute_de_cor_mat(
    de_ri,
    cell_types_use,
    regions_use,
    non_neuron_types,
    fdr_cutoff = fdr_cutoff
  )

  ## Write matrix cache as a reviewer-friendly TSV.
  cor_dt <- as.data.frame(cor_mat, check.names = FALSE)
  cor_dt <- cbind(cell_type = rownames(cor_mat), cor_dt)

  utils::write.table(
    cor_dt,
    file = cache_file,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  cor_mat
}


resolve_de_cor_paths <- function(de_region_interaction_dir = NULL,
                                 gene_to_chr_file = NULL,
                                 ct_file = NULL,
                                 outDir = NULL,
                                 data_cache_dir = NULL) {
  root <- .resolve_data_root_dir(NULL)

  rel <- list(
    de_region_interaction_dir =
      "differential_expression/results/LEVEL_6/sex_age/cell_type_region_interaction_absolute_effects",
    gene_to_chr_file =
      "metadata/gene_to_chromosome.txt",
    ct_file =
      "differential_expression/metadata/cell_types_for_de_filtering_plot.txt"
  )

  pick_in <- function(x, key) {
    if (is.null(x)) {
      return(file.path(root, rel[[key]]))
    }
    .resolve_under_root(root, x)
  }

  out <- .resolve_out_dir(outDir)
  cache <- .resolve_cache_dir(data_cache_dir)

  # if a cache wasn't set, then use the differential_expression subdirectiory.
  if (is.null(data_cache_dir)) {
    cache <- file.path(cache, "differential_expression")
  }

  .ensure_dir(out)
  .ensure_dir(cache)

  list(
    data_root_dir             = root,
    de_region_interaction_dir = pick_in(de_region_interaction_dir, "de_region_interaction_dir"),
    gene_to_chr_file          = pick_in(gene_to_chr_file, "gene_to_chr_file"),
    ct_file                   = pick_in(ct_file, "ct_file"),
    outDir                    = out,
    data_cache_dir            = cache
  )
}
