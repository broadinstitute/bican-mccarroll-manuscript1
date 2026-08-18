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
  grDevices::svg(out_file, width = 6, height = 6)

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
  grDevices::svg(out_file, width = 10, height = 10)

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
#' with some other PMID cell type, or with neither. Row/column labels are
#' distinguished only by a dataset text suffix (no forced grouping/split),
#' matching the existing \code{cell_type__region} labeling convention.
#'
#' The correlation matrix is cached to a plain-text TSV file to speed up
#' subsequent runs and to provide reviewer-inspectable intermediate data.
#' Cell type/dataset combinations with fewer than \code{min_num_genes} of
#' their own significant genes are dropped from the matrix (see
#' \code{bican.mccarroll.de.analysis::compute_de_cor_mat_datasets()}), so a
#' cache built at one \code{min_num_genes} value must be regenerated
#' (\code{force_recompute = TRUE}) if that value changes.
#'
#' @param min_num_genes Minimum number of a cell type/dataset combination's
#'   own significant genes required to keep it in the matrix. Defaults to
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

  cor_mat <- clean_cor_mat_names(cor_mat)

  out_file <- file.path(
    bican_paths$outDir,
    sprintf("de_cor_heatmap_bican_sea_ad_vs_pmid_39402379_%s.svg", contrast)
  )
  grDevices::svg(out_file, width = 12, height = 12)

  ht <- bican.mccarroll.de.analysis::plot_de_cor_heatmap_complex(
    cor_mat,
    clustering_method = "complete",
    breaks = seq(-1, 1, length.out = 101),
    palette_colors = c("steelblue", "white", "darkorange"),
    legend_title = NULL,
    show_dendrograms = TRUE
  )

  ComplexHeatmap::draw(ht)
  grDevices::dev.off()

  logger::log_info("DONE plotting DE correlation heatmap: {out_file}")

  invisible(cor_mat)
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
