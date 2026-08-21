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
# plot_de_hits_age_bmi_sex_smoking()

#' Plot the DE-gene-hit heatmap comparing age, sex, BMI, and smoking
#'
#' Wires the BICAN LEVEL_6 age/BMI/smoking-complete-cases differential
#' expression results into
#' \code{bican.mccarroll.differentialexpression::compute_de_region_contrast_counts()}
#' and \code{bican.mccarroll.differentialexpression::plot_de_region_contrast_heatmap()},
#' comparing the number of significant DE genes for the \code{age},
#' \code{female_vs_male}, \code{bmi}, and \code{smoker_vs_non_smoker} contrasts
#' (displayed in that order: Age, Sex, BMI, Smoking) across regions and cell
#' types. Writes two SVGs to the configured figure output directory (see
#' \code{\link{get_out_dir}}) — a color (viridis) version and a no-color
#' version (white cells, black borders, numbers only).
#'
#' The underlying gene-count table is cached as a plain-text TSV file to speed
#' up subsequent runs and to provide reviewer-inspectable intermediate data;
#' it is not a primary output.
#'
#' @param force_recompute Logical scalar. If \code{TRUE}, ignore an existing
#'   cache file, recompute the gene-count table, and overwrite the cache.
#'   Defaults to \code{FALSE}.
#'
#' @return Invisibly returns the data.table produced by
#'   \code{bican.mccarroll.differentialexpression::compute_de_region_contrast_counts()}.
#'
#' @export
plot_de_hits_age_bmi_sex_smoking <- function(force_recompute = FALSE) {
  # Make R CMD CHECK happy
  region <- contrast <- cell_type <- NULL

  paths <- resolve_de_hits_age_bmi_sex_smoking_paths()

  cache_file <- file.path(paths$data_cache_dir, "de_hits_age_bmi_sex_smoking.tsv")

  if (!force_recompute && file.exists(cache_file)) {
    de_counts <- data.table::fread(cache_file)

    # The cache is written in the exact row order produced by
    # compute_de_region_contrast_counts()'s setorder(region, contrast,
    # cell_type), so region/contrast factor levels can be recovered directly
    # from row-appearance order without re-reading any DE files. This is
    # exact for "complete cases" data where every contrast is tested across
    # the same full set of regions (the point of using complete cases here).
    de_counts[, region := factor(region, levels = unique(region))]
    de_counts[, contrast := factor(contrast, levels = unique(contrast))]

    ct_order <- scan(paths$ct_file, what = character(), quiet = TRUE)
    present <- unique(as.character(de_counts$cell_type))
    de_counts[, cell_type := factor(cell_type, levels = intersect(ct_order, present))]
  } else {
    de_counts <- bican.mccarroll.differentialexpression::compute_de_region_contrast_counts(
      paths = rep(paths$in_dir, 4),
      contrasts = c("age", "female_vs_male", "bmi", "smoker_vs_non_smoker"),
      cellTypeListFile = paths$ct_file
    )

    utils::write.table(
      de_counts,
      file = cache_file,
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE
    )
  }

  contrast_labels <- c(
    age = "Age",
    bmi = "BMI",
    female_vs_male = "Sex",
    smoker_vs_non_smoker = "Smoking"
  )

  p_color <- bican.mccarroll.differentialexpression::plot_de_region_contrast_heatmap(
    de_counts,
    label_threshold = 10000,
    contrast_labels = contrast_labels
  )
  p_no_color <- bican.mccarroll.differentialexpression::plot_de_region_contrast_heatmap(
    de_counts,
    label_threshold = 10000,
    contrast_labels = contrast_labels,
    no_color = TRUE
  )

  save_plot_svg(
    p_color,
    out_file = "de_hits_age_bmi_sex_smoking.svg",
    out_dir = paths$outDir,
    width = 10,
    height = 7.5
  )
  save_plot_svg(
    p_no_color,
    out_file = "de_hits_age_bmi_sex_smoking_no_color.svg",
    out_dir = paths$outDir,
    width = 10,
    height = 7.5
  )

  logger::log_info("DONE plotting DE hits heatmap for age/BMI/sex/smoking")

  invisible(de_counts)
}

resolve_de_hits_age_bmi_sex_smoking_paths <- function() {
  root <- .resolve_data_root_dir(NULL)
  out <- .resolve_out_dir(NULL)
  cache <- file.path(.resolve_cache_dir(NULL), "differential_expression")

  .ensure_dir(out)
  .ensure_dir(cache)

  list(
    data_root_dir = root,
    in_dir = file.path(
      root,
      "differential_expression/results/LEVEL_6/age_smoking_bmi_complete_cases/age_smoking_bmi"
    ),
    ct_file = file.path(root, "differential_expression/metadata/cell_types_for_de_filtering_plot.txt"),
    outDir = out,
    data_cache_dir = cache
  )
}
