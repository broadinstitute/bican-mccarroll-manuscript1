# Example call (explicit paths) ----
# estimate_donor_superpopulation(
#     kgp_donors_file = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/ancestry_pca/inputs/1kg_donors.txt",
#     kgp_ancestry_file = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/ancestry_pca/inputs/kgp_ancestry.txt",
#     pca_files = c(
#         "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/ancestry_pca/inputs/smartpca_2025-01-06.evec",
#         "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/ancestry_pca/inputs/smartpca_2025-03-31.evec"
#     ),
#     output_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/ancestry_pca",
#     bican_donors_file = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/metadata/donor_list_178.txt",
#     n_pcs = 5,
#     k_neighbors = 15
# )

# Example call (configured paths) ----
# source("R/paths.R")
# options(
#     bican.mccarroll.figures.data_root_dir =
#         "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis",
#     bican.mccarroll.figures.out_dir =
#         "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/figure_repository",
#     bican.mccarroll.figures.cache_dir =
#         "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/figure_repository/data_cache"
# )
# plot_ancestry_pca()

#' Estimate donor superpopulation labels from PCA coordinates
#'
#' Assigns BICAN donors to 1000 Genomes superpopulations using
#' distance-weighted k-nearest neighbors in principal component space.
#' Reference samples are labeled as AFR, AMR, EAS, EUR, or SAS, and each
#' BICAN donor is assigned to the population receiving the largest weighted
#' vote among its nearest reference samples.
#'
#' The function writes a table containing the predicted superpopulation and
#' weighted population scores for each BICAN donor. It also writes an SVG
#' showing the BICAN donors and 1000 Genomes reference samples in PC1-PC2
#' space. If \code{bican_donors_file} is provided, only the listed BICAN
#' donors are included in the predictions and plot.
#'
#' @param kgp_donors_file Path to a file containing the 1000 Genomes donor
#'   IDs.
#' @param kgp_ancestry_file Path to a file containing 1000 Genomes donor IDs
#'   and their superpopulation labels.
#' @param pca_files Character vector of paths to SmartPCA \code{.evec} files.
#'   The files must contain donor IDs followed by PC1-PC10 and a final
#'   SmartPCA label column.
#' @param output_dir Directory where the prediction table and PCA SVG are
#'   written.
#' @param bican_donors_file Optional path to a file containing BICAN donor IDs
#'   to retain. If \code{NULL}, all BICAN donors in the PCA data are used.
#' @param n_pcs Number of principal components to use when calculating nearest
#'   neighbors.
#' @param k_neighbors Number of nearest 1000 Genomes reference samples used
#'   for classification.
#'
#' @return Invisibly returns a data frame containing the BICAN donor
#'   predictions and population scores.
#'
#' @export
estimate_donor_superpopulation <- function(
  kgp_donors_file,
  kgp_ancestry_file,
  pca_files,
  output_dir,
  bican_donors_file = NULL,
  n_pcs = 5,
  k_neighbors = 15
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  pca_df <- .build_superpopulation_pca_df(
    kgp_donors_file,
    kgp_ancestry_file,
    pca_files,
    bican_donors_file,
    n_pcs,
    k_neighbors
  )

  bican_predictions <- .bican_superpopulation_predictions(pca_df, n_pcs)

  utils::write.table(
    bican_predictions,
    file.path(output_dir, "bican_superpopulation_predictions.txt"),
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
  )

  plot <- .plot_superpopulations(
    pca_df,
    n_pcs,
    k_neighbors
  )

  ggplot2::ggsave(
    file.path(output_dir, "pca_superpopulation_light_dark.svg"),
    plot = plot,
    device = svglite_manuscript,
    width = 8,
    height = 6
  )

  print(table(bican_predictions$predicted_superpop))

  invisible(bican_predictions)
}


#' Generate the ancestry PCA superpopulation figure using configured paths
#'
#' Options-driven wrapper around \code{estimate_donor_superpopulation()}. All
#' path arguments default to \code{NULL}, in which case they are resolved
#' using \code{resolve_ancestry_pca_paths()} relative to the configured data
#' root, output, and cache directory options (see \code{\link{get_data_root_dir}},
#' \code{\link{get_out_dir}}, \code{\link{get_cache_dir}}).
#'
#' The classified PCA table (1000 Genomes reference samples plus BICAN
#' donors) is cached as a TSV file under the resolved cache directory so
#' that repeated calls skip the k-nearest-neighbor classification. Changing
#' \code{n_pcs} or \code{k_neighbors} does not invalidate this cache; delete
#' \code{ancestry_pca_classified.tsv} from the cache directory to force
#' recomputation.
#'
#' @param kgp_donors_file Path to a file containing the 1000 Genomes donor
#'   IDs. If \code{NULL}, resolved under the data root.
#' @param kgp_ancestry_file Path to a file containing 1000 Genomes donor IDs
#'   and their superpopulation labels. If \code{NULL}, resolved under the
#'   data root.
#' @param pca_files Character vector of paths to SmartPCA \code{.evec} files.
#'   If \code{NULL}, resolved under the data root.
#' @param bican_donors_file Optional path to a file containing BICAN donor IDs
#'   to retain. If \code{NULL}, resolved under the data root.
#' @param outDir Output directory for the generated SVG. If \code{NULL},
#'   resolved via configured output directory options.
#' @param data_cache_dir Directory used to cache the classified PCA table and
#'   BICAN predictions table. If \code{NULL}, resolved via configured cache
#'   directory options.
#' @param n_pcs Number of principal components to use when calculating nearest
#'   neighbors.
#' @param k_neighbors Number of nearest 1000 Genomes reference samples used
#'   for classification.
#'
#' @return Invisibly returns a data frame containing the BICAN donor
#'   predictions and population scores.
#'
#' @export
plot_ancestry_pca <- function(
  kgp_donors_file = NULL,
  kgp_ancestry_file = NULL,
  pca_files = NULL,
  bican_donors_file = NULL,
  outDir = NULL,
  data_cache_dir = NULL,
  n_pcs = 5,
  k_neighbors = 15
) {
  paths <- resolve_ancestry_pca_paths(
    kgp_donors_file = kgp_donors_file,
    kgp_ancestry_file = kgp_ancestry_file,
    pca_files = pca_files,
    bican_donors_file = bican_donors_file,
    outDir = outDir,
    data_cache_dir = data_cache_dir
  )

  pca_df <- get_or_build_ancestry_pca_cache(paths, n_pcs, k_neighbors)

  bican_predictions <- .bican_superpopulation_predictions(pca_df, n_pcs)

  utils::write.table(
    bican_predictions,
    file.path(paths$data_cache_dir, "bican_superpopulation_predictions.txt"),
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
  )

  plot <- .plot_superpopulations(
    pca_df,
    n_pcs,
    k_neighbors
  )

  ggplot2::ggsave(
    file.path(paths$outDir, "ancestry_pca.svg"),
    plot = plot,
    device = svglite_manuscript,
    width = 8,
    height = 6
  )

  print(table(bican_predictions$predicted_superpop))

  invisible(bican_predictions)
}


.build_superpopulation_pca_df <- function(
  kgp_donors_file,
  kgp_ancestry_file,
  pca_files,
  bican_donors_file,
  n_pcs,
  k_neighbors
) {
  # Make R CMD CHECK Happy
  project <- donor_external_id <- NULL

  pca_df <- .gather_ancestry_data(
    kgp_donors_file,
    kgp_ancestry_file,
    pca_files
  )

  if (!is.null(bican_donors_file)) {
    bican_donors <- utils::read.table(bican_donors_file)$V1

    pca_df <- pca_df |>
      dplyr::filter(
        project == "1000 Genomes" |
          donor_external_id %in% bican_donors
      )
  }

  .classify_superpopulations(
    pca_df,
    n_pcs,
    k_neighbors
  )
}


.bican_superpopulation_predictions <- function(pca_df, n_pcs) {
  # Make R CMD CHECK Happy
  project <- donor_external_id <- predicted_superpop <- NULL
  superpop_score <- superpop_margin <- NULL
  AFR_score <- AMR_score <- EAS_score <- EUR_score <- SAS_score <- NULL

  pc_cols <- paste0("PC", seq_len(n_pcs))

  pca_df |>
    dplyr::filter(project == "BICAN") |>
    dplyr::select(
      donor_external_id,
      dplyr::all_of(pc_cols),
      predicted_superpop,
      superpop_score,
      superpop_margin,
      AFR_score,
      AMR_score,
      EAS_score,
      EUR_score,
      SAS_score
    )
}


resolve_ancestry_pca_paths <- function(
  kgp_donors_file = NULL,
  kgp_ancestry_file = NULL,
  pca_files = NULL,
  bican_donors_file = NULL,
  outDir = NULL,
  data_cache_dir = NULL
) {
  root <- .resolve_data_root_dir(NULL)

  rel <- list(
    kgp_donors_file =
      "ancestry_pca/inputs/1kg_donors.txt",
    kgp_ancestry_file =
      "ancestry_pca/inputs/kgp_ancestry.txt",
    pca_files = c(
      "ancestry_pca/inputs/smartpca_2025-01-06.evec",
      "ancestry_pca/inputs/smartpca_2025-03-31.evec"
    ),
    bican_donors_file =
      "metadata/donor_list_178.txt"
  )

  pick_in <- function(x, key) {
    if (is.null(x)) {
      return(file.path(root, rel[[key]]))
    }
    .resolve_under_root(root, x)
  }

  if (is.null(pca_files)) {
    pca_files <- file.path(root, rel[["pca_files"]])
  } else {
    pca_files <- vapply(
      pca_files,
      function(p) .resolve_under_root(root, p),
      character(1),
      USE.NAMES = FALSE
    )
  }

  out <- .resolve_out_dir(outDir)
  cache <- .resolve_cache_dir(data_cache_dir)

  # if a cache wasn't set, then use the ancestry_pca subdirectory.
  if (is.null(data_cache_dir)) {
    cache <- file.path(cache, "ancestry_pca")
  }

  .ensure_dir(out)
  .ensure_dir(cache)

  list(
    data_root_dir      = root,
    kgp_donors_file    = pick_in(kgp_donors_file, "kgp_donors_file"),
    kgp_ancestry_file  = pick_in(kgp_ancestry_file, "kgp_ancestry_file"),
    pca_files          = pca_files,
    bican_donors_file  = pick_in(bican_donors_file, "bican_donors_file"),
    outDir             = out,
    data_cache_dir     = cache
  )
}


get_or_build_ancestry_pca_cache <- function(paths, n_pcs, k_neighbors) {
  # Make R CMD CHECK Happy
  predicted_superpop <- NULL

  superpop_levels <- c("AFR", "AMR", "EAS", "EUR", "SAS")
  cache_file <- file.path(
    paths$data_cache_dir,
    "ancestry_pca_classified.tsv"
  )

  if (file.exists(cache_file)) {
    logger::log_info("Using cached data from {cache_file}")
    pca_df <- as.data.frame(
      data.table::fread(
        cache_file,
        header = TRUE,
        sep = "\t",
        na.strings = "NA"
      )
    )
    pca_df$predicted_superpop <- factor(
      pca_df$predicted_superpop,
      levels = superpop_levels
    )
    return(pca_df)
  }

  logger::log_info(
    "No cached data from {cache_file} regenerating data from sources.  This can take a few minutes"
  )

  pca_df <- .build_superpopulation_pca_df(
    paths$kgp_donors_file,
    paths$kgp_ancestry_file,
    paths$pca_files,
    paths$bican_donors_file,
    n_pcs,
    k_neighbors
  )

  utils::write.table(
    pca_df,
    file = cache_file,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    na = "NA"
  )

  pca_df
}


.gather_ancestry_data <- function(
  kgp_donors_file,
  kgp_ancestry_file,
  pca_files
) {
  # Make R CMD CHECK Happy
  V12 <- donor_external_id <- IID <- superpop <- NULL

  kgp_donors <- utils::read.table(kgp_donors_file)

  kgp_donor_ancestry <- utils::read.table(
    kgp_ancestry_file,
    header = TRUE
  )

  pca_list <- lapply(pca_files, utils::read.table)

  pca_df <- dplyr::bind_rows(pca_list) |>
    dplyr::distinct() |>
    dplyr::select(-V12)

  colnames(pca_df) <- c(
    "donor_external_id",
    "PC1",
    "PC2",
    "PC3",
    "PC4",
    "PC5",
    "PC6",
    "PC7",
    "PC8",
    "PC9",
    "PC10"
  )

  pca_df |>
    dplyr::mutate(
      project = ifelse(
        donor_external_id %in% kgp_donors$V1,
        "1000 Genomes",
        "BICAN"
      )
    ) |>
    dplyr::left_join(
      kgp_donor_ancestry |>
        dplyr::select(IID, superpop),
      by = c("donor_external_id" = "IID")
    )
}


.classify_superpopulations <- function(
  pca_df,
  n_pcs,
  k_neighbors
) {
  superpop_levels <- c("AFR", "AMR", "EAS", "EUR", "SAS")
  pc_cols <- paste0("PC", seq_len(n_pcs))

  reference_idx <- pca_df$project == "1000 Genomes" &
    pca_df$superpop %in% superpop_levels

  query_idx <- pca_df$project == "BICAN"

  reference_pcs <- as.matrix(
    pca_df[reference_idx, pc_cols, drop = FALSE]
  )

  query_pcs <- as.matrix(
    pca_df[query_idx, pc_cols, drop = FALSE]
  )

  if (anyNA(reference_pcs)) {
    stop("Missing PC values were found among 1000 Genomes reference samples.")
  }

  if (anyNA(query_pcs)) {
    stop("Missing PC values were found among BICAN samples.")
  }

  reference_labels <- factor(
    pca_df$superpop[reference_idx],
    levels = superpop_levels
  )

  k_neighbors <- min(k_neighbors, nrow(reference_pcs))

  nearest <- FNN::get.knnx(
    reference_pcs,
    query_pcs,
    k = k_neighbors
  )

  classify_one_donor <- function(i) {
    .classify_donor(
      nearest$nn.index[i, ],
      nearest$nn.dist[i, ],
      reference_labels,
      superpop_levels
    )
  }

  predictions <- do.call(
    rbind,
    lapply(seq_len(nrow(query_pcs)), classify_one_donor)
  )

  pca_df$predicted_superpop <- as.character(pca_df$superpop)
  pca_df$superpop_score <- NA_real_
  pca_df$superpop_margin <- NA_real_

  score_cols <- paste0(superpop_levels, "_score")

  for (score_col in score_cols) {
    pca_df[[score_col]] <- NA_real_
  }

  pca_df$predicted_superpop[query_idx] <- predictions$predicted_superpop
  pca_df$superpop_score[query_idx] <- predictions$superpop_score
  pca_df$superpop_margin[query_idx] <- predictions$superpop_margin

  for (score_col in score_cols) {
    pca_df[[score_col]][query_idx] <- predictions[[score_col]]
  }

  pca_df$predicted_superpop <- factor(
    pca_df$predicted_superpop,
    levels = superpop_levels
  )

  pca_df
}


.classify_donor <- function(
  neighbor_indices,
  neighbor_distances,
  reference_labels,
  superpop_levels
) {
  neighbor_labels <- as.character(
    reference_labels[neighbor_indices]
  )

  weights <- 1 / pmax(neighbor_distances, 1e-12)

  vote_totals <- stats::setNames(
    numeric(length(superpop_levels)),
    superpop_levels
  )

  weighted_votes <- tapply(
    weights,
    neighbor_labels,
    sum
  )

  vote_totals[names(weighted_votes)] <- weighted_votes
  vote_proportions <- vote_totals / sum(vote_totals)
  ordered_scores <- sort(vote_proportions, decreasing = TRUE)

  data.frame(
    predicted_superpop = names(ordered_scores)[1],
    superpop_score = unname(ordered_scores[1]),
    superpop_margin = unname(
      ordered_scores[1] - ordered_scores[2]
    ),
    AFR_score = unname(vote_proportions["AFR"]),
    AMR_score = unname(vote_proportions["AMR"]),
    EAS_score = unname(vote_proportions["EAS"]),
    EUR_score = unname(vote_proportions["EUR"]),
    SAS_score = unname(vote_proportions["SAS"]),
    row.names = NULL
  )
}


.plot_superpopulations <- function(
  pca_df,
  n_pcs,
  k_neighbors
) {
  superpop_levels <- c("AFR", "AMR", "EAS", "EUR", "SAS")
  bican_display_label <- "Burger et al."

  reference_colors <- c(
    AFR = "#A1D99B",
    AMR = "#B39DDB",
    EAS = "#FCAE91",
    EUR = "#9ECAE1",
    SAS = "#FDD0A2"
  )

  bican_colors <- c(
    AFR = "#238B45",
    AMR = "#6A51A3",
    EAS = "#CB181D",
    EUR = "#08519C",
    SAS = "#D95F02"
  )

  pca_plot_df <- pca_df |>
    dplyr::filter(!is.na(predicted_superpop)) |>
    dplyr::mutate(
      project_label = ifelse(project == "BICAN", bican_display_label, project),
      plot_group = paste(
        predicted_superpop,
        project_label,
        sep = " - "
      )
    )

  plot_group_levels <- c(
    paste(superpop_levels, "1000 Genomes", sep = " - "),
    paste(superpop_levels, bican_display_label, sep = " - ")
  )

  pca_plot_df$plot_group <- factor(
    pca_plot_df$plot_group,
    levels = plot_group_levels
  )

  plot_group_colors <- c(
    stats::setNames(
      reference_colors,
      paste(superpop_levels, "1000 Genomes", sep = " - ")
    ),
    stats::setNames(
      bican_colors,
      paste(superpop_levels, bican_display_label, sep = " - ")
    )
  )

  # Make R CMD CHECK Happy
  PC1 <- PC2 <- plot_group <- predicted_superpop <- project <- NULL

  ggplot2::ggplot() +
    ggplot2::geom_point(
      data = pca_plot_df[
        pca_plot_df$project == "1000 Genomes",
      ],
      ggplot2::aes(
        x = PC1,
        y = PC2,
        color = plot_group
      ),
      alpha = 0.65,
      size = 1.3
    ) +
    ggplot2::geom_point(
      data = pca_plot_df[
        pca_plot_df$project == "BICAN",
      ],
      ggplot2::aes(
        x = PC1,
        y = PC2,
        color = plot_group
      ),
      alpha = 0.9,
      size = 1.8
    ) +
    ggplot2::scale_color_manual(
      values = plot_group_colors,
      drop = FALSE,
      na.translate = FALSE
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(
          size = 2,
          alpha = 1
        )
      )
    ) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Donor PCA embedding with 1000 Genomes reference",
      subtitle = paste0(
        bican_display_label,
        " labels assigned by distance-weighted ",
        k_neighbors,
        "-nearest neighbors using PCs 1-",
        n_pcs
      ),
      color = "Population and dataset",
      x = "PC1",
      y = "PC2"
    )
}
