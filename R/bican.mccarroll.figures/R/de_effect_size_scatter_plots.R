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

# ct_file <- gene_to_chr_file <- de_results_dir<- de_region_interaction_dir<- data_cache_dir<- outDir<- NULL

#' Generate sex and age differential expression scatter plots for manuscript figures
#'
#' @param ct_file Path to the cell type list file. If \code{NULL}, resolved
#'   under the data root.
#' @param gene_to_chr_file Path to the gene-to-chromosome mapping file. If
#'   \code{NULL}, resolved under the data root.
#' @param de_results_dir Base directory for differential expression results.
#'   If \code{NULL}, resolved under the data root.
#' @param de_region_interaction_dir Directory containing region-interaction
#'   differential expression results. If \code{NULL}, resolved under
#'   \code{de_results_dir}.
#' @param data_cache_dir Directory used to store cached DE tables as TSV
#'   files. If \code{NULL}, resolved via configured cache directory options.
#' @param outDir Output directory for generated SVG plots. If \code{NULL},
#'   resolved via configured output directory options.
#'
#' @export
de_sex_age_scatter_plots <- function(ct_file = NULL,
                                     gene_to_chr_file = NULL,
                                     de_results_dir = NULL,
                                     de_region_interaction_dir = NULL,
                                     data_cache_dir = NULL,
                                     outDir = NULL) {
  paths <- .resolve_sex_age_scatter_plot_paths(
    ct_file = ct_file,
    gene_to_chr_file = gene_to_chr_file,
    de_results_dir = de_results_dir,
    de_region_interaction_dir = de_region_interaction_dir,
    data_cache_dir = data_cache_dir,
    outDir = outDir
  )

  cache_file_age <- file.path(paths$data_cache_dir, "de_age.tsv")
  cache_file_sex <- file.path(paths$data_cache_dir, "de_sex.tsv")

  age_df <- get_or_build_de_cache(
    test = "age",
    cache_file = cache_file_age,
    de_region_interaction_dir = paths$de_region_interaction_dir,
    ct_file = paths$ct_file,
    gene_to_chr_file = paths$gene_to_chr_file
  )

  sex_df <- get_or_build_de_cache(
    test = "female_vs_male",
    cache_file = cache_file_sex,
    de_region_interaction_dir = paths$de_region_interaction_dir,
    ct_file = paths$ct_file,
    gene_to_chr_file = paths$gene_to_chr_file
  )

  #########################
  # Main figure plot group
  #########################

  # plot_de_scatter_svg(
  #     df = sex_df,
  #     test = "sex",
  #     cell_type1 = "OPC",
  #     cell_type2 = "astrocyte",
  #     region1 = "CaH",
  #     region2 = "CaH",
  #     xlab_prefix="Sex DE, ",
  #     outDir = paths$outDir,
  #     fdr_cutoff = 0.05,
  #     add_fit = TRUE)

  plot_de_scatter_svg(
    df = age_df,
    test = "age",
    cell_type1 = "SPN_D1_matrix",
    cell_type2 = "SPN_D2_matrix",
    region1 = "CaH",
    region2 = "NAC",
    xlab_prefix = "Age DE, ",
    outDir = paths$outDir,
    fdr_cutoff = 0.05,
    add_fit = TRUE
  )

  plot_de_scatter_svg(
    df = age_df,
    test = "age",
    cell_type1 = "SPN_D1_matrix",
    cell_type2 = "SPN_D2_matrix",
    region1 = "CaH",
    region2 = "CaH",
    xlab_prefix = "Age DE, ",
    outDir = paths$outDir,
    fdr_cutoff = 0.05,
    add_fit = TRUE
  )

  plot_de_scatter_svg(
    df = age_df,
    test = "age",
    cell_type1 = "astrocyte",
    cell_type2 = "astrocyte",
    region1 = "CaH",
    region2 = "DFC",
    xlab_prefix = "Age DE, ",
    outDir = paths$outDir,
    fdr_cutoff = 0.05,
    add_fit = TRUE
  )

  plot_de_scatter_svg(
    df = age_df,
    test = "age",
    cell_type1 = "OPC",
    cell_type2 = "astrocyte",
    region1 = "CaH",
    region2 = "CaH",
    xlab_prefix = "Age DE, ",
    outDir = paths$outDir,
    fdr_cutoff = 0.05,
    add_fit = TRUE
  )

  ####################
  # Supplemental 1:
  ####################

  plot_de_scatter_svg(
    df = age_df,
    test = "age",
    cell_type1 = "SPN_D1_matrix",
    cell_type2 = "SPN_D1_matrix",
    region1 = "CaH",
    region2 = "Pu",
    xlab_prefix = "Age DE, ",
    outDir = paths$outDir,
    fdr_cutoff = 0.05,
    add_fit = TRUE
  )

  plot_de_scatter_svg(
    df = age_df,
    test = "age",
    cell_type1 = "SPN_D1_matrix",
    cell_type2 = "SPN_D1_matrix",
    region1 = "CaH",
    region2 = "NAC",
    xlab_prefix = "Age DE, ",
    outDir = paths$outDir,
    fdr_cutoff = 0.05,
    add_fit = TRUE
  )

  plot_de_scatter_svg(
    df = age_df,
    test = "age",
    cell_type1 = "SPN_D1_matrix",
    cell_type2 = "SPN_D1_matrix",
    region1 = "Pu",
    region2 = "NAC",
    xlab_prefix = "Age DE, ",
    outDir = paths$outDir,
    fdr_cutoff = 0.05,
    add_fit = TRUE
  )

  # Supplemental 2:
  plot_de_scatter_svg(
    df = age_df,
    test = "age",
    cell_type1 = "astrocyte",
    cell_type2 = "astrocyte",
    region1 = "CaH",
    region2 = "Pu",
    xlab_prefix = "Age DE, ",
    outDir = paths$outDir,
    fdr_cutoff = 0.05,
    add_fit = TRUE
  )

  plot_de_scatter_svg(
    df = age_df,
    test = "age",
    cell_type1 = "astrocyte",
    cell_type2 = "astrocyte",
    region1 = "CaH",
    region2 = "NAC",
    xlab_prefix = "Age DE, ",
    outDir = paths$outDir,
    fdr_cutoff = 0.05,
    add_fit = TRUE
  )

  plot_de_scatter_svg(
    df = age_df,
    test = "age",
    cell_type1 = "astrocyte",
    cell_type2 = "astrocyte",
    region1 = "Pu",
    region2 = "NAC",
    xlab_prefix = "Age DE, ",
    outDir = paths$outDir,
    fdr_cutoff = 0.05,
    add_fit = TRUE
  )
}


#' Per-cell-type scatter plots of sex DE effect versus age DE effect
#'
#' For every cell type that has data in \code{region_use} (default "CaH"),
#' writes one SVG scatter plot of the sex DE log2 fold change (x) against the
#' age DE log2 fold change (y), restricted to the \code{top_n_expressed} most
#' highly expressed genes for that cell type and colored by chromosome using
#' the same X / Y / autosome color scheme as the volcano figures.
#'
#' Genes with no chromosome annotation and mitochondrial genes are dropped
#' before gene selection, so each panel plots \code{top_n_expressed} genes that
#' all carry a chromosome assignment.
#'
#' The list of cell types is derived from the data rather than hardcoded,
#' because several cell types are DFC-only by design and simply have no rows in
#' the CaH region. By default the "SPN_D1" and "SPN_D2" parent cell types are
#' also included (see \code{include_spn_d1_d2}), in addition to the
#' "_matrix"/"_striosome" subtypes already covered by the shared cell type
#' list.
#'
#' SVGs are written to a \code{de_sex_vs_age_scatter} subdirectory of the
#' configured figure output directory.
#'
#' @param ct_file Path to the cell type list file. If \code{NULL}, resolved
#'   under the data root.
#' @param gene_to_chr_file Path to the gene-to-chromosome mapping file. If
#'   \code{NULL}, resolved under the data root.
#' @param de_results_dir Base directory for differential expression results.
#'   If \code{NULL}, resolved under the data root.
#' @param de_region_interaction_dir Directory containing region-interaction
#'   differential expression results. If \code{NULL}, resolved under
#'   \code{de_results_dir}.
#' @param data_cache_dir Directory holding the cached \code{de_age.tsv} and
#'   \code{de_sex.tsv} tables. If \code{NULL}, resolved via configured cache
#'   directory options.
#' @param outDir Output directory for generated SVG plots. If \code{NULL},
#'   resolved via configured output directory options.
#' @param region_use Region to plot. Defaults to "CaH".
#' @param top_n_expressed Number of most highly expressed genes per cell type.
#'   Defaults to 4000.
#' @param autosomes_only If \code{TRUE}, chrX and chrY genes are dropped
#'   BEFORE gene selection, so each panel plots \code{top_n_expressed} of the
#'   most highly expressed AUTOSOMAL genes only (rather than the top
#'   \code{top_n_expressed} overall genes, which are dominated on the x axis by
#'   a handful of chrX/chrY outliers). All points are plotted black and the
#'   chromosome legend is omitted, since there is only one group. Output
#'   filenames get an \code{_autosomes_only} suffix so they do not collide with
#'   the default (all-chromosome) output. Defaults to \code{FALSE}.
#' @param sex_df Optional pre-loaded sex DE table, to avoid re-reading the
#'   cache during interactive work. If \code{NULL} it is loaded from the cache.
#' @param age_df Optional pre-loaded age DE table. See \code{sex_df}.
#' @param include_spn_d1_d2 If \code{TRUE} (the default), also include the
#'   parent "SPN_D1" and "SPN_D2" cell types (distinct from "SPN_D1_matrix" /
#'   "SPN_D1_striosome" / "SPN_D2_matrix" / "SPN_D2_striosome", which are
#'   always included via \code{ct_file}). These two are not part of the shared
#'   \code{cell_types_for_de_filtering_plot.txt} list used elsewhere in the
#'   package, so they are read from the raw region-interaction results and
#'   cached separately (\code{de_age_spn_d1_d2.tsv} /
#'   \code{de_sex_spn_d1_d2.tsv} in \code{data_cache_dir}), rather than
#'   expanding the shared cell type list that many other figures depend on.
#'
#' @return Invisibly returns a character vector of the SVG files written.
#' @export
de_sex_vs_age_effect_scatter_plots <- function(ct_file = NULL,
                                               gene_to_chr_file = NULL,
                                               de_results_dir = NULL,
                                               de_region_interaction_dir = NULL,
                                               data_cache_dir = NULL,
                                               outDir = NULL,
                                               region_use = "CaH",
                                               top_n_expressed = 4000,
                                               autosomes_only = FALSE,
                                               sex_df = NULL,
                                               age_df = NULL,
                                               include_spn_d1_d2 = TRUE) {
  # Make R CMD CHECK happy
  cell_type <- region <- chr <- NULL

  paths <- .resolve_sex_age_scatter_plot_paths(
    ct_file = ct_file,
    gene_to_chr_file = gene_to_chr_file,
    de_results_dir = de_results_dir,
    de_region_interaction_dir = de_region_interaction_dir,
    data_cache_dir = data_cache_dir,
    outDir = outDir
  )

  ## -----------------------
  ## Manuscript parameters
  ## -----------------------

  # Same scheme (and the same trailing-space legend-spacing hack) as the
  # volcano figures, but with "autosome" listed FIRST so it is drawn first and
  # the much rarer X/Y points are drawn on top of the dense autosomal cloud.
  chr_color_map <- c(
    "autosome   " = "black",
    "X   " = "cornflowerblue",
    "Y   " = "tomato"
  )

  xlab_string <- "Sex DE log2 fold change"
  ylab_string <- "Age DE log2 fold change per decade"

  plot_width <- 6
  plot_height <- 6
  rasterize_points <- TRUE

  ## -----------------------
  ## Data
  ## -----------------------

  if (is.null(age_df)) {
    age_df <- get_or_build_de_cache(
      test = "age",
      cache_file = file.path(paths$data_cache_dir, "de_age.tsv"),
      de_region_interaction_dir = paths$de_region_interaction_dir,
      ct_file = paths$ct_file,
      gene_to_chr_file = paths$gene_to_chr_file
    )
  }

  if (is.null(sex_df)) {
    sex_df <- get_or_build_de_cache(
      test = "female_vs_male",
      cache_file = file.path(paths$data_cache_dir, "de_sex.tsv"),
      de_region_interaction_dir = paths$de_region_interaction_dir,
      ct_file = paths$ct_file,
      gene_to_chr_file = paths$gene_to_chr_file
    )
  }

  age_df <- data.table::as.data.table(age_df)
  sex_df <- data.table::as.data.table(sex_df)

  if (isTRUE(include_spn_d1_d2)) {
    spn_cell_types <- c("SPN_D1", "SPN_D2")

    # Only fetch what's actually missing, so passing pre-loaded age_df/sex_df
    # that already contain SPN_D1/SPN_D2 doesn't duplicate rows.
    missing_age <- setdiff(spn_cell_types, unique(age_df[["cell_type"]]))
    missing_sex <- setdiff(spn_cell_types, unique(sex_df[["cell_type"]]))

    if (length(missing_age) > 0L || length(missing_sex) > 0L) {
      # SPN_D1/SPN_D2 are not part of the shared cell_types_for_de_filtering_plot.txt
      # list, so they are read from the raw region-interaction results and
      # cached under their own private ct_file/cache files rather than
      # expanding the shared list many other figures depend on.
      extra_ct_file <- file.path(paths$data_cache_dir, "ct_extra_spn_d1_d2.txt")
      if (!file.exists(extra_ct_file)) {
        writeLines(spn_cell_types, extra_ct_file)
      }

      if (length(missing_age) > 0L) {
        extra_age_df <- get_or_build_de_cache(
          test = "age",
          cache_file = file.path(paths$data_cache_dir, "de_age_spn_d1_d2.tsv"),
          de_region_interaction_dir = paths$de_region_interaction_dir,
          ct_file = extra_ct_file,
          gene_to_chr_file = paths$gene_to_chr_file
        )
        age_df <- data.table::rbindlist(
          list(age_df, data.table::as.data.table(extra_age_df)),
          use.names = TRUE
        )
      }

      if (length(missing_sex) > 0L) {
        extra_sex_df <- get_or_build_de_cache(
          test = "female_vs_male",
          cache_file = file.path(paths$data_cache_dir, "de_sex_spn_d1_d2.tsv"),
          de_region_interaction_dir = paths$de_region_interaction_dir,
          ct_file = extra_ct_file,
          gene_to_chr_file = paths$gene_to_chr_file
        )
        sex_df <- data.table::rbindlist(
          list(sex_df, data.table::as.data.table(extra_sex_df)),
          use.names = TRUE
        )
      }
    }
  }

  # Restrict to the region of interest up front: this shrinks each table ~5x
  # and makes the per-cell-type loop cheap.
  age_df <- age_df[region == region_use, ]
  sex_df <- sex_df[region == region_use, ]

  # Drop genes with no chromosome assignment BEFORE gene selection, so each
  # panel plots a full complement of colorable genes. Then collapse to
  # X / Y / autosome exactly as the volcano figures do (modify_chromosome_labels
  # also drops chr == "M").
  n_before <- nrow(age_df) + nrow(sex_df)
  age_df <- age_df[!is.na(chr), ]
  sex_df <- sex_df[!is.na(chr), ]
  logger::log_info(
    "Dropped {n_before - nrow(age_df) - nrow(sex_df)} rows with no chromosome annotation"
  )

  age_df <- modify_chromosome_labels(age_df)
  sex_df <- modify_chromosome_labels(sex_df)
  age_df$chr <- paste0(age_df$chr, "   ")
  sex_df$chr <- paste0(sex_df$chr, "   ")

  if (isTRUE(autosomes_only)) {
    # Drop chrX/chrY BEFORE gene selection, so top_n_expressed is selected from
    # autosomal genes only. All plotted points fall in a single color group,
    # so use a single-entry map and hide the (uninformative) legend below.
    age_df <- age_df[chr == "autosome   ", ]
    sex_df <- sex_df[chr == "autosome   ", ]
    chr_color_map <- c("autosome   " = "black")
  }

  # Derive the cell type list from the data: several cell types are DFC-only by
  # design and have no rows in CaH.
  cell_types_use <- sort(intersect(
    unique(sex_df[["cell_type"]]),
    unique(age_df[["cell_type"]])
  ))

  if (length(cell_types_use) == 0L) {
    stop(
      "No cell types have both sex and age DE results in region '",
      region_use, "'.",
      call. = FALSE
    )
  }

  logger::log_info(
    "Plotting sex vs age DE scatter for {length(cell_types_use)} cell types in region {region_use}"
  )

  ## -----------------------
  ## Plots
  ## -----------------------

  dataset_out_dir <- file.path(paths$outDir, "de_sex_vs_age_scatter")
  .ensure_dir(dataset_out_dir)

  out_files <- character(0)

  for (ct in cell_types_use) {
    logger::log_info("  {ct}")

    p <- bican.mccarroll.de.analysis::plot_de_effect_pair_scatter_gg(
      de_dt_x = sex_df,
      de_dt_y = age_df,
      cell_type_use = ct,
      region_use = region_use,
      test_x = "female_vs_male",
      test_y = "age",
      top_n_expressed = top_n_expressed,
      chr_color_map = chr_color_map,
      show_title = FALSE
    )

    if (is.null(p)) {
      logger::log_warn("  skipping {ct}: too few genes to plot")
      next
    }

    # add_style_volcano applies theme_classic, the enlarged axis/legend text,
    # legend.position = "top" with a single 3-column row, and rasterizes the
    # point layers at 600 dpi.
    p <- add_style_volcano(p, rasterize_points = rasterize_points)

    p <- p +
      ggplot2::labs(
        x = xlab_string,
        y = ylab_string,
        # Display label only. The value used for filtering is never mutated.
        title = gsub("_", " ", ct, fixed = TRUE)
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = ggplot2::rel(2), face = "bold", hjust = 0
        )
      )

    if (isTRUE(autosomes_only)) {
      # Single color group: a legend would just show one black dot.
      p <- p + ggplot2::theme(legend.position = "none")
    }

    fileStr <- paste0(
      "de_sex_vs_age_scatter_", region_use, "_", ct,
      if (isTRUE(autosomes_only)) "_autosomes_only" else "",
      ".svg"
    )

    out_files <- c(
      out_files,
      save_plot_svg(
        p,
        out_file = fileStr,
        out_dir = dataset_out_dir,
        width = plot_width,
        height = plot_height
      )
    )
  }

  logger::log_info("DONE plotting sex vs age DE scatter plots")

  invisible(out_files)
}


plot_de_scatter_svg <- function(df,
                                test,
                                cell_type1,
                                cell_type2,
                                region1,
                                region2,
                                xlab_prefix = NULL,
                                outDir,
                                fdr_cutoff = 0.05,
                                add_fit = TRUE,
                                width = 6,
                                height = 6,
                                rasterize_points = TRUE) {
  fileStr <- paste(
    "de_scatter_plot_",
    test,
    "_", cell_type1, "_", cell_type2,
    "_", region1, "_vs_", region2,
    ".svg",
    sep = ""
  )

  out_file <- file.path(outDir, fileStr)

  format_cell_type <- function(x) {
    x <- gsub("_", " ", x)
    paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  }

  name1 <- format_cell_type(cell_type1)
  name2 <- format_cell_type(cell_type2)

  if (identical(test, "age")) {
    effect_label <- "Age effect fold change per decade"
  } else if (identical(test, "sex")) {
    effect_label <- "Sex effect fold change"
  } else {
    effect_label <- "Effect size, log2"
  }

  xlab_string <- paste(
    effect_label,
    paste0(name1, " [", region1, "]"),
    sep = "\n"
  )

  ylab_string <- paste(
    effect_label,
    paste0(name2, " [", region2, "]"),
    sep = "\n"
  )

  p <- bican.mccarroll.de.analysis::plot_de_scatter_gg(
    df,
    cell_type1,
    cell_type2,
    region1,
    region2,
    fdr_cutoff = fdr_cutoff
  ) +
    ggplot2::labs(
      x = xlab_string,
      y = ylab_string
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = ggplot2::rel(2)),
      axis.title.y = ggplot2::element_text(size = ggplot2::rel(2)),
      axis.text = ggplot2::element_text(size = ggplot2::rel(1.75))
    )

  if (rasterize_points) {
    for (i in seq_along(p$layers)) {
      layer <- p$layers[[i]]
      if (inherits(layer$geom, "GeomPoint")) {
        p$layers[[i]] <- ggrastr::rasterise(layer, dpi = 600)
      }
    }
  }

  ggplot2::ggsave(
    filename = out_file,
    plot = p,
    device = svglite_manuscript,
    width = width,
    height = height,
    units = "in"
  )

  invisible()
}


get_or_build_de_cache <- function(test,
                                  cache_file,
                                  de_region_interaction_dir,
                                  ct_file,
                                  gene_to_chr_file) {
  if (file.exists(cache_file)) {
    logger::log_info("Using cached data from {cache_file}")
    de_ri <- utils::read.table(
      cache_file,
      header = TRUE, sep = "\t",
      stringsAsFactors = FALSE, check.names = FALSE
    )
    data.table::setDT(de_ri)
  } else {
    logger::log_info(
      "No cached data from {cache_file} regenerating data from sources.  This can take a few minutes"
    )
    gene_to_chr <- bican.mccarroll.de.analysis::read_gene_to_chr(
      gene_to_chr_file
    )

    de_ri <- bican.mccarroll.de.analysis::read_de_results(
      de_region_interaction_dir,
      test,
      ct_file,
      gene_to_chr
    )

    utils::write.table(
      de_ri,
      file = cache_file,
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE
    )
  }

  de_ri
}

.resolve_sex_age_scatter_plot_paths <- function(de_results_dir = NULL,
                                                de_region_interaction_dir = NULL,
                                                gene_to_chr_file = NULL,
                                                ct_file = NULL,
                                                outDir = NULL,
                                                data_cache_dir = NULL) {
  root <- .resolve_data_root_dir(NULL)

  rel <- list(
    de_results_dir =
      "differential_expression/results",
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
    de_results_dir            = pick_in(de_results_dir, "de_results_dir"),
    de_region_interaction_dir = pick_in(de_region_interaction_dir, "de_region_interaction_dir"),
    gene_to_chr_file          = pick_in(gene_to_chr_file, "gene_to_chr_file"),
    ct_file                   = pick_in(ct_file, "ct_file"),
    outDir                    = out,
    data_cache_dir            = cache
  )
}
