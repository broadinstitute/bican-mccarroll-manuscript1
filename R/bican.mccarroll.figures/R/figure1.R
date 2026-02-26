## ------------------------------------------------------------------
## Set configuration (development only; comment out in package build)
## ------------------------------------------------------------------
source("R/paths.R")
source ("R/age_pred_figures.R") #for path resolution

options(
    bican.mccarroll.figures.data_root_dir =
        "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis",

    bican.mccarroll.figures.out_dir =
        "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/figure_repository",

    bican.mccarroll.figures.cache_dir =
        "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/figure_repository/data_cache"
)

#' Correlation plot of sample covariates used in the model
#'
#' Compute pairwise correlations among sample-level covariates used for the
#' analysis model and draw a publication-formatted correlation plot. This is
#' intended as a figure-generation wrapper that delegates data preparation and
#' correlation computation to the analysis package, then applies consistent
#' labeling, ordering, and output settings for the manuscript.
#'
#' The correlation matrix is computed from the `samples` table of a prepared
#' `DGEList` using `bican.mccarroll.differentialexpression::getVariableCorrelation()`.
#' Covariate names are mapped to human-readable labels, then reordered into
#' conceptual groups (for example donor, sample, and cell QC) before plotting
#' with `corrplot::corrplot()`. When `outDir` is not `NULL`, the plot is written
#' to an SVG file.
#'
#' @param metacell_dir Directory containing the input objects used by the analysis
#'   package to prepare plotting data.
#' @param data_name Basename of the input object (as expected by the analysis
#'   package) used to load and prepare the `DGEList`.
#' @param randVars Character vector of random-effect covariate names passed to
#'   the analysis data-prep step.
#' @param fixedVars Character vector of fixed-effect covariate names passed to
#'   the analysis data-prep step.
#' @param outDir Output directory for the SVG file. If `NULL`, the plot is drawn
#'   to the active graphics device and no file is written.
#'
#' @return This function is called for its side effects and returns `NULL`
#'   invisibly. The primary outputs are the plotted figure and, optionally, an
#'   SVG file saved to `outDir`.
#'
#' @seealso
#'   \code{\link[corrplot]{corrplot}}
#'   \code{\link[bican.mccarroll.differentialexpression]{prepareMDSPlotData}}
#'   \code{\link[bican.mccarroll.differentialexpression]{getVariableCorrelation}}
#' @export
#'
plot_sample_covariate_correlations<-function (
    metacell_dir = NULL,
    data_name="donor_rxn_DGEList",
    randVars=c("donor", "imputed_sex", "biobank", "single_cell_assay", "region", "hbcac_status", "toxicology_group"),
    fixedVars=c("age", "PC1", "PC2", "PC3", "PC4", "PC5", "pmi_hr", "pct_intronic", "frac_contamination"),
    #randVars=c("donor", "imputed_sex", "biobank", "single_cell_assay", "region", "hbcac_status"),
    #fixedVars=c("age", "PC1", "PC2", "PC3", "PC4", "PC5", "pct_intronic", "frac_contamination"),
    outDir=NULL) {

    paths <- .resolve_age_pred_paths(
        metacell_dir = metacell_dir,
        outDir = outDir
    )

    cache_file <- file.path(paths$data_cache_dir, "plot_sample_covariate_correlations_cache.txt")

    if (file.exists(cache_file)) {
        logger::log_info("Using cached data from {cache_file}")
        df <- utils::read.table(cache_file, header = TRUE, sep = "\t",
                                stringsAsFactors = FALSE)
    } else {
        logger::log_info(
            "No cached data from {cache_file} regenerating data from sources.  This can take a few minutes"
        )

        prep <- bican.mccarroll.differentialexpression::prepareMDSPlotData(
            data_dir = paths$metacell_dir,
            data_name = data_name,
            additionalDonorMetadata = c(),
            randVars = randVars,
            fixedVars = fixedVars
        )
        #get the samples table from the DGEList and write it to the cache file for future use
        df <- prep$dge$samples
        utils::write.table(df, file = cache_file, sep = "\t",
                           row.names = FALSE, quote = FALSE)
    }

    #reassemble the variables to plot.
    required_vars <- unique(c(randVars, fixedVars))
    required_vars<-intersect(required_vars, colnames(df))

    correlation_vars <- setdiff(required_vars, "donor")
    corr_matrix <- bican.mccarroll.differentialexpression::getVariableCorrelation(df, cols_to_test = correlation_vars)

    #make feature names nice!
    pretty_map <- get_pretty_feature_names(correlation_vars)

    colnames(corr_matrix) <- as.vector(pretty_map[colnames(corr_matrix)])
    rownames(corr_matrix) <- as.vector(pretty_map[rownames(corr_matrix)])

    groups <- list(
        donor  = c("Imputed sex", "Biobank", "Age", "Genetic PC1", "Genetic PC2", "Genetic PC3", "Genetic PC4", "Genetic PC5", "Toxicology group", "HBCAC status", "PMI (hours)"),
        sample = c("Region", "Single-cell assay"),
        cell   = c("Percent intronic", "Fraction contamination")
    )

    # groups <- list(
    #     donor  = c("Imputed sex", "Biobank", "Age", "Genetic PC1", "Genetic PC2", "Genetic PC3", "Genetic PC4", "Genetic PC5", "HBCAC status"),
    #     sample = c("Region", "Single-cell assay"),
    #     cell   = c("Percent intronic", "Fraction contamination")
    # )
    setdiff(colnames(corr_matrix), as.vector(unlist (groups)))

    #corrplot draws to the active device, so need to open a device and capture, then close.
    output_svg <- file.path(paths$outDir, "figure1_feature_correlation.svg")
    svglite::svglite(output_svg, width = 8, height = 8)
    on.exit(grDevices::dev.off(), add = TRUE)
    cm <- plot_corrplot_grouped(corr_matrix, groups)
    invisible(cm)

}

plot_corrplot_grouped <- function(corr_matrix,
                                  groups,
                                  method = "circle",
                                  type = "upper",
                                  tl.col = "black",
                                  tl.srt = 45,
                                  na.label = "NA") {
    stopifnot(is.matrix(corr_matrix))
    stopifnot(!is.null(colnames(corr_matrix)), !is.null(rownames(corr_matrix)))
    stopifnot(identical(colnames(corr_matrix), rownames(corr_matrix)))
    stopifnot(is.list(groups), length(groups) > 0)

    # Flatten variables in the desired order (donor -> sample -> cell)
    ord <- unlist(groups, use.names = FALSE)

    # Basic checks
    if (anyDuplicated(ord)) {
        dup <- unique(ord[duplicated(ord)])
        stop("Variables appear in multiple groups: ", paste(dup, collapse = ", "))
    }
    missing <- setdiff(ord, colnames(corr_matrix))
    if (length(missing) > 0) {
        stop("These group variables are not in corr_matrix: ", paste(missing, collapse = ", "))
    }
    extra <- setdiff(colnames(corr_matrix), ord)
    if (length(extra) > 0) {
        stop("corr_matrix contains variables not present in groups: ", paste(extra, collapse = ", "))
    }

    cm <- corr_matrix[ord, ord, drop = FALSE]

    corrplot::corrplot(cm,
                       method = method,
                       type = type,
                       tl.col = tl.col,
                       tl.srt = tl.srt,
                       na.label = na.label
    )

    invisible(cm)
}


# Hard coded pretty names for features.
get_pretty_feature_names <- function(vars) {
    manual_map <- c(
        donor               = "Donor",
        imputed_sex         = "Imputed sex",
        biobank             = "Biobank",
        single_cell_assay   = "Single-cell assay",
        region              = "Region",
        hbcac_status        = "HBCAC status",
        toxicology_group    = "Toxicology group",
        age                 = "Age",
        pmi_hr              = "PMI (hours)",
        pct_intronic        = "Percent intronic",
        frac_contamination  = "Fraction contamination",
        PC1                 = "Genetic PC1",
        PC2                 = "Genetic PC2",
        PC3                 = "Genetic PC3",
        PC4                 = "Genetic PC4",
        PC5                 = "Genetic PC5"
    )

    pretty <- manual_map[vars]

    missing <- is.na(pretty)
    if (any(missing)) {
        auto <- gsub("_", " ", vars[missing])
        auto <- tools::toTitleCase(auto)
        pretty[missing] <- auto
    }

    names(pretty) <- vars
    pretty
}
