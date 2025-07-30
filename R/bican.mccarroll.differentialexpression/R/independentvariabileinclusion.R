#BiocManager::install("Glimma")

# https://github.com/mritchielab/GlimmaV2

# Test which variables should be included in the differential expression model
# examine the effects of all potentital independent variables via
# 1. MDS plot
# 2. Variance Partition
# 3. Correlation plot

# library (vcd)
# library (Glimma)
# library (edgeR)
# library (corrplot)
# library (logger)
# library(ggplot2)


# data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metacells"
# cellTypeGroupFile="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metadata/cell_type_groups.txt"
# cellTypeGroupFile=NULL
#
# data_name="donor_rxn_DGEList"
#
# # Variance Partition variables
# randVars=c("donor", "imputed_sex", "biobank", "single_cell_assay", "region", "hbcac_status", "toxicology_group")
# fixedVars=c("age", "PC1", "PC2", "PC3", "PC4", "PC5", "pmi_hr", "pct_intronic", "frac_contamination")
# max_num_samples=4000
#
# outMDSPlotRoot="/Volumes/nemesh/private_html/BICAN/MDS"
# outPDF= "/Volumes/nemesh/private_html/BICAN/MDS/mds_qc_plots.pdf"


# runMDSPlots(data_dir = data_dir, data_name = data_name, randVars = randVars, fixedVars = fixedVars, max_num_samples = max_num_samples, cellTypeGroupFile = cellTypeGroupFile, outMDSPlotRoot = outMDSPlotRoot, outPDF = outPDF)


#' Run MDS Plots and QC Report for a DGEList
#'
#' Load a precomputed DGEList, validate sample variables, scale PCs,
#' filter samples by library size, generate per–cell‑type and group MDS plots,
#' and compile a QC PDF report including histograms and variable correlations.

#' @param data_dir                Character. Directory containing the precomputed DGEList.
#' @param data_name               Character. Prefix or name used to load the DGEList object.
#' @param randVars               Character vector. Names of “random” metadata variables for MDS coloring.
#' @param fixedVars              Character vector. Names of “fixed” metadata variables for MDS grouping.
#' @param max_num_samples        Integer. Maximum number of samples to include in each MDS plot (default 2500).
#' @param cellTypeGroupFile      Character or NULL. Path to a two‑column TSV/CSV with columns `cell_type` and `group_label`; if NULL, no group plots are made.
#' @param outMDSPlotRoot         Character. Directory in which to save the Glimma MDS HTML files.
#' @param outPDF                 Character. File path for the output QC PDF report.
#' @return
#' Invisibly returns NULL. Side effects:
#' - Writes HTML MDS plots to `outMDSPlotRoot`.
#' - Writes a multi‑page PDF of QC plots to `outPDF`.
#'
#' @importFrom logger log_info
#' @importFrom utils read.table
#' @importFrom grDevices pdf dev.off
#' @importFrom Glimma glimmaMDS
#' @importFrom vcd   assocstats
#' @importFrom corrplot corrplot
#' @importFrom edgeR DGEList
#' @export
runMDSPlots<-function (data_dir, data_name, randVars, fixedVars, max_num_samples=2500, cellTypeGroupFile=NULL, outMDSPlotRoot, outPDF) {
    # load the pre-computed DGEList object
    dge=bican.mccarroll.differentialexpression::loadDGEList(data_dir, prefix = data_name)

    #validate the variables are present in the data set.
    required_vars=c(randVars, fixedVars)
    validateSampleVars(dge, required_vars)

    #scale genetic PCs to unit variance.
    dge$samples=scale_PC_cols(dge$samples)

    # I don't want to have repeated features.
    # Turns out you need norm.factors for logCPM / MDS plot, so can't drop it despite it always being 1.
    dropFeatures=c("sample_name")
    dge$samples <- dge$samples[, !colnames(dge$samples) %in% dropFeatures, drop = FALSE]

    #add a log10 library size column
    dge$samples$lib_log10 <- log10(dge$samples$lib.size)

    # MDS plots by cell type
    cell_type_list=unique(dge$samples$cell_type)
    plotList=list()
    if (length(cell_type_list) > 0) {
        for (cellType in cell_type_list) {
            logger::log_info(paste("Creating MDS plot for cell type:", cellType))
            dge_cell <- dge[, dge$samples$cell_type == cellType, keep.lib.sizes = TRUE]
            r<- filter_by_libsize(dge_cell, threshold_sd = 1.96, bins = 50, strTitlePrefix = cellType)
            dge_cell<- r$dge
            plotList[[cellType]]=r$plot
            mdsPlot(dge_cell, required_vars, num_samples = max_num_samples, outMDSPlotRoot=outMDSPlotRoot, data_name= cellType)
        }
    } else {
        logger::log_info("No cell types found in the DGEList samples.")
    }

    #MDS Plots by cell type groups
    if (!is.null(cellTypeGroupFile)) {
        logger::log_info(paste("Loading cell type groups from:", cellTypeGroupFile))
        cellTypeGroups <- utils::read.table(cellTypeGroupFile, header = TRUE, stringsAsFactors = FALSE)

        # loop over each group
        for (cellTypeGroup in unique(cellTypeGroups$group_label)) {
            logger::log_info(paste("Creating MDS plot for cell type group:", cellTypeGroup))
            cell_type_list= cellTypeGroups$cell_type[cellTypeGroups$group_label == cellTypeGroup]
            idx=which(dge$samples$cell_type %in% cell_type_list)
            dge_cell_group <- dge[, idx, keep.lib.sizes = TRUE]
            r<- filter_by_libsize(dge_cell_group, threshold_sd = 1.96, bins = 50, strTitlePrefix = cellTypeGroup)
            dge_cell_group<- r$dge
            plotList[[cellTypeGroup]]=r$plot
            mdsPlot(dge_cell_group, required_vars, num_samples = max_num_samples, outMDSPlotRoot=outMDSPlotRoot, data_name= cellTypeGroup)
        }
    }

    # MDS plot for all nuclei
    data_name=paste("All Nuclei", max_num_samples, "samples")
    logger::log_info(paste("Creating MDS plot for:", data_name))
    r<- filter_by_libsize(dge, threshold_sd = 1.96, bins = 50, strTitlePrefix="All Nuclei")
    dge_cell<- r$dge
    plotList[[data_name]]=r$plot
    mdsPlot(dge, required_vars, num_samples = max_num_samples, outMDSPlotRoot=outMDSPlotRoot, data_name = data_name)

    #QC report plot
    grDevices::pdf(outPDF)

    #correlation plot
    # exclude donor!
    correlation_vars <- setdiff(required_vars, "donor")
    corr_matrix <- getVariableCorrelation(dge$samples, cols_to_test = correlation_vars)

    #loop over individual cell type library size filtering plots.
    if (length(plotList) > 0) {
        for (cellType in names(plotList)) {
            p <- plotList[[cellType]]
            if (!is.null(p)) {
                print(p)
            } else {
                logger::log_info(paste("No plot available for cell type:", cellType))
            }
        }
    } else {
        log_info("No plots were generated for any cell types.")
    }

    grDevices::dev.off()


}

#' @title Generate and Save Glimma MDS Plots
#' @description
#' Optionally subsample a DGEList to a fixed number of libraries, round metadata,
#' and produce continuous‑ and discrete‑colour MDS plots via Glimma. The HTML widgets
#' are saved to disk with titles indicating the data subset.
#'
#' @param dge            DGEList. An edgeR DGEList object containing counts and sample metadata.
#' @param required_vars  Character vector. Sample‑metadata columns to include in the MDS plot.
#' @param num_samples    Integer or NULL. Maximum number of samples to randomly select; if NULL, all samples are used.
#' @param outMDSPlotRoot Character. Directory in which to save the resulting HTML MDS plots.
#' @param data_name      Character. Prefix for output filenames and titles (e.g., “All Nuclei”).
#'
#' @return
#' Invisibly returns NULL. Side effects:
#' - Saves two HTML Glimma MDS plots (`<data_name>_continuous.html` and `<data_name>_discrete.html`) in `outMDSPlotRoot`.
#' @import htmlwidgets Glimma logger
mdsPlot<-function (dge, required_vars, num_samples=NULL, outMDSPlotRoot, data_name="All nuclei") {
    #optionally limit the number of samples used in analysis - mostly to reduce compute costs for high N.
    if (!is.null(num_samples)) {
        s=min(num_samples, dim(dge)[2])
        idx=sample(1:dim(dge)[2], size=s)
        dgeThis=dge[,idx,keep.lib.sizes=TRUE]
    } else {
        dgeThis=dge
    }

    # format the sample metadata for glimma - this data is only for display purposes.
    dgeThis$samples=round_df_sig(dgeThis$samples, digits = 3L)
    num_samples_retained=dim(dgeThis$samples )[1]
    # Version 2 of the framework.
    logger::log_info(paste("Creating glimma MDS plot - continuous colour with num samples [", num_samples_retained, "]"))
    r=Glimma::glimmaMDS(dgeThis, continuous.colour=TRUE, launch=F, main="TEST V2", width=1600, height=900, var.explained=TRUE, title="All Nuclei")
    logger::log_info(paste("Creating glimma MDS plot - discrete colour with num samples [", num_samples_retained, "]"))
    r2=Glimma::glimmaMDS(dgeThis, continuous.colour=FALSE, launch=F, main="TEST V2", width=1600, height=900, var.explained=TRUE, title="All Nuclei")
    logger::log_info(paste("Finished creating glimma MDS plots"))

    # https://stackoverflow.com/questions/74379298/argument-selfcontained-deprecated-in-htmlwidgetssavewidget
    # work around for selfcontained not working?
    myoriginalwd=getwd()
    setwd(outMDSPlotRoot)
    outFileC= paste0(data_name, "_continuous.html")
    htmlwidgets::saveWidget(r, file=outFileC, selfcontained = TRUE, title=paste(data_name, "continuous"))
    outFileD= paste0(data_name, "_discrete.html")
    htmlwidgets::saveWidget(r2,outFileD, selfcontained = TRUE, title=paste(data_name, "discrete"))
    setwd(myoriginalwd)
}



#' Round numeric columns in a data.frame to a specified number of significant digits
#'
#' @param df A data.frame containing numeric columns to round
#' @param digits Integer; number of significant digits to round to (default 3)
#' @return The original data.frame with numeric columns rounded to the specified significant digits
#' @export
round_df_sig <- function(df, digits = 3L) {
    # sanity check
    if (!is.data.frame(df)) {
        stop("Input must be a data.frame", call. = FALSE)
    }
    # find numeric columns
    num_cols <- vapply(df, is.numeric, logical(1))
    # round each numeric column to 'digits' significant digits
    df[num_cols] <- lapply(
        df[num_cols],
        function(x) signif(x, digits = digits)
    )
    df
}





#' Compute Pairwise Variable Associations
#'
#' Calculates a matrix of pairwise association measures for the specified columns
#' in a data frame. Numeric–numeric pairs use Pearson correlation; factor–factor
#' pairs use Cramér's V
#' @param df             A data.frame containing the variables to test.
#' @param cols_to_test   Character vector of column names in df to include.
#' @return A square numeric matrix (length(cols_to_test) × length(cols_to_test))
#'   with row/column names = cols_to_test, containing:
#'   \itemize{
#'     \item Pearson correlation for numeric–numeric
#'     \item Cramér's V for factor–factor
#'     \item √(SS_between/SS_total) for numeric–factor or factor–numeric
#'   }
#'
#' @details
#' Rows with any missing values in the selected columns are dropped. Character
#' columns are coerced to factors before testing.
#'
#' @importFrom vcd assocstats
#' @import corrplot
#' @export
getVariableCorrelation <- function(df, cols_to_test=c()) {
    # Combine base and candidate predictors
    all_predictors <- cols_to_test
    df_subset <- df[, all_predictors, drop = FALSE]
    df_complete <- df_subset[stats::complete.cases(df_subset), ]

    # Coerce character columns to factors
    df_complete[] <- lapply(df_complete, function(x) if (is.character(x)) as.factor(x) else x)

    # Prepare result matrix
    n <- length(all_predictors)
    corr_matrix <- matrix(NA, nrow = n, ncol = n, dimnames = list(all_predictors, all_predictors))

    # Compute pairwise correlations
    for (i in seq_along(all_predictors)) {
        for (j in seq_along(all_predictors)) {
            var1 <- df_complete[[all_predictors[i]]]
            var2 <- df_complete[[all_predictors[j]]]

            if (is.numeric(var1) && is.numeric(var2)) {
                corr_matrix[i, j] <- stats::cor(var1, var2, use = "pairwise.complete.obs")
            } else if (is.factor(var1) && is.factor(var2)) {
                tbl <- table(var1, var2)
                suppressWarnings({
                    stat <- vcd::assocstats(tbl)$cramer
                })
                corr_matrix[i, j] <- stat
            } else if (is.factor(var1) && is.numeric(var2)) {
                corr_matrix[i, j] <- sqrt(summary(stats::aov(var2 ~ var1))[[1]]$`Sum Sq`[1] /
                                              sum(summary(stats::aov(var2 ~ var1))[[1]]$`Sum Sq`))
            } else if (is.numeric(var1) && is.factor(var2)) {
                corr_matrix[i, j] <- sqrt(summary(stats::aov(var1 ~ var2))[[1]]$`Sum Sq`[1] /
                                              sum(summary(stats::aov(var1 ~ var2))[[1]]$`Sum Sq`))
            }
        }
    }

    # Plot the matrix
    corrplot::corrplot(corr_matrix, method = "circle", type = "upper", tl.col = "black", tl.srt = 45, na.label = "NA")

    return(corr_matrix)
}

#' Validate that required sample variables are present
#'
#' @param dge A DGEList object
#' @param required_vars Character vector of variable names that must exist in dge$samples
#' @return TRUE if all are present, otherwise throws an error listing missing variables
#' @export
validateSampleVars <- function(dge, required_vars) {
    if (is.null(dge$samples) || nrow(dge$samples) == 0) {
        stop("DGEList contains no samples metadata to validate against.")
    }
    missing <- setdiff(required_vars, colnames(dge$samples))
    if (length(missing) > 0) {
        stop(sprintf(
            "Missing sample variables in DGEList: %s",
            paste(missing, collapse = ", ")
        ))
    }
    return(TRUE)
}

#' Scale all PC columns in a data.frame to unit variance
#'
#' @param df A data.frame containing PC columns (named "PC1", "PC2", etc.)
#' @return The original data.frame with PC columns scaled to unit variance
scale_PC_cols <- function(df) {
    # find columns starting with "PC"
    pc.cols <- grep("^PC", names(df))

    # loop over those columns and replace with scaled values
    for (j in pc.cols) {
        # scale() returns a 1‐column matrix, so pull out the vector
        df[[j]] <- as.numeric(scale(df[[j]], center = TRUE, scale = TRUE))
    }
    return (df)
}


#' Filter a DGEList by library size, plotting the distribution and returning filtered DGEList
#'
#' @param dge A DGEList object containing sample metadata with library sizes
#' @param threshold_sd Numeric. Number of standard deviations below the mean to set as the filtering threshold (default 1.96)
#' @param bins Integer. Number of bins for the histogram (default 50)
#' @param strTitlePrefix Character. Prefix for the plot title (default NULL)
#' @return A list containing:
#' - `plot`: A ggplot histogram of log10 library sizes with a vertical line at the threshold
#' - `dge`: The filtered DGEList object retaining only samples with library size above the threshold
#' @export
#' @import ggplot2
#' @import logger
filter_by_libsize <- function(dge, threshold_sd = 1.96, bins = 50, strTitlePrefix=NULL) {
    if (!inherits(dge, "DGEList")) {
        stop("`dge` must be a DGEList object", call. = FALSE)
    }

    # compute log10 library sizes
    lib_log <- log10(dge$samples$lib.size)
    n_total <- length(lib_log)
    m       <- mean(lib_log)
    s       <- stats::sd(lib_log)
    threshold <- m - threshold_sd * s

    # determine which samples to keep
    keep <- lib_log > threshold
    n_retained <- base::sum(keep)

    # build plot title
    title <- base::sprintf(
        "%d samples, retained %d with lib.size > %.2f",
        n_total, n_retained, threshold
    )
    if (!is.null(strTitlePrefix)) {
        title <- paste(strTitlePrefix, title)
    }

    # make the ggplot
    df <- base::data.frame(lib_log = lib_log)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = lib_log)) +
        ggplot2::geom_histogram(bins = bins) +
        ggplot2::geom_vline(xintercept = threshold, colour = "red", linewidth = 1) +
        ggplot2::labs(
            x     = "library size [log10]",
            title = title
        ) +
        ggplot2::theme_minimal()

    # filter the DGEList
    dge_filt <- dge[, keep, keep.lib.sizes = TRUE]

    # log the filtering summary
    logger::log_info(
        "Filtered DGEList to retain {n_retained} of {n_total} samples with library size > {round(threshold, 3)} [log10 UMIs]"
    )

    list(plot = p, dge = dge_filt)
}

# filter_by_libsize <- function(dge, threshold_sd = 1.96, bins = 50) {
#     if (!inherits(dge, "DGEList")) {
#         stop("`dge` must be a DGEList object", call. = FALSE)
#     }
#
#     # 1. compute log10 lib sizes
#     lib_log <- log10(dge$samples$lib.size)
#     n_total <- length(lib_log)
#     m <- mean(lib_log)
#     s <- sd(lib_log)
#     threshold <- m - threshold_sd * s
#
#     # 2. how many survive?
#     keep <- lib_log > threshold
#     n_retained <- sum(keep)
#
#     # 3. build title
#     title <- sprintf(
#         "MDS plot for %d samples, retained %d with lib.size > %.2f",
#         n_total, n_retained, threshold
#     )
#
#     # 4. make the ggplot
#     df <- data.frame(lib_log = lib_log)
#     p <- ggplot(df, aes(x = lib_log)) +
#         geom_histogram(bins = bins) +
#         geom_vline(xintercept = threshold, colour = "red", linewidth = 1) +
#         labs(
#             x     = "library size [log10]",
#             title = title
#         ) +
#         theme_minimal()
#
#     # 5. filter the DGEList
#     dge_filt <- dge[, keep, keep.lib.sizes = TRUE]
#
#     log_info("Filtered DGEList to retain ", dim (dge)[2], " of ", dim (dge)[2]," samples with library size > ", round(threshold,3), " [log 10 UMIs]")
#     # 6. return both
#     list(plot = p, dge = dge_filt)
# }

