# library(bican.mccarroll.differentialexpression)
# library(variancePartition)
# library(Glimma)
# library(ggplot2)
# library(ggrepel)

###################################
# CELL TYPE TESTS MERGED REGIONS
###################################

# Dropping PMI, HBCAC from model. Not all donors have PMI, PMI does not contribute very much to variance explained.
# data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metacells"
# data_name="donor_rxn_DGEList"
# randVars=c("donor", "village")
# #note: "imputed_sex" moves to a fixed effect!
# fixedVars=c("age", "PC1", "PC2", "PC3", "PC4", "PC5", "pct_intronic", "frac_contamination", "imputed_sex", "toxicology_group", "single_cell_assay", "region", "biobank")
# contrast_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metadata/differential_expression_contrasts_all.txt"
# result_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age_toxiciology/cell_type"
# cellTypeListFile=NULL
# outPDF=paste(result_dir, "volcano_plots.pdf", sep="/")
# interaction_var=NULL; absolute_effects=FALSE

# Dropping toxicology from model to include additional donors
# Only running on sex and age to be more donor inclusive.
# data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metacells"
# data_name="donor_rxn_DGEList"
# randVars=c("donor", "village")
# fixedVars=c("age", "PC1", "PC2", "PC3", "PC4", "PC5", "pct_intronic", "frac_contamination", "imputed_sex", "single_cell_assay", "region", "biobank")
# contrast_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metadata/differential_expression_contrasts_sex_age.txt"
# result_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age/cell_type"
# outPDF=paste(result_dir, "volcano_plots.pdf", sep="/")
# cellTypeListFile=NULL


# Example run - no interaction or absolute effects.
# bican.mccarroll.differentialexpression::differential_expression(data_dir, data_name, randVars, fixedVars, contrast_file, interaction_var=NULL, absolute_effects=FALSE, cellTypeListFile, outPDF, result_dir)

###################################
# CELL TYPE PER REGION TESTS - data partitioned by region and fit.
###################################

# data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metacells"
# data_name="donor_rxn_DGEList"
# randVars=c("donor", "village")
# #note: "imputed_sex" moves to a fixed effect!  region will be automatically removed.
# fixedVars=c("age", "PC1", "PC2", "PC3", "PC4", "PC5", "pct_intronic", "frac_contamination", "imputed_sex", "toxicology_group", "single_cell_assay", "region", "biobank")
# contrast_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metadata/differential_expression_contrasts_all.txt"
# result_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age_toxiciology/cell_type_subset_region"
# result_dir="/downloads/differential_expression/sex_age_toxiciology/cell_type_subset_region"
# cellTypeListFile=NULL
# outPDF=paste(result_dir, "volcano_plots.pdf", sep="/")

# Example run
# bican.mccarroll.differentialexpression::differential_expression_region(data_dir, data_name, randVars, fixedVars, contrast_file, cellTypeListFile, outPDF, result_dir)

###################################
# CELL TYPE REGION INTERACTION TESTS
###################################

# data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metacells"
# data_name="donor_rxn_DGEList"
# randVars=c("donor", "village")
# #note: "imputed_sex" moves to a fixed effect!  region will be automatically removed.
# fixedVars=c("age", "PC1", "PC2", "PC3", "PC4", "PC5", "pct_intronic", "frac_contamination", "imputed_sex", "toxicology_group", "single_cell_assay", "region", "biobank")
# contrast_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metadata/differential_expression_contrasts_all.txt"
# result_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/cell_type_region_interaction_CaH_baseline"
# result_dir="/downloads/differential_expression/sex_age_toxiciology/cell_type_region_interaction_CaH_baseline"
# cellTypeListFile=NULL
# outPDF=paste(result_dir, "volcano_plots.pdf", sep="/")
# interaction_var="region" #set to null to not compute interactions.
# absolute_effects = FALSE #set to TRUE to compute absolute effects per region (only when interaction_var is not NULL)

# Example run
# bican.mccarroll.differentialexpression::differential_expression(data_dir, data_name, randVars, fixedVars, contrast_file, interaction_var, absolute_effects, cellTypeListFile, outPDF, result_dir)

###################################
# CELL TYPE REGION INTERACTION WITH ABSOLUTE EFFECTS
###################################

# Toxicology
# data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metacells"
# data_name="donor_rxn_DGEList"
# randVars=c("donor", "village")
# #note: "imputed_sex" moves to a fixed effect!  region will be automatically removed.
# fixedVars=c("age", "PC1", "PC2", "PC3", "PC4", "PC5", "pct_intronic", "frac_contamination", "imputed_sex", "toxicology_group", "single_cell_assay", "region", "biobank")
# contrast_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metadata/differential_expression_contrasts_all.txt"
# result_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/cell_type_region_interaction_absolute_effects"
# result_dir="/downloads/differential_expression/sex_age_toxiciology/cell_type_region_interaction_absolute_effects"
# cellTypeListFile=NULL
# outPDF=paste(result_dir, "volcano_plots.pdf", sep="/")
# interaction_var="region" #set to null to not compute interactions.
# absolute_effects = TRUE #set to TRUE to compute absolute effects per region (only when interaction_var is not NULL)

# No toxicology
# data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metacells"
# data_name="donor_rxn_DGEList"
# randVars=c("donor", "village")
# #note: "imputed_sex" moves to a fixed effect!  region will be automatically removed.
# fixedVars=c("age", "PC1", "PC2", "PC3", "PC4", "PC5", "pct_intronic", "frac_contamination", "imputed_sex", "single_cell_assay", "region", "biobank")
# contrast_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metadata/differential_expression_contrasts_sex_age.txt"
# result_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/cell_type_results_sex_age_region_interaction_absolute_effects"
# result_dir="/downloads/differential_expression/sex_age/cell_type_region_interaction_absolute_effects"
# cellTypeListFile=NULL
# outPDF=paste(result_dir, "volcano_plots.pdf", sep="/")
# interaction_var="region" #set to null to not compute interactions.
# absolute_effects = TRUE #set to TRUE to compute absolute effects per region (only when interaction_var is not NULL)

# Example run with interactions.
# bican.mccarroll.differentialexpression::differential_expression(data_dir, data_name, randVars, fixedVars, contrast_file, interaction_var, absolute_effects, cellTypeListFile, outPDF, result_dir)



#' Run differential expression analysis for each cell type in the DGEList.
#'
#' This uses a means model (~0 + fixedVars) and a random effects model for the specified random variables.
#' For each contrast group, the fixed effects are reordered so the contrast group is first, which
#' makes all levels of the contrast available for comparison.
#' @param data_dir Directory containing the DGEList data.
#' @param data_name Name of the DGEList data file (without extension).
#' @param randVars Vector of random effect variables.
#' @param fixedVars Vector of fixed effect variables.
#' @param contrast_file Path to the file containing contrast definitions.
#' @param interaction_var Optional name of a variable to test for interactions with the contrast variable. If NULL, no interaction terms are added.
#' @param absolute_effects If TRUE, computes absolute effects per level for continuous-by-categorical interactions. Only used if interaction_var is not NULL.
#' @param cellTypeListFile A file containing an explicit list of cell types to test.  If NULL, all cell types in the DGEList will be tested.
#' @param outPDF Optional path to output PDF file for plots.
#' @param result_dir Directory to save the differential expression results.
#' @param n_cores Integer. Number of cores for parallel processing.
#' @export
differential_expression <- function(data_dir, data_name, randVars, fixedVars, contrast_file, interaction_var=NULL, absolute_effects=FALSE, cellTypeListFile=NULL, outPDF=NULL, result_dir, n_cores = parallel::detectCores() - 2) {
    #validate the output directory exists
    if (!dir.exists(result_dir)) {
        logger::log_info(paste("Creating result directory:", result_dir))
        dir.create(result_dir, recursive=TRUE)
    }
    if (!dir.exists(result_dir)) {
        stop("Result directory does not exist: ", result_dir)
    }

    #load the DGEList and prepare the data
    d=bican.mccarroll.differentialexpression::prepare_data_for_differential_expression(data_dir, data_name, randVars, fixedVars)
    dge=d$dge; fixedVars=d$fixedVars; randVars=d$randVars

    dge=filter_dgelist_by_celltype_list(dge, cellTypeListFile)

    contrast_defs <- read.table(contrast_file, stringsAsFactors = FALSE, sep="\t", header=TRUE)

    # Variance Partition by cell type
    cell_type_list=unique(dge$samples$cell_type)
    if (length(cell_type_list) == 0) {
        logger::log_info("No cell types found in the DGEList samples.")
        return(NULL)
    }

    #cellType=cell_type_list[1]
    #cellType="GABA_MGE_DFC"
    lineStr <- strrep("=", 80)

    plot_list= list()

    for (cellType in cell_type_list) {
        logger::log_info(lineStr)
        logger::log_info(paste("Creating differential expression analysis for cell type:", cellType))
        logger::log_info(lineStr)

        dge_cell <- dge[, dge$samples$cell_type == cellType, keep.lib.sizes = TRUE]
        #filtering samples by library size
        r<- filter_by_libsize(dge_cell, threshold_sd = 1.96, bins = 50, strTitlePrefix = cellType)
        dge_cell<- r$dge

        #filter to the top 75% of highly expressed genes as a first pass.
        dge_cell<-filter_top_expressed_genes(dge_cell, gene_filter_frac = 0.75, verbose = TRUE)
        #filter to cpm cutoff of 1.
        r2=plot_logCPM_density_quantiles(dge_cell, cpm_cutoff = 1, logCPM_xlim = c(-5, 15), lower_quantile = 0.05, upper_quantile = 0.95, quantile_steps = 5, min_samples=1, fraction_samples=0.1)
        dge_cell=r2$filtered_dge

        #run differential expression
        #this produces one list per contrast comparison.
        z<-differential_expression_one_cell_type(dge_cell, fixedVars, randVars, contrast_defs,
                                                 interaction_var=interaction_var, absolute_effects=absolute_effects,
                                                 verbose = TRUE, n_cores = n_cores)

        # flatten the results for summary and plotting
        z_flat <- flatten_de_results(z)

        if (length(z_flat) == 0) {
            logger::log_warn("No non-empty DE result tables to write for this cellType/region.")
            next
        }

        #save the results
        for (contrast in names(z_flat)) {
            out=z_flat[[contrast]]
            # replace ":" with "_" in contrast name for filenames
            contrast_name_clean <- gsub(":", "_", contrast, fixed = TRUE)
            n=paste(cellType, contrast_name_clean, sep="_")
            outFile <- file.path(result_dir, paste0(n, "_DE_results.txt"))
            logger::log_info(paste("Saving results to:", outFile))
            write.table(out, file = outFile, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
        }

        #make a volcano plot for each contrast
        for (contrast in names(z_flat)) {
            # replace ":" with "_" in contrast name for filenames
            contrast_name_clean <- gsub(":", "_", contrast, fixed = TRUE)
            n=paste(cellType, contrast, sep="_")
            df <- z_flat[[contrast]]
            if (nrow(df) > 0) {
                p <- make_volcano(df, fdr_thresh = 0.05, lfc_thresh = 0,
                                  top_n_each = 10, title = paste(cellType, contrast_name_clean))
                plot_list[[n]] <- p
            }
        }

    }

    if (!is.null(outPDF)) {
        logger::log_info(paste("Saving all plots to PDF:", outPDF))
        grDevices::pdf(outPDF)
        pages=paginate_plots(plot_list, plots_per_page = 2)
        for (i in 1:length(pages)) {
            print(pages[[i]])
        }
        grDevices::dev.off()
    }

}

#' Run differential expression analysis for each cell type and region in the DGEList.
#'
#' This uses a means model (~0 + fixedVars) and a random effects model for the specified random variables.
#' For each contrast group, the fixed effects are reordered so the contrast group is first, which
#' makes all levels of the contrast available for comparison.
#' @param data_dir Directory containing the DGEList data.
#' @param data_name Name of the DGEList data file (without extension).
#' @param randVars Vector of random effect variables.
#' @param fixedVars Vector of fixed effect variables.
#' @param contrast_file Path to the file containing contrast definitions.
#' @param cellTypeListFile A file containing an explicit list of cell types to test.  If NULL, all cell types in the DGEList will be tested.
#' @param outPDF Optional path to output PDF file for plots.
#' @param result_dir Directory to save the differential expression results.
#' @param n_cores Integer. Number of cores for parallel processing.
#' @export
differential_expression_region <- function(data_dir, data_name, randVars, fixedVars, contrast_file, cellTypeListFile=NULL, outPDF=NULL, result_dir, n_cores = parallel::detectCores() - 2) {
    #load the DGEList and prepare the data
    d=bican.mccarroll.differentialexpression::prepare_data_for_differential_expression(data_dir, data_name, randVars, fixedVars)
    dge=d$dge; fixedVars=d$fixedVars; randVars=d$randVars

    dge=filter_dgelist_by_celltype_list(dge, cellTypeListFile)

    #if region is listed in the fixedVars, remove it.
    if ("region" %in% fixedVars) {
        fixedVars=setdiff(fixedVars, "region")
        logger::log_info("region found in fixedVars, removing it for region-specific differential expression analysis.")
    }

    contrast_defs <- read.table(contrast_file, stringsAsFactors = FALSE, sep="\t", header=TRUE)

    # Variance Partition by cell type
    cell_type_list=unique(dge$samples$cell_type)

    if (length(cell_type_list) == 0) {
        logger::log_info("No cell types found in the DGEList samples.")
        return(NULL)
    }

    #cellType="astrocyte";
    line <- strrep("=", 80)

    plot_list= list()

    for (cellType in cell_type_list) {

        logger::log_info(line)
        logger::log_info(paste("Creating differential expression analysis for cell type:", cellType))
        logger::log_info(line)

        dge_cell <- dge[, dge$samples$cell_type == cellType, keep.lib.sizes = TRUE]

        region_list<-unique(dge_cell$samples$region)

        for (region in region_list) {
            dge_cell_region <- dge_cell[, dge_cell$samples$region == region, keep.lib.sizes = TRUE]
            logger::log_info(paste("  Analyzing region:", region, "with", dim(dge_cell_region$samples)[1], "samples."))

            #filtering samples by library size
            r<- filter_by_libsize(dge_cell_region, threshold_sd = 1.96, bins = 50, strTitlePrefix = cellType)
            dge_cell_region<- r$dge

            #filter to the top 75% of highly expressed genes as a first pass.
            dge_cell_region<-filter_top_expressed_genes(dge_cell_region, gene_filter_frac = 0.75, verbose = TRUE)
            #filter to cpm cutoff of 1.
            r2=plot_logCPM_density_quantiles(dge_cell_region, cpm_cutoff = 1, logCPM_xlim = c(-5, 15), lower_quantile = 0.05, upper_quantile = 0.95, quantile_steps = 5, min_samples=1, fraction_samples=0.1)
            dge_cell_region=r2$filtered_dge

            #run differential expression
            #this produces one list per contrast comparison.
            # no interaction or absolute effects for region-specific tests.
            z<-differential_expression_one_cell_type(dge_cell_region, fixedVars, randVars, contrast_defs,
                                                     interaction_var=NULL, absolute_effects=FALSE,
                                                     verbose = TRUE, n_cores = n_cores)

            # flatten the results for summary and plotting
            # keep only data frames, keep ONLY inner names, preserve order
            z_flat <- flatten_de_results(z)

            if (length(z_flat) == 0) {
                logger::log_warn("No non-empty DE result tables to write for this cellType/region.")
                next
            }

            #save the results
            for (contrast in names(z_flat)) {
                out=z_flat[[contrast]]
                n=paste(cellType, region, contrast, sep="_")
                outFile <- file.path(result_dir, paste0(n, "_DE_results.txt"))
                logger::log_info(paste("Saving results to:", outFile))
                write.table(out, file = outFile, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
            }

            #make a volcano plot for each contrast
            for (contrast in names(z_flat)) {
                n=paste(cellType, region, contrast, sep="_")
                df <- z_flat[[contrast]]
                if (nrow(df) > 0) {
                    p <- make_volcano(df, fdr_thresh = 0.05, lfc_thresh = 0,
                                      top_n_each = 10, title = paste(cellType, contrast))
                    plot_list[[n]] <- p
                }
            }

        }
    }

    if (!is.null(outPDF)) {
        logger::log_info(paste("Saving all plots to PDF:", outPDF))
        grDevices::pdf(outPDF)
        pages=paginate_plots(plot_list, plots_per_page = 2)
        for (i in 1:length(pages)) {
            print(pages[[i]])
        }
        grDevices::dev.off()
    }

}

flatten_de_results <- function(z) {
    z_flat <- do.call(c, lapply(unname(z), function(x) {
        if (!is.list(x)) return(list())
        x[vapply(x, is.data.frame, logical(1))]
    }))

    if (length(z_flat) == 0) return(list())

    keep <- vapply(z_flat, function(df) nrow(df) > 0, logical(1))
    z_flat[keep]
}

########################
# DELEGATION FUNCTION
# Decides which of the 4 modes to use based on contrast_def, interaction_var and absolute_effects.
########################

# Dispatcher: pooled DE modes per contrast_group with validation for interactions
differential_expression_one_cell_type <- function(
        dge_cell,
        fixedVars,
        randVars,
        contrast_defs,
        interaction_var = NULL,          # NULL or factor name
        absolute_effects = FALSE,        # FALSE: use differential_expression_one_cell_type_contrast_group for Mode 1/2
        verbose = TRUE,
        n_cores = parallel::detectCores() - 2
){
    msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))

    stopifnot(is.list(dge_cell), is.data.frame(contrast_defs))
    dge <- dge_cell
    if (!is.null(dge$samples) && exists("sanitize_levels", mode = "function")) {
        dge$samples <- sanitize_levels(dge$samples)
    }

    # Continuous if contrast_defs rows for var have both levels NA; else fall back to sample type
    is_continuous_var <- function(var, df, samples){
        rows <- df[df$variable == var, , drop = FALSE]
        if (nrow(rows)) {
            all(is.na(rows$reference_level) & is.na(rows$comparison_level))
        } else {
            is.numeric(samples[[var]])
        }
    }

    contrast_groups <- unique(as.character(contrast_defs$variable))
    out <- vector("list", length(contrast_groups))
    names(out) <- contrast_groups

    for (cg in contrast_groups) {
        cg_is_cont <- is_continuous_var(cg, contrast_defs, dge$samples)

        msg("DE for contrast_group='%s' (interaction_var=%s, absolute_effects=%s, continuous=%s)",
            cg, ifelse(is.null(interaction_var), "NULL", interaction_var),
            absolute_effects, cg_is_cont)

        if (!absolute_effects) {
            # Validate: if an interaction_var is requested for relative interactions,
            # the tested contrast must be continuous. Otherwise skip.
            if (!is.null(interaction_var) && !cg_is_cont) {
                msg("SKIP: '%s' has categorical contrasts in contrast_defs; relative interactions require continuous.", cg)
                out[[cg]] <- list()
                next
            }

            #Validate: if an interaction_var is requested, the baseline must be present in the data
            #for some data sets they only have a single region, and the interaction is inappropriate.
            if (!is.null(interaction_var)) {
                baseline_level <- contrast_defs[contrast_defs$variable == cg & is.na(contrast_defs$comparison_level), "baseline_region"]
                levs=unique(dge$samples[[interaction_var]])
                if (!(baseline_level %in% unique(dge$samples[[interaction_var]]))) {
                    msg("SKIP: baseline level '%s' for interaction_var '%s' not found in data for contrast_group '%s'.",
                        baseline_level, interaction_var, cg)
                    out[[cg]] <- list()
                    next
                }
                if (length(levs) < 2) {
                    msg("SKIP: interaction_var '%s' has only one level in data for contrast_group '%s'.",
                        interaction_var, cg)
                    out[[cg]] <- list()
                    next
                }
            }

            # Unified path for Mode 1 (no interaction) and Mode 2 (relative interaction vs baseline)
            out[[cg]] <- differential_expression_one_cell_type_contrast_group(
                dge_cell       = dge,
                fixedVars      = fixedVars,
                randVars       = randVars,
                contrast_defs  = contrast_defs,
                contrast_group = cg,
                interaction_var= interaction_var,   # NULL => average effect; factor name => relative interactions
                verbose        = verbose,
                n_cores        = n_cores
            )
            next
        }

        # absolute_effects = TRUE
        if (cg_is_cont) {
            # Mode 3: absolute per-level slopes for continuous covariate
            if (is.null(interaction_var)) stop("interaction_var required for absolute_effects=TRUE with continuous '", cg, "'.")
            out[[cg]] <- continuous_by_factor_differential_expression(
                dge_cell       = dge,
                fixedVars      = fixedVars,
                randVars       = randVars,
                interaction_var= interaction_var,
                continuous_var = cg,
                verbose        = verbose,
                n_cores        = n_cores
            )
        } else {
            # Mode 4: categorical-by-categorical within-stratum contrasts (names from contrast_defs)
            if (is.null(interaction_var)) stop("interaction_var required for absolute_effects=TRUE with categorical '", cg, "'.")
            out[[cg]] <- categorical_by_categorical_differential_expression(
                dge_cell       = dge,
                fixedVars      = fixedVars,
                randVars       = randVars,
                contrast_defs  = contrast_defs,
                factor_var     = cg,
                interaction_var= interaction_var,
                verbose        = verbose,
                n_cores        = n_cores
            )
        }
    }

    out
}

#########################################
# CELL TYPE + REGION ABSOLUTE EFFECS
#########################################

#' Fit region-specific slopes for a continuous variable using limma/voom + dream
#'
#' This function estimates absolute slopes of a continuous covariate
#' (e.g. age) within each level of a categorical factor (e.g. brain region).
#' It constructs explicit per-level slope terms (`continuous_var * I(interaction_var==level)`)
#' in the design matrix, fits the model with `voomWithDreamWeights` and `dream`,
#' and returns one `topTable` per level with statistics for the slope.
#'
#' @param dge_cell A DGEList-like object (counts with \code{samples} metadata).
#' @param fixedVars Character vector of fixed effect variables to include
#'   (do not include the continuous variable or its interaction term).
#' @param randVars Character vector of random effect variables to include
#'   (random intercepts).
#' @param interaction_var Character string. The name of the categorical
#'   factor variable in \code{dge_cell$samples} (e.g. "region").
#' @param continuous_var Character string. The name of the continuous
#'   covariate variable in \code{dge_cell$samples} (e.g. "age").
#' @param verbose Logical. If \code{TRUE}, print progress messages.
#' @param n_cores Integer. Number of cores for parallel processing.
#'
#' @return A named list of \code{data.frame}s (one per level of the factor).
#'   Each element is the result of \code{variancePartition::topTable} for the corresponding
#'   per-level slope, with genes in rows and standard limma statistics.
#'
#' @details
#' This strategy avoids contrast matrices by creating explicit regressors
#' for each level of the factor multiplied by the continuous variable.
#' For example, with \code{continuous_var="age"} and \code{interaction_var="region"}
#' having levels \code{"CaH","NAC"}, the design will include columns
#' \code{age_regionCaH} and \code{age_regionNAC}. Their coefficients are the
#' estimated slopes of age in those regions.
#'
#' @import variancePartition
#' @import BiocParallel
#' @import stats
#' @export
continuous_by_factor_differential_expression <- function(
        dge_cell,
        fixedVars,                          # include interaction_var and covariates; do NOT include continuous_var or its interaction
        randVars,
        interaction_var = "region",         # factor
        continuous_var = "age",             # numeric
        verbose = TRUE,
        n_cores = parallel::detectCores() - 2
){
    msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))

    # data
    dge_this <- dge_cell
    dge_this$samples <- droplevels(dge_this$samples)
    samp <- dge_this$samples

    if (!(interaction_var %in% names(samp))) stop("interaction_var '", interaction_var, "' not found.")
    if (!(continuous_var %in% names(samp))) stop("continuous_var '", continuous_var, "' not found.")
    if (!is.factor(samp[[interaction_var]])) samp[[interaction_var]] <- factor(samp[[interaction_var]])
    if (!is.numeric(samp[[continuous_var]])) stop("continuous_var must be numeric.")

    levs <- levels(samp[[interaction_var]])

    # explicit per-level slope columns
    cont_cols <- paste0(continuous_var, "_", interaction_var, levs)
    for (i in seq_along(levs)) {
        lev <- levs[i]
        col <- cont_cols[i]
        samp[[col]] <- as.numeric(samp[[interaction_var]] == lev) * samp[[continuous_var]]
    }

    # fixed/random effects
    fv <- unique(fixedVars)
    # if there's only one level for interaction_var, drop the baseline term,
    # we'll still evaluate the cont_col for the level that remains - it should be the same as the global result.
    if (length(levs) < 2) {
        message("Only one level for ", interaction_var, "; dropping term [", interaction_var, "] from fixed effects.")
        fv <- setdiff(fv, interaction_var)
    }

    #TODO: is this reasonable here? - get rid of fixed effects with only 1 level.
    fv <- drop_single_level_rand_effects(fv, metadata = samp, verbose = verbose)

    fv <- setdiff(fv, c(continuous_var, paste0(continuous_var, ":", interaction_var)))  # ensure no global cont or interaction
    fv <- unique(c(fv, cont_cols))                                                      # add explicit slope cols
    rv <- prune_random_effects_insufficient_replication(randVars, data = samp)

    # formulas
    rhs_fixed <- paste(fv, collapse = " + ")
    fixed_form <- stats::as.formula(paste("~ 0 +", rhs_fixed))
    rand_part  <- if (length(rv)) paste0("(1|", rv, ")", collapse = " + ") else NULL
    full_form  <- if (!is.null(rand_part)) stats::as.formula(paste("~ 0 +", rhs_fixed, "+", rand_part)) else fixed_form

    # design checks
    X <- stats::model.matrix(fixed_form, data = samp)
    if (qr(X)$rank < ncol(X)) stop("Design not full rank. Check fixed effects.")
    miss <- setdiff(cont_cols, colnames(X))
    if (length(miss)) stop("Missing per-level slope columns in design: ", paste(miss, collapse = ", "))

    # Check for any other issues with the fit, and return early if they are detected.
    chk <- should_skip_dream_subset(fixed_form, samp, min_n = 50)

    if (chk$skip) {
        logger::log_warn(paste("Skipping dream fit:", chk$reason))
        return(list())
    }

    # voom + dream + eBayes
    #param <- BiocParallel::MulticoreParam(workers = n_cores)
    param <- make_bpparam(n_cores=n_cores)
    v1 <- variancePartition::voomWithDreamWeights(dge_this, full_form, data = samp, span=0.3, BPPARAM = param)
    keep <- filter_high_weight_genes(v1, dge_this, quantile_threshold = 0.999)
    dge_this  <- dge_this[keep, ]
    v2   <- variancePartition::voomWithDreamWeights(dge_this, full_form, data = samp, span=0.3, BPPARAM = param, plot = FALSE)

    fit <- capture_dream_warnings({
        variancePartition::dream(v2, full_form, data = samp, BPPARAM = param)
    })
    fit <- variancePartition::eBayes(fit, trend = TRUE, robust = TRUE)

    # outputs: one topTable per level's slope
    nice_names <- paste0(continuous_var, "_", levs)
    tabs <- stats::setNames(
        lapply(seq_along(cont_cols), function(i) variancePartition::topTable(fit, coef = cont_cols[i], number = Inf)),
        nice_names
    )

    tabs
}

#' Categorical x categorical DE via a combined factor (named by contrast_defs)
#'
#' Builds a combined factor \code{combo = interaction(factor_var, interaction_var)} so each
#' coefficient is one cell of the grid. Fits \code{~ 0 + combo + other_fixed + (1|rand...)}
#' with voom+dream, constructs within-\code{interaction_var} contrasts for the pairs listed
#' in \code{contrast_defs} where \code{variable == factor_var}, and returns a named list of
#' \code{topTable} results. Names are \code{<contrast_name>_<region>}, e.g. \code{female_vs_male_CaH}.
#'
#' @param dge_cell DGEList-like with counts and \code{samples}.
#' @param fixedVars Character vector of extra fixed effects to include. Do NOT include
#'   \code{factor_var}, \code{interaction_var}, or their interaction; this function
#'   constructs the combined factor internally.
#' @param randVars Character vector of random-effect grouping variables (random intercepts).
#' @param contrast_defs data.frame with columns \code{contrast_name}, \code{variable},
#'   \code{reference_level}, \code{comparison_level}. Rows for \code{variable == factor_var}
#'   define which pairs to test, and provide human-friendly names.
#' @param factor_var Character scalar. Primary categorical variable (e.g. "imputed_sex").
#' @param interaction_var Character scalar. Stratifying categorical variable (e.g. "region").
#' @param verbose Logical.
#' @param n_cores Integer workers for \code{BiocParallel::MulticoreParam()}.
#'
#' @return Named list of \code{data.frame}s. Each is a \code{variancePartition::topTable} for a contrast
#'   \code{<contrast_name>_<region>} testing mean(comparison, region) - mean(reference, region) = 0.
#'
#' @import variancePartition
#' @import BiocParallel
#' @export
categorical_by_categorical_differential_expression <- function(
        dge_cell,
        fixedVars,
        randVars,
        contrast_defs,
        factor_var      = "imputed_sex",
        interaction_var = "region",
        verbose         = TRUE,
        n_cores         = parallel::detectCores() - 2
){
    msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))

    # --- data prep ---
    dge <- dge_cell
    dge$samples <- droplevels(dge$samples)
    samp <- dge$samples

    if (!(factor_var %in% names(samp))) stop("factor_var '", factor_var, "' not found.")
    if (!(interaction_var %in% names(samp))) stop("interaction_var '", interaction_var, "' not found.")
    if (!is.factor(samp[[factor_var]]))      samp[[factor_var]]      <- factor(samp[[factor_var]])
    if (!is.factor(samp[[interaction_var]])) samp[[interaction_var]] <- factor(samp[[interaction_var]])

    levF <- levels(samp[[factor_var]])
    levR <- levels(samp[[interaction_var]])
    if (length(levF) < 2) stop("Need >= 2 levels for ", factor_var)

    # Pairs and names from contrast_defs
    cd <- contrast_defs[contrast_defs$variable == factor_var, c("contrast_name","reference_level","comparison_level"), drop = FALSE]
    if (!nrow(cd)) stop("contrast_defs has no rows for variable=='", factor_var, "'.")
    if (any(is.na(cd$reference_level) | is.na(cd$comparison_level)))
        stop("contrast_defs rows for ", factor_var, " must have reference_level and comparison_level.")

    # map a token from contrast_defs to an actual level of factor_var
    map_token_to_level <- function(tok, levels_vec){
        tok <- as.character(tok)
        if (tok %in% levels_vec) return(tok)
        sani <- make.names(tok, allow_ = TRUE)
        if (sani %in% levels_vec) return(sani)
        # numeric like "1","2": levels might also be numeric-coded as characters
        if (suppressWarnings(!is.na(as.numeric(tok))) && (tok %in% levels_vec)) return(tok)
        stop("Token '", tok, "' not found among levels: ", paste(levels_vec, collapse = ", "))
    }
    cd$reference_level  <- vapply(cd$reference_level,  map_token_to_level, character(1), levels_vec = levF)
    cd$comparison_level <- vapply(cd$comparison_level, map_token_to_level, character(1), levels_vec = levF)

    # Combined factor with explicit labels
    # we're going to use samp as the design matrix, and are just appending a single column
    # that is the interaction of factor_var and interaction_var so we can
    # construct contrasts within each level of interaction_var.
    # for example, this might contain all combinations of sex and region that are tested as explicit levels.
    samp$combo <- interaction(samp[[factor_var]], samp[[interaction_var]], sep="__", drop=TRUE)
    levC <- levels(samp$combo)

    # --- fixed/random effects ---
    fv <- setdiff(unique(fixedVars), c(factor_var, interaction_var, paste0(factor_var, ":", interaction_var)))
    fv <- c("combo", fv)

    rv <- prune_random_effects_insufficient_replication(randVars, data = samp)

    # --- formulas ---
    rhs_fixed <- paste(fv, collapse = " + ")
    fixed_form <- stats::as.formula(paste("~ 0 +", rhs_fixed))
    rand_part  <- if (length(rv)) paste0("(1|", rv, ")", collapse = " + ") else NULL
    full_form  <- if (!is.null(rand_part)) stats::as.formula(paste("~ 0 +", rhs_fixed, "+", rand_part)) else fixed_form

    # --- design ---
    X <- stats::model.matrix(fixed_form, data = samp)
    if (qr(X)$rank < ncol(X)) stop("Design not full rank.")

    # Check for any other issues with the fit, and return early if they are detected.
    chk <- should_skip_dream_subset(fixed_form, samp, min_n = 50)

    if (chk$skip) {
        logger::log_warn(paste("Skipping dream fit:", chk$reason))
        return(list())
    }

    # Map "A__R" labels to the actual design columns created by model.matrix
    combocols <- grep("^combo", colnames(X), value = TRUE)
    if (!length(combocols)) stop("No 'combo' columns found in design.")
    find_combo_col <- function(level_label) {
        cand <- c(
            paste0("combo", level_label),
            paste0("combo", make.names(level_label, allow_ = TRUE)),
            paste0("combo", make.names(level_label, allow_ = FALSE))
        )
        hit <- cand[cand %in% combocols]
        if (length(hit)) return(hit[1])
        # fallback on suffix matching
        suff <- sub("^combo", "", combocols)
        m1 <- make.names(level_label, allow_ = TRUE)
        m2 <- make.names(level_label, allow_ = FALSE)
        if (level_label %in% suff) return(paste0("combo", level_label))
        if (m1 %in% suff)         return(paste0("combo", m1))
        if (m2 %in% suff)         return(paste0("combo", m2))
        NA_character_
    }
    combo_map <- stats::setNames(vapply(levC, find_combo_col, character(1)), levC)
    if (anyNA(combo_map)) {
        bad <- names(combo_map)[is.na(combo_map)]
        stop("Could not map combo level(s) to design columns: ", paste(bad, collapse = ", "))
    }

    # --- build contrast matrix L (all specified pairs, in every region) ---
    # columns named "<contrast_name>_<region>"
    pairs <- do.call(rbind, lapply(levR, function(r) {
        data.frame(region = r,
                   contrast_name    = cd$contrast_name,
                   reference_level  = cd$reference_level,
                   comparison_level = cd$comparison_level,
                   stringsAsFactors = FALSE)
    }))

    coef_names <- colnames(X)
    L <- matrix(0, nrow = length(coef_names), ncol = nrow(pairs),
                dimnames = list(coef_names, paste0(pairs$contrast_name, "_", pairs$region)))

    for (i in seq_len(nrow(pairs))) {
        A <- pairs$reference_level[i]
        B <- pairs$comparison_level[i]
        R <- pairs$region[i]
        lvlB <- paste0(B, "__", R)
        lvlA <- paste0(A, "__", R)
        colB <- combo_map[[lvlB]]
        colA <- combo_map[[lvlA]]
        if (is.na(colB) || is.na(colA)) stop("Missing combo cell(s) for ", B, " or ", A, " at ", R)

        L[colB, i] <-  1
        L[colA, i] <- -1
    }

    # --- voom + dream with L ---
    param <- make_bpparam(n_cores=n_cores)
    v1 <- variancePartition::voomWithDreamWeights(dge, full_form, data = samp, span=0.3, BPPARAM = param)
    keep <- filter_high_weight_genes(v1, dge, quantile_threshold = 0.999)
    dge2 <- dge[keep, ]
    v2 <- variancePartition::voomWithDreamWeights(dge2, full_form, data = samp, span=0.3, BPPARAM = param, plot = FALSE)

    fit <- capture_dream_warnings({
        variancePartition::dream(v2, full_form, data = samp, BPPARAM = param, L = L)
    })
    fit <- variancePartition::eBayes(fit, trend = TRUE, robust = TRUE)

    # --- results ---
    tabs <- stats::setNames(
        lapply(colnames(L), function(nm) variancePartition::topTable(fit, coef = nm, number = Inf)),
        colnames(L)
    )

    tabs
}


#####################################
# CELL TYPE + Optional INTERACTION
#####################################

#' Differential expression for one cell type with optional continuous-by-factor interactions
#'
#' Fit a `variancePartition::dream` model for one cell type, optionally adding an
#' interaction between a **continuous** contrast variable (e.g., `age`) and a
#' multi-level factor (default `region`). If an interaction is requested and a
#' baseline level is provided in `contrast_defs`, the function:
#' (1) enforces treatment coding for the interaction factor with that baseline,
#' (2) adds the interaction term to the fixed effects,
#' (3) returns topTables for the main continuous effect (renamed to
#' `contrast_group:interaction_var<baseline>`), the per-level interaction
#' differences (`contrast_group:interaction_var<level>`), and an additional
#' data frame of absolute per-level effects (`contrast_group_absolute_effects`).
#'
#' @note In this version, interaction terms are supported **only for continuous
#'   variables**. For categorical contrast variables (e.g., multi-level treatment),
#'   interaction testing is not constructed here.
#'
#' @param dge_cell A `DGEList` with counts and sample metadata in
#'   `dge_cell$samples`.
#' @param fixedVars `character()`. Fixed effects to include. The function will
#'   move `contrast_group` to the front and may append `contrast_group:interaction_var`.
#' @param randVars `character()`. Random-effect grouping variables. Terms with
#'   insufficient replication are pruned.
#' @param contrast_defs `data.frame` describing contrasts. Must contain columns:
#'   `contrast_name`, `variable`, `reference_level`, `comparison_level`, and
#'   `baseline_region`. The `baseline_region` entry for `variable == contrast_group`
#'   is used as the interaction baseline when `interaction_var` is not `NULL`.
#' @param contrast_group `character(1)`. contrast_name to test (e.g., `"age"`).
#' @param interaction_var `character(1)` or `NULL`. Factor to interact with
#'   `contrast_group` (default `"region"`). If `NULL`, no interaction is added.
#'   If non-`NULL`, a baseline level must be supplied via
#'   `contrast_defs$baseline_region` where `variable == contrast_group`.
#' @param verbose `logical(1)`. Verbose logging.
#' @param n_cores `integer(1)`. Worker count for `BiocParallel::MulticoreParam`.
#'
#' @importFrom variancePartition voomWithDreamWeights dream eBayes topTable
#' @importFrom BiocParallel MulticoreParam
#' @importFrom stats model.matrix as.formula relevel contr.treatment
differential_expression_one_cell_type_contrast_group <- function(
        dge_cell, fixedVars, randVars, contrast_defs,
        contrast_group = "age",
        interaction_var = "region",          # set NULL to disable interactions
        verbose = TRUE,
        n_cores = parallel::detectCores() - 2
){
    # ---- helpers -----------------------------------------------------------
    .get_baseline <- function(df, var) {
        x <- unique(na.omit(df$baseline_region[df$variable == var]))
        if (length(x) == 1) x else NULL
    }

    .ensure_treatment_coding <- function(df, var, baseline) {
        # Make sure the column is a factor
        if (!is.factor(df[[var]])) df[[var]] <- factor(df[[var]])

        # Relevel so that 'baseline' is the reference level
        df[[var]] <- stats::relevel(df[[var]], ref = baseline)

        # Apply treatment coding contrasts with that baseline
        contrasts(df[[var]]) <- stats::contr.treatment(
            nlevels(df[[var]]),
            base = which(levels(df[[var]]) == baseline)
        )
        df
    }

    #
    .pick_int_name <- function(cols, a, b, lev) {
        x <- paste0(a, ":", b, lev)
        y <- paste0(b, lev, ":", a)
        if (x %in% cols) x else if (y %in% cols) y else NA_character_
    }

    # ---- data --------------------------------------------------------------
    dge_cell_this <- dge_cell
    dge_cell_this$samples <- droplevels(dge_cell_this$samples)

    # random effects
    rv <- prune_random_effects_insufficient_replication(randVars, data = dge_cell_this$samples)

    # fixed effects
    fv <- fixedVars
    fv <- move_to_front(fv, contrast_group)
    fv <- drop_single_level_rand_effects(fv, metadata = dge_cell_this$samples, verbose = verbose)

    # interaction toggle + baseline
    baseline <- NULL
    add_interaction <- FALSE
    if (!is.null(interaction_var)) {
        baseline <- .get_baseline(contrast_defs, contrast_group)
        add_interaction <- !is.null(baseline)
    }

    # if requested, enforce treatment coding for the interaction factor on DATA
    if (add_interaction) {
        stopifnot(interaction_var %in% names(dge_cell_this$samples))
        dge_cell_this$samples <- .ensure_treatment_coding(dge_cell_this$samples, interaction_var, baseline)
        inter_term <- paste0(contrast_group, ":", interaction_var)
        if (!(inter_term %in% fv)) fv <- c(fv, inter_term)
    }

    # ---- formulas ----------------------------------------------------------
    rand_part  <- if (length(rv)) paste0("(1|", rv, ")", collapse = " + ") else NULL
    fixed_part <- paste(fv, collapse = " + ")
    fixed_form <- stats::as.formula(paste("~ 0 +", fixed_part))
    full_form  <- stats::as.formula(paste("~ 0 +", paste(c(fixed_part, rand_part), collapse = " + ")))

    design <- stats::model.matrix(fixed_form, data = dge_cell_this$samples)
    if (qr(design)$rank < ncol(design)) stop("Design matrix not full rank.")

    # contrasts for factor main-effects of the contrast_group (if any)
    contrast_defs_this <- contrast_defs[contrast_defs$variable == contrast_group, , drop = FALSE]
    contrast_defs_this <- sanitize_contrast_levels(contrast_defs_this, design, verbose = verbose)
    contrast_matrix    <- generate_contrasts_from_defs(contrast_defs_this, design)
    has_contrasts_groups <- !all(is.na(contrast_defs_this$reference_level) & is.na(contrast_defs_this$comparison_level))
    L <- if (has_contrasts_groups) contrast_matrix else NULL

    # Check for any other issues with the fit, and return early if they are detected.
    chk <- should_skip_dream_subset(fixed_form, dge_cell_this$samples, min_n = 50)

    if (chk$skip) {
        logger::log_warn(paste("Skipping dream fit:", chk$reason))
        return(list())
    }

    # ---- fit ---------------------------------------------------------------
    #param <- BiocParallel::MulticoreParam(workers = n_cores)
    param <- make_bpparam(n_cores=n_cores)
    vobj <- variancePartition::voomWithDreamWeights(
        counts = dge_cell_this, formula = full_form, data = dge_cell_this$samples, span=0.3, BPPARAM = param
    )
    keep <- filter_high_weight_genes(vobj, dge_cell_this, quantile_threshold = 0.999)
    dge_cell_this <- dge_cell_this[keep, ]
    vobj <- variancePartition::voomWithDreamWeights(
        dge_cell_this, full_form, data = dge_cell_this$samples, span=0.3, BPPARAM = param, plot = FALSE
    )

    #keep the pre ebayes fit for absolute effects
    fit <- capture_dream_warnings({
        variancePartition::dream(exprObj = vobj, formula = full_form,
                                 data = dge_cell_this$samples, BPPARAM = param, L = L)
    })
    fitmm <- variancePartition::eBayes(fit, trend = TRUE, robust = TRUE)

    log_decide_tests_summary(fitmm, L = L, label = paste("DREAM DE summary for", contrast_group))

    # ---- collect results ---------------------------------------------------
    tt <- list()

    # 1) factor contrasts (if any)
    if (!is.null(L)) {
        have <- intersect(colnames(L), colnames(coef(fitmm)))
        for (cn in have) tt[[cn]] <- variancePartition::topTable(fitmm, coef = cn, number = Inf)
    }

    coef_names <- colnames(coef(fitmm))

    # 2) main continuous effect always, but rename if interaction is active
    if (contrast_group %in% coef_names) {
        main_tbl <- variancePartition::topTable(fitmm, coef = contrast_group, number = Inf)
        main_name <- contrast_group
        if (add_interaction) {
            main_name <- paste0(contrast_group, ":", interaction_var, baseline)
        }
        tt[[main_name]] <- main_tbl
    }

    # 3) interaction terms, if requested
    if (add_interaction) {
        levs <- levels(dge_cell_this$samples[[interaction_var]])
        # for each level including baseline, create a name contrast_group:interaction_var<lev>
        # baseline uses the renamed main effect; others use explicit interaction coefs
        for (lev in levs) {
            nm <- paste0(contrast_group, ":", interaction_var, lev)
            if (lev == baseline) {
                # already stored as main_name
                next
            } else {
                ic <- .pick_int_name(coef_names, contrast_group, paste0(interaction_var), lev)
                if (is.na(ic)) next
                tt[[nm]] <- variancePartition::topTable(fitmm, coef = ic, number = Inf)
            }
        }
    }

    tt
}

generate_contrasts_from_defs <- function(contrast_defs, design_matrix) {
    # escape any regex metacharacters (incl. hyphen)
    .rex_escape <- function(x) gsub("([][{}()+*^$|\\.?<>\\-])", "\\\\\\1", x)

    stopifnot(is.data.frame(contrast_defs))
    stopifnot(is.matrix(design_matrix) || is.data.frame(design_matrix))

    # Drop interaction columns entirely
    design_cols_raw  <- colnames(design_matrix)
    design_cols_raw  <- design_cols_raw[!grepl(":", design_cols_raw, fixed = TRUE)]

    # Safe names for makeContrasts
    design_cols_safe <- make.names(design_cols_raw, unique = TRUE)
    raw2safe <- stats::setNames(design_cols_safe, design_cols_raw)
    safe2raw <- stats::setNames(design_cols_raw,  design_cols_safe)

    # --- translators work only on NON-interaction columns ---
    translate_side <- function(expr, var, design_cols) {
        if (is.na(expr) || is.null(expr) || nchar(trimws(expr)) == 0) return("0")
        s <- gsub("\\s+", "", as.character(expr))

        var_pat  <- paste0("^", .rex_escape(var))
        var_cols <- grep(var_pat, design_cols, value = TRUE)
        if (length(var_cols) == 0)
            stop("No design columns found for factor '", var, "'. Did you use '~ 0 + ", var, " + ...'?")

        levels_available <- sub(var_pat, "", var_cols)

        # remap numeric tokens like "1" -> "X1" if present
        m <- gregexpr("[A-Za-z0-9_.-]+", s, perl = TRUE)
        toks <- regmatches(s, m)[[1]]
        if (length(toks)) {
            mapped <- vapply(toks, function(tok) {
                if (grepl("^[0-9]+(\\.[0-9]+)?$", tok)) {
                    sani <- make.names(tok)
                    if (sani %in% levels_available) sani else tok
                } else tok
            }, character(1))
            regmatches(s, m)[[1]] <- mapped
        }

        # replace level tokens with full column names (var + level)
        levels_available <- levels_available[order(nchar(levels_available), decreasing = TRUE)]
        for (lev in levels_available) {
            s <- gsub(paste0("(?<![A-Za-z0-9_.])", .rex_escape(lev), "(?![A-Za-z0-9_.])"),
                      paste0(var, lev), s, perl = TRUE)
        }
        s
    }

    contrast_list <- list()
    for (i in seq_len(nrow(contrast_defs))) {
        row <- contrast_defs[i, ]
        cname <- as.character(row$contrast_name)
        var   <- as.character(row$variable)
        ref   <- row$reference_level
        comp  <- row$comparison_level

        # Continuous: expect a single column named exactly <var> (no interactions)
        if ((is.na(ref) || length(ref) == 0) && (is.na(comp) || length(comp) == 0)) {
            if (!(var %in% design_cols_raw)) stop("Continuous term '", var, "' not found in design.")
            contrast_list[[cname]] <- var
            next
        }

        comp_str <- translate_side(comp, var, design_cols_raw)
        ref_str  <- translate_side(ref,  var, design_cols_raw)

        contrast_list[[cname]] <-
            if (identical(ref_str, "0")) comp_str else
                if (identical(comp_str, "0")) paste0("0 - (", ref_str, ")") else
                    paste0("(", comp_str, ") - (", ref_str, ")")
    }

    # safeify expressions for limma
    safeify_expr <- function(expr) {
        s <- expr
        raws <- names(raw2safe)[order(nchar(names(raw2safe)), decreasing = TRUE)]
        for (r in raws) {
            pat <- paste0("(?<![A-Za-z0-9_.])", .rex_escape(r), "(?![A-Za-z0-9_.])")
            s <- gsub(pat, raw2safe[[r]], s, perl = TRUE)
        }
        s
    }
    contrast_list_safe <- lapply(contrast_list, safeify_expr)

    CM_safe <- do.call(limma::makeContrasts,
                       args = c(contrast_list_safe, list(levels = design_cols_safe)))

    #TODO: I'd like to switch to variancePartition::makeContrastsDream for consistency, but I would need to pass in data.
    # CM_safe <- do.call(variancePartition::makeContrastsDream,
    #                     args = c(contrast_list_safe, list(levels = design_cols_safe)))

    # map rows back so contrasts.fit aligns with the fit
    rownames(CM_safe) <- unname(safe2raw[rownames(CM_safe)])
    CM_safe
}




sanitize_levels <- function(df, exclude = character()) {
    for (col in setdiff(names(df), exclude)) {
        if (is.factor(df[[col]])) {
            levels(df[[col]]) <- make.names(levels(df[[col]]))
        } else if (is.character(df[[col]])) {
            df[[col]] <- make.names(df[[col]])
        }
    }
    df
}



# Simple, design-aware sanitizer:
# - Maps tokens to actual factor levels present in the design (~ 0 + variable).
# - If any unknown token appears in an expression, that side becomes NA.
# - After sanitizing both sides: drop rows where exactly one side is NA.
# - Keep rows where both sides are NA (continuous) or both valid (factor contrast).
sanitize_contrast_levels <- function(contrast_defs, design_matrix, verbose = TRUE) {
    stopifnot(is.data.frame(contrast_defs))
    stopifnot(is.matrix(design_matrix) || is.data.frame(design_matrix))

    design_cols <- colnames(design_matrix)

    sanitize_expr_for_var <- function(expr, var) {
        # NA means "no expression" (allowed for continuous rows or 0-side)
        if (length(expr) == 0 || is.na(expr)) return(NA_character_)
        s <- as.character(expr)

        # pull suffixes for this var from design (means coding: ~ 0 + var)
        var_pat  <- paste0("^", var)
        var_cols <- grep(var_pat, design_cols, value = TRUE)
        if (!length(var_cols)) return(NA_character_)  # no columns -> treat as invalid for this var
        lvl_suffixes <- sub(var_pat, "", var_cols)

        # tokenize: contiguous level-like chunks; leave operators/parens
        m <- gregexpr("[A-Za-z0-9_.-]+", s, perl = TRUE)
        toks <- regmatches(s, m)[[1]]
        if (!length(toks)) return(NA_character_)  # nothing meaningful -> invalid

        had_unknown <- FALSE
        mapped <- vapply(toks, function(tok) {
            # numeric literal? could be a level like "1"/"2" *or* a true number (e.g., "/2")
            if (grepl("^[0-9]+(\\.[0-9]+)?$", tok)) {
                # prefer direct suffix match, else try make.names("1")->"X1"
                if (tok %in% lvl_suffixes) return(tok)
                tok_sani <- make.names(tok)
                if (tok_sani %in% lvl_suffixes) return(tok_sani)
                return(tok)  # true numeric constant
            }
            # text token already a suffix?
            if (tok %in% lvl_suffixes) return(tok)
            # try sanitized text
            tok_sani <- make.names(tok)
            if (tok_sani %in% lvl_suffixes) return(tok_sani)

            had_unknown <<- TRUE
            tok  # keep to reinsert (we'll null out the whole side below)
        }, character(1))

        if (had_unknown) return(NA_character_)  # this side invalid -> caller will drop the row

        # put mapped tokens back (still suffixes, not full var+suffix)
        regmatches(s, m)[[1]] <- mapped
        s
    }

    out <- contrast_defs
    # force character columns and sanitize per-row
    out$reference_level  <- NA_character_
    out$comparison_level <- NA_character_
    for (i in seq_len(nrow(out))) {
        var <- as.character(contrast_defs$variable[i])
        out$reference_level[i]  <- sanitize_expr_for_var(contrast_defs$reference_level[i],  var)
        out$comparison_level[i] <- sanitize_expr_for_var(contrast_defs$comparison_level[i], var)
    }

    # Decide which rows to keep:
    # - keep if both sides NA (continuous)
    # - keep if both sides non-NA (valid factor contrast)
    # - drop if exactly one side is NA
    both_na     <- is.na(out$reference_level) & is.na(out$comparison_level)
    both_non_na <- !is.na(out$reference_level) & !is.na(out$comparison_level)
    keep <- both_na | both_non_na

    dropped <- out$contrast_name[!keep]
    if (verbose && length(dropped)) {
        message(sprintf("Dropping %d contrast(s) due to missing/unknown levels: %s",
                        length(dropped), paste(dropped, collapse = ", ")))
    }

    out[keep, , drop = FALSE]
}


move_to_front <- function(vec, name_to_move) {
    if (!(name_to_move %in% vec)) {
        stop(sprintf("'%s' not found in vector", name_to_move))
    }
    c(name_to_move, setdiff(vec, name_to_move))
}

log_decide_tests_summary <- function(fit, L, label = "DE summary") {
    contrast_names <- colnames(L)
    out <- summary(limma::decideTests(fit)[, contrast_names, drop = FALSE])
    logger::log_info("{label}:\n{paste(capture.output(print(out)), collapse = '\n')}")
}


# Filter genes with extreme voom weights
#TODO: remove dge argument, no longer needed.
filter_high_weight_genes <- function(vobj, dge, quantile_threshold = 0.999, max_threshold=1e10, verbose = TRUE) {
    stopifnot(!is.null(vobj), !is.null(vobj$E), !is.null(vobj$weights))

    tested_genes <- rownames(vobj$E)
    weights <- vobj$weights

    # If weights have rownames, enforce ordering to tested_genes
    if (!is.null(rownames(weights))) {
        weights <- weights[tested_genes, , drop = FALSE]
    } else {
        if (nrow(weights) != length(tested_genes)) {
            stop("weights rows do not match tested genes and have no rownames to align.")
        }
    }

    # row-wise max, tolerate all-NA rows
    max_w <- apply(weights, 1, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
    finite_idx <- is.finite(max_w)
    if (!any(finite_idx)) stop("No finite weights to compute threshold.")

    thr <- stats::quantile(max_w[finite_idx], quantile_threshold, na.rm = TRUE)
    #if the threshold is still above some very high threshold, reduce it.
    if(thr > max_threshold)
        thr <- max_threshold
    keep_idx <- finite_idx & (max_w < thr)

    if (verbose) {
        message(sprintf(
            "Filtering %d of %d genes (%.2f%%) with extreme weights (quantile %.3f -> %.2e)",
            sum(!keep_idx), length(keep_idx), 100 * sum(!keep_idx) / length(keep_idx),
            quantile_threshold, thr
        ))
    }

    tested_genes[keep_idx]
}


#' Decide whether to skip a dream fit for a small or unstable subset
#'
#' Uses the fixed-effects design matrix to decide whether to skip a fit because
#' there are too few samples or the design is numerically unstable.
#'
#' Assumes upstream code has already enforced complete cases and a full-rank
#' design for the same \code{fixed_form} and \code{samp}.
#'
#' @param fixed_form Fixed-effects formula used to build the design matrix.
#' @param samp Sample metadata (\code{data.frame}), typically \code{dge$samples}.
#' @param min_n Minimum number of samples required.
#' @param min_df_resid Minimum required residual degrees of freedom (\code{n - p}).
#' @param max_kappa Maximum allowed condition number for the design matrix.
#'
#' @return A list with \code{skip}, \code{reason}, \code{n}, \code{p},
#'   \code{df_resid}, and \code{kappa}.
#'
should_skip_dream_subset <- function(fixed_form, samp,
                                           min_n = 25,
                                           min_df_resid = 10,
                                           max_kappa = 1e10) {

    n <- nrow(samp)
    if (n < min_n) {
        return(list(skip = TRUE,
                    reason = sprintf("Too few samples: n=%d < %d", n, min_n),
                    n = n, p = NA_integer_, df_resid = NA_integer_, kappa = NA_real_))
    }

    X <- stats::model.matrix(fixed_form, data = samp)
    p <- ncol(X)

    df_resid <- n - p
    if (df_resid < min_df_resid) {
        return(list(skip = TRUE,
                    reason = sprintf("Too few residual degrees of freedom: df_resid=%d < %d", df_resid, min_df_resid),
                    n = n, p = p, df_resid = df_resid, kappa = NA_real_))
    }

    sv <- base::svd(X, nu = 0, nv = 0)$d
    kappa_est <- max(sv) / min(sv)

    if (!is.finite(kappa_est) || kappa_est > max_kappa) {
        return(list(skip = TRUE,
                    reason = sprintf("Ill-conditioned design: kappa~%.2e > %.2e", kappa_est, max_kappa),
                    n = n, p = p, df_resid = df_resid, kappa = kappa_est))
    }

    list(skip = FALSE,
         reason = "",
         n = n,
         p = p,
         df_resid = df_resid,
         kappa = kappa_est)
}



capture_dream_warnings <- function(expr) {
    warning_msgs <- character()

    result <- withCallingHandlers(
        expr,
        warning = function(w) {
            msg <- conditionMessage(w)
            if (grepl("Model failed to converge with .*negative eigenvalue", msg)) {
                warning_msgs <<- c(warning_msgs, msg)
                invokeRestart("muffleWarning")  # Suppress the warning
            }
        }
    )

    failed_count <- length(warning_msgs)
    if (failed_count > 0) {
        message(sprintf("%d genes failed to converge", failed_count))
    }

    return(result)
}

summarize_top_tables_for_celltype <- function(topTables_list,
                                              lfc_threshold = 0.5,
                                              pval_threshold = 0.05) {
    # Initialize result holder
    summary_list <- list()

    for (contrast_group in names(topTables_list)) {
        group_results <- topTables_list[[contrast_group]]
        for (contrast_name in names(group_results)) {
            df <- group_results[[contrast_name]]

            # Define status per gene
            status <- rep("NotSig", nrow(df))
            status[df$logFC >= lfc_threshold & df$adj.P.Val < pval_threshold] <- "Up"
            status[df$logFC <= -lfc_threshold & df$adj.P.Val < pval_threshold] <- "Down"

            # Tabulate and align
            tab <- table(factor(status, levels = c("Down", "NotSig", "Up")))
            summary_list[[contrast_name]] <- tab
        }
    }

    # Combine to a matrix
    summary_matrix <- do.call(cbind, summary_list)
    return(summary_matrix)
}



make_volcano <- function(df,
                             fdr_thresh = 0.05,
                             lfc_thresh = 0,
                             top_n_each = 10,
                             title = NULL,
                             point_alpha = 0.6,
                             add_counts_inset = TRUE) {

    if (!all(c("logFC", "P.Value", "adj.P.Val") %in% names(df)))
        stop("df must contain: logFC, P.Value, adj.P.Val")
    if (is.null(rownames(df)))
        stop("Row names must contain gene symbols for labeling.")

    d <- df
    d$neglog10FDR <- -log10(pmax(d$adj.P.Val, .Machine$double.xmin))
    d$sig <- d$adj.P.Val <= fdr_thresh & abs(d$logFC) >= lfc_thresh

    d$dir <- ifelse(d$sig & d$logFC < 0, "down",
                    ifelse(d$sig & d$logFC > 0, "up", "ns"))
    # Fix display order to match your legend preference
    d$dir <- factor(d$dir, levels = c("down", "ns", "up"))

    # labels: top |logFC| among FDR-significant hits, split by direction
    sig_up   <- d[d$dir == "up",   , drop = FALSE]
    sig_down <- d[d$dir == "down", , drop = FALSE]
    ord_up   <- if (nrow(sig_up))   order(-abs(sig_up$logFC), sig_up$adj.P.Val)   else integer(0)
    ord_down <- if (nrow(sig_down)) order(-abs(sig_down$logFC), sig_down$adj.P.Val) else integer(0)
    lab_up   <- if (length(ord_up))   utils::head(sig_up[ord_up,   , drop = FALSE], top_n_each) else d[0,]
    lab_down <- if (length(ord_down)) utils::head(sig_down[ord_down, , drop = FALSE], top_n_each) else d[0,]
    labs <- rbind(lab_up, lab_down)
    if (nrow(labs)) labs$label <- rownames(labs)

    # counts (ordered down/ns/up)
    n_down <- sum(d$dir == "down", na.rm = TRUE)
    n_ns   <- sum(d$dir == "ns",   na.rm = TRUE)
    n_up   <- sum(d$dir == "up",   na.rm = TRUE)

    col_map <- c("down" = "steelblue3", "ns" = "gray70", "up" = "firebrick2")

    legend_labels <- c(
        down = sprintf("down (%s)", format(n_down, big.mark = ",")),
        ns   = sprintf("ns (%s)",   format(n_ns,   big.mark = ",")),
        up   = sprintf("up (%s)",   format(n_up,   big.mark = ","))
    )

    # Make R CMD HAPPY
    logFC<- neglog10FDR <- label <- NULL

    p <- ggplot2::ggplot(d, ggplot2::aes(x = logFC, y = neglog10FDR, color = dir)) +
        ggplot2::geom_point(alpha = point_alpha, size = 1.2) +
        ggplot2::scale_color_manual(
            values = c("down" = "steelblue3", "ns" = "gray70", "up" = "firebrick2"),
            breaks = c("down", "ns", "up"),
            labels = legend_labels,
            drop   = FALSE
        ) +
        ggplot2::geom_hline(yintercept = -log10(fdr_thresh), linetype = "dashed", linewidth = 0.4) +
        ggplot2::labs(
            title = title,
            x = "log2 fold change",
            y = expression(-log[10]("FDR")),
            color = "Direction"
        ) +
        ggplot2::theme_classic(base_size = 12)

    # Symmetric x-axis
    xmax <- max(abs(d$logFC), lfc_thresh, na.rm = TRUE)
    if (is.finite(xmax) && xmax > 0) {
        xlim <- c(-1, 1) * (1.1 * xmax)
        p <- p + ggplot2::scale_x_continuous(limits = xlim,
                                             expand = ggplot2::expansion(mult = 0))
    }


    if (lfc_thresh > 0) {
        p <- p + ggplot2::geom_vline(xintercept = c(-lfc_thresh, lfc_thresh),
                                     linetype = "dotted", linewidth = 0.4)
    }

    if (nrow(labs) > 0) {
        p <- p + ggrepel::geom_text_repel(
            data = labs,
            ggplot2::aes(x = logFC, y = neglog10FDR, label = label),
            size = 3,
            color = "black",
            min.segment.length = 0,
            max.overlaps = 10000,
            box.padding = 0.3,
            point.padding = 0.2,
            show.legend = FALSE  # keep legend showing points
        )
    }


    return (p)
}

paginate_plots <- function(plots, plots_per_page = 2) {
    if (!length(plots)) return(invisible(list()))
    ppp <- as.integer(plots_per_page)
    if (is.na(ppp) || ppp < 1) stop("plots_per_page must be a positive integer.")

    # pad with blanks so length is a multiple of ppp
    n_pad <- (ppp - (length(plots) %% ppp)) %% ppp
    if (n_pad > 0) plots <- c(plots, rep(list(cowplot::ggdraw()), n_pad))

    idx <- split(seq_along(plots), ceiling(seq_along(plots) / ppp))

    pages <- lapply(idx, function(ii) {
        cowplot::plot_grid(
            plotlist = plots[ii],
            ncol = 1, nrow = ppp,
            align = "v",
            rel_heights = rep(1, ppp)
        )
    })

    return (pages)
}


#quickly regenerate the PDF from the files in the result directory.
generate_pdf_from_files<-function (result_dir, outPDF) {
    files=list.files(result_dir, pattern="_DE_results.txt", full.names = TRUE, recursive = TRUE)
    plot_list= list()
    i=1
    for (f in files) {
        df=read.table(f, header=T, stringsAsFactors=F, sep="\t")
        n=gsub("_DE_results.txt", "", basename(f))
        p <- make_volcano(df, fdr_thresh = 0.05, lfc_thresh = 0,
                          top_n_each = 10, title = n)
        plot_list[[i]] <- p
        i=i+1
    }

    if (!is.null(outPDF)) {
        logger::log_info(paste("Saving all plots to PDF:", outPDF))
        grDevices::pdf(outPDF)
        pages=paginate_plots(plot_list, plots_per_page = 2)
        for (i in 1:length(pages)) {
            print(pages[[i]])
        }
        grDevices::dev.off()
    }
}

filter_dgelist_by_celltype_list<-function (dge, cellTypeListFile=NULL) {
    if (is.null(cellTypeListFile)) {
        return(dge)
    }

    size_prefilter <- dim(dge)[2]

    cell_type_list=read.table(cellTypeListFile, stringsAsFactors = FALSE, sep="\t", header=FALSE)$V1
    idx=dge$samples$cell_type %in% cell_type_list
    if (any(is.na(idx))) {
        stop("Some cell types in the list are not present in the DGEList samples.")
    }
    dge_filtered <- dge[, idx, keep.lib.sizes = TRUE]
    size_postfilter <- dim(dge_filtered)[2]
    logger::log_info(paste("Filtered DGEList from", size_prefilter, "to", size_postfilter, "metacells based on cell type list."))
    dge_filtered$samples$cell_type <- factor(dge_filtered$samples$cell_type, levels = cell_type_list)
    return(dge_filtered)
}

# for parallel processing on macOS and UGER/Linux.
make_bpparam <- function(n_cores) {
    # If we only want one core, or running on UGER (SGE/UGE),
    # fall back to SerialParam for safety.

    if (n_cores <= 1) {
        param=BiocParallel::SerialParam()
        logger::log_info("Using SerialParam for single-core execution.")
    } else {
        param=BiocParallel::MulticoreParam(
            workers       = n_cores,
            stop.on.error = TRUE,
            progressbar   = FALSE
        )
        logger::log_info(paste("Using MulticoreParam with", n_cores, "workers."))
    }
    return (param)
}


`%||%` <- function(a, b) if (is.null(a)) b else a
