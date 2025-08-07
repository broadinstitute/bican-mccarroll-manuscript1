# library (vcd)
# library (Glimma)
# library (edgeR)
# library (corrplot)
# library (logger)
# library(ggplot2)

#library(variancePartition)
#library(edgeR)
#library(limma)

#library(BiocParallel)
#library (variancePartition)
#library(lme4)

# data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metacells"
# cellTypeGroupFile="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metadata/cell_type_groups.txt"
# data_name="donor_rxn_DGEList"

# randVars=c("donor", "imputed_sex", "biobank", "single_cell_assay", "region", "hbcac_status", "toxicology_group", "village")
# fixedVars=c("age", "PC1", "PC2", "PC3", "PC4", "PC5", "pmi_hr", "pct_intronic", "frac_contamination")
# variance_partition_result_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/variance_partition_results"
# outPDF=paste(variance_partition_result_dir, "/variance_partition_plots.pdf", sep="")


# randVars=c("donor", "imputed_sex", "biobank", "single_cell_assay", "region", "hbcac_status", "toxicology_group")
# fixedVars=c("age", "PC1", "PC2", "PC3", "PC4", "PC5", "pmi_hr", "pct_intronic", "frac_contamination")
# variance_partition_result_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/variance_partition_results_no_village"
# outPDF=paste(variance_partition_result_dir, "/variance_partition_plots-No_village.pdf", sep="")

# randVars=c("donor", "imputed_sex", "single_cell_assay", "region", "toxicology_group", "village")
# fixedVars=c("age", "PC1", "PC2", "PC3", "PC4", "PC5", "pmi_hr", "pct_intronic", "frac_contamination")
# variance_partition_result_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/variance_partition_results_no_hbcac_biobank"
# outPDF=paste(variance_partition_result_dir, "/variance_partition_plots-no_hbcac_biobank.pdf", sep="")

# randVars=c("donor", "imputed_sex", "single_cell_assay", "region", "toxicology_group", "village", "biobank")
# fixedVars=c("age", "PC1", "PC2", "PC3", "PC4", "PC5", "pmi_hr", "pct_intronic", "frac_contamination")
# variance_partition_result_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/variance_partition_results_no_hbcac"
# outPDF=paste(variance_partition_result_dir, "/variance_partition_plots-no_hbcac.pdf", sep="")

# library (bican.mccarroll.differentialexpression)
# bican.mccarroll.differentialexpression::runVariancePartition(data_dir = data_dir, data_name = data_name, randVars = randVars, fixedVars = fixedVars, outPDF = outPDF, variance_partition_result_dir = variance_partition_result_dir)

#' Run variance partition and generate QC plots
#'
#' Load a precomputed DGEList, validate sample variables, scale PCs,
#' filter samples by library size and runs variance partition both on the full
#' data set and by cell type.
#'
#' A special field in the DGEList samples metadata is num_nuclei, which is the total number of
#' nuclei captured in the single cell assay. If present, this is used in the variance partition model
#' to control for technical variation.  The value is converted to a Z-score of log10(num_nuclei) and
#' added to the fixed effects.
#'
#' @param data_dir                Character. Directory containing the precomputed DGEList.
#' @param data_name               Character. Prefix or name used to load the DGEList object.
#' @param randVars               Character vector. Names of random metadata variables for MDS coloring.
#' @param fixedVars              Character vector. Names of fixed metadata variables for MDS grouping.
#' @param outPDF                 QC and variance partition plots are saved to this single PDF file.
#' @param variance_partition_result_dir Character. Directory to save variance partition RDS objects to inspect later.
#' @importFrom logger log_info
#' @importFrom utils read.table
#' @importFrom grDevices pdf dev.off
#' @importFrom Glimma glimmaMDS
#' @importFrom vcd   assocstats
#' @importFrom corrplot corrplot
#' @importFrom edgeR DGEList
#' @export
runVariancePartition<-function (data_dir, data_name, randVars, fixedVars, outPDF, variance_partition_result_dir) {

    #load the DGEList and prepare the data
    d=prepare_data_for_differential_expression(data_dir, data_name, randVars, fixedVars)
    dge=d$dge;fixedVars=d$fixedVars;randVars=d$randVars

    # Variance Partition by cell type
    cell_type_list=unique(dge$samples$cell_type)
    #cell_type_list=cell_type_list[2]
    plotList=list()
    line <- strrep("=", 80)
    #cellType="GABA_CGE"
    if (length(cell_type_list) > 0) {
        for (cellType in cell_type_list) {
            logger::log_info(line)
            logger::log_info(paste("Creating variance partition analysis for cell type:", cellType))
            logger::log_info(line)

            dge_cell <- dge[, dge$samples$cell_type == cellType, keep.lib.sizes = TRUE]
            #filtering samples by library size
            r<- filter_by_libsize(dge_cell, threshold_sd = 1.96, bins = 50, strTitlePrefix = cellType)
            dge_cell<- r$dge

            #filter to the top 75% of highly expressed genes as a first pass.
            dge_cell<-filter_top_expressed_genes(dge_cell, gene_filter_frac = 0.75, verbose = TRUE)
            #filter to cpm cutoff of 1.
            r2=plot_logCPM_density_quantiles(dge_cell, cpm_cutoff = 1, logCPM_xlim = c(-5, 15), lower_quantile = 0.05, upper_quantile = 0.95, quantile_steps = 5)
            dge_cell=r2$filtered_dge

            p1=r$plot
            p1<-p1 + ggplot2::theme(plot.title = ggplot2::element_text(size = 10))

            p2=r2$plot
            p2<- p2 + ggplot2::theme(plot.title = ggplot2::element_text(size = 10))

            filtering_qc_plots=cowplot::plot_grid(p1, p2, ncol=1, label_size = 12)

            # batch to age correlation using canCorPairs(formula, data)
            # fString=paste(required_vars, collapse = " + ")
            # f<-formula(paste("~", fString, sep=" "))
            # C<- canCorPairs(formula=f, data=dge$samples)
            # plotCorrMatrix(C)


            varPart=run_variance_partition(dge_cell, fixedVars, randVars, verbose = TRUE, allCellTypes = FALSE)
            p3=variancePartition::plotVarPart(varPart, main = paste("Variance Partition", cellType), label.angle = 30)
            p3<- p3 + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 9))

            one_page<-cowplot::plot_grid(filtering_qc_plots, p3, ncol=1, label_size = 12)
            plotList[[cellType]]=one_page

            outRDS=paste(variance_partition_result_dir, "/", cellType, "_variance_partition.rds", sep="")
            saveRDS(varPart, file=outRDS)

        }
    } else {
        logger::log_info("No cell types found in the DGEList samples.")
    }

    if (!is.null(outPDF)) {
        logger::log_info(paste("Saving all plots to PDF:", outPDF))
        grDevices::pdf(outPDF)
        for (i in 1:length(plotList)) {
            print(plotList[[i]])
        }
        grDevices::dev.off()
    }

    #some QC plots up front.
    # z1=plot_variable_by_donor(dge,  variable = "frac_contamination", return_data = TRUE)
    # z2=plot_variable_by_region(dge, variable = "frac_contamination", return_data = TRUE)
    # z3=plot_variable_by_donor(dge, variable = "pct_intronic", return_data = TRUE)
    # z4=plot_variable_by_region(dge, variable = "pct_intronic", return_data = TRUE)

}

#' Prepare data for variance partition or differential expression
#'
#' Load a precomputed DGEList, validate sample variables, scale PCs,
#' filter samples by library size and runs variance partition both on the full
#' data set and by cell type.
#'
#' A special field in the DGEList samples metadata is num_nuclei, which is the total number of
#' nuclei captured in the single cell assay. If present, this is used in the variance partition model
#' to control for technical variation.  The value is converted to a Z-score of log10(num_nuclei) and
#' added to the fixed effects.
#'
#' @param data_dir                Character. Directory containing the precomputed DGEList.
#' @param data_name               Character. Prefix or name used to load the DGEList object.
#' @param randVars               Character vector. Names of random metadata variables for MDS coloring.
#' @param fixedVars              Character vector. Names of fixed metadata variables for MDS grouping.
#' @return A list containing the DGEList object, fixed variables, and random variables.
prepare_data_for_differential_expression<-function (data_dir, data_name, randVars, fixedVars) {
    # load the pre-computed DGEList object
    logger::log_info(paste("Loading DGEList from:", data_dir, "with prefix:", data_name))
    dge=bican.mccarroll.differentialexpression::loadDGEList(data_dir, prefix = data_name)

    #restrict to the groups that should be used for variance partition (which is the same as differential expression).
    dge=dge[,dge$samples$differential_expression, keep.lib.sizes = TRUE]

    #if the num_nuclei variable is present, convert it to log10 for regressions
    if ("num_nuclei" %in% colnames(dge$samples)) {
        logger::log_info("Converting num_nuclei to Z-score of log10(num_nuceli) and adding to fixed effects")
        dge$samples$z_log10_nuclei <- scale(log10(dge$samples$num_nuclei))
        #remove the original num_nuclei column
        dge$samples$num_nuclei <- NULL
        fixedVars=c("z_log10_nuclei", fixedVars)
    }

    # Validate the variables are present in the data set.
    required_vars=c(randVars, fixedVars)
    validateSampleVars(dge, required_vars)

    #scale genetic PCs to unit variance.
    dge$samples=bican.mccarroll.differentialexpression::scale_PC_cols(dge$samples)

    #convert to factors
    #don't convert norm.factors!
    dge$samples<-auto_factorize_df(dge$samples, exclude= c("norm.factors", "lib.size"), max_unique_numeric = 10)

    #scale age to decades for numeric stability.
    if ("age" %in% colnames(dge$samples)) {
        logger::log_info("Scaling age to decades for numeric stability")
        dge$samples$age <- dge$samples$age / 10
    }

    #Subset to complete cases
    dge=subset_dge_to_complete_cases(dge, required_vars)

    result<- list(
        dge = dge,
        fixedVars = fixedVars,
        randVars = randVars
    )
    return(result)
}

# if allCellTypes, add cell_type to the random variables
run_variance_partition <- function(dge_subset, fixedVars, randVars, verbose = TRUE, allCellTypes=FALSE, n_cores = parallel::detectCores() - 2) {
    #dge_subset<-dge_cell
    if (allCellTypes==TRUE) {
        logger::log_info(paste("Running variance partition for all cell types"))
        rv=c("cell_type", randVars)
    } else {
        logger::log_info(paste("Running variance partition for a single cell type"))
        rv=randVars
    }

    #It's possible after subsetting that some random effects only have one obseration per level
    #for example, in a cell type of a single region most donors might only have 1 observation,
    #so could not be estimated.
    #this prunes out those random effects.
    rv <- prune_random_effects_insufficient_replication(rv, data=dge_subset$samples)

    #alternative pruning where only factors with more than one observation per level are kept.
    #this is the bare minimum pruning that needs to take place to prevent errors in variancePartition.
    #This can result in numeric instability warnings from lme4 if the random effects are not well estimated.
    #rv <- drop_single_level_rand_effects(rv, metadata=dge_subset$samples, verbose = TRUE)

    # Build formula string
    rand_part <- paste0("(1|", rv, ")", collapse = " + ")
    fixed_part <- paste(fixedVars, collapse = " + ")
    formula_str <- paste(fixed_part, rand_part, sep = " + ")
    form <- stats::as.formula(paste(" ~ ", formula_str))

    if (verbose) {
        logger::log_info("Model formula:\n{paste(deparse(form), collapse = '\n')}")
    }

    # Voom transformation
    v <- limma::voom(dge_subset, plot = FALSE)

    # Set up parallel backend
    param <- BiocParallel::MulticoreParam(workers = n_cores)

    #running in batches with logging.
    varPart <- profile_variancePartition_runtime(exprObj = v, formula = form, data = dge_subset$samples, batch_size = 1000, BPPARAM = param, verbose=FALSE)

    # logger::log_info(paste("Running variance partition with", n_cores, "cores"))
    # varPart2 <- variancePartition::fitExtractVarPartModel(exprObj=v, formula=form, data=dge_subset$samples, BPPARAM = param)
    # varPartFull <- variancePartition::fitVarPartModel(exprObj=v[1:4], formula=form, data=dge_subset$samples[1:4], BPPARAM = param)
    # logger::log_info(paste("Variance partition completed with", nrow(varPart), "genes."))

    return(varPart)

    # r=as.data.frame(varPart)

    # plot_stratified_gene(dgeThis, varPart, variable="age", gene_name = NULL, main = NULL)
    # plot_stratified_gene(dgeThis, varPart, variable="pmi_hr", gene_name = NULL, main = NULL)
    # plot_stratified_gene(dgeThis, varPart, variable="toxicology_group", gene_name = NULL, main = NULL)
    # plot_stratified_gene(dgeThis, varPart, variable="single_cell_assay", gene_name = NULL, main = NULL)
    #
    # plotAdjustedResidualsByGroup(dgeThis, fixedVars, randVars, variable="age", gene = NULL)
    # plotAdjustedResidualsByGroup(dgeThis, fixedVars, randVars, variable="frac_contamination", gene = NULL)
    # plotAdjustedResidualsByGroup(dgeThis, fixedVars, randVars, variable="pmi_hr", gene = NULL)
    #
    # plotAdjustedResidualsByGroup(dgeThis, fixedVars, randVars, variable="toxicology_group", gene = NULL)
    # plotAdjustedResidualsByGroup(dgeThis, fixedVars, randVars, variable="single_cell_assay", gene = NULL)
    #
    # plotPercentBars(varPart[rownames(varPart)=="CD2AP", ])
    # plotPercentBars(varPart[rownames(varPart)=="FAM66C", ])

    #an interesting finding - the gene with the highest variance explained by toxicology is CD2AP
    #but that gene is better explained by cell type.



}

#' Build a Variance Partitioning Model Formula
#'
#' Constructs a model formula suitable for use with `variancePartition::fitExtractVarPartModel()`
#' by combining fixed and random effects into a single right-hand side formula.
#'
#' @param fixedVars Character vector of fixed effect variable names. These will be included
#'   in the formula as standard linear terms.
#' @param randVars Character vector of random effect variable names. These will be added to
#'   the formula in the format `(1|var)`, appropriate for use with mixed models.
#'
#' @return An object of class `formula` representing the combined fixed and random effects
#'   model, with a right-hand side suitable for use in variance partitioning.
#'
#' @export
buildVariancePartitionModelFormula<-function (fixedVars, randVars) {
    rand_part <- paste0("(1|", randVars, ")", collapse = " + ")
    fixed_part <- paste(fixedVars, collapse = " + ")
    formula_str <- paste(fixed_part, rand_part, sep = " + ")
    form <- stats::as.formula(paste("~", formula_str))
    return(form)
}

auto_factorize_df <- function(df, max_unique_numeric = 10, exclude = NULL) {
    stopifnot(is.data.frame(df))
    stopifnot(is.null(exclude) || all(exclude %in% colnames(df)))

    df[] <- lapply(names(df), function(col_name) {
        col <- df[[col_name]]

        # Skip excluded columns
        if (!is.null(exclude) && col_name %in% exclude) {
            return(col)
        }

        # Factorize logic
        if (is.character(col)) {
            factor(col)
        } else if (is.numeric(col) && length(unique(col)) < max_unique_numeric) {
            factor(col)
        } else {
            col
        }
    })

    return(as.data.frame(df))
}


subset_dge_to_complete_cases <- function(dge, vars, verbose=TRUE) {
    stopifnot("DGEList" %in% class(dge))
    stopifnot(all(vars %in% colnames(dge$samples)))

    # Subset to relevant metadata columns and find complete rows
    metadata_subset <- dge$samples[, vars, drop = FALSE]
    keep_rows <- stats::complete.cases(metadata_subset)

    # Subset DGEList
    dge_subset <- dge[, keep_rows, keep.lib.sizes = FALSE]

    # Also subset the samples metadata
    # Drop unused factor levels in all metadata columns after subsetting
    dge_subset$samples <- droplevels(dge$samples[keep_rows, , drop = FALSE])

    return(dge_subset)
}


plot_stratified_gene <- function(dge, varPart, variable, gene_name = NULL, main = NULL) {
    # Validate inputs
    stopifnot("DGEList" %in% class(dge))

    # Get gene expression matrix (logCPM)
    logCPM <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

    # Identify gene to plot
    if (is.null(gene_name)) {
        i <- which.max(varPart[[variable]])
        gene_name <- rownames(varPart)[i]
        message(sprintf("Plotting top gene for '%s': %s", variable, gene_name))
    } else {
        i <- match(gene_name, rownames(logCPM))
        if (is.na(i)) stop(sprintf("Gene '%s' not found in DGEList.", gene_name))
    }

    variance_explained <- varPart[gene_name, variable]
    var_pct <- round(variance_explained * 100, 1)

    # Prepare data frame
    GE <- data.frame(
        Expression = logCPM[i, ],
        Group = dge$samples[[variable]]
    )

    # Set plot title
    if (is.null(main)) {
        main <- sprintf("%s (%.1f%% variance explained by %s)", gene_name, var_pct, variable)
    }

    # Plot using variancePartition's stratified boxplot
    variancePartition::plotStratify(Expression ~ Group, GE, main = main, xlab = variable)
}

plotAdjustedResidualsByGroup <- function(varPart, dge, fixedVars, randVars, variable, gene = NULL) {
    stopifnot("DGEList" %in% class(dge))
    stopifnot(variable %in% colnames(dge$samples))

    # Drop variable of interest from both fixed and random
    fixedVars_adj <- setdiff(fixedVars, variable)
    randVars_adj <- setdiff(randVars, variable)

    # Build adjusted formula
    form <- buildVariancePartitionModelFormula(fixedVars_adj, randVars_adj)

    # Compute logCPM
    logCPM <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

    # Pick gene with highest variance explained by variable (if gene is NULL)
    if (is.null(gene)) {
        i <- which.max(varPart[[variable]])
        gene <- rownames(varPart)[i]
        message(sprintf("Plotting gene with highest variance explained by '%s': %s", variable, gene))
    } else {
        i <- match(gene, rownames(logCPM))
        if (is.na(i)) stop(sprintf("Gene '%s' not found in DGEList.", gene))
    }

    # Expression and metadata
    expr <- logCPM[i, ]
    metadata <- dge$samples
    metadata$expr <- expr

    # Fit mixed model
    fit <- lme4::lmer(stats::update(form, expr ~ .), data = metadata)

    # Prepare data for plotting
    df_plot <- data.frame(
        Residual = stats::residuals(fit),
        Group = metadata[[variable]]
    )

    variance_explained <- varPart[gene, variable]
    var_pct <- round(variance_explained * 100, 1)

    # Set plot title
    main_label <- sprintf("%s (%.1f%% variance explained by %s)", gene, var_pct, variable)

    # Plot
    # make R CMD CHECK happy
    Group<-Residual<-NULL

    if (is.numeric(df_plot$Group)) {
        ggplot2::ggplot(df_plot, aes(x = Group, y = Residual)) +
            ggplot2::geom_point(alpha = 0.7, size = 1.5) +
            ggplot2::geom_smooth(method = "loess", se = TRUE, color = "steelblue", formula = y ~ x) +
            ggplot2::labs(
                title = main_label,
                x = variable,
                y = "Residual expression"
            ) +
            ggplot2::theme_minimal()
    } else {
        variancePartition::plotStratify(Residual ~ Group, df_plot, main = main_label, xlab = variable)
    }
}


#################
# ADHOC STUFFS
#################


plot_variable_by_donor <- function(dge, variable = "frac_contamination", return_data = TRUE) {
    stopifnot("DGEList" %in% class(dge))
    df <- dge$samples
    stopifnot("donor" %in% colnames(df))
    stopifnot(variable %in% colnames(df))

    donors <- unique(df$donor)
    donor_stats <- data.frame(
        donor = donors,
        median = numeric(length(donors)),
        IQR = numeric(length(donors)),
        n = integer(length(donors)),
        stringsAsFactors = FALSE
    )

    for (i in seq_along(donors)) {
        vals <- df[df$donor == donors[i], variable]
        donor_stats$median[i] <- stats::median(vals, na.rm = TRUE)
        donor_stats$IQR[i] <- stats::IQR(vals, na.rm = TRUE)
        donor_stats$n[i] <- length(vals)
    }

    donor_order <- donor_stats$donor[order(donor_stats$median)]
    df$donor <- factor(df$donor, levels = donor_order)

    # make R CMD CHECK happy
    donor <- NULL

    p <- ggplot(df, aes(x = donor, y = .data[[variable]])) +
        geom_boxplot(outlier.alpha = 0.4, fill = "skyblue") +
        labs(
            title = paste(variable, "by Donor"),
            x = "Donor",
            y = variable
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

    print(p)

    if (return_data) {
        return(donor_stats[order(donor_stats$median, decreasing = TRUE), ])
    } else {
        invisible(NULL)
    }
}




plot_variable_by_region <- function(dge, variable = "frac_contamination", return_data = TRUE) {
    stopifnot("DGEList" %in% class(dge))
    df <- dge$samples
    stopifnot("region" %in% colnames(df))
    stopifnot(variable %in% colnames(df))

    regions <- unique(df$region)
    region_stats <- data.frame(
        region = regions,
        median = numeric(length(regions)),
        IQR = numeric(length(regions)),
        n = integer(length(regions)),
        stringsAsFactors = FALSE
    )

    for (i in seq_along(regions)) {
        vals <- df[df$region == regions[i], variable]
        region_stats$median[i] <- stats::median(vals, na.rm = TRUE)
        region_stats$IQR[i] <- stats::IQR(vals, na.rm = TRUE)
        region_stats$n[i] <- length(vals)
    }

    region_order <- region_stats$region[order(region_stats$median)]
    df$region <- factor(df$region, levels = region_order)

    # make R CMD CHECK happy
    region <- NULL

    p <- ggplot(df, aes(x = region, y = .data[[variable]])) +
        geom_boxplot(outlier.shape = NA, fill = "lightgreen", alpha = 0.6) +
        geom_jitter(width = 0.2, alpha = 0.5, size = 1.2, color = "black") +
        labs(
            title = paste(variable, "by Region"),
            x = "Region",
            y = variable
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

    print(p)

    if (return_data) {
        return(region_stats[order(region_stats$median, decreasing = TRUE), ])
    } else {
        invisible(NULL)
    }
}

plot_variable_by_group <- function(dge, variable = "frac_contamination", group_col = "region", return_data = TRUE) {
    stopifnot("DGEList" %in% class(dge))
    df <- dge$samples
    stopifnot(group_col %in% colnames(df))
    stopifnot(variable %in% colnames(df))

    groups <- unique(df[[group_col]])
    group_stats <- data.frame(
        group = groups,
        median = numeric(length(groups)),
        IQR = numeric(length(groups)),
        n = integer(length(groups)),
        stringsAsFactors = FALSE
    )

    for (i in seq_along(groups)) {
        vals <- df[df[[group_col]] == groups[i], variable]
        group_stats$median[i] <- stats::median(vals, na.rm = TRUE)
        group_stats$IQR[i] <- stats::IQR(vals, na.rm = TRUE)
        group_stats$n[i] <- length(vals)
    }

    group_order <- group_stats$group[order(group_stats$median)]
    df[[group_col]] <- factor(df[[group_col]], levels = group_order)

    # make R CMD CHECK happy
    group <- NULL

    p <- ggplot2::ggplot(df, ggplot2::aes_string(x = group_col, y = variable)) +
        ggplot2::geom_boxplot(outlier.shape = NA, fill = "lightgreen", alpha = 0.6) +
        ggplot2::geom_jitter(width = 0.2, alpha = 0.5, size = 1.2, color = "black") +
        ggplot2::labs(
            title = paste(variable, "by", group_col),
            x = group_col,
            y = variable
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    print(p)

    if (return_data) {
        colnames(group_stats)[1] <- group_col  # Rename 'group' back to original column name
        return(group_stats[order(group_stats$median, decreasing = TRUE), ])
    } else {
        invisible(NULL)
    }
}


plot_logCPM_density_quantiles <- function(dge, cpm_cutoff = 0.5, logCPM_xlim = c(-5, 15), lower_quantile = 0.05, upper_quantile = 0.95, quantile_steps = 5) {
    #Make R CMD CHECK HAPPY
    y<- x <- ymin <- ymax <- band <- grid_x <- mean_density <- NULL

    stopifnot("DGEList" %in% class(dge))

    nsamples <- ncol(dge$counts)
    orig_gene_count <- nrow(dge$counts)

    # Compute log-CPM matrix
    lcpm <- edgeR::cpm(dge, log = TRUE)

    # Filter genes by CPM cutoff
    keep_genes <- rowSums(edgeR::cpm(dge) >= cpm_cutoff) > 0
    filtered_dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
    filtered_gene_count <- sum(keep_genes)

    # Recompute log-CPM matrix on filtered genes
    lcpm_filtered <- edgeR::cpm(filtered_dge, log = TRUE)

    # Create grid for density estimation
    grid_x <- seq(logCPM_xlim[1], logCPM_xlim[2], length.out = 1000)
    density_matrix <- matrix(0, nrow = length(grid_x), ncol = nsamples)

    # Compute density per sample
    for (i in seq_len(nsamples)) {
        den <- stats::density(lcpm_filtered[, i], from = logCPM_xlim[1], to = logCPM_xlim[2], n = length(grid_x))
        density_matrix[, i] <- den$y
    }

    # Compute mean density
    mean_density <- rowMeans(density_matrix)

    # Prepare quantile bands
    quantile_levels <- seq(0.5, upper_quantile, length.out = quantile_steps)

    band_data <- data.frame()
    for (i in seq_along(quantile_levels)) {
        lower_prob <- 1 - quantile_levels[i]
        upper_prob <- quantile_levels[i]
        lower_band <- apply(density_matrix, 1, stats::quantile, probs = lower_prob)
        upper_band <- apply(density_matrix, 1, stats::quantile, probs = upper_prob)
        band_data <- rbind(band_data, data.frame(
            x = grid_x,
            ymin = lower_band,
            ymax = upper_band,
            band = i
        ))
    }

    # Blues color palette (lighter for outer quantiles)
    band_colors <- rev(grDevices::colorRampPalette(c("lightblue", "blue"))(length(quantile_levels)))

    # Plot exclusive bands (non-stacking)
    p <- ggplot2::ggplot()


    for (i in seq_along(quantile_levels)) {
        band_df <- subset(band_data, band == i)
        p <- p + ggplot2::geom_ribbon(
            data = band_df,
            ggplot2::aes(x = x, ymin = ymin, ymax = ymax),
            fill = band_colors[i],
            alpha = 0.6
        )
    }

    # Add mean line and cutoff line
    strTitle=paste0("Expression filtering by CPM (All Genes: ", orig_gene_count, " filtered to ", filtered_gene_count, ")")



    p <- p +
        ggplot2::geom_line(
            data = data.frame(x = grid_x, y = mean_density),
            ggplot2::aes(x = x, y = y),
            color = "black",
            linewidth = 1.2
        ) +
        ggplot2::geom_vline(xintercept = log2(cpm_cutoff), linetype = "dashed", color = "red") +
        ggplot2::labs(
            title =strTitle,
            x = "Log-CPM",
            y = "Density"
        ) +
        ggplot2::theme_minimal()



    return(list(filtered_dge = filtered_dge, plot = p))
}

filter_top_expressed_genes <- function(dge, gene_filter_frac = 0.5, verbose = TRUE) {
    stopifnot("DGEList" %in% class(dge))
    stopifnot(gene_filter_frac > 0 && gene_filter_frac <= 1)

    total_counts <- rowSums(dge$counts)
    n_keep <- floor(length(total_counts) * gene_filter_frac)

    # Rank genes by total counts and select top fraction
    keep <- rank(-total_counts, ties.method = "first") <= n_keep

    # Subset DGEList and recalculate normalization factors
    dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]
    dge_filtered <- edgeR::calcNormFactors(dge_filtered)

    if (verbose) {
        logger::log_info(paste("Number of genes after filtering:", nrow(dge_filtered)))
    }

    return(dge_filtered)
}


prune_random_effects_insufficient_replication <- function(randVars, data, min_replicated_levels = 2, min_fraction_obs = 0.1, verbose = TRUE) {
    pruned_randVars <- randVars

    for (var in randVars) {
        n_obs_per_level <- table(data[[var]])
        n_levels_with_replication <- sum(n_obs_per_level >= 2)
        n_obs_in_replicated_levels <- sum(n_obs_per_level[n_obs_per_level >= 2])
        total_obs <- nrow(data)
        fraction_obs_in_replicated_levels <- n_obs_in_replicated_levels / total_obs

        if (n_levels_with_replication < min_replicated_levels || fraction_obs_in_replicated_levels < min_fraction_obs) {
            if (verbose) {
                logger::log_info(paste0("Dropping random effect ", var, " - insufficient replication (",
                                        n_levels_with_replication, " levels with >=2 observations, covering ",
                                        round(100 * fraction_obs_in_replicated_levels, 1), "% of observations)."))
            }
            pruned_randVars <- setdiff(pruned_randVars, var)
        }
    }

    return(pruned_randVars)
}

drop_single_level_rand_effects <- function(randVars, metadata, verbose = TRUE) {
    stopifnot(all(randVars %in% colnames(metadata)))

    retainedVars <- randVars
    for (var in randVars) {
        unique_levels <- unique(na.omit(metadata[[var]]))
        if (length(unique_levels) <= 1) {
            retainedVars <- setdiff(retainedVars, var)
            if (verbose) {
                message(paste0("Dropping random effect '", var, "' - only ", length(unique_levels), " level(s) present in data."))
            }
        }
    }

    return(retainedVars)
}


diagnose_random_effects_replication <- function(randVars, data) {
    stopifnot(is.data.frame(data))
    stopifnot(all(randVars %in% colnames(data)))

    diagnostic_report <- data.frame(
        Variable = character(),
        Levels = integer(),
        Singleton_Levels = integer(),
        Fraction_Singleton_Levels = numeric(),
        Levels_with_2plus_Obs = integer(),
        Replicated_Obs_Count = integer(),
        stringsAsFactors = FALSE
    )

    for (var in randVars) {
        n_obs_per_level <- table(data[[var]])
        levels <- length(n_obs_per_level)
        singleton_levels <- sum(n_obs_per_level == 1)
        singleton_fraction <- singleton_levels / levels
        levels_with_replicates <- sum(n_obs_per_level >= 2)
        replicated_obs_count <- sum(n_obs_per_level[n_obs_per_level >= 2])

        diagnostic_report <- rbind(diagnostic_report, data.frame(
            Variable = var,
            Levels = levels,
            Singleton_Levels = singleton_levels,
            Fraction_Singleton_Levels = round(singleton_fraction, 2),
            Levels_with_2plus_Obs = levels_with_replicates,
            Replicated_Obs_Count = replicated_obs_count,
            stringsAsFactors = FALSE
        ))
    }

    return(diagnostic_report)
}

#' Profile Runtime of Variance Partitioning by Gene Batches
#'
#' Runs `variancePartition::fitExtractVarPartModel()` in batches to profile runtime
#' across a large gene expression dataset, reporting per-batch and total runtimes.
#' Designed for benchmarking or runtime monitoring on large datasets.
#'
#' @param exprObj A gene expression object such as a matrix, `DGEList`, or any format
#'   compatible with `variancePartition::fitExtractVarPartModel()`. Rows represent genes.
#' @param formula A model formula specifying the variance partitioning design, typically
#'   involving covariates in `data` (e.g., `~ (1|Individual) + Age + Sex`).
#' @param data A data frame with metadata corresponding to the columns (samples) in `exprObj`.
#'   Must contain all variables referenced in `formula`.
#' @param batch_size Integer. Number of genes to include in each batch (default is 1000).
#' @param BPPARAM A BiocParallel parameter object used to control parallel processing.
#'   Defaults to `BiocParallel::SerialParam()` (serial execution).
#' @param verbose Logical. If `TRUE` (default), print progress messages and per-batch runtimes.
#'
#' @return A variancePartition object with variance estimates for all genes,
#'   as returned by `variancePartition::fitExtractVarPartModel()`, combined across batches.
#' @export

profile_variancePartition_runtime <- function(exprObj, formula, data, batch_size = 1000, BPPARAM = BiocParallel::SerialParam(), verbose = TRUE) {
    stopifnot(!is.null(exprObj))  # Accept any supported exprObj types (matrix, DGEList, etc.)

    n_genes <- nrow(exprObj)
    batches <- split(1:n_genes, ceiling(seq_along(1:n_genes) / batch_size))
    n_batches <- length(batches)

    varPart_results <- list()  # Collect results per batch

    logger::log_info(paste0("Running variance partition in ", n_batches, " batches of size ", batch_size, " (total genes: ", n_genes, ")"))

    pb <- progress::progress_bar$new(
        format = "Progress [:bar] :current/:total (:percent) eta: :eta",
        total = n_batches,
        clear = FALSE,
        width = 60
    )

    total_start_time <- Sys.time()

    for (i in seq_along(batches)) {
        batch_idx <- batches[[i]]

        if (verbose) {
            pb$message(paste0(
                "Processing batch ", i,
                " (genes ", min(batch_idx), "-", max(batch_idx),
                ") of ", n_genes, "..."
            ))
        }

        batch_start_time <- Sys.time()

        varPart_batch <- variancePartition::fitExtractVarPartModel(
            exprObj = exprObj[batch_idx, ],
            formula = formula,
            data = data,
            BPPARAM = BPPARAM
        )

        batch_end_time <- Sys.time()
        batch_runtime_sec <- as.numeric(difftime(batch_end_time, batch_start_time, units = "secs"))
        per_gene_runtime <- batch_runtime_sec / length(batch_idx)

        varPart_results[[i]] <- varPart_batch
        if (verbose) {
            pb$message(paste0(
                "Batch ", i, " processed in ",
                round(batch_runtime_sec, 2), " seconds (",
                round(per_gene_runtime, 3), " sec/gene)"
            ))
        }
        pb$tick()
    }

    total_end_time <- Sys.time()
    total_runtime_sec <- as.numeric(difftime(total_end_time, total_start_time, units = "secs"))

    logger::log_info(paste0("Total runtime: ", round(total_runtime_sec, 2), " seconds"))

    #Convert back to a standard variancePartition object.
    # Combine all varPart batch results into a single object
    combined_varPart <- bindVarPartResults(varPart_results)
    return(combined_varPart)
}

bindVarPartResults <- function(varPart_list) {
    stopifnot(is.list(varPart_list))
    stopifnot(all(sapply(varPart_list, function(x) inherits(x, "varPartResults"))))

    # Combine batches using base R rbind (since it's a plain data.frame subclass)
    combined <- do.call(rbind, varPart_list)

    # Restore attributes
    attr(combined, "class") <- structure("varPartResults", package = "variancePartition")
    attr(combined, "package") <- "variancePartition"
    attr(combined, "row.names") <- rownames(combined)
    attr(combined, "names") <- colnames(combined)
    attr(combined, ".S3Class") <- "data.frame"

    # These metadata attributes are static
    attr(combined, "type") <- "linear mixed model"
    attr(combined, "method") <- "Variance explained (%)"

    # Restore formula from first object
    attr(combined, "formula") <- attr(varPart_list[[1]], "formula")

    return(combined)
}


