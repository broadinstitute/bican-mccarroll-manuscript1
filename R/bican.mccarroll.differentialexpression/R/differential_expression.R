

# data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metacells"
# data_name="donor_rxn_DGEList"
# randVars=c("donor", "imputed_sex", "single_cell_assay", "region", "toxicology_group", "village", "biobank")
# fixedVars=c("age", "PC1", "PC2", "PC3", "PC4", "PC5", "pmi_hr", "pct_intronic", "frac_contamination")
# contrast_file


differential_expression <- function(data_dir, data_name, randVars, fixedVars, contrast_file, outPDF, result_dir) {
    #load the DGEList and prepare the data
    d=prepare_data_for_differential_expression(data_dir, data_name, randVars, fixedVars)
    dge=d$dge; fixedVars=d$fixedVars; randVars=d$randVars

    contrast_defs <- read.table(contrast_file, stringsAsFactors = FALSE, sep="\t", header=TRUE)

    # Variance Partition by cell type
    cell_type_list=unique(dge$samples$cell_type)
    #cell_type_list=cell_type_list[2]
    plotList=list()
    line <- strrep("=", 80)
    if (length(cell_type_list) > 0) {
        for (cellType in cell_type_list) {
            logger::log_info(line)
            logger::log_info(paste("Creating differential expression analysis for cell type:", cellType))
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

            #run differential expression
            z<-differential_expression_one_cell_type(dge_cell, fixedVars, randVars, verbose = TRUE)
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

}

differential_expression_one_cell_type<-function (dge_subset, fixedVars, randVars, verbose = TRUE, n_cores = parallel::detectCores() - 2) {

    # Drop random effects if they have insufficient replication
    rv <- prune_random_effects_insufficient_replication(randVars, data=dge_subset$samples)

    # Build formula string
    rand_part <- paste0("(1|", rv, ")", collapse = " + ")
    fixed_part <- paste(fixedVars, collapse = " + ")
    formula_str <- paste(fixed_part, rand_part, sep = " + ")
    form <- stats::as.formula(paste(" ~ ", formula_str))

    #Can't I just build the fixed effect formula directly?
    fixed_formula <- extract_fixed_effects_formula(form)
    design_matrix <- model.matrix(fixed_formula, data = dge_cell$samples)

    # 4. Generate contrast matrix from flat file & design
    contrast_matrix <- generate_contrasts_from_file(contrast_file, design_matrix)

    # 5. Run dream()
    param <- BiocParallel::MulticoreParam(workers = n_cores)
    fit <- variancePartition::dream(exprObj = dge_cell$counts, formula = form, data = dge_cell$samples, BPPARAM = param)

    # 6. Extract topTable results per contrast
    topTables_list <- list()
    for (contrast_name in colnames(contrast_matrix)) {
        topTables_list[[contrast_name]] <- limma::topTable(fit, coef = contrast_name, number = Inf)
    }


}

# Read contrast definitions
# contrast_defs <- read.delim(contrast_file, stringsAsFactors = FALSE)

generate_contrasts_from_file <- function(contrast_defs, design_matrix) {

    # Initialize list to hold contrast formulas
    contrast_list <- list()

    for (i in seq_len(nrow(contrast_defs))) {
        row <- contrast_defs[i, ]

        var <- row$variable
        ref <- row$reference_level
        comp <- row$comparison_level
        cname <- row$contrast_name

        if (row$type == "fixed") {
            # Continuous variables (simple contrasts)
            contrast_list[[cname]] <- var

        } else if (row$type == "random") {
            # Categorical contrasts: compare levels
            # This assumes model.matrix() uses factor level coding with level names embedded
            contrast_vector <- paste0(var, comp, " - ", var, ref)
            contrast_list[[cname]] <- contrast_vector
        }
    }

    # Build the contrast matrix
    contrast_formula <- paste(unlist(contrast_list), collapse = ", ")
    contrast_matrix <- limma::makeContrasts(contrasts = contrast_formula, levels = design_matrix)

    return(contrast_matrix)
}

generate_design_matrix <- function(dge_samples, formula_string) {
    # Extract only fixed effects part from the formula string
    # e.g., "~ age + PC1 + PC2 + imputed_sex + toxicology_group"
    form <- as.formula(formula_string)

    design_matrix <- model.matrix(form, data = dge_samples)
    return(design_matrix)
}
