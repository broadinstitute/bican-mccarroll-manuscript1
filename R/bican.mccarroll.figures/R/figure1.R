
# feature correlation plot
figure1_feature_correlation<-function (
    data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/metacells/LEVEL_3",
    data_name="donor_rxn_DGEList",
    randVars=c("donor", "imputed_sex", "biobank", "single_cell_assay", "region", "hbcac_status", "toxicology_group"),
    fixedVars=c("age", "PC1", "PC2", "PC3", "PC4", "PC5", "pmi_hr", "pct_intronic", "frac_contamination"),
    outDir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/figure_repository") {

    prep <- bican.mccarroll.differentialexpression::prepareMDSPlotData(
        data_dir = data_dir,
        data_name = data_name,
        additionalDonorMetadata = c(),
        randVars = randVars,
        fixedVars = fixedVars
    )
    dge <- prep$dge
    required_vars <- prep$required_vars

    correlation_vars <- setdiff(required_vars, "donor")
    corr_matrix <- bican.mccarroll.differentialexpression::getVariableCorrelation(dge$samples, cols_to_test = correlation_vars)

    #make feature names nice!
    pretty_map <- get_pretty_feature_names(correlation_vars)

    colnames(corr_matrix) <- as.vector(pretty_map[colnames(corr_matrix)])
    rownames(corr_matrix) <- as.vector(pretty_map[rownames(corr_matrix)])

    if (!is.null(outDir)) {
        output_svg <- file.path(outDir, "figure1_feature_correlation.svg")
        svg(output_svg, width = 8, height = 8)
        on.exit(dev.off(), add = TRUE)
    }

    #corrplot::corrplot(corr_matrix, method = "circle", type = "upper", tl.col = "black", tl.srt = 45, na.label = "NA")

    groups <- list(
        donor  = c("Imputed sex", "Biobank", "Age", "PC1", "PC2", "PC3", "PC4", "PC5", "Toxicology group", "HBCAC status", "PMI (hours)"),
        sample = c("Region", "Single-cell assay"),
        cell   = c("Percent intronic", "Fraction contamination")
    )

    setdiff(colnames(corr_matrix), as.vector(unlist (groups)))

    plot_corrplot_with_group_rects(corr_matrix, groups)

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
        PC1                 = "PC1",
        PC2                 = "PC2",
        PC3                 = "PC3",
        PC4                 = "PC4",
        PC5                 = "PC5"
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
