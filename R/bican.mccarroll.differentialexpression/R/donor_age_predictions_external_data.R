# library (ggplot2)
# library(dplyr)
# library(gganimate)
# library(cowplot)


# Test the donor age prediction model on external datasets using already trained models
# Test the donor age prediction model the region by region data to look at both prediction accuracy and the
# correlation of residuals within cell type across regions.

# data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metacells"
# data_name="donor_rxn_DGEList"
# model_file_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/age_prediction"
# result_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/age_prediction"
# #
# retained_features=c("donor", "age")
# donor_col = "donor"
# age_col = "age"
# out_region_pdf_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/age_prediction/age_prediction_residuals_by_region.pdf"
# out_celltype_pdf_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/age_prediction/age_prediction_residuals_by_celltype.pdf"
# out_region_predictions_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/age_prediction/age_prediction_region_results.txt"
#
# # Restrict the cell correlation analysis to a subset of cell types
# cellTypeListFile="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metadata/mash_cell_type_list_simple.txt"

# Test fits on each region.
# This function uses a joint across region average model to try and predict age in each region separately.
# I'm starting to think this is a not good idea.
# test_donor_age_predictions_region_data <- function(data_dir, data_name,
#                                                    model_file_dir,
#                                                    result_dir,
#                                                    contig_yaml_file = NULL,
#                                                    reduced_gtf_file = NULL,
#                                                    retained_features = c("donor", "age"),
#                                                    donor_col = "donor",
#                                                    age_col = "age",
#                                                    out_region_pdf_file = NULL,
#                                                    out_celltype_pdf_file = NULL,
#                                                    out_region_predictions_file = NULL) {
#     #validate the output directory exists
#     if (!dir.exists(result_dir)) {
#         stop("Result directory does not exist: ", result_dir)
#     }
#
#     #Load the model coefficients.
#     all_models = load_models(model_file_dir)
#     model_metrics = load_model_metrics(model_file_dir)
#
#     #load the DGEList and prepare the data
#     d = bican.mccarroll.differentialexpression::prepare_data_for_differential_expression(data_dir,
#                                                                                          data_name,
#                                                                                          randVars = c(),
#                                                                                          fixedVars = c())
#     dge = d$dge
#
#     #filter to autosomes, so we don't learn sex / age biases
#     dge = filter_dge_to_autosomes (dge, contig_yaml_file, reduced_gtf_file)
#
#     #validate all model genes in DGE
#     model_genes = unique(all_models$feature)
#     missing_genes = setdiff(model_genes, rownames(dge$counts))
#     if (length(missing_genes) > 0) {
#         stop(
#             "The following model genes are missing from the DGEList counts matrix: ",
#             paste(missing_genes, collapse = ", ")
#         )
#     }
#
#     cell_type_list = unique(dge$samples$cell_type)
#     #cell_type_list="microglia" #hard coded for now.
#     lineStr <- strrep("=", 80)
#
#     results = list()
#
#     #cellType="microglia"; region="CaH"
#     for (cellType in cell_type_list) {
#         logger::log_info(lineStr)
#         logger::log_info(paste(
#             "Learning donor age model from expression for cell type:",
#             cellType
#         ))
#         logger::log_info(lineStr)
#
#         regions = unique (dge$samples[dge$samples$cell_type == cellType, ]$region)
#         model = all_models[all_models$cell_type == cellType, ]
#
#         for (region in regions) {
#             logger::log_info(paste("  Testing region:", region))
#
#             #subset to the cell type and region
#             dge_sub = dge[, dge$samples$cell_type == cellType &
#                               dge$samples$region == region]
#
#             #collapse to one observation per donor
#             dge_sub <- collapse_by_donor(dge_sub,
#                                          donor_col = donor_col,
#                                          keep_cols = retained_features)
#
#             #filtering samples by library size - not predicting very small samples.
#             r <- bican.mccarroll.differentialexpression::filter_by_libsize(
#                 dge_sub,
#                 threshold_sd = 1.96,
#                 bins = 50,
#                 strTitlePrefix = cellType
#             )
#             dge_sub <- r$dge
#
#             #test the model
#             result <- predict_age_from_dge(dge_sub, model, prior.count = 1)
#             result <- cbind(cell_type = cellType, region = region, result)
#             results[[paste(cellType, region, sep = "_")]] = result
#         }
#     }
#
#     results = do.call(rbind, results)
#
#     if (!is.null(out_region_predictions_file))
#         write.table(
#             results,
#             file = out_region_predictions_file,
#             sep = "\t",
#             quote = FALSE,
#             row.names = FALSE
#         )
#
#     #compute metrics per region
#     metrics = get_age_prediction_metrics(results)
#
#     if (!is.null(out_region_pdf_file))
#         pdf(out_region_pdf_file)
#
#     for (cellType in cell_type_list) {
#         m = model_metrics[model_metrics$cell_type == cellType &
#                               model_metrics$set == "Held-out donors", ]
#         titleStr <- paste(
#             "Region specific predictions for ",
#             cellType,
#             " (on held-out donors)\n",
#             "Overall model metrics: r = ",
#             sprintf("%.3f", m$r),
#             ", MedAE = ",
#             sprintf("%.3f", m$median_abs_error),
#             ", MAE = ",
#             sprintf("%.3f", m$mean_abs_error),
#             sep = ""
#         )
#
#
#         p1 = plot_age_predictions_by_region(
#             results,
#             cellType = cellType,
#             titleStr = titleStr,
#             facet_col = "region",
#             metrics_df = metrics
#         )
#         print (p1)
#         res_mat <- compute_residual_matrix(results, cellType = cellType)
#         if (ncol(res_mat) >= 2) {
#             p2 <- plot_residual_pair_scatter(res_mat, cellType)
#             print (p2)
#             p3 <- plot_residual_corr_heatmap(res_mat, cellType)
#             ComplexHeatmap::draw(p3, heatmap_legend_side = "right")
#         }
#     }
#
#     if (!is.null(out_region_pdf_file))
#         dev.off()
# }

# Compare residuals across cell types
compare_age_residuals_celltype <- function (model_file_dir,
                                            cellTypeListFile = NULL,
                                            out_celltype_pdf_file = NULL) {
    model_predictions = load_model_predictions(model_file_dir)

    if (!is.null(cellTypeListFile)) {
        cell_types = read.table(cellTypeListFile,
                                header = FALSE,
                                stringsAsFactors = FALSE)$V1
        model_predictions = model_predictions[model_predictions$cell_type %in% cell_types, ]
    }

    res_mat <- compute_residual_matrix(
        model_predictions,
        group_cols = c("donor", "cell_type"),
        drop_empty = TRUE
    )
    p2 <- plot_residual_pair_scatter_paged(
        res_mat,
        cellType = NULL,
        per_page = 4,
        facet_font_size = 8
    )

    p3 <- plot_residual_corr_heatmap(res_mat, cellType = NULL, annotate_cells = TRUE)

    p4 <- compare_age_model_features(model_file_dir, cellTypeListFile= cellTypeListFile)

    if (!is.null(out_celltype_pdf_file))
        pdf (out_celltype_pdf_file)

    ComplexHeatmap::draw(p4, heatmap_legend_side = "right")
    ComplexHeatmap::draw(p3, heatmap_legend_side = "right")
    for (p in p2) {
        print (p)
    }

    if (!is.null(out_celltype_pdf_file))
        dev.off()

}

# Compare overlap of genes across models - this is also used in compare_age_residuals_celltype
compare_age_model_features <- function (model_file_dir, cellTypeListFile = NULL) {
    all_models = load_models(model_file_dir)

    if (!is.null(cellTypeListFile)) {
        cell_types = read.table(cellTypeListFile,
                                header = FALSE,
                                stringsAsFactors = FALSE)$V1
        all_models = all_models[all_models$cell_type %in% cell_types, ]
    }

    J <- jaccard_by_celltype(all_models, coef_thresh = 0)  # or small error, e.g., 1e-8
    ht <- plot_jaccard_heatmap(J, title = "Cell-type gene overlap in aging programs", annotate_cells = TRUE)

    return (ht)
}

# Use cell type model to predict external data.
# This is the "old" original data.
# training_data_samples_file = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metacells/donor_rxn_DGEList_samples.tsv.gz"
#
# covariate_file = "/broad/mccarroll/dropulation/analysis/latent_factor/metadata/BA46/BA46.n191/BA46.n191.knowncovars.txt"
# cell_type = "astrocyte"
# metacell_file = "/broad/mccarroll/dropulation/analysis/eQTL/BA46/BA46.n191.raw.dge.astrocyte.All/maf_0.05_cisDist_10kb/BA46.n191.raw.dge.astrocyte.All.noOutliers.metacells.txt"
# cell_type = "microglia"
# metacell_file = "/broad/mccarroll/dropulation/analysis/eQTL/BA46/BA46.n191.raw.dge.microglia.All/maf_0.05_cisDist_10kb/BA46.n191.raw.dge.microglia.All.noOutliers.metacells.txt"
# source_gtf_file = "/broad/mccarroll/software/metadata/individual_reference/GRCh38_ensembl_v43/GRCh38_ensembl_v43.reduced.gtf"
# target_gtf_file = "/broad/mccarroll/software/metadata/individual_reference/GRCh38_maskedAlt.89/GRCh38_maskedAlt.reduced.gtf"

# Reprocessed data on gencode v44
# covariate_file = "/broad/mccarroll/dropulation/analysis/cellarium_upload/SNAP200_freeze1/BA46.n180.knowncovars_simple.txt"
# cell_type = "astrocyte"
# metacell_file = "/broad/mccarroll/dropulation/analysis/cellarium_upload/SNAP200_freeze1/metacells/astrocyte__BA46.metacells.txt.gz"
# # cell_type ="microglia"; metacell_file="/broad/mccarroll/dropulation/analysis/cellarium_upload/SNAP200_freeze1/metacells/microglia__BA46.metacells.txt.gz"
# source_gtf_file = "/broad/mccarroll/software/metadata/individual_reference/GRCh38_ensembl_v43/GRCh38_ensembl_v43.reduced.gtf"
# target_gtf_file = "/broad/mccarroll/dropulation/analysis/cellarium_upload/SNAP200_freeze1/modified_v44.annotation.reduced.gtf"
# model_file_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/age_prediction"
# model_file_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/age_prediction_alpha_0"

##########################
# MINI-TASK: Optimize donor subset to better match mean expression
##########################

optimize_test_donor_subset <- function (model_file_dir,
                                       training_data_samples_file,
                                       covariate_file,
                                       metacell_file,
                                       cell_type,
                                       source_gtf_file,
                                       target_gtf_file) {

    all_models = load_models(model_file_dir)
    model = all_models[all_models$cell_type == cell_type, ]

    #read in the training data so we can get the age distribution for this cell type.
    training_samples = read.table(training_data_samples_file,
                                  header = TRUE,
                                  stringsAsFactors = FALSE)
    training_samples = training_samples[training_samples$cell_type == cell_type, ]
    age_df_train = unique(training_samples[, c("donor", "age")])
    age_df_train$age = as.numeric(age_df_train$age) / 10

    #load the data and make a DGEList object
    dge <- load_mccarroll_metacells(covariate_file, metacell_file)
    dge$samples$age <- as.numeric(dge$samples$age) / 10

    #only run the controls
    # dge <- dge[, dge$samples$schizophrenia==FALSE]

    #because the gene symbols may differ between the training and test data, map them via ENSG IDs
    gene_symbol_map = gene_symbol_mapper(
        gene_symbols_source = model$feature,
        source_gtf_file = source_gtf_file,
        target_gtf_file = target_gtf_file
    )
    #map the model features to the target gene symbols
    #line up the map and model - they should be by default
    idx = match(model$feature, gene_symbol_map$source_gene_symbol)
    #map in the updated gene symbols
    model$feature_new = gene_symbol_map$target_gene_symbol[idx]
    #update the model to use the new gene symbols, as long as the new symbol is not NA
    model$feature <- ifelse(!is.na(model$feature_new),
                            model$feature_new,
                            model$feature)
    model$feature_new = NULL

    #filtering samples by library size - not predicting very small samples.
    r <- bican.mccarroll.differentialexpression::filter_by_libsize(
        dge,
        threshold_sd = 1.96,
        bins = 50,
        strTitlePrefix = cell_type
    )
    dge <- r$dge

    #filter to cpm cutoff of 1, but keep genes in the model too!
    r2 = bican.mccarroll.differentialexpression::plot_logCPM_density_quantiles(
        dge,
        cpm_cutoff = 1,
        logCPM_xlim = c(-5, 15),
        lower_quantile = 0.05,
        upper_quantile = 0.95,
        quantile_steps = 5,
        min_samples = 1,
        fraction_samples = 0.1
    )

    #if there are any features that aren't in both the filtered dge and the model, need to drop.
    genesMissingFromUnfilteredData = setdiff(model$feature, rownames(dge$counts))
    if (length(genesMissingFromUnfilteredData) > 0) {
        warning(
            paste0(
                "The following model features are missing from the unfiltered data and will be removed from the model: ",
                paste(genesMissingFromUnfilteredData, collapse = ", ")
            )
        )
        model = model[!model$feature %in% genesMissingFromUnfilteredData, ]
    }

    #keep the genes that are either significant figures originally, or pass the filtering threshold.
    genesMissingFromModel = setdiff(model$feature, rownames(r2$filtered_dge$counts))
    genesFinal = union (genesMissingFromModel, rownames(r2$filtered_dge$counts))
    dge <- dge[genesFinal, ]

    #just use the model coefficients that are non-zero for this analysis.
    model = model[model$coef != 0, ]

    logger::log_info(paste("Number of features in model", nrow(model)))
    #run the prediction on the original data
    result1 <- predict_age_from_dge(dge, model, prior.count = 1)
    plot_age_pred_external(result1,
                           cell_type,
                           strTitle = paste("Emi BA46 [default mean/sd]", cell_type, "\n"))

    #run predictions on controls only
    dge_sub <- dge[, dge$samples$schizophrenia == FALSE]
    result1 <- predict_age_from_dge(dge_sub, model, prior.count = 1)
    plot_age_pred_external(
        result1,
        cell_type,
        strTitle = paste("Emi BA46 [default mean/sd] controls", cell_type, "\n")
    )

    #run the prediction on the relearned mean/sd using all donors
    result2 <- predict_age_from_dge(dge,
                                    model,
                                    prior.count = 1,
                                    override_model_params_dge = dge)
    plot_age_pred_external(result2,
                           cell_type,
                           strTitle = paste("Emi BA46 [retrained mean/sd]", cell_type, "\n"))

    #run the prediction on the relearned mean/sd using all control donors
    dge_sub <- dge[, dge$samples$schizophrenia == FALSE]
    result2 <- predict_age_from_dge(dge_sub,
                                    model,
                                    prior.count = 1,
                                    override_model_params_dge = dge)
    plot_age_pred_external(
        result2,
        cell_type,
        strTitle = paste("Emi BA46 [retrained mean/sd] controls", cell_type, "\n")
    )

    #Learn a subset of donors that better match the training mean expression by optimization
    lcpm <- edgeR::cpm(dge$counts, log = TRUE, prior.count = 1)
    res <- greedy_match_donors(lcpm, model, max_drop_frac = 0.75)

    filtered_lcpm = lcpm[, res$kept_donors]
    p1 = plot_mean_expression_comparison(model, lcpm, cell_type)
    p2 = plot_mean_expression_comparison(model, filtered_lcpm, cell_type, strTitle = " (after donor filtering)")

    age_df_test = dge$samples[, c("donor", "age")]

    #animate_mean_matching_panels_cowplot(model, lcpm, res, "microglia", age_df_test,
    #                                     out_file = "/downloads/mean_match_panels.mp4", fps = 10)

    #the age distributions before/after age matching
    plot_age_distribution(
        age_df_train,
        age_df_test,
        strTitle = paste("Donor Age Distribution - all data", cell_type),
        name_train = "BICAN",
        name_test = "SNAP200 - all donors"
    )

    #plot the age distributions of the original data and optimized subset
    optimized_age_df = age_df = age_df_test[age_df_test$donor %in% res$kept_donors, ]
    plot_age_distribution(
        age_df_train,
        optimized_age_df,
        strTitle = paste("Donor Age Distribution - all data", cell_type),
        name_train = "BICAN",
        name_test = "SNAP200 - optimized subset"
    )

    #run the prediction on the relearned mean/sd
    dge_sub <- dge[, colnames(dge) %in% res$kept_donors]
    #train using the new model, but the original data
    result2 <- predict_age_from_dge(dge,
                                    model,
                                    prior.count = 1,
                                    override_model_params_dge = dge_sub)
    plot_age_pred_external(result2,
                           cell_type,
                           strTitle = paste("Emi BA46 [retrained mean/sd]", cell_type, "\n"))

    #what if we just predict controls?
    dge_controls <- dge[, dge$samples$schizophrenia == FALSE]
    result2 <- predict_age_from_dge(
        dge_controls,
        model,
        prior.count = 1,
        override_model_params_dge = dge_sub
    )
    plot_age_pred_external(
        result2,
        cell_type,
        strTitle = paste("Emi BA46 Controls [retrained mean/sd] controls", cell_type)
    )

    #plot the age distribution of the original and optimized results
    plot_age_distribution(
        age_df_train,
        age_df_test,
        strTitle = paste("Donor Age Distribution - all data", cell_type),
        name_train = "BICAN",
        name_test = "SNAP200 - all donors"
    )

    #plot the age distributions of the original data and optimized subset
    age_df_test_subset = age_df_test[age_df_test$donor %in% res$kept_donors, ]
    plot_age_distribution(
        age_df_train,
        age_df_test_subset,
        strTitle = paste("Donor Age Distribution - all data", cell_type),
        name_train = "BICAN",
        name_test = "SNAP200 - optimized subset"
    )

    ##########################
    # MATCH AGE DISTRIBUTION
    ##########################

    # learn a subset of donors that better match the training age distribution
    age_matched <- sample_test_ages_to_match_train(age_df_train, age_df_test, max_drop_frac = 0.75)
    print (age_matched$prob_plot)
    print (age_matched$ks_plot)

    #the age distributions before/after age matching
    plot_age_distribution(
        age_df_train,
        age_df_test,
        strTitle = paste("Donor Age Distribution - all data", cell_type),
        name_train = "BICAN",
        name_test = "SNAP200 - all donors"
    )

    #plot the age distributions of the original data and optimized subset
    plot_age_distribution(
        age_df_train,
        age_matched$subset,
        strTitle = paste("Donor Age Distribution - all data", cell_type),
        name_train = "BICAN",
        name_test = "SNAP200 - age matched subset"
    )

    # Run the age predictions on the age-matched subset
    #run the prediction on the relearned mean/sd
    dge_sub <- dge[, colnames(dge) %in% age_matched$subset$donor]
    #train using the new model, but the original data
    result2 <- predict_age_from_dge(dge,
                                    model,
                                    prior.count = 1,
                                    override_model_params_dge = dge_sub)
    plot_age_pred_external(result2,
                           cell_type,
                           strTitle = paste("Emi BA46 [age matched mean/sd]", cell_type, "\n"))

    #what if we just predict controls?
    dge_sub <- dge[, colnames(dge) %in% age_matched$subset$donor]
    dge_controls <- dge[, dge$samples$schizophrenia == FALSE]
    result3 <- predict_age_from_dge(
        dge_controls,
        model,
        prior.count = 1,
        override_model_params_dge = dge_sub
    )
    plot_age_pred_external(
        result3,
        cell_type,
        strTitle = paste("Emi BA46 Controls [age matched mean/sd]", cell_type)
    )


}

plot_age_distribution <- function (age_df_train,
                                   age_df_test,
                                   strTitle = NULL,
                                   name_train = "BICAN",
                                   name_test = "SNAP200") {
    age_df_train$dataset <- name_train
    age_df_test$dataset  <- name_test

    age_df <- rbind(age_df_train, age_df_test)

    if (is.null(strTitle))
        strTitle <- "Donor Age Distribution"

    # Ensure dataset is a factor with clean levels
    age_df$dataset <- factor(age_df$dataset)

    # Build a shared colour mapping whose *names* match dataset levels
    # Assume 2 datasets; extend this vector if you add more.
    dataset_colors <- c("steelblue", "firebrick")
    names(dataset_colors) <- levels(age_df$dataset)

    #Make R CMD CHECK happy
    age<-dataset<-uniform<-NULL

    # Histogram
    p1 <- ggplot(age_df, aes(x = age, fill = dataset)) +
        geom_histogram(alpha = 0.4,
                       position = "identity",
                       bins = 50) +
        labs(x = "Age [Decades]", y = "Count") +
        scale_fill_manual(values = dataset_colors) +
        theme_bw() +
        ggtitle(strTitle)

    p2 <- ggplot(age_df, aes(x = age, colour = dataset)) +
        stat_ecdf(size = 1) +
        geom_line(
            data = data.frame(age     = c(min(age_df$age), max(age_df$age)), uniform = c(0, 1)),
            aes(x = age, y = uniform),
            inherit.aes = FALSE,
            linetype = 2
        ) +
        labs(x = "Age [Decades]", y = "ECDF") +
        scale_colour_manual(values = dataset_colors) +
        theme_bw()

    cowplot::plot_grid(p1, p2, ncol = 1)

}

sample_test_ages_to_match_train <- function(age_df_train,
                                            age_df_test,
                                            max_drop_frac = 0.75,
                                            max_ks = NULL,
                                            n_reps_per_size = 1000,
                                            seed = 1) {
    # age_df_train, age_df_test: data frames with an "age" column
    # max_drop_frac: maximum fraction of test donors allowed to be dropped
    # max_ks: optional maximum acceptable KS distance.
    #         If not NULL, we:
    #           1) restrict to sizes with D <= max_ks
    #           2) within that restricted set, pick the "end of the decline"
    #              (rule below). If none satisfy D <= max_ks, fall back to
    #              using all sizes.
    # n_reps_per_size: number of random resamples per subset size.

    set.seed(seed)

    if (max_drop_frac < 0 || max_drop_frac >= 1) {
        stop("max_drop_frac must be in [0, 1).")
    }
    if (!("age" %in% names(age_df_train)) ||
        !("age" %in% names(age_df_test))) {
        stop("Both data frames must have an 'age' column.")
    }

    n_test <- nrow(age_df_test)
    if (n_test == 0L) {
        stop("age_df_test has no rows.")
    }

    # Minimum allowed subset size
    n_keep_min <- max(1L, ceiling((1 - max_drop_frac) * n_test))

    #### 1. Importance weights via density ratio ####

    dens_train <- density(age_df_train$age)
    dens_test  <- density(age_df_test$age, bw = dens_train$bw)

    p_train <- approx(dens_train$x,
                      dens_train$y,
                      xout = age_df_test$age,
                      rule = 2)$y
    p_test  <- approx(dens_test$x,
                      dens_test$y,
                      xout = age_df_test$age,
                      rule = 2)$y

    positive_test <- p_test[p_test > 0]
    if (length(positive_test) == 0L) {
        prob <- rep(1 / n_test, n_test)
    } else {
        tiny <- min(positive_test)
        p_test[p_test <= 0] <- tiny

        weights <- p_train / p_test
        weights[!is.finite(weights) | weights < 0] <- 0

        if (all(weights == 0)) {
            prob <- rep(1 / n_test, n_test)
        } else {
            prob <- weights / sum(weights)
        }
    }

    #### 2. Candidate subset sizes ####

    if (n_test <= 1000L) {
        step <- 1L
    } else {
        step <- max(1L, floor((n_test - n_keep_min) / 100L))
    }

    sizes <- seq(from = n_test,
                 to = n_keep_min,
                 by = -step)
    if (sizes[length(sizes)] != n_keep_min) {
        sizes <- c(sizes, n_keep_min)
    }

    #### 3. For each size: repeat sampling, keep best KS ####

    n_sizes  <- length(sizes)
    D_vec    <- numeric(n_sizes)
    m_vec    <- integer(n_sizes)
    idx_list <- vector("list", n_sizes)

    for (i in seq_len(n_sizes)) {
        m <- sizes[i]
        m_vec[i] <- m

        best_D_for_m   <- Inf
        best_idx_for_m <- NULL

        for (r in seq_len(n_reps_per_size)) {
            idx <- sample.int(
                n_test,
                size = m,
                replace = FALSE,
                prob = prob
            )
            ks_res <- suppressWarnings(stats::ks.test(age_df_test$age[idx], age_df_train$age))
            D <- as.numeric(ks_res$statistic)

            if (D < best_D_for_m) {
                best_D_for_m   <- D
                best_idx_for_m <- idx
            }
        }

        D_vec[i]      <- best_D_for_m
        idx_list[[i]] <- best_idx_for_m
    }

    #### 4. Selection rule: first increase in KS when going to smaller sizes ####
    # m_vec is constructed so that i = 1 has the largest subset size (n_test),
    # and i increases as subset size decreases. So we scan i = 2..n and look
    # for the first i where D_vec[i] > D_vec[i - 1]; we then choose i - 1.

    choose_end_of_decline <- function(idx_seq) {
        # idx_seq must be ordered so that m_vec[idx_seq] is decreasing.
        if (length(idx_seq) == 1L)
            return(idx_seq)

        for (k in seq(2L, length(idx_seq))) {
            i_prev <- idx_seq[k - 1L]
            i_cur  <- idx_seq[k]
            if (D_vec[i_cur] > D_vec[i_prev]) {
                return(i_prev)
            }
        }
        # If we never see an increase, pick the smallest size (last index)
        idx_seq[length(idx_seq)]
    }

    chosen_i <- NA_integer_

    if (!is.null(max_ks)) {
        ok <- which(D_vec <= max_ks)
        if (length(ok) > 0L) {
            chosen_i <- choose_end_of_decline(ok)
        }
    }

    if (is.na(chosen_i)) {
        all_idx <- seq_len(n_sizes)
        chosen_i <- choose_end_of_decline(all_idx)
    }

    chosen_idx <- idx_list[[chosen_i]]
    chosen_m   <- m_vec[chosen_i]
    chosen_D   <- D_vec[chosen_i]

    #### 5a. ggplot2 KS-distance vs subset-size diagnostic ####

    df_ks <- data.frame(subset_size = m_vec, ks_distance = D_vec)

    strTitle = paste("KS distance vs retained subset size (optimal n=[",
                     chosen_m,
                     "])",
                     sep = "")

    #Make R CMD CHECK happy
    subset_size<-ks_distance<-label<-NULL

    ks_plot <- ggplot(df_ks, aes(x = subset_size, y = ks_distance)) +
        geom_line() +
        geom_point() +
        geom_point(
            data = data.frame(subset_size = chosen_m, ks_distance = chosen_D),
            aes(x = subset_size, y = ks_distance),
            color = "red",
            size = 3
        ) +
        geom_text(
            data = data.frame(subset_size = chosen_m, ks_distance = chosen_D),
            aes(label = "chosen"),
            hjust = -0.2,
            vjust = 0,
            size = 3
        ) +
        labs(x = "Number of test donors kept", y = "KS distance to train distribution", title = strTitle) +
        theme_bw()

    #### 5b. ggplot2 sampling probability vs age ####

    df_prob <- data.frame(age  = age_df_test$age, prob = prob)

    #Make R CMD CHECK happy
    age<-prob<-NULL

    prob_plot <- ggplot(df_prob, aes(x = age, y = prob)) +
        geom_point() +
        labs(x = "SNAP200 donor age [decades]", y = "sampling probability", title = "Sampling probability vs test donor age") +
        theme_bw()

    #### Return results ####
    list(
        subset       = age_df_test[chosen_idx, , drop = FALSE],
        ks_chosen    = chosen_D,
        size_chosen  = chosen_m,
        sizes        = m_vec,
        ks_values    = D_vec,
        prob         = prob,
        ks_plot      = ks_plot,
        prob_plot    = prob_plot
    )
}






# lcpm_ext: genes x donors log-CPM matrix for external dataset
# model: data.frame with at least columns 'feature' and 'mean_train'
# max_drop_frac: stop if more than this fraction of donors would be removed
# tol: minimum improvement in total absolute error required to drop a donor
# lcpm_ext=lcpm
greedy_match_donors <- function(lcpm_ext,
                                model,
                                max_drop_frac = 0.5,
                                tol = 1e-8) {

    stopifnot(is.matrix(lcpm_ext) || is.data.frame(lcpm_ext))
    lcpm_ext <- as.matrix(lcpm_ext)

    # intersect genes
    feats <- intersect(model$feature, rownames(lcpm_ext))
    if (length(feats) == 0L)
        stop("No overlapping genes between model and lcpm_ext")

    train_mu <- model$mean_train[match(feats, model$feature)]
    ext_mat  <- lcpm_ext[feats, , drop = FALSE]

    donors <- colnames(ext_mat)
    n_donors <- length(donors)
    if (is.null(donors))
        donors <- paste0("donor", seq_len(n_donors))

    # initial sums / means over all donors
    kept      <- rep(TRUE, n_donors)
    cur_sums  <- rowSums(ext_mat)
    n_kept    <- n_donors
    cur_mu    <- cur_sums / n_kept
    cur_err   <- sum(abs(cur_mu - train_mu))

    history <- data.frame(
        step            = 0L,
        n_kept          = n_kept,
        total_abs_error = cur_err,
        dropped_donor   = NA_character_,
        stringsAsFactors = FALSE
    )

    max_drops <- floor(n_donors * max_drop_frac)

    for (step in seq_len(max_drops)) {
        best_err <- Inf
        best_j   <- NA_integer_

        # try dropping each currently kept donor
        for (j in which(kept)) {
            new_sums <- cur_sums - ext_mat[, j]
            new_n    <- n_kept - 1L
            mu_new   <- new_sums / new_n
            err_new  <- sum(abs(mu_new - train_mu))

            if (err_new < best_err) {
                best_err <- err_new
                best_j   <- j
            }
        }

        # stop if no improvement
        if (is.na(best_j) || best_err > cur_err - tol)
            break

        # commit the best drop
        kept[best_j] <- FALSE
        cur_sums     <- cur_sums - ext_mat[, best_j]
        n_kept       <- n_kept - 1L
        cur_err      <- best_err

        history <- rbind(
            history,
            data.frame(
                step            = step,
                n_kept          = n_kept,
                total_abs_error = cur_err,
                dropped_donor   = donors[best_j],
                stringsAsFactors = FALSE
            )
        )
    }

    kept_donors    <- donors[kept]
    dropped_donors <- donors[!kept]

    list(
        kept_donors        = kept_donors,
        dropped_donors     = dropped_donors,
        history            = history,
        initial_error      = history$total_abs_error[1],
        final_error        = utils::tail(history$total_abs_error, 1),
        improvement        = history$total_abs_error[1] - utils::tail(history$total_abs_error, 1)
    )
}



# --------------------------------------------------------------------
# Helper: build one frame's scatter data
# --------------------------------------------------------------------
# .mean_frame <- function(model, lcpm, keep_donors, step_int, step_lab, n_kept) {
#     lcpm_sub <- lcpm[, keep_donors, drop = FALSE]
#     feats <- intersect(model$feature, rownames(lcpm_sub))
#     if (length(feats) == 0L) {
#         return(data.frame(
#             feature    = character(0),
#             train_mean = numeric(0),
#             ext_mean   = numeric(0),
#             coef_sign  = character(0),
#             step_int   = integer(0),
#             step_lab   = character(0),
#             n_kept     = integer(0),
#             stringsAsFactors = FALSE
#         ))
#     }
#     idx <- match(feats, model$feature)
#     train_mean <- model$mean_train[idx]
#     ext_mean   <- rowMeans(lcpm_sub[feats, , drop = FALSE])
#     coef_sign  <- ifelse(model$coef[idx] >= 0, "Positive", "Negative")
#
#     data.frame(
#         feature    = feats,
#         train_mean = train_mean,
#         ext_mean   = ext_mean,
#         coef_sign  = coef_sign,
#         step_int   = step_int,
#         step_lab   = step_lab,
#         n_kept     = n_kept,
#         stringsAsFactors = FALSE
#     )
# }
#
# # --------------------------------------------------------------------
# # Main animation function
# # --------------------------------------------------------------------
# animate_mean_matching_panels_cowplot <- function(
        #         model, lcpm, res, cell_type, age_df,
#         out_file = "mean_match_panels.mp4",
#         fps = 2,
#         width = 1900, height = 1080,
#         dpi = 300) {
#
#     history <- res$history
#
#     # Enforce ascending step order for consistent playback
#     ord <- order(history$step)
#     history <- history[ord, , drop = FALSE]
#
#     donors_all <- colnames(lcpm)
#
#     # Age table must cover all donors referenced by lcpm
#     if (!"donor" %in% colnames(age_df)) age_df$donor <- rownames(age_df)
#     if (!"age" %in% colnames(age_df)) stop("age_df must have an 'age' column")
#     age_df$age <- as.numeric(age_df$age)
#     missing_age <- setdiff(donors_all, age_df$donor)
#     if (length(missing_age)) {
#         stop(sprintf("age_df missing %d donor(s), e.g. %s",
#                      length(missing_age), paste(utils::head(missing_age, 5L), collapse = ", ")))
#     }
#
#     # Reconstruct donor inclusion list per (sorted) step
#     keep_list <- vector("list", nrow(history))
#     keep <- donors_all
#     keep_list[[1]] <- keep
#     if (nrow(history) > 1) {
#         for (i in 2:nrow(history)) {
#             d_drop <- history$dropped_donor[i]
#             if (!is.na(d_drop)) keep <- setdiff(keep, d_drop)
#             keep_list[[i]] <- keep
#         }
#     }
#
#     # Stable state labels used across all panels
#     state_levels <- sprintf("Step %d | kept %d", history$step, history$n_kept)
#
#     # Assemble per-step scatter + summary (align with .mean_frame signature: step_int, step_lab)
#     pieces <- lapply(seq_along(keep_list), function(i) {
#         if (length(keep_list[[i]]) < 1L) {
#             stop(sprintf("No donors kept at step %d; animation requires >=1 donor per step.", history$step[i]))
#         }
#         sc <- .mean_frame(
#             model       = model,
#             lcpm        = lcpm,
#             keep_donors = keep_list[[i]],
#             step_int    = history$step[i],
#             step_lab    = state_levels[i],
#             n_kept      = history$n_kept[i]
#         )
#         err <- if (nrow(sc)) mean(abs(sc$ext_mean - sc$train_mean), na.rm = TRUE) else NA_real_
#         list(
#             scatter = sc,
#             summary = data.frame(
#                 step_int = history$step[i],
#                 step_lab = state_levels[i],
#                 n_kept   = history$n_kept[i],
#                 frac     = history$n_kept[i]/ncol(lcpm),
#                 error_ma = err,
#                 stringsAsFactors = FALSE
#             )
#         )
#     })
#
#     scatters <- if (length(pieces)) do.call(rbind, lapply(pieces, `[[`, "scatter")) else NULL
#     sums     <- if (length(pieces)) do.call(rbind, lapply(pieces, `[[`, "summary")) else NULL
#
#     # Ordered state factor for strict play order in transition_manual()
#     scatters$state <- factor(scatters$step_lab, levels = state_levels)
#     sums$state     <- factor(sums$step_lab,     levels = state_levels)
#
#     # Panel B data: ages per step for kept donors (must exist for every step)
#     ages_list <- lapply(seq_along(keep_list), function(i) {
#         donors_i <- keep_list[[i]]
#         m <- match(donors_i, age_df$donor)
#         ages <- age_df$age[m]
#         if (anyNA(ages)) {
#             miss <- donors_i[is.na(ages)]
#             stop(sprintf("Missing ages for %d donor(s) at step %d, e.g. %s",
#                          length(miss), history$step[i], paste(utils::head(miss, 5L), collapse = ", ")))
#         }
#         data.frame(
#             state    = factor(rep(state_levels[i], length(ages)), levels = state_levels),
#             step_int = rep(history$step[i], length(ages)),
#             age      = ages,
#             stringsAsFactors = FALSE
#         )
#     })
#     ages_long <- do.call(rbind, ages_list)
#
#     # Global axis limits for scatter
#     lim_all <- range(c(scatters$train_mean, scatters$ext_mean), na.rm = TRUE)
#
#     # -------------------------------
#     # Panel A: mean expression scatter
#     # -------------------------------
#     p_scatter <- ggplot2::ggplot(
#         scatters, ggplot2::aes(x = train_mean, y = ext_mean, color = coef_sign, group = feature)
#     ) +
#         ggplot2::geom_point(alpha = 0.7, size = 2) +
#         ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
#         ggplot2::scale_x_continuous("Model mean expression (log2)", limits = lim_all, expand = ggplot2::expansion(mult = 0.02)) +
#         ggplot2::scale_y_continuous("External mean expression (log2)", limits = lim_all, expand = ggplot2::expansion(mult = 0.02)) +
#         ggplot2::scale_color_manual(values = c("Negative" = "orange", "Positive" = "steelblue")) +
#         ggplot2::labs(
#             title = sprintf("%s mean expression comparison", cell_type),
#             subtitle = "{current_frame}"
#         ) +
#         ggplot2::theme_classic(base_size = 8) +
#         ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = "top") +
#         gganimate::transition_manual(scatters$state)
#
#     # ----------------------------------------------
#     # Panel B: histogram of donor ages (replaces frac)
#     # ----------------------------------------------
#     age_rng <- range(ages_long$age, na.rm = TRUE)
#     bw <- max(0.25, (age_rng[2] - age_rng[1]) / 20)
#     boundary <- floor(age_rng[1])
#
#     p_age <- ggplot2::ggplot(ages_long, ggplot2::aes(x = age)) +
#         ggplot2::geom_histogram(binwidth = bw, boundary = boundary, closed = "right", fill = "grey40") +
#         ggplot2::ggtitle("Donor Age Distribution") +
#         ggplot2::theme_classic(base_size = 10) +
#         ggplot2::theme(
#             axis.title.x = ggplot2::element_blank(),
#             axis.title.y = ggplot2::element_blank(),
#             axis.text.x  = ggplot2::element_text(size = 8),
#             axis.text.y  = ggplot2::element_text(size = 8),
#             plot.title   = ggplot2::element_text(hjust = 0.5, size = 10, face = "bold")
#         ) +
#         gganimate::transition_manual(ages_long$state)
#
#     # -----------------------
#     # Panel C: error trace
#     # -----------------------
#     sums$step_f <- factor(sums$step_int, levels = history$step)
#     p_err <- ggplot2::ggplot(sums, ggplot2::aes(x = as.integer(sums$step_f), y = error_ma)) +
#         ggplot2::geom_line(linewidth = 0.6, na.rm = TRUE) +
#         ggplot2::geom_point(size = 3, na.rm = TRUE) +
#         ggplot2::ggtitle("Error Trace") +
#         ggplot2::coord_cartesian(ylim = c(0, max(sums$error_ma, na.rm = TRUE))) +  # force baseline = 0
#         ggplot2::theme_classic(base_size = 12) +
#         ggplot2::theme(
#             axis.title.x = ggplot2::element_blank(),
#             axis.text.x  = ggplot2::element_blank(),
#             axis.ticks.x = ggplot2::element_blank(),
#             axis.title.y = ggplot2::element_text(size = 10),
#             axis.text.y  = ggplot2::element_text(size = 8),
#             plot.title   = ggplot2::element_text(hjust = 0.5, size = 11, face = "bold"),
#             plot.margin  = ggplot2::margin(2, 4, 2, 4)
#         ) +
#         gganimate::transition_manual(sums$state)
#
#
#     # Render each panel
#     render_panel <- function(p) {
#         gganimate::animate(
#             p,
#             fps = fps,
#             renderer = gganimate::magick_renderer(),
#             width = width, height = height, res = dpi
#         )
#     }
#     a_frames <- render_panel(p_scatter)
#     b_frames <- render_panel(p_age)
#     c_frames <- render_panel(p_err)
#
#     # Combine frames using cowplot
#     tdir <- tempfile("frames_")
#     dir.create(tdir)
#     nF <- length(a_frames)
#     png_paths <- character(nF)
#
#     for (i in seq_len(nF)) {
#         combined <- cowplot::plot_grid(
#             cowplot::ggdraw() +
#                 ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", color = NA)) +
#                 cowplot::draw_image(a_frames[i], interpolate = TRUE),
#             cowplot::plot_grid(
#                 cowplot::ggdraw() +
#                     ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", color = NA)) +
#                     cowplot::draw_image(b_frames[i], interpolate = TRUE),
#                 cowplot::ggdraw() +
#                     ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", color = NA)) +
#                     cowplot::draw_image(c_frames[i], interpolate = TRUE),
#                 ncol = 1,
#                 rel_heights = c(1, 1)  # was (1, 1)
#             ),
#             ncol = 2,
#             rel_widths = c(3.0, 1.0)      # was (2.4, 1)
#         )
#
#         ggplot2::ggsave(
#             filename = (png_paths[i] <- file.path(tdir, sprintf("frame_%04d.png", i))),
#             plot = combined, width = width/dpi, height = height/dpi, dpi = dpi, bg = "white"
#         )
#     }
#
#
#     # Encode MP4
#     av::av_encode_video(png_paths, output = out_file, framerate = fps)
#
#
#     invisible(out_file)
# }


##########################
# END MINI-TASK : Optimize donor subset to better match mean expression
##########################


#This has a bunch of junk for mapping between ENSG builds, which is imperfect at best.
predict_external_data <- function (model_file_dir, cell_type, covariate_file,
                                   metacell_file,
                                   source_gtf_file,
                                   target_gtf_file) {

    all_models = load_models(model_file_dir)
    model = all_models[all_models$cell_type == cell_type, ]
    model = model[model$coef != 0, ]
    #load the data and make a DGEList object
    dge <- load_mccarroll_metacells(covariate_file, metacell_file)
    dge$samples$age <- as.numeric(dge$samples$age) / 10

    #only run the controls
    dge <- dge[, dge$samples$schizophrenia == FALSE]

    #because the gene symbols may differ between the training and test data, map them via ENSG IDs
    gene_symbol_map = gene_symbol_mapper(
        gene_symbols_source = model$feature,
        source_gtf_file = source_gtf_file,
        target_gtf_file = target_gtf_file
    )
    #map the model features to the target gene symbols
    #line up the map and model - they should be by default
    idx = match(model$feature, gene_symbol_map$source_gene_symbol)
    #map in the updated gene symbols
    model$feature_new = gene_symbol_map$target_gene_symbol[idx]
    #update the model to use the new gene symbols, as long as the new symbol is not NA
    model$feature <- ifelse(!is.na(model$feature_new),
                            model$feature_new,
                            model$feature)

    #what features are missing?
    df <- data.frame(coef = model$coef,
                     feature_missing = is.na(model$feature_new))

    strTitle = paste(
        "Model Coefficients\n",
        cell_type,
        " missing features ",
        length(which(df$feature_missing)),
        " / ",
        dim(df)[1],
        sep = ""
    )

    #Make R CMD CHECK happy
    coef <- feature_missing <- NULL

    ggplot(df, aes(
        x = seq_along(coef),
        y = coef,
        fill = feature_missing
    )) +
        geom_col() +
        labs(x = NULL, y = "Coefficient") +
        ggtitle(strTitle) +
        theme_classic()


    #filtering samples by library size - not predicting very small samples.
    r <- bican.mccarroll.differentialexpression::filter_by_libsize(
        dge,
        threshold_sd = 1.96,
        bins = 50,
        strTitlePrefix = cell_type
    )
    dge <- r$dge

    #run the model
    # result <- predict_age_from_dge(dge, model, prior.count = 1)
    # plot_age_pred_external(result,
    #                        cell_type,
    #                        strTitle = paste("Emi BA46 Controls", cell_type, "\n"))
    # plot_age_pred_external(result, cell_type, strTitle=paste("Emi BA46", cell_type, "\n"))

    #what about a dumb hack?  Replace the model mean/sd with the observed data.
    #probably terrible idea.
    # result2 <- predict_age_from_dge(dge,
    #                                 model,
    #                                 prior.count = 1,
    #                                 override_model_params_dge = dge)
    # plot_age_pred_external(
    #     result3,
    #     cell_type,
    #     strTitle = paste("Emi BA46 [retrained mean/sd] Controls", cell_type, "\n")
    # )

    #compare the model coefficients and means to the test expression data.
    plots = compare_external_data_to_model (dge, model, cell_type)
    for (p in plots)
        print (p)

}



compare_external_data_to_model <- function (dge, model, cell_type) {
    #are some of the features not correlated with age?
    dge_new <- dge
    lcpm <- edgeR::cpm(dge_new, log = TRUE, prior.count = 1)  # genes x samples
    cor_list <- lapply(model$feature, function(gene) {
        if (gene %in% rownames(lcpm)) {
            expr <- lcpm[gene, ]
            age <- dge_new$samples$age
            r <- cor(expr, age, use = "pairwise.complete.obs")
            return (data.frame(
                gene = gene,
                cor = r,
                stringsAsFactors = FALSE
            ))
        } else {
            return (data.frame(
                gene = gene,
                cor = NA,
                stringsAsFactors = FALSE
            ))
        }
    })
    cor_df = do.call(rbind, cor_list)
    idx = match(model$feature, cor_df$gene)
    model$feature_cor = cor_df$cor[idx]


    #compare the correlation and sign of the external data to the model.
    p1 = plot_coef_vs_feature_cor(model, cell_type = cell_type, nonzero_only = FALSE)

    #check the mean and SD of each gene in the model vs training data
    p2 = plot_mean_expression_comparison(model, lcpm, cell_type)

    #what if we split the model genes into test/train many times, and compute the difference in the mean expression for each feature?
    #how often is the observed shift larger than the sampling shift?

    result = list(p1, p2)
    return (result)


}

plot_coef_vs_feature_cor <- function(model,
                                     cell_type = NULL,
                                     nonzero_only = FALSE) {
    stopifnot(is.data.frame(model), all(c("coef", "feature_cor") %in% colnames(model)))

    df <- model
    if (nonzero_only)
        df <- df[is.finite(df$coef) & df$coef != 0, , drop = FALSE]

    idx <- is.finite(df$coef) & is.finite(df$feature_cor)
    if (!any(idx))
        stop("No finite pairs to plot.")

    r_val <- suppressWarnings(stats::cor(df$coef[idx], df$feature_cor[idx], use = "pairwise.complete.obs"))
    agree <- sign(df$coef[idx]) == sign(df$feature_cor[idx])
    frac_same <- mean(agree)

    title_str <- sprintf(
        "%s%s r = %.3f, sign agreement = %.1f%%",
        if (!is.null(cell_type))
            paste0(cell_type, ": ")
        else
            "",
        "",
        r_val,
        100 * frac_same
    )

    df$sign_agreement <- factor(
        ifelse(
            is.finite(df$coef) & is.finite(df$feature_cor) &
                sign(df$coef) == sign(df$feature_cor),
            "Same",
            "Different"
        ),
        levels = c("Same", "Different")
    )

    #Make R CMD CHECK happy
    coef <- feature_cor <- sign_agreement <- NULL

    ggplot2::ggplot(df[idx, , drop = FALSE],
                    ggplot2::aes(x = coef, y = feature_cor, color = sign_agreement)) +
        ggplot2::geom_point(alpha = 0.7, size = 2) +
        ggplot2::geom_abline(intercept = 0,
                             slope = 1,
                             color = "red") +
        ggplot2::geom_hline(
            yintercept = 0,
            linetype = "dashed",
            color = "grey50"
        ) +
        ggplot2::geom_vline(
            xintercept = 0,
            linetype = "dashed",
            color = "grey50"
        ) +
        ggplot2::scale_color_manual(values = c(
            "Same" = "steelblue",
            "Different" = "orange"
        ),
        name = "Sign match") +
        ggplot2::labs(title = title_str, x = "Model coefficient", y = "Gene correlation with age") +
        ggplot2::theme_classic(base_size = 12) +
        ggplot2::theme(legend.position = "top")
}


plot_mean_expression_comparison <- function(model, lcpm, cell_type, strTitle =
                                                NULL) {
    stopifnot(is.data.frame(model), is.matrix(lcpm))

    # Compute external mean/sd
    model$ext_mean <- rowMeans(lcpm[model$feature, , drop = FALSE], na.rm = TRUE)
    model$ext_sd   <- apply(lcpm[model$feature, , drop = FALSE], 1, sd)
    model$mean_shift <- model$ext_mean - model$mean_train
    model$sd_ratio   <- model$ext_sd / model$sd_train

    # Prep for plotting
    df <- subset(model, is.finite(mean_train) & is.finite(ext_mean))
    df$coef_sign <- factor(
        sign(df$coef),
        levels = c(-1, 0, 1),
        labels = c("Negative", "Zero", "Positive")
    )

    xy_lim <- range(c(df$mean_train, df$ext_mean), na.rm = TRUE)
    r_val  <- suppressWarnings(stats::cor(df$mean_train, df$ext_mean))

    # Test: mean shift vs. coefficient sign
    df_test <- subset(df, coef_sign %in% c("Negative", "Positive"))
    wtest <- suppressWarnings(wilcox.test(mean_shift ~ coef_sign, data = df_test))
    p_val <- signif(wtest$p.value, 3)

    if (!is.null(strTitle)) {
        plot_title <- strTitle
    } else {
        plot_title <- paste0(cell_type, " mean expression comparison")
    }

    #Make R CMD CHECK happy
    mean_train <- ext_mean <- coef_sign <- NULL

    ggplot2::ggplot(df,
                    ggplot2::aes(x = mean_train, y = ext_mean, color = coef_sign)) +
        ggplot2::geom_abline(
            intercept = 0,
            slope = 1,
            linetype = "dashed",
            color = "red"
        ) +
        ggplot2::geom_point(alpha = 0.7, size = 1.8) +
        ggplot2::scale_x_continuous(limits = xy_lim,
                                    expand = ggplot2::expansion(mult = 0.02)) +
        ggplot2::scale_y_continuous(limits = xy_lim,
                                    expand = ggplot2::expansion(mult = 0.02)) +
        ggplot2::scale_color_manual(values = c(
            "Negative" = "#D55E00",
            "Zero" = "grey70",
            "Positive" = "#0072B2"
        )) +
        ggplot2::labs(
            title = plot_title,
            subtitle = sprintf(
                "Pearson r = %.3f | Wilcoxon p = %.3g (external - model means partitioned by coefficient sign)",
                r_val,
                p_val
            ),
            x = "Model mean expression (log2)",
            y = "External mean expression (log2)",
            color = "Coefficient sign"
        ) +
        ggplot2::theme_classic(base_size = 12) +
        ggplot2::theme(
            legend.position = "top",
            legend.title = ggplot2::element_text(size = 9),
            legend.text  = ggplot2::element_text(size = 9)
        )
}



plot_age_pred_external <- function(result,
                                   cell_type,
                                   strTitle = paste("Emi BA46", cell_type, "\n")) {
    metrics = compute_age_metrics(result$pred, result$age)
    label <- paste0(
        "r = ",
        sprintf("%.3f", metrics$r),
        "\nMedAE = ",
        sprintf("%.3f", metrics$median_abs_error),
        "\nMAE = ",
        sprintf("%.3f", metrics$mean_abs_error)
    )

    p <- plot_age_predictions(result, cell_type, titleStr = strTitle)

    # Determine plot limits for corner placement
    plot_limits <- layer_scales(p)
    x_pos <- plot_limits$x$range$range[1]
    y_pos <- plot_limits$y$range$range[2]

    p <- p +
        annotate(
            "text",
            x = x_pos,
            y = y_pos,
            label = label,
            hjust = 0,
            vjust = 1,
            size = 6,
            color = "black"
        ) +
        guides(color = "none")

    print (p)
}




load_mccarroll_metacells <- function (covariate_file, metacell_file) {
    md = read.table(
        covariate_file,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE
    )
    e = read.table(
        metacell_file,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE
    )
    # gene column can be either gene_symbol or GENE
    if ("gene_symbol" %in% colnames (e)) {
        rownames (e) = e$gene_symbol
        e$gene_symbol = NULL
    } else if ("GENE" %in% colnames (e)) {
        rownames (e) = e$GENE
        e$GENE = NULL
    } else {
        stop("No gene symbol column found in metacell file.")
    }

    counts = as.matrix(e)
    #intersect
    common_donors = intersect(colnames(counts), md$donor)
    counts = counts[, common_donors]
    md = md[match(common_donors, md$donor), ]
    dge <- edgeR::DGEList(counts = counts, samples = md)
    return (dge)
}


# Build Jaccard matrix of selected features per cell type
jaccard_by_celltype <- function(all_models, coef_thresh = 0) {
    stopifnot(all(c("cell_type", "feature", "coef") %in% names(all_models)))

    # binary incidence matrix: features x cell_types
    sel <- abs(all_models$coef) > coef_thresh
    df  <- all_models[sel, c("feature", "cell_type")]
    df  <- unique(df)  # one row per (feature, cell_type)

    feats <- sort(unique(df$feature))
    cts   <- sort(unique(df$cell_type))
    A <- matrix(
        FALSE,
        nrow = length(feats),
        ncol = length(cts),
        dimnames = list(feats, cts)
    )
    A[cbind(match(df$feature, feats), match(df$cell_type, cts))] <- TRUE

    # Jaccard: J = (A' A) / (n_i + n_j - A' A)
    XtX <- crossprod(A)                           # |intersection|
    n   <- matrix(colSums(A), ncol = ncol(A), nrow = ncol(A))
    J   <- as.matrix(XtX / (n + t(n) - XtX))
    diag(J) <- 1
    J
}

# Heatmap with optional numeric annotation
# plot_jaccard_heatmap <- function(J,
#                                  title = "Feature overlap (Jaccard Index)",
#                                  annotate_cells = TRUE,
#                                  cluster = TRUE) {
#     stopifnot(
#         requireNamespace("ComplexHeatmap", quietly = TRUE),
#         requireNamespace("circlize", quietly = TRUE),
#         requireNamespace("grid", quietly = TRUE)
#     )
#
#     col_fun <- circlize::colorRamp2(c(0, 0.5, 1), c("#f7fbff", "#6baed6", "#08306b"))
#
#     cf <- if (isTRUE(annotate_cells)) {
#         function(j, i, x, y, w, h, fill) {
#             grid::grid.text(sprintf("%.2f", J[i, j]), x, y, gp = grid::gpar(fontsize = 10))
#         }
#     } else
#         NULL
#
#     ComplexHeatmap::Heatmap(
#         J,
#         name = "Jaccard",
#         col = col_fun,
#         cluster_rows = cluster,
#         cluster_columns = cluster,
#         show_row_dend = cluster,
#         show_column_dend = cluster,
#         row_names_gp = grid::gpar(fontsize = 10),
#         column_names_gp = grid::gpar(fontsize = 10),
#         heatmap_legend_param = list(at = c(0, 0.5, 1)),
#         cell_fun = cf,
#         column_title = title
#     )
# }

plot_age_predictions_by_region <- function(results,
                                           cellType,
                                           titleStr = NULL,
                                           facet_col = "region",
                                           metrics_df = NULL) {
    if (is.null(titleStr))
        titleStr <- paste("Region specific predictions for", cellType)

    # facet variable
    df <- results[results$cell_type == cellType, ]
    facet_var <- df[[facet_col]]
    df2 <- cbind(df, facet_var = facet_var)

    # symmetric limits
    rng <- range(c(df2$age, df2$pred), na.rm = TRUE)
    rng <- c(rng[1] * 0.9, rng[2] * 1.1)

    #Make R CMD CHECK happy
    age <- pred <- x <- y <- label <- NULL

    p <- ggplot(df2, aes(x = age, y = pred)) +
        geom_abline(
            intercept = 0,
            slope = 1,
            color = "black",
            linetype = "dashed"
        ) +
        geom_point(aes(color = facet_var), size = 2, alpha = 0.6) +
        geom_smooth(
            method = "lm",
            formula = y ~ x,
            se = FALSE,
            linewidth = 0.6,
            color = "red"
        ) +
        facet_wrap( ~ facet_var) +
        guides(color = "none") +
        labs(title = titleStr, x = "Chronological Age", y = "Predicted Age") +
        scale_x_continuous(limits = rng) +
        scale_y_continuous(limits = rng) +
        theme_classic(base_size = 12)

    # optional per-facet metrics annotation
    if (!is.null(metrics_df)) {
        keep <- metrics_df$cell_type == cellType &
            metrics_df[[facet_col]] %in% unique(df2$facet_var)
        if (any(keep)) {
            ann <- metrics_df[keep, c(facet_col,
                                      "r",
                                      "median_abs_error",
                                      "mean_abs_error")]
            ann$label <- paste0(
                "r = ",
                sprintf("%.3f", ann$r),
                "\nMedAE = ",
                sprintf("%.3f", ann$median_abs_error),
                "\nMAE = ",
                sprintf("%.3f", ann$mean_abs_error)
            )
            ann_df <- data.frame(
                facet_var = ann[[facet_col]],
                label = ann$label,
                x = -Inf,
                y = Inf,
                stringsAsFactors = FALSE
            )
            p <- p +
                geom_text(
                    data = ann_df,
                    aes(
                        x = x,
                        y = y,
                        label = label,
                        color = facet_var
                    ),
                    hjust = -0.1,
                    vjust = 1.1,
                    size = 3.2,
                    inherit.aes = FALSE
                ) +
                guides(color = "none")
        }
    }

    p
}




# if there's more than one observation for the grouping, aggregate using mean or median
compute_residual_matrix <- function(results,
                                    group_cols = c("donor", "region"),
                                    cellType = NULL,
                                    agg = c("mean", "median"),
                                    drop_empty = TRUE) {
    agg <- match.arg(agg)

    stopifnot(all(c("age", "pred") %in% names(results)))
    stopifnot(length(group_cols) == 2, all(group_cols %in% names(results)))

    df <- if (is.null(cellType))
        results
    else
        results[results$cell_type == cellType, , drop = FALSE]
    if (!nrow(df))
        stop("No rows after filtering")

    # residuals
    resid <- df$pred - df$age

    # aggregator
    fun <- if (agg == "mean")
        function(x)
            mean(x, na.rm = TRUE)
    else
        function(x)
            median(x, na.rm = TRUE)

    # donors x second group (e.g., region or cell_type)
    M <- tapply(resid, list(df[[group_cols[1]]], df[[group_cols[2]]]), fun)
    M <- as.matrix(M)

    if (drop_empty) {
        keep_r <- rowSums(is.finite(M)) > 0
        keep_c <- colSums(is.finite(M)) > 0
        M <- M[keep_r, keep_c, drop = FALSE]
    }
    if (!nrow(M) ||
        !ncol(M))
        stop("Residual matrix empty after dropping empty groups")

    M
}

# 2. Scatter plots with annotated correlations
plot_residual_pair_scatter <- function(res_mat, cellType) {
    regs <- colnames(res_mat)
    pairs <- utils::combn(regs, 2, simplify = FALSE)
    df_list <- lapply(pairs, function(p) {
        x <- res_mat[, p[1]]
        y <- res_mat[, p[2]]
        keep <- is.finite(x) & is.finite(y)
        r <- cor(x[keep], y[keep])
        data.frame(
            x = x[keep],
            y = y[keep],
            pair = paste(p[1], "vs", p[2]),
            corr = r,
            stringsAsFactors = FALSE
        )
    })
    D <- do.call(rbind, df_list)

    rng <- range(c(D$x, D$y), na.rm = TRUE)
    rng <- c(rng[1] * 0.9, rng[2] * 1.1)

    title = "Age Prediction residuals (predicted - actual)"
    title = paste(cellType, title, sep = "\n")

    #Make R CMD CHECK happy
    x <- y <- pair <- corr <- NULL

    ggplot(D, aes(x = x, y = y)) +
        geom_abline(
            intercept = 0,
            slope = 1,
            color = "black",
            linetype = "dashed"
        ) +
        geom_point(size = 2,
                   alpha = 0.7,
                   color = "steelblue") +
        geom_smooth(
            method = "lm",
            formula = y ~ x,
            se = FALSE,
            linewidth = 0.6,
            color = "red"
        ) +
        facet_wrap( ~ pair) +
        geom_text(
            data = D |> dplyr::distinct(pair, corr),
            aes(
                x = -Inf,
                y = Inf,
                label = sprintf("r = %.2f", corr)
            ),
            hjust = -0.1,
            vjust = 1.2,
            size = 3.2,
            inherit.aes = FALSE
        ) +
        labs(title = title, x = "Residual in region A", y = "Residual in region B") +
        scale_x_continuous(limits = rng) +
        scale_y_continuous(limits = rng) +
        theme_classic(base_size = 10)
}

# 2. Scatter plots with annotated correlations, paged if too many pairs
# plot_residual_pair_scatter_paged <- function(res_mat,
#                                              cellType = NULL,
#                                              per_page = 12,
#                                              facet_font_size = 10,
#                                              ncol = NULL) {
#     stopifnot(is.matrix(res_mat))
#     regs  <- colnames(res_mat)
#     prs   <- if (length(regs) >= 2)
#         utils::combn(regs, 2, simplify = FALSE)
#     else
#         list()
#     if (!length(prs))
#         stop("Need >=2 regions (columns) in res_mat")
#
#     # build subplots
#     make_panel <- function(name, x, y) {
#         keep <- is.finite(x) & is.finite(y)
#         if (sum(keep) < 2)
#             return(NULL)
#         x <- x[keep]
#         y <- y[keep]
#         r <- cor(x, y)
#         rng <- range(c(x, y))
#         pad <- diff(rng) * 0.1
#         lim <- c(rng[1] - pad, rng[2] + pad)
#
#         ggplot(data.frame(x, y), aes(x, y)) +
#             geom_abline(
#                 intercept = 0,
#                 slope = 1,
#                 color = "black",
#                 linetype = "dashed"
#             ) +
#             geom_point(size = 2,
#                        alpha = 0.7,
#                        color = "steelblue") +
#             geom_smooth(
#                 method = "lm",
#                 formula = y ~ x,
#                 se = FALSE,
#                 linewidth = 0.6,
#                 color = "red"
#             ) +
#             annotate(
#                 "text",
#                 x = -Inf,
#                 y = Inf,
#                 label = sprintf("r = %.2f", r),
#                 hjust = -0.1,
#                 vjust = 1.2,
#                 size = 3.2
#             ) +
#             ggtitle(name) +
#             coord_cartesian(xlim = lim, ylim = lim) +
#             labs(x = "Residual in cell type A", y = "Residual in cell type B") +
#             theme_classic(base_size = 12) +
#             theme(plot.title = element_text(size = facet_font_size, hjust = 0.5))
#     }
#
#     plots <- list()
#     for (p in prs) {
#         name <- paste(p[1], "vs", p[2])
#         pan  <- make_panel(name, res_mat[, p[1]], res_mat[, p[2]])
#         if (!is.null(pan))
#             plots[[length(plots) + 1]] <- pan
#     }
#     if (!length(plots))
#         stop("No valid pairs after filtering")
#
#     # paging with cowplot
#     if (is.null(ncol))
#         ncol <- ceiling(sqrt(per_page))
#     nrow <- ceiling(per_page / ncol)
#     blanks <- function(n)
#         replicate(n, ggplot() + theme_void(), simplify = FALSE)
#
#     page_title <- "Age Prediction residuals (predicted - actual)"
#     if (!is.null(cellType))
#         page_title <- paste(cellType, page_title, sep = "\n")
#
#     pages <- list()
#     for (s in seq(1, length(plots), by = per_page)) {
#         page_plots <- plots[s:min(s + per_page - 1, length(plots))]
#         if (length(page_plots) < per_page)
#             page_plots <- c(page_plots, blanks(per_page - length(page_plots)))
#
#         grid <- cowplot::plot_grid(plotlist = page_plots,
#                                    ncol = ncol,
#                                    nrow = nrow)
#         pg <- cowplot::ggdraw() +
#             cowplot::draw_label(
#                 page_title,
#                 x = 0,
#                 y = 1,
#                 hjust = 0,
#                 vjust = 1,
#                 size = 14
#             ) +
#             cowplot::draw_plot(grid, y = 0, height = 0.94)
#         pages[[length(pages) + 1]] <- pg
#     }
#
#     pages
# }



# m: donors x regions residual matrix (from compute_residual_matrix)
# method: "pearson" or "spearman"
# cluster: TRUE = hierarchical clustering rows/cols; FALSE = keep input order
# annotate_cells: TRUE = show correlation values in heatmap cells; FALSE = no annotation
# plot_residual_corr_heatmap <- function(res_mat,
#                                        cellType = NULL,
#                                        annotate_cells = TRUE) {
#     stopifnot(is.matrix(res_mat))
#     C <- cor(res_mat, use = "pairwise.complete.obs")  # regions x regions
#
#     col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("#3b4cc0", "white", "#b40426"))
#
#     column_title <- "Age Prediction residuals (predicted - actual)"
#     if (!is.null(cellType))
#         column_title <- paste(cellType, column_title, sep = "\n")
#
#     # optional cell annotation
#     cf <- if (isTRUE(annotate_cells)) {
#         function(j, i, x, y, width, height, fill) {
#             grid::grid.text(sprintf("%.2f", C[i, j]), x, y, gp = grid::gpar(fontsize = 12))
#         }
#     } else
#         NULL
#
#     ComplexHeatmap::Heatmap(
#         C,
#         name = "r",
#         col = col_fun,
#         cluster_rows = TRUE,
#         cluster_columns = TRUE,
#         row_title = NULL,
#         column_title = column_title,
#         show_row_dend = TRUE,
#         show_column_dend = TRUE,
#         row_names_gp = grid::gpar(fontsize = 9),
#         column_names_gp = grid::gpar(fontsize = 9),
#         heatmap_legend_param = list(at = c(-1, -0.5, 0, 0.5, 1), title = "Correlation"),
#         cell_fun = cf
#     )
# }




# get_age_prediction_metrics <- function(results) {
#     groups <- split(results, list(results$cell_type, results$region), drop = TRUE)
#
#     rows <- lapply(groups, function(df) {
#         m <- compute_age_metrics(df$pred, df$age)  # named numeric
#         data.frame(
#             cell_type = df$cell_type[1],
#             region    = df$region[1],
#             as.list(m),
#             row.names = NULL,
#             check.names = FALSE
#         )
#     })
#
#     r = do.call(rbind, rows)
#     rownames (r) <- NULL
#     r = r[order(r$cell_type, r$region), ]
#
#     return (r)
# }

# load_model_metrics <- function (model_file_dir) {
#     model_files = list.files(model_file_dir,
#                              pattern = "donor_age_model_metrics*",
#                              full.names = TRUE)
#
#     #x=model_files[1]
#     parseOne <- function (x) {
#         a = read.table(
#             x,
#             header = TRUE,
#             sep = "\t",
#             stringsAsFactors = FALSE
#         )
#
#         return (a)
#     }
#
#     all_models = lapply(model_files, parseOne)
#     logger::log_info(paste("Loaded model metrics for "),
#                      length(all_models),
#                      " cell types")
#     all_models = do.call(rbind, all_models)
#     return(all_models)
#
#     # Some ad-hoc plots.
#     #plot (all_models[all_models$set=="Cross-validation",]$median_abs_error, all_models[all_models$set=="Held-out donors",]$median_abs_error, xlab="Cross Validation", ylab="Held-out donors", main="median absolute error", )
#     #abline(0,1, col="red")
#
#     #plot (all_models[all_models$set=="Cross-validation",]$r, all_models[all_models$set=="Held-out donors",]$r, xlab="Cross Validation", ylab="Held-out donors", main="correlation")
#     #abline(0,1, col="red")
#
#     #z=data.frame(cell_type=all_models[all_models$set=="Cross-validation",]$cell_type, cv=all_models[all_models$set=="Cross-validation",]$r, ho=all_models[all_models$set=="Held-out donors",]$r)
#     #z$abs=abs(z$cv-z$ho)
#     #z=z[order(z$abs, decreasing = T),]
#
# }

# load_models <- function (model_file_dir) {
#     model_files = list.files(model_file_dir,
#                              pattern = "donor_age_model_coefficients_*",
#                              full.names = TRUE)
#
#     #x=model_files[1]
#     parseOne <- function (x) {
#         a = read.table(
#             x,
#             header = TRUE,
#             sep = "\t",
#             stringsAsFactors = FALSE
#         )
#         return (a)
#     }
#
#     all_models = lapply(model_files, parseOne)
#     logger::log_info(paste("Loaded model coefficients for "),
#                      length(all_models),
#                      " cell types")
#     all_models = do.call(rbind, all_models)
#     return(all_models)
# }
#
# load_model_predictions <- function (model_file_dir) {
#     model_files = list.files(model_file_dir,
#                              pattern = "donor_age_predictions*",
#                              full.names = TRUE)
#
#     #x=model_files[1]
#     parseOne <- function (x) {
#         a = read.table(
#             x,
#             header = TRUE,
#             sep = "\t",
#             stringsAsFactors = FALSE
#         )
#         return (a)
#     }
#
#     all_models = lapply(model_files, parseOne)
#     logger::log_info(paste("Loaded model predictions for "),
#                      length(all_models),
#                      " cell types")
#     all_models = do.call(rbind, all_models)
#     return(all_models)
# }
