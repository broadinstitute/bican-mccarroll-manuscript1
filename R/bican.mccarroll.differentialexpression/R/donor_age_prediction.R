## ------------------------------------------------------------------
## DEBUG ONLY — interactive development block
## Uncomment for manual experimentation in RStudio
## ------------------------------------------------------------------

# library(edgeR)
# library(glmnet)
# library(ggplot2)
# library(logger)
# library (cowplot)
# library (dplyr)
# source("R/age_model_fit_core.R", local=TRUE)


#' Train donor age prediction models across all cell types and write results
#'
#' Loads a donor-level RNA-seq dataset (as a \code{DGEList}), restricts to autosomal genes,
#' and then fits a donor age prediction model per cell type . For each successfully fit
#' cell type, this function generates QC plots and writes
#' concatenated tab-delimited outputs (donor predictions, model coefficients, and metrics)
#' across all cell types.
#'
#' Cell types may be skipped due to insufficient genes after filtering.
#'
#' @param data_dir Directory containing the serialized input object used by
#'   \code{bican.mccarroll.differentialexpression::prepare_data_for_differential_expression()}.
#' @param data_name Name of the dataset object to load (passed through to
#'   \code{prepare_data_for_differential_expression()}).
#' @param result_dir Output directory for covariates, metrics, and predictions/residuals
#' @param age_de_results_dir Directory containing per-cell-type age differential expression
#'   results used by \code{\link{get_age_de_results}} and \code{\link{predict_age_celltype}}.
#' @param outPDFFile Optional path to a PDF file. If provided, QC plots generated during
#'   processing are printed to this PDF; otherwise plots are printed to the active device.
#' @param contig_yaml_file Path to the contig-groups YAML file used by
#'   \code{\link{filter_dge_to_autosomes}}.
#' @param reduced_gtf_file Path to the reduced GTF used by \code{\link{filter_dge_to_autosomes}}.
#' @param donor_col Column name in \code{dge$samples} identifying the donor. Default \code{"donor"}.
#' @param age_col Column name in \code{dge$samples} giving chronological age. Default \code{"age"}.
#' @param seed Integer random seed used for splitting/training inside \code{\link{predict_age_celltype}}.
#' @param fdr_threshold FDR cutoff used when filtering to age DE genes inside
#'   \code{\link{predict_age_celltype}}. Default \code{0.05}.
#' @param optimize_alpha Logical; if \code{TRUE}, select elastic net \code{alpha} by grid search
#'   (via \code{\link{train_enet_cv_optimize_alpha}}). If \code{FALSE}, use \code{alpha_fixed}.
#' @param alpha_fixed Numeric elastic net \code{alpha} used when \code{optimize_alpha = FALSE}
#'   (\code{0} ridge, \code{1} lasso). Default \code{0.5}.
#' @param mc_repeats Integer; number of repeated stratified 80/20 splits used for the final
#'   Monte Carlo residual summaries inside \code{\link{predict_age_celltype}}.
#' @param min_donors Minimum number of donors required to fit a model for a given cell type.
#' @param donor_age_range Optional numeric vector of length 2 giving the minimum and maximum age in decades.
#' @param n_cores Integer; number of CPU cores to use for Monte Carlo Cross Validation
#' @importFrom logger log_info
#' @importFrom grDevices pdf dev.off
#'
#' @export
predict_age_by_celltype<-function (data_dir, data_name, result_dir, age_de_results_dir, outPDFFile=NULL,
                                  contig_yaml_file, reduced_gtf_file,
                                  donor_col = "donor", age_col = "age",
                                  seed =42, fdr_threshold=0.05,
                                  optimize_alpha=FALSE, alpha_fixed=0.5,
                                  mc_repeats=200, min_donors=50,
                                  donor_age_range=c(), n_cores=1) {

    #validate the output directory exists
    if (!dir.exists(result_dir)) {
        dir.create(result_dir, recursive = TRUE)
    }

    #load the DGEList and prepare the data
    d=bican.mccarroll.differentialexpression::prepare_data_for_differential_expression(data_dir, data_name, randVars=c(), fixedVars=c())
    dge=d$dge

    #filter the genes to autosomal only
    dge=filter_dge_to_autosomes (dge, contig_yaml_file, reduced_gtf_file)

    #filter to age if requested.
    dge=filter_dge_to_donor_age(dge, donor_age_range)

    cell_type_list=unique(dge$samples$cell_type)
    #cell_type_list="microglia" #hard coded for now.
    lineStr <- strrep("=", 80)

    if (!is.null(outPDFFile)) {
        pdf(outPDFFile)
    }

    #get the output base name for text
    output_basename <- get_output_basename(outPDFFile)

    #infer the columns you'll need to retain
    retained_features=c(donor_col, age_col)

    test_set_metrics_list=list()
    outputs_list <- list()

    for (cellType in cell_type_list) {
        logger::log_info(lineStr)
        logger::log_info(paste("Learning donor age model from expression for cell type:", cellType))
        logger::log_info(lineStr)
        r=predict_age_celltype(cellType, dge, retained_features=retained_features,
                               donor_col = donor_col, age_de_results_dir=age_de_results_dir,
                               region=NULL, fdr_threshold=fdr_threshold,
                               optimize_alpha=optimize_alpha, alpha_fixed=alpha_fixed,
                               mc_repeats=mc_repeats, min_donors=min_donors,
                               seed=seed, n_cores=n_cores)

        # gracefully skip cell types that short-circuit (e.g. too few genes)
        if (is.null(r)) {
            next
        }

        test_set_metrics_list[[cellType]]=r$test_set_metrics

        # Keep only what will be written
        outputs_list[[cellType]] <- extract_age_outputs(r, cellType = cellType, region = NA_character_)

        # QC plots
        plot_list=age_prediction_qc_plots(r, cellType)
        plot_list$mc_donor_pred_plot<-plot_mc_donor_predictions(outputs_list[[cellType]]$donor_predictions)

        for (p in plot_list) {
            if (!is.null(p))
                print(p)
        }

    }

    # Write one file per artifact type
    write_age_outputs_all(outputs_list, result_dir, output_basename)

    test_set_metrics_df <- do.call(rbind, test_set_metrics_list)
    mean_metrics <- apply (test_set_metrics_df, 2, mean)
    logger::log_info(
        "Mean metrics of held out test sets across all cell types: {paste(names(mean_metrics), sprintf('[%.3f]', mean_metrics), collapse = ', ')}"
    )

    if (!is.null(outPDFFile)) {
        dev.off()
    }

}

#' Train donor age prediction models across cell types and regions and write results
#'
#' Loads a donor-level RNA-seq dataset (as a \code{DGEList}), restricts to autosomal genes,
#' and then fits a donor age prediction model for each (cell type, region) combination.
#' Regions are discovered per cell type from \code{dge$samples$region}. For each
#' successfully fit combination, this function generates QC plots and writes concatenated
#' tab-delimited outputs (donor predictions, model coefficients, and metrics) across all
#' combinations.
#'
#' Combinations may be skipped due to insufficient genes after filtering.
#'
#' @param data_dir Directory containing the serialized input object used by
#'   \code{bican.mccarroll.differentialexpression::prepare_data_for_differential_expression()}.
#' @param data_name Name of the dataset object to load (passed through to
#'   \code{prepare_data_for_differential_expression()}).
#' @param result_dir Output directory for covariates, metrics, and predictions/residuals.
#' @param age_de_results_dir Directory containing per-cell-type or per-(cell type, region)
#'   age differential expression results used by \code{\link{get_age_de_results}} and
#'   \code{\link{predict_age_celltype}}.
#' @param outPDFFile Optional path to a PDF file. If provided, QC plots generated during
#'   processing are printed to this PDF; otherwise plots are printed to the active device.
#' @param contig_yaml_file Path to the contig-groups YAML file used by
#'   \code{\link{filter_dge_to_autosomes}}.
#' @param reduced_gtf_file Path to the reduced GTF used by \code{\link{filter_dge_to_autosomes}}.
#' @param donor_col Column name in \code{dge$samples} identifying the donor. Default \code{"donor"}.
#' @param age_col Column name in \code{dge$samples} giving chronological age. Default \code{"age"}.
#' @param seed Integer random seed used for splitting/training inside \code{\link{predict_age_celltype}}.
#' @param fdr_threshold FDR cutoff used when filtering to age DE genes inside
#'   \code{\link{predict_age_celltype}}. Default \code{0.05}.
#' @param optimize_alpha Logical; if \code{TRUE}, select elastic net \code{alpha} by grid search
#'   (via \code{\link{train_enet_cv_optimize_alpha}}). If \code{FALSE}, use \code{alpha_fixed}.
#' @param alpha_fixed Numeric elastic net \code{alpha} used when \code{optimize_alpha = FALSE}
#'   (\code{0} ridge, \code{1} lasso). Default \code{0.5}.
#' @param mc_repeats Integer; number of repeated stratified 80/20 splits used for the final
#'   Monte Carlo residual summaries inside \code{\link{predict_age_celltype}}.
#' @param min_donors Minimum number of donors required to fit a model for a given cell type.
#' @param donor_age_range Optional numeric vector of length 2 giving the minimum and maximum age in decades.
#' @param n_cores Integer; number of CPU cores to use for Monte Carlo Cross Validation
#' @importFrom logger log_info
#' @importFrom grDevices pdf dev.off
#'
#' @export
predict_age_by_celltype_region <- function(data_dir, data_name, result_dir, age_de_results_dir,
                                           outPDFFile = NULL,
                                           contig_yaml_file, reduced_gtf_file,
                                           donor_col = "donor", age_col = "age",
                                           seed = 42, fdr_threshold = 0.05,
                                           optimize_alpha = FALSE, alpha_fixed = 0.5,
                                           mc_repeats = 200, min_donors=50,
                                           donor_age_range=c(), n_cores=1) {

    if (!dir.exists(result_dir)) {
        dir.create(result_dir, recursive = TRUE)
    }

    #load the DGEList and prepare the data
    d=bican.mccarroll.differentialexpression::prepare_data_for_differential_expression(data_dir, data_name, randVars=c(), fixedVars=c())
    dge=d$dge

    #filter the genes to autosomal only
    dge=filter_dge_to_autosomes (dge, contig_yaml_file, reduced_gtf_file)

    #filter to age if requested.
    dge=filter_dge_to_donor_age(dge, donor_age_range)

    cell_type_list <- unique(dge$samples$cell_type)
    lineStr <- strrep("=", 80)

    if (!is.null(outPDFFile)) {
        grDevices::pdf(outPDFFile)
    }

    output_basename <- get_output_basename(outPDFFile)
    retained_features <- c(donor_col, age_col)

    test_set_metrics_list <- list()
    outputs_list <- list()

    for (cellType in cell_type_list) {

        logger::log_info(lineStr)
        logger::log_info(paste("Learning donor age model from expression for cell type:", cellType))
        logger::log_info(lineStr)

        regions <- unique(dge$samples$region[dge$samples$cell_type == cellType])
        regions <- regions[!is.na(regions)]

        # If region is missing entirely for this cell type, skip (or treat as one "NA" bin).
        if (length(regions) == 0) {
            next
        }

        for (region in regions) {

            logger::log_info(paste0("Region: [", region, "]"))

            r <- predict_age_celltype(
                cellType,
                dge,
                retained_features = retained_features,
                donor_col = donor_col,
                age_de_results_dir = age_de_results_dir,
                region = region,
                fdr_threshold = fdr_threshold,
                optimize_alpha = optimize_alpha,
                alpha_fixed = alpha_fixed,
                mc_repeats = mc_repeats,
                min_donors=min_donors,
                seed = seed,
                n_cores=n_cores
            )

            if (is.null(r)) {
                next
            }

            key <- paste(cellType, region, sep = "__")
            test_set_metrics_list[[key]] <- r$test_set_metrics

            outputs_list[[key]] <- extract_age_outputs(
                r,
                cellType = cellType,
                region = as.character(region)
            )

            plot_list <- age_prediction_qc_plots(r, cellType, region=region)

            plot_list$mc_donor_pred_plot <- plot_mc_donor_predictions_report(
                outputs_list[[key]]$donor_predictions,
                gam_fit_df = outputs_list[[key]]$gam_fit_df)$combined_plot

            for (p in plot_list) {
                if (!is.null(p) & !is.null(outPDFFile)) {
                    print(p)
                }
            }
        }
    }

    write_age_outputs_all(outputs_list, result_dir, output_basename)

    if (length(test_set_metrics_list) > 0) {
        test_set_metrics_df <- do.call(rbind, test_set_metrics_list)
        mean_metrics <- apply(test_set_metrics_df, 2, mean)
        logger::log_info(
            "Mean metrics of held out test sets across all fits: {paste(names(mean_metrics), sprintf('[%.3f]', mean_metrics), collapse = ', ')}"
        )
    } else {
        logger::log_info("No successful (cell type, region) fits; no mean metrics to report.")
    }

    if (!is.null(outPDFFile)) {
        grDevices::dev.off()
    }

    invisible(NULL)
}




age_prediction_qc_plots<-function (r, cellType, region=NULL) {
    age_dist_plot <- plot_age_hist_stacked(r$dge_train, r$dge_test, r$cellType)

    #Region aware title
    title_celltype_str<-cellType
    if (!is.null(region)) {
        title_celltype_str<-paste0(cellType, " [", region, "]")
    }

    pred_age_plot_train <- plot_age_predictions(r$cv_model$donor_age_predictions, cellType, titleStr=paste("TRAIN", title_celltype_str, "\n"))
    pred_age_plot_test <- plot_age_predictions(r$test_set_predictions, cellType, titleStr=paste("TEST", title_celltype_str, "\n [20% held out data]"))

    #Region aware title
    perf_plot_title_str=NULL
    if (!is.null(region))
        perf_plot_title_str = paste("Performance", paste0(cellType, " [", region, "]"))


    perf_plot<-plot_model_performance(
        cv_metrics = r$cv_model$overall_metrics,
        per_fold_metrics = r$cv_model$per_fold_metrics,
        test_metrics = r$test_set_metrics,
        cellType = cellType,
        str_title=perf_plot_title_str
    )

    p=cowplot::plot_grid(pred_age_plot_train, pred_age_plot_test, age_dist_plot, perf_plot, ncol=2)

    coef1=plot_fold_coefficients(r$cv_model$fold_models, positive_only = FALSE, min_nonzero = 1, top_n=NA, cellType=cellType)

    #plot the age DE results if available against the model average coefficients
    model_vs_de_plot=NULL
    if (!is.null(r$age_de_results))
        model_vs_de_plot=plot_model_vs_de(r$cv_model$final_model, r$age_de_results, cellType = NULL)

    return (list(main_qc_plot=p, coefficient_plot=coef1, model_vs_de_plot=model_vs_de_plot))

}

#' Fit a donor-level age prediction model for a single cell type (optionally within a region)
#'
#' Trains an elastic net (or ridge/lasso depending on \code{alpha_fixed}) age
#' prediction model from gene expression for one cell type. If \code{region} is
#' provided, the input \code{DGEList} is first subset to the requested cell type
#' and region before collapsing multiple observations per donor into a single
#' donor-level profile.
#'
#' The function:
#' \enumerate{
#'   \item Subsets the input \code{DGEList} to the requested cell type (and region, if set).
#'   \item Collapses multiple observations per donor into a single donor-level profile via
#'         \code{collapse_by_donor()}.
#'   \item Applies standard expression QC and gene filtering steps.
#'   \item Optionally filters genes to a set significant in external age DE results.
#'   \item Performs an outer 80/20 donor-level holdout split for validation.
#'   \item Trains the model on the 80\% training donors using inner stratified K-fold CV to
#'         choose the regularization strength (and optionally alpha).
#'   \item Evaluates the model on the held-out 20\% donors.
#'   \item Computes "final" donor residuals by running repeated stratified 80/20 Monte Carlo
#'         splits on the full donor-level dataset and summarizing the per-donor out-of-fold
#'         residual distribution.
#'   \item Fits a final model on all donors for external prediction.
#' }
#'
#' @param cellType Character scalar; the cell type label used to subset
#'   \code{dge$samples$cell_type}.
#' @param dge A \code{DGEList} containing counts and sample metadata. Must contain
#'   \code{dge$samples$cell_type}, \code{dge$samples$region} (if using \code{region}),
#'   donor IDs, and ages.
#' @param retained_features Character vector of sample-level columns to retain
#'   during donor collapsing (passed to \code{collapse_by_donor()}).
#' @param donor_col Column name in \code{dge$samples} containing donor IDs.
#' @param age_de_results_dir Directory containing age differential expression
#'   results used by \code{get_age_de_results()}.
#' @param region Optional character scalar; if not \code{NULL}, subset
#'   \code{dge$samples$region == region} before collapsing across samples. When set,
#'   \code{get_age_de_results()} is called with the same \code{region} so that region-
#'   specific DE results can be used (if present).
#' @param fdr_threshold Numeric; FDR cutoff for selecting genes from age DE
#'   results. If \code{get_age_de_results()} returns \code{NULL}, no DE-based
#'   feature filtering is performed.
#' @param optimize_alpha Logical; if \code{TRUE}, optimize elastic net alpha via
#'   \code{train_enet_cv_optimize_alpha()}. If \code{FALSE}, use \code{alpha_fixed}.
#' @param alpha_fixed Numeric in [0, 1]; elastic net alpha used when
#'   \code{optimize_alpha = FALSE}. \code{0} corresponds to ridge and \code{1}
#'   corresponds to lasso.
#' @param mc_repeats Integer; number of repeated stratified 80/20 splits to run in
#'   \code{compute_final_residuals_repeated_splits()} when computing final donor
#'   residual summaries on the full dataset.
#' @param min_donors Integer; minimum number of donors required to fit the model.
#'   If the donor-level dataset has fewer donors, the function returns \code{NULL}.
#' @param seed Integer random seed used for the outer holdout split, inner CV fold
#'   construction, Monte Carlo repeated splits, and the final refit for external prediction.
#' @param n_cores Integer; number of parallel cores to use for Monte Carlo cross fold validation.
#' @return A list with components:
#' \itemize{
#'   \item \code{cell_type}: cell type name.
#'   \item \code{region}: region used for training (\code{NA_character_} if not provided).
#'   \item \code{dge_train}: donor-level training \code{DGEList} (80\% donors).
#'   \item \code{dge_test}: donor-level held-out \code{DGEList} (20\% donors).
#'   \item \code{cv_model}: result of \code{train_enet_cv()} or \code{train_enet_cv_optimize_alpha()}
#'         on \code{dge_train}.
#'   \item \code{test_set_predictions}: data.frame of predictions on held-out donors with columns
#'         \code{donor}, \code{age}, and \code{pred}.
#'   \item \code{test_set_metrics}: performance metrics on held-out donors as returned by
#'         \code{compute_age_metrics()}.
#'   \item \code{age_de_results}: age DE results returned by \code{get_age_de_results()}
#'         (or \code{NULL} if unavailable).
#'   \item \code{final_oof_predictions}: per-donor residual summary returned by
#'         \code{compute_final_residuals_repeated_splits()} (one row per donor).
#'   \item \code{final_oof_metrics}: overall metrics computed from the per-donor average prediction
#'         across repeats (\code{final_repeat$avg_pred_metrics}).
#'   \item \code{final_cv_model}: final model fit on all donors for external prediction
#'         (computed with \code{compute_oof = FALSE}).
#' }
#'
#' @details
#' The "final" donor residuals are computed using repeated stratified 80/20 splits on the full
#' donor-level dataset. For each repeat, donors in the test split are predicted by a model
#' trained on the corresponding training split, and only those out-of-fold predictions contribute
#' to the donor's residual distribution.
#'
#' @export
predict_age_celltype <- function(cellType, dge, retained_features = c("donor", "age"),
                                 donor_col = "donor", age_de_results_dir,
                                 region = NULL,
                                 fdr_threshold = 0.05,
                                 optimize_alpha = FALSE, alpha_fixed = 0.5,
                                 mc_repeats = 100, min_donors=50,
                                 seed = 1, n_cores=1) {

    if (optimize_alpha) {
        logger::log_info("Optimizing alpha via grid search.")
    } else {
        logger::log_info(paste("Using fixed alpha =", alpha_fixed))
    }

    # Subset to cell type (and region if provided) BEFORE collapsing
    if (is.null(region)) {
        dge_cell <- dge[, dge$samples$cell_type == cellType, keep.lib.sizes = TRUE]
    } else {
        dge_cell <- dge[, dge$samples$cell_type == cellType & dge$samples$region == region,
                        keep.lib.sizes = TRUE]
    }

    # Collapse multiple samples per donor
    dge_cell <- collapse_by_donor(dge_cell, donor_col = donor_col, keep_cols = retained_features, sum_cols="num_nuclei")

    # Filtering samples by library size
    r <- bican.mccarroll.differentialexpression::filter_by_libsize(
        dge_cell, threshold_sd = 1.96, bins = 50, strTitlePrefix = cellType
    )
    dge_cell <- r$dge

    # Filter to the top 75% of highly expressed genes as a first pass
    dge_cell <- bican.mccarroll.differentialexpression::filter_top_expressed_genes(
        dge_cell, gene_filter_frac = 0.75, verbose = TRUE
    )

    # Filter to cpm cutoff of 1
    r2 <- bican.mccarroll.differentialexpression::plot_logCPM_density_quantiles(
        dge_cell, cpm_cutoff = 1, logCPM_xlim = c(-5, 15),
        lower_quantile = 0.05, upper_quantile = 0.95, quantile_steps = 5,
        min_samples = 1, fraction_samples = 0.1
    )
    dge_cell <- r2$filtered_dge

    # Get age DE results (optional), region-aware
    age_de_results <- get_age_de_results(cellType, age_de_results_dir,
                                         region = region,
                                         fdr_threshold = fdr_threshold)
    if (!is.null(age_de_results)) {
        genes_to_keep <- intersect(rownames(dge_cell), rownames(age_de_results))
        logger::log_info(paste("Filtering to", length(genes_to_keep),
                               "genes significant in age DE results at FDR <=", fdr_threshold))
        dge_cell <- dge_cell[genes_to_keep, ]
    }

    # check that there are at least 2 genes to fit the model
    if (!has_min_features(dge_cell, min_features = 2L, context = cellType)) {
        return(NULL)
    }

    #check that there are a minimum number of donors to fit the model
    if (!has_min_donors(dge_cell, min_donors = min_donors, context = cellType)) {
        return(NULL)
    }

    # Outer holdout split (80/20)
    n_bins <- 5
    d <- make_holdout_stratified(dge_cell, age_col = "age", donor_col = "donor",
                                 test_prop = 0.20, n_bins = n_bins, seed = seed)
    dge_train <- d$dge_train
    dge_test  <- d$dge_test

    # Validate train/test split are completely disjoint on donors
    if (length(intersect(dge_train$samples$donor, dge_test$samples$donor)) > 0) {
        stop("Train and test donors are not disjoint!")
    }

    logger::log_info(paste("Training data has", ncol(dge_train), "donors across",
                           nrow(dge_train$counts), "genes"))
    logger::log_info(paste("Final holdout test data has", ncol(dge_test), "donors across",
                           nrow(dge_test$counts), "genes"))

    # Inner CV folds on the 80% training data
    k_fold_index <- make_cv_folds_stratified(dge_train, age_col = "age", donor_col = "donor",
                                             K = 5, n_bins = n_bins, seed = seed)

    if (optimize_alpha) {
        cv_model <- train_enet_cv_optimize_alpha(dge_train, k_fold_index,
                                                 age_col = "age", donor_col = "donor",
                                                 seed = seed)
    } else {
        cv_model <- train_enet_cv(dge_train, k_fold_index,
                                  age_col = "age", donor_col = "donor",
                                  alpha_fixed = alpha_fixed, seed = seed)
    }

    num_nonzero_coefs <- sum(cv_model$final_model$coef != 0)
    logger::log_info(paste("Final (80% train) model has", num_nonzero_coefs, "non-zero coefficients."))

    # Predict on 20% held-out data using the 80% train final model
    test_set_predictions <- predict_age_from_dge(dge_test, cv_model$final_model, prior.count = 1)
    metrics_test <- compute_age_metrics(test_set_predictions$pred, test_set_predictions$age)

    # Final repeated-split OOF residual computation on full dataset
    str=paste0("Computing final repeated-split OOF residuals on all donors [",ncol(dge_cell),"].")
    logger::log_info(str)

    final_repeat <- compute_final_residuals_repeated_splits_parallel(
        dge_cell,
        age_col = "age",
        donor_col = "donor",
        n_bins = n_bins,
        test_prop = 0.20,
        n_repeats = mc_repeats,
        optimize_alpha = optimize_alpha,
        alpha_fixed = alpha_fixed,
        seed = seed,
        verbose_every = mc_repeats / 10,
        n_cores=n_cores,
    )



    # Final model for external prediction (fit on ALL donors; no OOF)
    logger::log_info("Fitting final model on all donors for external prediction.")

    if (optimize_alpha) {
        final_cv_model <- train_enet_cv_optimize_alpha(
            dge_cell,
            k_fold_index = NULL,
            age_col = "age",
            donor_col = "donor",
            alpha_grid = seq(0, 1, by = 0.1),
            compute_oof = FALSE,
            seed = seed
        )
    } else {
        final_cv_model <- train_enet_cv(
            dge_cell,
            k_fold_index = NULL,
            age_col = "age",
            donor_col = "donor",
            alpha_fixed = alpha_fixed,
            compute_oof = FALSE,
            seed = seed
        )
    }

    final_oof_predictions <- final_repeat$donor_residual_summary
    final_oof_metrics <- final_repeat$avg_pred_metrics
    gam_fit_df <- final_repeat$gam_fit_df

    #add the cell type and region to the gam_fit_df for later use, and add them as the first columns for easier parsing
    region_value <- if (is.null(region)) NA_character_ else as.character(region)
    gam_fit_df <- data.frame(cell_type = cellType, region = region_value, gam_fit_df)

    list(
        cell_type = cellType,
        region = if (is.null(region)) NA_character_ else as.character(region),
        dge_train = dge_train,
        dge_test = dge_test,
        cv_model = cv_model,
        test_set_predictions = test_set_predictions,
        test_set_metrics = metrics_test,
        age_de_results = age_de_results,
        final_cv_model = final_cv_model,
        final_oof_predictions = final_oof_predictions,
        final_oof_metrics = final_oof_metrics,
        gam_fit_df = gam_fit_df
    )
}


#' Check that a DGEList has enough features for model fitting
#'
#' Internal guard to ensure that downstream modeling steps (e.g. glmnet)
#' have a sufficient number of features. Intended to be used to short-circuit
#' model fitting for cell types with too few genes after filtering.
#'
#' @param dge A \code{DGEList}.
#' @param min_features Integer; minimum number of features (rows) required.
#' @param context Optional character string used in log messages
#'   (e.g. cell type name).
#'
#' @return Logical; \code{TRUE} if \code{nrow(dge) >= min_features}, otherwise
#'   \code{FALSE}.
#'
#' @keywords internal
has_min_features <- function(dge, min_features = 2L, context = NULL) {
    stopifnot("DGEList" %in% class(dge))

    if (nrow(dge) < min_features) {
        msg <- paste0(
            "Skipping model fit: only ", nrow(dge),
            " feature(s) available (need >= ", min_features, ")."
        )
        if (!is.null(context)) {
            msg <- paste0(context, ": ", msg)
        }
        logger::log_warn(msg)
        return(FALSE)
    }

    TRUE
}

#' Check that a DGEList has enough features for model fitting
#'
#' Internal guard to ensure that downstream modeling steps (e.g. glmnet)
#' have a sufficient number of features. Intended to be used to short-circuit
#' model fitting for cell types with too few genes after filtering.
#'
#' @param dge A \code{DGEList}.
#' @param min_donors Integer; minimum number of features (rows) required.
#' @param context Optional character string used in log messages
#'   (e.g. cell type name).
#'
#' @return Logical; \code{TRUE} if \code{nrow(dge) >= min_features}, otherwise
#'   \code{FALSE}.
#'
#' @keywords internal
has_min_donors <- function(dge, min_donors = 50L, context = NULL) {
    stopifnot("DGEList" %in% class(dge))

    if (ncol(dge) < min_donors) {
        msg <- paste0(
            "Skipping model fit: only ", ncol(dge),
            " donors available (need >= ", min_donors, ")."
        )
        if (!is.null(context)) {
            msg <- paste0(context, ": ", msg)
        }
        logger::log_warn(msg)
        return(FALSE)
    }
    TRUE
}



########################
# PLOTS
########################

plot_age_predictions <- function(df, cellType, titleStr = NULL) {
    if (is.null(titleStr))
        titleStr <- paste("Test set predictions for", cellType)

    # compute symmetric range
    range_all <- range(c(df$age, df$pred), na.rm = TRUE)
    range_all[1]=range_all[1]*0.9
    range_all[2]=range_all[2]*1.1

    xlim_val <- range_all
    ylim_val <- range_all

    #Make R CMD CHECK happy
    age <- pred <- NULL

    ggplot(df, aes(x = age, y = pred)) +
        geom_point(size = 2, alpha = 0.7, color = "steelblue") +
        geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
        geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.6, color = "red") +
        labs(
            title = titleStr,
            x = "Chronological Age [decades]",
            y = "Predicted Age [decades]"
        ) +
        scale_x_continuous(limits = xlim_val) +
        scale_y_continuous(limits = ylim_val) +
        theme_classic(base_size = 12)
}




plot_age_hist_stacked <- function(dge_train, dge_test, cellType = NULL) {
    stopifnot("DGEList" %in% class(dge_train), "DGEList" %in% class(dge_test))

    df <- rbind(
        data.frame(age = dge_train$samples$age, set = "Train"),
        data.frame(age = dge_test$samples$age, set = "Held-out donors")
    )
    df$set <- factor(df$set, levels = c("Train", "Held-out donors"))
    x_range <- range(df$age, na.rm = TRUE)

    #Make R CMD CHECK happy
    age <- set <- NULL

    ggplot(df, aes(x = age, fill = set)) +
        geom_histogram(color = "white", bins = 30, na.rm = TRUE) +
        facet_grid(rows = vars(set)) +
        scale_x_continuous(expand = expansion(mult = 0)) +
        coord_cartesian(xlim = x_range) +
        scale_fill_manual(values = c("Train" = "steelblue", "Held-out donors" = "orange")) +
        labs(
            title = paste0("Age distribution",
                           if (!is.null(cellType)) paste0(" (", cellType, ")")),
            x = "Age (years)", y = "Count"
        ) +
        theme_classic(base_size = 12) +
        theme(
            strip.text.y = element_blank(),
            legend.position = "top",
            legend.title = element_blank()
        )
}


plot_model_performance <- function(cv_metrics, per_fold_metrics, test_metrics, cellType = NULL, str_title = NULL) {
    stopifnot(is.data.frame(cv_metrics), is.data.frame(per_fold_metrics), is.data.frame(test_metrics))

    #format the title string
    if (is.null(str_title)) {
        str_title = paste0(
            "Performance",
            if (!is.null(cellType)) paste0(" (", cellType, ")")
        )
    }

    # Long per-fold metrics
    df_folds <- tidyr::pivot_longer(per_fold_metrics,
                             cols = c("r","mean_abs_error","median_abs_error"),
                             names_to = "metric", values_to = "value")
    df_folds$metric[df_folds$metric=="mean_abs_error"]   <- "Mean AE"
    df_folds$metric[df_folds$metric=="median_abs_error"] <- "Median AE"
    df_folds$source <- "Per-fold"
    df_folds$x <- "All"  # single column per facet

    # Points for CV mean and held-out
    df_points <- rbind(
        data.frame(source="Cross-validation mean", metric="r",         value=cv_metrics$r),
        data.frame(source="Cross-validation mean", metric="Mean AE",       value=cv_metrics$mean_abs_error),
        data.frame(source="Cross-validation mean", metric="Median AE", value=cv_metrics$median_abs_error),
        data.frame(source="Held-out donors",         metric="r",         value=test_metrics$r),
        data.frame(source="Held-out donors",         metric="Mean AE",       value=test_metrics$mean_abs_error),
        data.frame(source="Held-out donors",         metric="Median AE", value=test_metrics$median_abs_error)
    )
    df_points$x <- "All"

    #Make R CMD CHECK happy
    x <- value <- metric <- source <- NULL

    ggplot() +
        # box + jitter for per-fold metrics
        geom_boxplot(
            data = df_folds,
            aes(x = x, y = value),
            width = 0.5,
            alpha = 0.25,
            fill = "lightblue",
            color = "lightblue",
            outlier.shape = NA
        ) +
        geom_jitter(
            data = df_folds,
            aes(x = x, y = value),
            width = 0.06,
            alpha = 0.6,
            size = 1.8,
            color = "lightblue"
        ) +
        # points for CV mean (dark blue filled circle) and held-out (orange)
        geom_point(
            data = df_points,
            aes(x = x, y = value, color = source, shape = source, fill = source),
            size = 3, stroke = 0.5
        ) +
        facet_wrap(~ metric, scales = "free_y", nrow = 1) +
        scale_color_manual(values = c(
            "Cross-validation mean" = "darkblue",
            "Held-out donors" = "orange"
        )) +
        scale_fill_manual(values = c(
            "Cross-validation mean" = "darkblue",
            "Held-out donors" = "orange"
        )) +
        scale_shape_manual(values = c(
            "Cross-validation mean" = 19,  # filled circle
            "Held-out donors" = 21           # circle with border
        )) +
        labs(
            title = str_title,
            x = NULL, y = NULL
        ) +
        theme_classic(base_size = 12) +
        theme(
            legend.position = "top",
            legend.title = element_blank(),
            legend.text = element_text(size = 8),
            legend.key.size = unit(0.4, "lines"),
            legend.margin = margin(t = -3, b = -3),
            strip.background = element_blank(),
            strip.text = element_text(size = 9),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
        )

}

# Plot per-feature coefficient distributions across K folds
# - fold_models: list of data.frames from CV (each with columns: feature, coef, fold, mean_train, sd_train, intercept)
# - positive_only: keep features with coef > 0 in at least one fold (default TRUE); if FALSE, keep any nonzero
# - min_nonzero: keep features with at least this many nonzero coefficients across folds
# - top_n: limit to top-N features by median coefficient (after filtering); use NA to keep all
# - cellType: optional title suffix
# Plot per-feature coefficient values across CV folds (nonzero coefficients only)
plot_fold_coefficients <- function(fold_models, positive_only = FALSE,
                                   min_nonzero = 1, top_n = NA, cellType = NULL) {
    stopifnot(is.list(fold_models), length(fold_models) >= 1)

    # Make R CMD check happy (non-standard eval + ggplot2)
    coef <- feature <- fold <- med <- any_pos <- nnz <- NULL

    df <- dplyr::bind_rows(fold_models)

    if (!("fold" %in% names(df))) stop("fold_models must include a 'fold' column.")
    if (!("feature" %in% names(df))) stop("fold_models must include a 'feature' column.")
    if (!("coef" %in% names(df))) stop("fold_models must include a 'coef' column.")

    df$fold <- factor(df$fold, levels = sort(unique(df$fold)))

    # Remove zero coefficients
    df <- dplyr::filter(df, coef != 0)

    # Summarize feature-level medians
    feature_summary <- dplyr::arrange(
        dplyr::filter(
            dplyr::summarise(
                dplyr::group_by(df, feature),
                any_pos = any(coef > 0),
                nnz = sum(coef != 0),
                med = stats::median(coef),
                .groups = "drop"
            ),
            (!positive_only | any_pos),
            nnz >= min_nonzero
        ),
        dplyr::desc(med)
    )

    if (!is.na(top_n) && is.finite(top_n)) {
        feature_summary <- utils::head(feature_summary, top_n)
    }

    # Keep only selected features (semi_join filters rows only; no columns added)
    df_plot <- dplyr::semi_join(df, feature_summary, by = "feature")

    # Set plotting order for features (largest medians at top)
    df_plot$feature <- factor(df_plot$feature, levels = rev(feature_summary$feature))

    # Median ticks: use feature_summary (contains 'med'), restricted to plotted features
    median_lines <- dplyr::filter(feature_summary, feature %in% levels(df_plot$feature))
    median_lines$feature <- factor(median_lines$feature, levels = levels(df_plot$feature))

    ggplot2::ggplot(df_plot, ggplot2::aes(x = coef, y = feature, color = fold)) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
        ggplot2::geom_point(
            size = 2.2,
            alpha = 0.8,
            position = ggplot2::position_jitter(height = 0.15, width = 0)
        ) +
        # Draw medians as horizontal tick marks at x = med for each feature row
        ggplot2::geom_point(
            data = median_lines,
            ggplot2::aes(x = med, y = feature),
            inherit.aes = FALSE,
            color = "black",
            shape = 95,   # horizontal line glyph (robust on discrete y)
            size = 1.5
        ) +
        ggplot2::scale_color_brewer(palette = "Dark2") +
        ggplot2::labs(
            title = paste0(
                "Elastic Net coefficients across folds",
                if (!is.null(cellType)) paste0(" (", cellType, ")")
            ),
            x = "Coefficient",
            y = NULL,
            color = "Fold"
        ) +
        ggplot2::theme_classic(base_size = 12) +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 9),
            axis.ticks.y = ggplot2::element_blank(),
            axis.ticks.length.y = grid::unit(0, "pt"),
            legend.position = "top",
            panel.grid.major.x = ggplot2::element_line(color = "grey90", linewidth = 0.3)
        )
}


# Compare elastic net model coefficients with DE logFC (nonzero only)
#final_model=r$cv_model$final_model; age_de_results=r$age_de_results; cellType=cellType
plot_model_vs_de <- function(final_model, age_de_results, cellType = NULL, feature="logFC") {

    df <- dplyr::inner_join(
        dplyr::select(
            dplyr::filter(final_model, coef != 0),
            feature,
            model_coef = coef
        ),
        dplyr::select(
            tibble::rownames_to_column(age_de_results, "feature"),
            feature,
            logFC
        ),
        by = "feature"
    )

    # fraction with same direction
    same_direction <- sign(df$model_coef) == sign(df$logFC)
    n_total <- nrow(df)
    n_total_de<-length(which(age_de_results$adj.P.Val<=0.05))
    frac_same <- if (n_total > 0) mean(same_direction) else NA_real_

    strTitle <- if (!is.null(cellType)) cellType else ""
    strTitle <- paste0(strTitle, sprintf("sign %.1f%% agreement - genes with non-zero coefficients [%d]", 100 * frac_same, n_total))
    strTitle <- paste0(strTitle, sprintf("\ngenes significant in age DE: [%d]", n_total_de))

    #Make R CMD CHECK happy
    logFC <- model_coef <- NULL

    ggplot(df, aes(x = logFC, y = model_coef)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
        geom_point(alpha = 0.8, color = "steelblue", size = 2) +
        geom_smooth(method = "lm", formula = y ~ x, color = "black", linewidth = 0.6, se = FALSE) +
        labs(
            title = strTitle,
            x = "Differential Expression logFC (age effect)",
            y = "Elastic Net model coefficient"
        ) +
        theme_classic(base_size = 12)
}

#' Plot donor-level Monte Carlo age predictions with optional bias curve overlay
#'
#' Plots donor chronological age (x) against either:
#' \itemize{
#'   \item \code{pred_mean}: donor mean out-of-fold prediction (default), or
#'   \item \code{pred_mean_corrected}: bias-corrected prediction (post hoc)
#' }
#'
#' When \code{gam_fit_df} is provided and \code{y_var = "pred_mean"}, the bias
#' curve (evaluated on an age grid) is overlaid as a red line to visualize the
#' age-dependent regression-to-the-mean bias. The curve itself is not intended
#' for interpretation.
#'
#' Error bars reflect Monte Carlo variability of out-of-fold predictions across
#' repeated splits (\code{resid_sd}) and remain valid after deterministic mean
#' shifts. Therefore, error bars are drawn when \code{resid_sd} is available for
#' both supported y-variables.
#'
#' @param donor_predictions data.frame with one row per donor. Must contain
#'   \code{age}, \code{pred_mean}. If \code{y_var = "pred_mean_corrected"}, must
#'   contain \code{pred_mean_corrected}. If \code{resid_sd} is present, it is used
#'   for error bars. Also expects \code{cell_type} and \code{region} (single
#'   unique values) for the plot title.
#' @param gam_fit_df Optional data.frame with columns \code{age} and
#'   \code{gam_pred}. Only used to overlay the bias curve when
#'   \code{y_var = "pred_mean"}.
#' @param y_var Either \code{"pred_mean"} (default) or \code{"pred_mean_corrected"}.
#' @param color_var Column in \code{donor_predictions} to map to point color.
#'   Defaults to \code{"resid_mean"}.
#' @param legend_title Legend title for the color scale. Defaults to \code{"Residuals"}.
#' @param alpha_points Alpha for points.
#' @param errorbar_width Width for \code{geom_errorbar()}.
#'
#' @return A ggplot object.
#' @export
plot_mc_donor_predictions <- function(donor_predictions,
                                      gam_fit_df = NULL,
                                      y_var = "pred_mean",
                                      color_var = "resid_mean",
                                      legend_title = "Residuals",
                                      alpha_points = 0.8,
                                      errorbar_width = 0.1) {

    stopifnot(is.data.frame(donor_predictions))

    allowed_y <- c("pred_mean", "pred_mean_corrected")
    if (!(y_var %in% allowed_y)) {
        stop("y_var must be one of: ", paste(allowed_y, collapse = ", "))
    }

    if (!("age" %in% colnames(donor_predictions))) {
        stop("donor_predictions must contain column 'age'")
    }
    if (!("pred_mean" %in% colnames(donor_predictions))) {
        stop("donor_predictions must contain column 'pred_mean'")
    }
    if (!(y_var %in% colnames(donor_predictions))) {
        stop("y_var not found in donor_predictions: ", y_var)
    }
    if (!(color_var %in% colnames(donor_predictions))) {
        stop("color_var not found in donor_predictions: ", color_var)
    }

    ct <- unique(donor_predictions$cell_type)
    rg <- unique(donor_predictions$region)
    if (length(ct) != 1 || length(rg) != 1) {
        stop("donor_predictions must have exactly one unique cell_type and one unique region")
    }

    title_str <- paste0(
        "Monte Carlo CV donor age predictions\n",
        "Cell type: ", ct,
        if (!is.na(rg)) paste0(" | Region: ", rg) else ""
    )

    y_lab <- switch(
        y_var,
        pred_mean = "Predicted age (MC mean)",
        pred_mean_corrected = "Predicted age (MC mean corrected)",
        y_var
    )

    rng <- range(c(donor_predictions$age, donor_predictions[[y_var]]), na.rm = TRUE)
    pad <- 0.05 * diff(rng)
    lims <- c(rng[1] - pad, rng[2] + pad)

    # MAKE R CMD CHECK happy (ggplot NSE)
    age <- ymin <- ymax <- gam_pred <- NULL

    p <- ggplot2::ggplot(
        donor_predictions,
        ggplot2::aes(x = age, y = .data[[y_var]], color = .data[[color_var]])
    )

    if ("resid_sd" %in% colnames(donor_predictions)) {
        df_err <- donor_predictions
        df_err$ymin <- df_err[[y_var]] - df_err$resid_sd
        df_err$ymax <- df_err[[y_var]] + df_err$resid_sd

        p <- p + ggplot2::geom_errorbar(
            data = df_err,
            ggplot2::aes(x = age, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            width = errorbar_width,
            alpha = 0.4
        )
    }

    p <- p +
        ggplot2::geom_point(size = 2, alpha = alpha_points) +
        ggplot2::scale_color_gradient2(
            low = "steelblue",
            mid = "grey80",
            high = "firebrick",
            midpoint = 0,
            name = legend_title
        ) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = 2, linewidth = 0.8) +
        ggplot2::coord_equal(xlim = lims, ylim = lims) +
        ggplot2::labs(
            x = "Chronological age",
            y = y_lab,
            title = title_str
        ) +
        ggplot2::theme_bw(base_size = 10)

    if (y_var == "pred_mean" && !is.null(gam_fit_df)) {
        stopifnot(is.data.frame(gam_fit_df))
        stopifnot(all(c("age", "gam_pred") %in% colnames(gam_fit_df)))

        df_line <- gam_fit_df[order(gam_fit_df$age), c("age", "gam_pred")]
        p <- p + ggplot2::geom_line(
            data = df_line,
            ggplot2::aes(x = age, y = gam_pred),
            inherit.aes = FALSE,
            color = "red",
            linewidth = 0.9
        )
    }

    p
}


#' Generate a donor-level prediction report with raw and bias-corrected panels
#'
#' This function produces a two-panel report summarizing donor-level model
#' predictions before and after age-bias correction. The top panel displays
#' raw predicted values (which may exhibit age-dependent bias), and the
#' bottom panel displays bias-corrected predictions. In both panels,
#' points are colored by residuals to visualize model deviation patterns.
#'
#' The input data should be for a single grouping the model was trained on,
#' for example cell type and region.
#'
#' Internally, the function calls \code{plot_mc_donor_predictions()} twice:
#' once using the original predicted values and GAM fit (for raw predictions),
#' and once using the bias-corrected predictions without overlaying the GAM fit.
#' The two plots are then vertically stacked using \code{cowplot::plot_grid()}.
#'
#' @param donor_predictions A data.frame containing donor-level predictions
#'   and residuals. Must contain columns \code{pred_mean},
#'   \code{pred_mean_corrected}, \code{resid_mean}, and
#'   \code{resid_mean_corrected}, along with any additional columns required
#'   by \code{plot_mc_donor_predictions()}.
#' @param gam_fit_df A data.frame containing fitted GAM values for overlay
#'   in the raw prediction panel. This is passed directly to
#'   \code{plot_mc_donor_predictions()}. Set to \code{NULL} if no GAM fit
#'   should be shown.
#'
#' @return A named list containing:
#'   \itemize{
#'     \item \code{raw_plot}: The ggplot object showing raw predictions.
#'     \item \code{corrected_plot}: The ggplot object showing bias-corrected predictions.
#'     \item \code{combined_plot}: A vertically stacked ggplot object containing both panels.
#'   }
#'
#'
#' @export
plot_mc_donor_predictions_report <- function(donor_predictions,
                                             gam_fit_df) {

    # Silence R CMD CHECK notes
    resid_mean <- resid_mean_corrected <- NULL

    p_raw <- plot_mc_donor_predictions(
        donor_predictions,
        gam_fit_df = gam_fit_df,
        y_var = "pred_mean",
        color_var = "resid_mean",
        legend_title = "Residuals"
    ) +
        ggplot2::labs(subtitle = "Raw predictions (with age-dependent bias)")

    p_corr <- plot_mc_donor_predictions(
        donor_predictions,
        gam_fit_df = NULL,
        y_var = "pred_mean_corrected",
        color_var = "resid_mean_corrected",
        legend_title = "Residuals"
    ) +
        ggplot2::labs(subtitle = "Bias-corrected predictions")

    p_both=cowplot::plot_grid(
        p_raw,
        p_corr,
        ncol = 1,
        nrow = 2,
        align = "v"
    )

    result=list(raw_plot=p_raw,
                corrected_plot=p_corr,
                combined_plot=p_both)
    return (result)
}




#' Create a custom ggplot2 theme with smaller plot titles
#'
#' @param base_size Base font size for the theme (default 12).
#' @param title_size Relative size for the plot title (default 10).
#' @param axis_text_size Size of axis text (default 9).
#' @param legend_text_size Size of legend text (default 9).
#'
#' @return A ggplot2 theme object.
#' @export
custom_theme_small_title <- function(base_size = 12, title_size = 10,
                                     axis_text_size = 9, legend_text_size = 9) {
    ggplot2::theme_classic(base_size = base_size) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(size = title_size, face = "bold", hjust = 0.5),
            axis.title = ggplot2::element_text(size = base_size * 0.9),
            axis.text  = ggplot2::element_text(size = axis_text_size),
            legend.title = ggplot2::element_text(size = legend_text_size, face = "bold"),
            legend.text  = ggplot2::element_text(size = legend_text_size),
            strip.text = ggplot2::element_text(size = base_size * 0.9, face = "bold")
        )
}



#######################
# WRITE OUTPUTS
#######################

.round_columns <- function(df, digits_map) {
    for (nm in names(digits_map)) {
        if (nm %in% colnames(df)) {
            df[[nm]] <- round(df[[nm]], digits_map[[nm]])
        }
    }
    df
}

get_output_basename <- function(pdf_file, default = "age_prediction_results") {

    if (is.null(pdf_file) || is.na(pdf_file) || !nzchar(pdf_file)) {
        return(default)
    }

    sub("\\.pdf$", "", basename(pdf_file))
}


extract_age_outputs <- function(r, cellType, region = NA_character_) {

    meta_cols <- c("cell_type", "region")

    # -----------------------------
    # Donor predictions (MC OOF)
    # -----------------------------
    base_cols <- c(
        "donor", "age", "num_features", "num_nuclei", "num_umis", "n_oof", "pred_mean",
        "resid_mean", "resid_median", "resid_sd"
    )

    # New/optional correction columns (preferred names)
    corr_cols <- intersect(
        c("pred_mean_corrected", "resid_mean_corrected"),
        names(r$final_oof_predictions)
    )

    preds <- r$final_oof_predictions[
        , c(base_cols, corr_cols),
        drop = FALSE
    ]

    preds$cell_type <- cellType
    preds$region <- region
    preds <- preds[, c(meta_cols, setdiff(colnames(preds), meta_cols)), drop = FALSE]

    # -----------------------------
    # Model coefficients
    # -----------------------------
    coefs <- r$final_cv_model$final_model
    coefs$cell_type <- cellType
    coefs$region <- region
    coefs <- coefs[, c(meta_cols, setdiff(colnames(coefs), meta_cols)), drop = FALSE]

    # -----------------------------
    # Summary metrics
    # -----------------------------
    mets <- rbind(
        cbind(set = "Final MC OOF (all donors; mean pred)", r$final_oof_metrics),
        cbind(set = "Inner CV OOF (80% train)", r$cv_model$overall_metrics),
        cbind(set = "Outer holdout (20% donors)", r$test_set_metrics)
    )

    mets$cell_type <- cellType
    mets$region <- region
    mets <- mets[, c(meta_cols, setdiff(colnames(mets), meta_cols)), drop = FALSE]

    # -----------------------------
    # Per-fold metrics (inner CV)
    # -----------------------------
    fold_mets <- r$cv_model$per_fold_metrics
    fold_mets$cell_type <- cellType
    fold_mets$region <- region
    fold_mets <- fold_mets[, c(meta_cols, setdiff(colnames(fold_mets), meta_cols)), drop = FALSE]

    # -----------------------------
    # GAM fit curve df (age grid) for plotting
    # -----------------------------
    gam_fit_df <- NULL
    if (!is.null(r$gam_fit_df)) {
        gam_fit_df <- r$gam_fit_df
    } else if (!is.null(r$final_oof_gam_fit_df)) {
        gam_fit_df <- r$final_oof_gam_fit_df
    } else if (!is.null(r$final_repeat$gam_fit_df)) {
        gam_fit_df <- r$final_repeat$gam_fit_df
    }

    list(
        donor_predictions  = preds,
        model_coefficients = coefs,
        model_metrics      = mets,
        per_fold_metrics   = fold_mets,
        gam_fit_df         = gam_fit_df
    )
}


write_age_outputs_all <- function(outputs, result_dir, output_basename) {



    donor_predictions <- do.call(rbind, lapply(outputs, function(x) x$donor_predictions))
    model_coefficients <- do.call(rbind, lapply(outputs, function(x) x$model_coefficients))
    model_metrics <- do.call(rbind, lapply(outputs, function(x) x$model_metrics))
    per_fold_metrics <- do.call(rbind, lapply(outputs, function(x) x$per_fold_metrics))
    gam_fit_dfs <- do.call(rbind, lapply(outputs, function(x) x$gam_fit_df))
    rownames (gam_fit_dfs)=NULL

    #round the outputs.
    #the number of digits to round to
    donor_pred_digits <- list(pred_mean = 3, resid_mean = 3, resid_median = 3, resid_sd = 3)
    metrics_digits <- list(r = 3, median_abs_error = 3, mean_abs_error = 3)

    #round the outputs
    donor_predictions <- .round_columns(donor_predictions, donor_pred_digits)
    model_metrics <- .round_columns(model_metrics, metrics_digits)
    per_fold_metrics <- .round_columns(per_fold_metrics, metrics_digits)

    out_pred <- file.path(result_dir,
                          paste0(output_basename, "_donor_predictions.txt"))
    out_coef <- file.path(result_dir,
                          paste0(output_basename, "_model_coefficients.txt"))
    out_met  <- file.path(result_dir,
                          paste0(output_basename, "_model_metrics.txt"))
    out_fmet <- file.path(result_dir,
                          paste0(output_basename, "_model_per_fold_metrics.txt"))

    out_gam_fit <- file.path(result_dir,
                          paste0(output_basename, "_gam_fit.txt"))

    write.table(donor_predictions, file = out_pred,
                sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = TRUE)

    write.table(model_coefficients, file = out_coef,
                sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = TRUE)

    write.table(model_metrics, file = out_met,
                sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = TRUE)

    write.table(per_fold_metrics, file = out_fmet,
                sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = TRUE)

    write.table(gam_fit_dfs, file = out_gam_fit,
                sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = TRUE)


    invisible(list(
        donor_predictions = out_pred,
        model_coefficients = out_coef,
        model_metrics = out_met,
        per_fold_metrics = out_fmet
    ))
}



####################
# READ IN DATA AND COLLAPSE TO DONOR OBSERVATIONS
######################

#' Collapse counts in a DGEList to one column per donor
#'
#' Aggregates raw counts in an edgeR \code{DGEList} so that each donor
#' contributes a single column whose counts are the sum of all samples
#' belonging to that donor. This is useful when a donor has multiple
#' sequencing libraries or replicates and you need one observation per donor
#' for downstream modeling.
#'
#' @param dge A \code{DGEList} containing raw counts and a \code{samples} data.frame.
#' @param donor_col Character scalar giving the column in \code{dge$samples}
#'   that identifies the donor for each sample. Each donor must appear in one
#'   or more rows.
#' @param keep_cols Character vector of column names in \code{dge$samples}
#'   to retain in the output. For each donor, the function checks that these
#'   columns have exactly one unique (non-NA) value across all donor samples.
#'   An error is thrown if a column has multiple distinct values within a donor.
#' @param sum_cols Character vector of numeric columns in \code{dge$samples} to
#' sum within donor (e.g. \code{"num_nuclei"}).
#'
#' @return A new \code{DGEList} where:
#'   \itemize{
#'     \item \code{$counts} is a matrix of gene-by-donor integer counts
#'       summed across all samples of each donor.
#'     \item \code{$samples} contains one row per donor, restricted to
#'       \code{keep_cols}, plus standard edgeR fields \code{group},
#'       \code{lib.size}, and \code{norm.factors}.
#'   }
#'   Library sizes are recomputed as column sums of the collapsed counts,
#'   and normalization factors are set to 1.
collapse_by_donor <- function(dge, donor_col = "donor",
                              keep_cols = character(0),
                              sum_cols  = character(0)) {
    stopifnot(is.list(dge), !is.null(dge$counts), !is.null(dge$samples),
              donor_col %in% colnames(dge$samples))

    smp <- dge$samples
    donors <- as.character(smp[[donor_col]])
    if (anyNA(donors)) stop("donor_col has NA values.")
    f <- factor(donors)

    # donor indicator (samples x donors)
    X <- model.matrix(~ 0 + f)
    colnames(X) <- levels(f)

    # sum counts across samples of each donor (genes x donors)
    C <- as.matrix(dge$counts) %*% X
    storage.mode(C) <- "integer"

    # helper: initialize an NA vector matching the class of x
    init_na_like <- function(x, n) {
        if (is.integer(x))            return(rep(NA_integer_, n))
        if (is.numeric(x))            return(rep(NA_real_, n))
        if (is.logical(x))            return(rep(NA, n))
        return(rep(NA_character_, n)) # character or factor handled as character
    }

    # build samples data.frame with keep_cols and validate uniqueness per donor
    if (length(keep_cols)) {
        missing <- setdiff(keep_cols, colnames(smp))
        if (length(missing)) stop("keep_cols not found in dge$samples: ", paste(missing, collapse = ", "))

        nD <- nlevels(f)
        out_list <- vector("list", length(keep_cols))
        names(out_list) <- keep_cols

        for (j in seq_along(keep_cols)) {
            colname <- keep_cols[j]
            colj <- smp[[colname]]
            if (is.factor(colj)) colj <- as.character(colj)

            v <- init_na_like(colj, nD)
            names(v) <- levels(f)

            for (d in levels(f)) {
                vals <- unique(colj[f == d])
                vals <- vals[!is.na(vals)]
                if (length(vals) != 1) {
                    stop(sprintf("Column '%s' not unique within donor '%s': values = {%s}",
                                 colname, d, paste(utils::head(vals, 10), collapse = ", ")))
                }
                v[d] <- vals
            }
            out_list[[j]] <- v
        }
        samples_keep <- as.data.frame(out_list, stringsAsFactors = FALSE, row.names = levels(f))
    } else {
        samples_keep <- data.frame(row.names = levels(f))
    }

    # add donor-wise sums for selected numeric columns
    if (length(sum_cols)) {
        missing <- setdiff(sum_cols, colnames(smp))
        if (length(missing)) stop("sum_cols not found in dge$samples: ", paste(missing, collapse = ", "))

        for (colname in sum_cols) {
            x <- smp[[colname]]
            if (!is.numeric(x)) {
                stop("sum_cols must be numeric; column '", colname, "' has class ", class(x)[1], ".")
            }
            if (anyNA(x)) {
                stop("sum_cols column '", colname, "' contains NA values.")
            }

            summed <- as.numeric(crossprod(x, X))
            names(summed) <- colnames(X)

            samples_keep[[colname]] <- summed[rownames(samples_keep)]
        }
    }

    # required edgeR-like fields
    lib.size <- colSums(C)
    samples_out <- data.frame(
        group = 1L,
        lib.size = lib.size,
        norm.factors = rep(1, length(lib.size)),
        samples_keep,
        row.names = colnames(C),
        check.names = FALSE
    )

    dge_out <- dge
    dge_out$counts  <- C
    dge_out$samples <- samples_out
    dge_out
}


#' Load age DE results for a given cell type and region
#'
#' Loads differential expression results from a specified directory
#' for a given cell type and region. Filters results to those with adjusted
#' p-value below a specified FDR threshold.
#'
#' The expected encoding of DE results filenames is:
#' - General age DE results: \code{<cellType>_age_DE_results.txt}
#' - Region-specific age DE results: \code{<cellType>_age_<region>_DE_results.txt}
#'
#' @param cellType Character scalar; cell type name.
#' @param age_de_results_dir Character scalar; path to directory containing DE results files.
#' @param region Optional character scalar; region name. If NULL, looks for
#'  general age DE results file for the cell type. If provided, looks for
#'  region-specific DE results file.
#' @param fdr_threshold Numeric scalar; FDR threshold for filtering DE results
#'  (default 0.05).
#' @return A data.frame of DE results filtered to adjusted p-value <= fdr_threshold.
#'  If age_de_results_dir is NULL, returns NULL.
#'
get_age_de_results <- function(cellType,
                               age_de_results_dir,
                               region = NULL,
                               fdr_threshold = 0.05) {

    if (is.null(age_de_results_dir))
        return(NULL)

    # Construct expected filename
    if (is.null(region)) {
        # e.g. astrocyte__age_DE_results.txt
        expectedFileName <- paste(cellType, "__age", "_DE_results.txt", sep = "")
    } else {
        # e.g. astrocyte__CaH__age_DE_results.txt
        expectedFileName <- paste(cellType, "__", region, "__age", "_DE_results.txt", sep = "")
    }

    files <- list.files(path = age_de_results_dir, full.names = TRUE)

    match <- grep(paste0("^", expectedFileName, "$"), basename(files), value = TRUE)

    if (length(match) != 1) {
        stop(
            "Did not find exactly one DE results file for cell type ",
            cellType,
            if (!is.null(region)) paste0(" and region ", region),
            " in directory ", age_de_results_dir, ". Expected file: ", expectedFileName
        )
    }

    f <- file.path(age_de_results_dir, match)
    de_results <- utils::read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    de_results <- de_results[de_results$adj.P.Val <= fdr_threshold, , drop = FALSE]

    de_results
}


#' Filter an edgeR DGEList to autosomal genes
#'
#' Subsets an \code{edgeR::DGEList} to genes annotated as autosomal. Autosomes are
#' determined from a contig YAML file (contig \eqn{\rightarrow} class mapping),
#' and genes are identified from a reduced GTF-like TSV containing gene records.
#'
#' If either \code{contig_yaml_file} or \code{reduced_gtf_file} is \code{NULL},
#' the input \code{dge} is returned unchanged.
#'
#' @param dge An \code{edgeR::DGEList}. Gene identifiers are expected to be in
#'   \code{rownames(dge)}.
#' @param contig_yaml_file Path to a YAML file mapping contig names to contig
#'   classes. Contigs with value \code{"autosome"} are retained.
#' @param reduced_gtf_file Path to a tab-delimited file (reduced GTF-like) with
#'   at least the columns \code{chr}, \code{annotationType}, and \code{gene_name}.
#'
#' @return An \code{edgeR::DGEList} containing only rows whose gene names match
#' autosomal genes in \code{reduced_gtf_file}. Library sizes are preserved
#' (\code{keep.lib.sizes = TRUE}).
#' @export
filter_dge_to_autosomes<-function (dge, contig_yaml_file, reduced_gtf_file) {
    if (is.null(contig_yaml_file) | is.null(reduced_gtf_file)) {
        return (dge)
    }
    logger::log_info("Filtering DGEList to autosomal genes")
    z=unlist (yaml::yaml.load_file(contig_yaml_file))
    autosomes=names (z[z=="autosome"])

    gtf=data.table::fread(reduced_gtf_file, header=T, sep="\t", stringsAsFactors = FALSE)
    gtf=gtf[gtf$chr %in% autosomes & gtf$annotationType=="gene", ]

    autosomal_genes=unique (gtf$gene_name)
    dge_filtered=dge[rownames(dge)%in%autosomal_genes, , keep.lib.sizes = TRUE]
    logger::log_info("Filtered from {nrow(dge)} to {nrow(dge_filtered)} genes")
    return (dge_filtered)
}

filter_dge_to_donor_age<-function (dge, donor_age_range) {
    if (is.null(donor_age_range)) {
        return (dge)
    }

    logger::log_info("Filtering DGEList to donors with ages in range [{donor_age_range[1]}, {donor_age_range[2]}] decades")
    idx=dge$samples$age>=donor_age_range[1] & dge$samples$age<=donor_age_range[2]
    dge_filtered=dge[, idx, keep.lib.sizes = TRUE]
    logger::log_info("Filtered from {ncol(dge)} to {ncol(dge_filtered)} pseudobulked donor samples")
    return (dge_filtered)
}

filter_dge_to_gene_functions<-function (dge, gene_functions=c("protein_coding"), gtf_path) {
    if (is.null(gene_functions) | is.null(gtf_path)) {
        return (dge)
    }

    logger::log_info("Filtering DGEList to genes with functions: {paste(gene_functions, collapse=', ')}")
    gtf=parse_gtf_genes_df(gtf_path)
    gtf=gtf[gtf$gene_type %in% gene_functions, ]

    selected_genes=unique (gtf$gene_name)
    dge_filtered=dge[rownames(dge)%in%selected_genes, , keep.lib.sizes = TRUE]
    logger::log_info("Filtered from {nrow(dge)} to {nrow(dge_filtered)} genes")
    return (dge_filtered)
}

parse_gtf_genes_df <- function(gtf_path) {
    gr <- import(gtf_path)
    df <- as.data.frame(gr)

    # keep only gene-level rows
    df <- df[df$type == "gene", ]

    # normalize ENSG IDs (drop version suffix)
    df$ENSG <- sub("\\.\\d+$", "", as.character(df$gene_id))

    # handle type naming differences
    if ("gene_type" %in% names(df)) {
        df$gene_type <- as.character(df$gene_type)
    } else if ("gene_biotype" %in% names(df)) {
        df$gene_type <- as.character(df$gene_biotype)
    } else {
        df$gene_type <- NA_character_
    }

    df$gene_name <- as.character(df$gene_name)

    df=df[, c("seqnames", "gene_type", "gene_name", "ENSG")]
    return (df)
}

