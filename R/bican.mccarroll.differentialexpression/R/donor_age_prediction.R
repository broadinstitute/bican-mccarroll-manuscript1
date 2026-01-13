library(edgeR)
library(glmnet)
library(ggplot2)
library(logger)
library (cowplot)
library (dplyr)
#
data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/metacells/LEVEL_2"
data_name="donor_rxn_DGEList"

age_de_results_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results/LEVEL_2/sex_age/cell_type"
result_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/age_prediction/age_prediction_alpha_0"
outPDFFile="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/age_prediction/age_prediction_alpha_0/age_prediction_results_alpha0.pdf"
# #
# #filtering to autosomal genes
contig_yaml_file="/broad/mccarroll/software/metadata/individual_reference/GRCh38_ensembl_v43/GRCh38_ensembl_v43.contig_groups.yaml"
reduced_gtf_file="/broad/mccarroll/software/metadata/individual_reference/GRCh38_ensembl_v43/GRCh38_ensembl_v43.reduced.gtf"
#
# #filtering to gene functional types.  #TODO collapse this and autosomal processes.
# gtf_path="/broad/mccarroll/software/metadata/individual_reference/GRCh38_ensembl_v43/GRCh38_ensembl_v43.gtf"
# #

donor_col = "donor"
age_col = "age"
seed =12345; fdr_threshold=0.05; optimize_alpha=FALSE; alpha_fixed=0

#run Emi's data:
# data_dir="/broad/mccarroll/dropulation/analysis/cellarium_upload/SNAP200_freeze1/metacells"
# data_name="donor_rxn_DGEList"

#still use the same set of overall features from BICAN.
#I don't want to have to revisit this unless neccesary.
# age_de_results_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age/cell_type"
# age_de_results_dir="/broad/mccarroll/dropulation/analysis/SNAP200/differential_expression/sex_age"
# result_dir="/broad/mccarroll/dropulation/analysis/SNAP200/differential_expression/age_prediction"
# outPDFFile="/broad/mccarroll/dropulation/analysis/SNAP200/differential_expression/age_prediction/age_prediction_results_snap200_DE_genes.pdf"

predict_age_by_celltype<-function (data_dir, data_name, result_dir, age_de_results_dir, outPDFFile=NULL,
                                  contig_yaml_file, reduced_gtf_file,
                                  donor_col = "donor", age_col = "age",
                                  seed =42, fdr_threshold=0.05,
                                  optimize_alpha=FALSE, alpha_fixed=0.5) {

    #validate the output directory exists
    if (!dir.exists(result_dir)) {
        stop("Result directory does not exist: ", result_dir)
    }

    #load the DGEList and prepare the data
    d=bican.mccarroll.differentialexpression::prepare_data_for_differential_expression(data_dir, data_name, randVars=c(), fixedVars=c())
    dge=d$dge

    #filter the genes to autosomal only
    dge=filter_dge_to_autosomes (dge, contig_yaml_file, reduced_gtf_file)

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

        r=predict_age_celltype(cellType, dge, retained_features=retained_features, donor_col = donor_col, age_de_results_dir=age_de_results_dir, fdr_threshold=fdr_threshold, optimize_alpha=optimize_alpha, alpha_fixed=alpha_fixed, seed=seed)
        test_set_metrics_list[[cellType]]=r$test_set_metrics

        # QC plots
        plot_list=age_prediction_qc_plots(r, cellType)

        for (p in plot_list) {
            if (!is.null(p))
                print(p)
        }

        # Keep only what will be written
        outputs_list[[cellType]] <- extract_age_outputs(r, cellType = cellType, region = NA_character_)

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



age_prediction_qc_plots<-function (r, cellType) {
    age_dist_plot <- plot_age_hist_stacked(r$dge_train, r$dge_test, r$cellType)
    pred_age_plot_train <- plot_age_predictions(r$cv_model$donor_age_predictions, cellType, titleStr=paste("TRAIN", cellType, "\n"))
    pred_age_plot_test <- plot_age_predictions(r$test_set_predictions, cellType, titleStr=paste("TEST", cellType, "\n [20% held out data]"))
    perf_plot<-plot_model_performance(
        cv_metrics = r$cv_model$overall_metrics,
        per_fold_metrics = r$cv_model$per_fold_metrics,
        test_metrics = r$test_set_metrics,
        cellType = cellType
    )

    p=cowplot::plot_grid(pred_age_plot_train, pred_age_plot_test, age_dist_plot, perf_plot, ncol=2)

    coef1=plot_fold_coefficients(r$cv_model$fold_models, positive_only = FALSE, min_nonzero = 1, top_n=NA, cellType=cellType)

    #plot the age DE results if available against the model average coefficients
    model_vs_de_plot=NULL
    if (!is.null(r$age_de_results))
        model_vs_de_plot=plot_model_vs_de(r$cv_model$final_model, r$age_de_results, cellType = NULL)

    return (list(main_qc_plot=p, coefficient_plot=coef1, model_vs_de_plot=model_vs_de_plot))

}

compare_age_residuals_to_donor_coefficients<-function () {
    #Load the DGEList, get the sample metadata
    #simplify to per-donor covariates
    #compare covariates to each cell type donor age residuals (prediction - chronological age)

    d=bican.mccarroll.differentialexpression::prepare_data_for_differential_expression(data_dir, data_name, randVars=c(), fixedVars=c())
    dge=d$dge

}

#' Fit a donor-level age prediction model for a single cell type
#'
#' Trains an elastic net (or ridge/lasso depending on \code{alpha_fixed}) age
#' prediction model from gene expression for one cell type. The function:
#' \enumerate{
#'   \item Subsets the input \code{DGEList} to the requested cell type.
#'   \item Collapses multiple observations per donor (e.g., across regions) into a
#'         single donor-level profile via \code{collapse_by_donor()}.
#'   \item Applies standard expression QC and gene filtering steps.
#'   \item Optionally filters genes to a set significant in external age DE results.
#'   \item Performs an outer 80/20 donor-level holdout split for validation.
#'   \item Trains the model on the 80\% training donors using inner stratified
#'         K-fold CV to choose the regularization strength (and optionally alpha).
#'   \item Evaluates the model on the held-out 20\% donors.
#'   \item Computes "final" donor residuals by running repeated stratified 80/20
#'         Monte Carlo splits on the full donor-level dataset and summarizing the
#'         per-donor out-of-fold residual distribution.
#' }
#'
#' @param cellType Character scalar; the cell type label used to subset
#'   \code{dge$samples$cell_type}.
#' @param dge A \code{DGEList} containing counts and sample metadata. Must contain
#'   \code{dge$samples$cell_type}, donor IDs, and ages.
#' @param retained_features Character vector of sample-level columns to retain
#'   during donor collapsing (passed to \code{collapse_by_donor()}).
#' @param donor_col Column name in \code{dge$samples} containing donor IDs.
#' @param age_de_results_dir Directory containing age differential expression
#'   results used by \code{get_age_de_results()}.
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
#' @param seed Integer random seed used for the outer holdout split, inner CV fold
#'   construction, and Monte Carlo repeated splits.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{cell_type}: cell type name.
#'   \item \code{dge_train}: donor-level training \code{DGEList} (80\% donors).
#'   \item \code{dge_test}: donor-level held-out \code{DGEList} (20\% donors).
#'   \item \code{cv_model}: result of \code{train_enet_cv()} or
#'         \code{train_enet_cv_optimize_alpha()} on \code{dge_train}.
#'   \item \code{test_set_predictions}: data.frame of predictions on held-out
#'         donors with columns \code{donor}, \code{age}, and \code{pred}.
#'   \item \code{test_set_metrics}: performance metrics on held-out donors as
#'         returned by \code{compute_age_metrics()}.
#'   \item \code{age_de_results}: age DE results returned by \code{get_age_de_results()}
#'         (or \code{NULL} if unavailable).
#'   \item \code{final_oof_predictions}: per-donor residual summary returned by
#'         \code{compute_final_residuals_repeated_splits()} (one row per donor).
#'   \item \code{final_oof_metrics}: overall metrics computed from the per-donor
#'         average prediction across repeats (\code{final_repeat$avg_pred_metrics}).
#' }
#'
#' @details
#' The "final" donor residuals are computed using repeated stratified 80/20 splits
#' on the full donor-level dataset. For each repeat, donors in the test split are
#' predicted by a model trained on the corresponding training split, and only
#' those out-of-fold predictions contribute to the donor's residual distribution.
#' This reduces dependence of residual estimates on a single fold assignment and
#' can be used to assess residual stability with respect to training-set
#' composition.
#'
#' @export
predict_age_celltype <- function(cellType, dge, retained_features = c("donor", "age"),
                                 donor_col = "donor", age_de_results_dir,
                                 fdr_threshold = 0.05,
                                 optimize_alpha = FALSE, alpha_fixed = 0.5,
                                 mc_repeats=100,
                                 seed = 1) {

    if (optimize_alpha) {
        logger::log_info("Optimizing alpha via grid search.")
    } else {
        logger::log_info(paste("Using fixed alpha =", alpha_fixed))
    }

    dge_cell <- dge[, dge$samples$cell_type == cellType, keep.lib.sizes = TRUE]

    # Collapse across regions and multiple samples per donor
    dge_cell <- collapse_by_donor(dge_cell, donor_col = donor_col, keep_cols = retained_features)

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

    # Get age DE results (optional)
    age_de_results <- get_age_de_results(cellType, age_de_results_dir, fdr_threshold = fdr_threshold)
    if (!is.null(age_de_results)) {
        genes_to_keep <- intersect(rownames(dge_cell), rownames(age_de_results))
        logger::log_info(paste("Filtering to", length(genes_to_keep),
                               "genes significant in age DE results at FDR <=", fdr_threshold))
        dge_cell <- dge_cell[genes_to_keep, ]
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

    # ------------------------------------------------------------------
    # NEW: Many fold repeated-split OOF residual computation on full dataset
    # ------------------------------------------------------------------
    logger::log_info("Computing final repeated-split OOF residuals on all donors (full dataset).")

    final_repeat <- compute_final_residuals_repeated_splits(
        dge_cell,
        age_col = "age",
        donor_col = "donor",
        n_bins = n_bins,
        test_prop = 0.20,
        n_repeats = mc_repeats,
        optimize_alpha = optimize_alpha,
        alpha_fixed = alpha_fixed,
        seed = seed
    )

    final_oof_predictions <- final_repeat$donor_residual_summary
    final_oof_metrics <- final_repeat$avg_pred_metrics

    result <- list(
        cell_type = cellType,
        dge_train = dge_train,
        dge_test = dge_test,

        # Existing outputs
        cv_model = cv_model,
        test_set_predictions = test_set_predictions,
        test_set_metrics = metrics_test,
        age_de_results = age_de_results,

        # New outputs
        final_cv_model = final_cv_model,
        final_oof_predictions = final_oof_predictions,
        final_oof_metrics = final_oof_metrics
    )

    return(result)
}




#########################
# TRAIN AND PREDICT AGE
#########################

#' Train Elastic Net age model with K-fold cross-validation
#'
#' @param dge_train A DGEList object containing log-CPM data for training.
#' @param k_fold_index Named integer vector: donor IDs mapped to fold number.
#' @param age_col Column in dge$samples giving donor age.
#' @param donor_col Column in dge$samples giving donor ID.
#' @param alpha_fixed Elastic net alpha (0=L2 ridge, 1=L1 lasso).
#' @param refit_model_all_data If set to true, refit the final model using all of the data.
#' If set to false, the final model is obtained by averaging the per-fold model coefficients.
#'
#' @param seed Random seed for reproducibility.
#'
#' @return A list with:
#' \itemize{
#'   \item per_fold — metrics for each CV fold.
#'   \item oof — per-donor predicted and true ages.
#'   \item overall — combined metrics across all folds.
#'   \item coef_matrix — feature × fold coefficient table.
#'   \item intercepts — intercept per fold.
#'   \item scaling — mean and sd of each feature per fold.
#' }
#' @export
# Donor-level K-fold elastic-net with per-fold data frames and reusable metrics.
train_enet_cv <- function(dge_train, k_fold_index,
                          age_col = "age", donor_col = "donor",
                          alpha_fixed = 0.5,
                          refit_model_all_data=T,
                          seed = 1) {

    stopifnot("DGEList" %in% class(dge_train))
    set.seed(seed)

    # helpers
    scale_train <- function(M) {
        mu <- colMeans(M); s <- apply(M, 2, sd); s[s == 0] <- 1
        list(X = sweep(sweep(M, 2, mu, "-"), 2, s, "/"), mu = mu, s = s)
    }

    scale_apply <- function(M, mu, s) sweep(sweep(M, 2, mu, "-"), 2, s, "/")

    smp   <- dge_train$samples
    donor <- as.character(smp[[donor_col]])
    age   <- as.numeric(smp[[age_col]])

    stopifnot(!is.null(names(k_fold_index)))
    fold_id <- as.integer(k_fold_index[donor])
    K <- max(fold_id, na.rm = TRUE)

    dge_train <- edgeR::calcNormFactors(dge_train)
    X_all <- t(edgeR::cpm(dge_train, log = TRUE, prior.count = 1))  # samples x genes

    oof_pred <- rep(NA_real_, length(age))
    per_fold <- data.frame(fold=integer(0), n_train=integer(0), n_test=integer(0),
                           lambda=double(0), r=double(0), median_abs_error=double(0), mean_abs_error=double(0))
    fold_predictions <- vector("list", K)
    fold_models <- vector("list", K)

    for (k in seq_len(K)) {
        te <- which(fold_id == k)
        tr <- which(fold_id != k)

        Xtr_s <- scale_train(X_all[tr, , drop = FALSE])
        Xte   <- scale_apply(X_all[te, , drop = FALSE], Xtr_s$mu, Xtr_s$s)
        ytr <- age[tr]; yte <- age[te]

        cv <- glmnet::cv.glmnet(Xtr_s$X, ytr, family = "gaussian",
                                alpha = alpha_fixed, standardize = FALSE,
                                type.measure = "mae")
        fit <- glmnet::glmnet(Xtr_s$X, ytr, family = "gaussian",
                              alpha = alpha_fixed, lambda = cv$lambda.min,
                              standardize = FALSE)

        pred <- as.numeric(predict(fit, Xte))
        oof_pred[te] <- pred

        m <- compute_age_metrics(pred, yte)
        per_fold[nrow(per_fold)+1, ] <- cbind(fold = k, n_train = length(tr), n_test = length(te), cv$lambda.min, m)

        fold_predictions[[k]] <- data.frame(
            donor = donor[te], age = yte, pred = pred, fold = k,
            row.names = rownames(smp)[te], check.names = FALSE
        )

        # per-fold model dataframe: coefficients + scaling + intercept
        feats_k <- rownames(fit$beta)
        coef_df_k <- data.frame(
            feature = feats_k,
            coef = as.numeric(fit$beta),
            mean_train = Xtr_s$mu[feats_k],
            sd_train = Xtr_s$s[feats_k],
            intercept = as.numeric(fit$a0),
            fold = k,
            row.names = NULL, check.names = FALSE
        )
        fold_models[[k]] <- coef_df_k
    }

    overall <- compute_age_metrics(oof_pred, age)
    oof_df <- do.call(rbind, fold_predictions)
    model_df <- do.call(rbind, fold_models)  # long form: one row per (feature, fold)

    #refit the model using all of the data.
    #this will find the final lambda
    if (refit_model_all_data) {
        logger::log_info("Refitting final model using all training data.")
        X_all_s <- scale_train(X_all)

        cv_final <- glmnet::cv.glmnet(X_all_s$X, age, family="gaussian",
                                      alpha=alpha_fixed, standardize=FALSE,
                                      type.measure="mae")

        fit_final <- glmnet::glmnet(X_all_s$X, age, family="gaussian",
                                    alpha=alpha_fixed, lambda=cv_final$lambda.min,
                                    standardize=FALSE)

        feats_final <- rownames(fit_final$beta)
        final_model <- data.frame(
            feature = feats_final,
            coef = as.numeric(fit_final$beta),
            mean_train = X_all_s$mu[feats_final],
            sd_train = X_all_s$s[feats_final],
            intercept = as.numeric(fit_final$a0),
            row.names = NULL, check.names = FALSE
        )
        final_model$cor<-compute_age_expression_correlation(X_all[, final_model$feature], age)
        logger::log_info("Final model lambda: {cv_final$lambda.min}")
    } else {
        #take the average of the model parameters across the folds.
        final_model <- aggregate(cbind(coef, mean_train, sd_train, intercept) ~ feature,
                                 data = model_df, FUN = mean)
    }

    result=list(
        per_fold_metrics = per_fold,
        overall_metrics = overall,          # OOF metrics
        donor_age_predictions = oof_df,               # donor-level predictions
        fold_models = fold_models,  # list of per-fold model data frames
        model_df = model_df,         # combined long data frame
        final_model = final_model
    )
}

#######################
# GRID SEARCH ALPHA
# In practice, this performs no better than fixing alpha=0.5
# This is the same finding as the paper: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0237006
#######################

scale_train <- function(M) {
    mu <- colMeans(M)
    s <- apply(M, 2, sd)
    s[s == 0] <- 1
    list(
        X = sweep(sweep(M, 2, mu, "-"), 2, s, "/"),
        mu = mu,
        s = s
    )
}

scale_apply <- function(M, mu, s) {
    sweep(sweep(M, 2, mu, "-"), 2, s, "/")
}

# Search for best (alpha, lambda)
select_alpha_lambda <- function(X, y, alpha_grid, foldid = NULL, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    best <- list(cvmin = Inf)
    for (a in alpha_grid) {
        cv <- glmnet::cv.glmnet(
            X, y, family = "gaussian",
            alpha = a, foldid = foldid,
            standardize = FALSE, type.measure = "mae"
        )
        cvmin <- min(cv$cvm)
        if (cvmin < best$cvmin) {
            best <- list(
                alpha = a,
                lambda_min = cv$lambda.min,
                lambda_1se = cv$lambda.1se,
                cvmin = cvmin
            )
        }
    }
    best
}

train_enet_cv_optimize_alpha <- function(dge_train, k_fold_index,
                          age_col = "age", donor_col = "donor",
                          alpha_grid = seq(0,1, by = 0.1),
                          refit_model_all_data = TRUE,
                          seed = 1) {

    stopifnot("DGEList" %in% class(dge_train))
    set.seed(seed)

    smp <- dge_train$samples
    donor <- as.character(smp[[donor_col]])
    age <- as.numeric(smp[[age_col]])

    stopifnot(!is.null(names(k_fold_index)))
    fold_id <- as.integer(k_fold_index[donor])
    if (any(is.na(fold_id))) stop("Missing fold assignment for some donors")
    K <- max(fold_id, na.rm = TRUE)

    dge_train <- edgeR::calcNormFactors(dge_train)
    X_all <- t(edgeR::cpm(dge_train, log = TRUE, prior.count = 1))  # samples x genes

    oof_pred <- rep(NA_real_, length(age))
    per_fold <- data.frame(
        fold = integer(0),
        n_train = integer(0),
        n_test = integer(0),
        alpha = double(0),
        lambda = double(0),
        r = double(0),
        median_abs_error = double(0),
        mean_abs_error = double(0)
    )
    fold_predictions <- vector("list", K)
    fold_models <- vector("list", K)

    for (k in seq_len(K)) {
        te <- which(fold_id == k)
        tr <- which(fold_id != k)

        Xtr_s <- scale_train(X_all[tr, , drop = FALSE])
        Xte <- scale_apply(X_all[te, , drop = FALSE], Xtr_s$mu, Xtr_s$s)
        ytr <- age[tr]; yte <- age[te]

        set.seed(seed + k)
        inner_foldid <- sample(rep(1:5, length.out = length(tr)))

        sel <- select_alpha_lambda(Xtr_s$X, ytr, alpha_grid, foldid = inner_foldid)
        fit <- glmnet::glmnet(
            Xtr_s$X, ytr, family = "gaussian",
            alpha = sel$alpha, lambda = sel$lambda_min,
            standardize = FALSE
        )

        pred <- as.numeric(predict(fit, newx = Xte))
        oof_pred[te] <- pred

        m <- compute_age_metrics(pred, yte)
        per_fold <- rbind(per_fold, data.frame(
            fold = k, n_train = length(tr), n_test = length(te),
            alpha = sel$alpha, lambda = sel$lambda_min,
            r = m[["r"]],
            median_abs_error = m[["median_abs_error"]],
            mean_abs_error = m[["mean_abs_error"]]
        ))

        fold_predictions[[k]] <- data.frame(
            donor = donor[te], age = yte, pred = pred, fold = k,
            row.names = rownames(smp)[te], check.names = FALSE
        )

        feats_k <- rownames(fit$beta)
        coef_df_k <- data.frame(
            feature = feats_k,
            coef = as.numeric(fit$beta),
            mean_train = Xtr_s$mu[feats_k],
            sd_train = Xtr_s$s[feats_k],
            intercept = as.numeric(fit$a0),
            fold = k,
            row.names = NULL, check.names = FALSE
        )
        fold_models[[k]] <- coef_df_k
    }

    overall <- compute_age_metrics(oof_pred, age)
    oof_df <- do.call(rbind, fold_predictions)
    model_df <- do.call(rbind, fold_models)

    if (isTRUE(refit_model_all_data)) {
        logger::log_info("Refitting final model using all training data.")
        X_all_s <- scale_train(X_all)
        sel_full <- select_alpha_lambda(X_all_s$X, age, alpha_grid, seed = seed)
        fit_final <- glmnet::glmnet(
            X_all_s$X, age, family = "gaussian",
            alpha = sel_full$alpha, lambda = sel_full$lambda_min,
            standardize = FALSE
        )

        feats_final <- rownames(fit_final$beta)
        final_model <- data.frame(
            feature = feats_final,
            coef = as.numeric(fit_final$beta),
            mean_train = X_all_s$mu[feats_final],
            sd_train = X_all_s$s[feats_final],
            intercept = as.numeric(fit_final$a0),
            row.names = NULL, check.names = FALSE
        )
        final_alpha <- sel_full$alpha
        final_lambda <- sel_full$lambda_min
        final_model$cor<-compute_age_expression_correlation(X_all[, final_model$feature], age)

        logger::log_info("Final model alpha: {final_alpha}, lambda: {final_lambda}")
    } else {
        final_model <- aggregate(cbind(coef, mean_train, sd_train, intercept) ~ feature,
                                 data = model_df, FUN = mean)
        final_alpha <- NA_real_
        final_lambda <- NA_real_
    }

    list(
        per_fold_metrics = per_fold,
        overall_metrics = overall,
        donor_age_predictions = oof_df,
        fold_models = fold_models,
        model_df = model_df,
        final_model = final_model,
        final_alpha = final_alpha,
        final_lambda = final_lambda
    )
}

compute_age_expression_correlation<-function (X_all, age) {
    cor_vec=apply (X_all, 2, function (x) {
        suppressWarnings(stats::cor(x, age))
    })
    return (cor_vec)

}


#' Compute correlation and error metrics for predicted vs. actual ages
#'
#' @param pred Numeric vector of predicted ages.
#' @param actual Numeric vector of actual (chronological) ages.
#'
#' @return A one-row data.frame with:
#' \itemize{
#'   \item r — Pearson correlation.
#'   \item median_abs_error — median absolute difference.
#'   \item mean_abs_error — mean absolute error.
#' }
#' @export
compute_age_metrics <- function(pred, actual) {
    stopifnot(length(pred) == length(actual))
    r     <- suppressWarnings(stats::cor(pred, actual))
    median_abs_error <- stats::median(abs(pred - actual))
    mean_abs_error   <- mean(abs(pred - actual))
    data.frame(r = r, median_abs_error = median_abs_error, mean_abs_error = mean_abs_error)
}


# Predict ages from a DGEList using the finalized model.
# Starts from RAW counts: calcNormFactors -> cpm(log=TRUE, prior.count=1) -> z-score using training means/SDs -> linear predictor.
# dge_new: DGEList with raw counts
# model_final: data.frame from finalize_enet_model()
# override_model_params (optional) Override the gene mean, sd, and model intercept with the expression values from this DGEList.
# Returns: dataframe with donor, age, pred columns.
predict_age_from_dge <- function(dge_new, model_final, prior.count = 1, override_model_params_dge=NULL, update_intercept=T) {
    stopifnot("DGEList" %in% class(dge_new))

    if (!is.null(override_model_params_dge)) {
        logger::log_info("Overriding model parameters with provided DGEList expression values.")
        #compute mean and sd from the override DGEList
        #dge_override <- edgeR::calcNormFactors(override_model_params_dge)
        dge_override <- override_model_params_dge
        lcpm_override <- edgeR::cpm(dge_override, log = TRUE, prior.count = prior.count)  # genes x samples

        common <- intersect(model_final$feature, rownames(lcpm_override))
        if (length(common) > 0) {
            idx_common <- match(common, model_final$feature)
            mu_override <- rowMeans(lcpm_override[common, , drop = FALSE])
            sd_override <- apply(lcpm_override[common, , drop = FALSE], 1, sd)
            model_final$mean_train[idx_common] <- mu_override
            model_final$sd_train[idx_common] <- sd_override
        }
        #the intercept has to be from the population to be predicted, not from the subset used to train the model.
        if (update_intercept)
            model_final$intercept=mean((dge_new$samples$age))
    }

    # 1) normalize + log-CPM
    if (is.null(dge_new$samples$norm.factors)) dge_new <- edgeR::calcNormFactors(dge_new)
    lcpm <- edgeR::cpm(dge_new, log = TRUE, prior.count = prior.count)  # genes x samples

    # 2) align to model features and build samples × features matrix
    feats <- model_final$feature
    samp_ids <- colnames(lcpm)
    X <- matrix(0.0, nrow = length(samp_ids), ncol = length(feats),
                dimnames = list(samp_ids, feats))

    common <- intersect(feats, rownames(lcpm))
    if (length(common) > 0)
        X[, common] <- t(lcpm[common, , drop = FALSE])

    # 3) impute missing genes at training mean (so z = 0 for those features)
    if (length(common) < length(feats)) {
        miss <- setdiff(feats, common)
        logger::log_warn(paste("Genes in model not found in new data:", length(miss), "of", nrow(model_final)))
        if (length(miss) > 0) {
            idx_miss <- match(miss, model_final$feature)
            mu_miss  <- model_final$mean_train[idx_miss]
            for (g in seq_along(miss)) X[, miss[g]] <- mu_miss[g]
        }
    }

    # Precompute index to pull params in 'feats' order
    idx <- match(feats, model_final$feature)

    # 4) z-score using training means/SDs
    mu <- model_final$mean_train[idx]
    sd <- model_final$sd_train[idx]
    sd[!is.finite(sd) | sd == 0] <- 1
    Xz <- sweep(sweep(X, 2, mu, "-"), 2, sd, "/")

    # 5) linear predictor
    beta <- model_final$coef[idx]
    intercept <- unique(model_final$intercept)[1]
    pred <- as.numeric(intercept + Xz %*% beta)
    stats::setNames(pred, samp_ids)

    data.frame(donor = dge_new$samples$donor,
               age   = dge_new$samples$age,
               pred  = pred,
               row.names = rownames(dge_new$samples),
               check.names = FALSE)
}



#' Compute donor age residuals using repeated stratified train/test splits
#'
#' Performs repeated Monte Carlo cross-validation by repeatedly splitting donors
#' into stratified 80/20 train/test sets, fitting an age prediction model on the
#' training donors, and predicting ages for the held-out donors. For each donor,
#' predictions and residuals are accumulated across all repeats in which that
#' donor was held out, allowing estimation of average residuals and their
#' stability with respect to training-set composition.
#'
#' This approach reduces dependence on any single train/test split and can
#' mitigate repeatable biases introduced by fixed fold assignments, at the cost
#' of increased computation. It is intended for within-cell-type residual
#' analysis rather than external model deployment.
#'
#' @param dge_cell A donor-collapsed \code{DGEList} for a single cell type, with
#'   one sample per donor and donor age stored in \code{dge_cell$samples}.
#' @param age_col Column name in \code{dge_cell$samples} containing donor age.
#' @param donor_col Column name in \code{dge_cell$samples} containing donor IDs.
#' @param n_bins Number of age quantile bins used for stratified splitting.
#' @param test_prop Fraction of donors assigned to the test set in each repeat
#'   (e.g., 0.20 for an 80/20 split).
#' @param n_repeats Number of random stratified splits to perform.
#' @param optimize_alpha Logical; if \code{TRUE}, optimize elastic net alpha
#'   within each repeat using inner cross-validation.
#' @param alpha_fixed Elastic net alpha value to use when
#'   \code{optimize_alpha = FALSE}.
#' @param seed Random seed used to initialize the sequence of random splits.
#'   Each repeat uses \code{seed + repeat - 1}.
#' @param refit_model_all_data Logical; passed to \code{train_enet_cv} to control
#'   whether the final model in each repeat is refit on all training donors before
#'   prediction.
#' @param verbose_every Integer; if non-NULL and positive, log progress every
#'   \code{verbose_every} repeats.
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item \code{donor_oof_predictions}: Long-format data.frame with one row per
#'     donor per repeat in which the donor was held out, containing donor ID,
#'     true age, predicted age, residual (\code{pred - age}), and repeat index.
#'   \item \code{donor_residual_summary}: Data.frame with one row per donor,
#'     summarizing residual behavior across repeats, including mean prediction,
#'     mean and median residuals, residual standard deviation, and the number of
#'     times the donor was held out (\code{n_oof}).
#'   \item \code{per_repeat_metrics}: Data.frame of performance metrics (e.g.,
#'     correlation, median absolute error, mean absolute error) computed
#'     separately for each repeat.
#'   \item \code{avg_pred_metrics}: Performance metrics computed using the
#'     per-donor average predicted age across repeats.
#'   \item \code{params}: List of parameters used to generate the repeated splits
#'     (test fraction, number of bins, number of repeats, alpha settings, seed).
#' }
#'
#' @details
#' Only out-of-fold predictions are used for residual computation; donors are
#' never evaluated on models trained using their own data. The number of
#' out-of-fold predictions per donor may vary slightly depending on random
#' sampling. The resulting residual summaries can be used to assess both donor
#' age deviation and the stability of that deviation across different training
#' backgrounds.
#'
#' @export
compute_final_residuals_repeated_splits <- function(dge_cell,
                                                    age_col = "age",
                                                    donor_col = "donor",
                                                    n_bins = 5,
                                                    test_prop = 0.20,
                                                    n_repeats = 100,
                                                    optimize_alpha = FALSE,
                                                    alpha_fixed = 0.5,
                                                    seed = 1,
                                                    refit_model_all_data = TRUE,
                                                    verbose_every = 10) {

    stopifnot("DGEList" %in% class(dge_cell))
    stopifnot(test_prop > 0, test_prop < 1)
    stopifnot(n_repeats >= 1)

    oof_list <- vector("list", n_repeats)
    metrics_list <- vector("list", n_repeats)

    fit_predict_once <- function(dge_train, dge_test, seed_fit) {

        k_fold_index <- make_cv_folds_stratified(dge_train,
                                                 age_col = age_col,
                                                 donor_col = donor_col,
                                                 K = 5,
                                                 n_bins = n_bins,
                                                 seed = seed_fit)

        if (optimize_alpha) {
            r <- train_enet_cv_optimize_alpha(dge_train, k_fold_index,
                                              age_col = age_col,
                                              donor_col = donor_col,
                                              seed = seed_fit)
        } else {
            r <- train_enet_cv(dge_train, k_fold_index,
                               age_col = age_col,
                               donor_col = donor_col,
                               alpha_fixed = alpha_fixed,
                               refit_model_all_data = refit_model_all_data,
                               seed = seed_fit)
        }

        pred_df <- predict_age_from_dge(dge_test, r$final_model, prior.count = 1)
        list(pred_df = pred_df, model = r)
    }

    set.seed(seed)

    for (rep_i in seq_len(n_repeats)) {

        seed_i <- seed + rep_i - 1L

        if (!is.null(verbose_every) && verbose_every > 0 && (rep_i %% verbose_every == 0)) {
            logger::log_info(paste("Repeated split", rep_i, "of", n_repeats))
        }

        h <- make_holdout_stratified(dge_cell,
                                     age_col = age_col,
                                     donor_col = donor_col,
                                     test_prop = test_prop,
                                     n_bins = n_bins,
                                     seed = seed_i)

        dge_train <- h$dge_train
        dge_test  <- h$dge_test

        fp <- fit_predict_once(dge_train, dge_test, seed_fit = seed_i)

        pred_df <- fp$pred_df
        pred_df[["mc_iter"]] <- rep_i
        pred_df[["resid"]] <- pred_df[["pred"]] - pred_df[["age"]]

        oof_list[[rep_i]] <- pred_df[, c("donor", "age", "pred", "resid", "mc_iter")]

        m <- compute_age_metrics(pred_df$pred, pred_df$age)
        metrics_list[[rep_i]] <- cbind(mc_iter = rep_i, m)
    }

    oof_all <- do.call(rbind, oof_list)
    per_repeat_metrics <- do.call(rbind, metrics_list)

    # Aggregate donor residual stability across repeats where donor was OOF
    age_unique <- unique(dge_cell$samples[, c(donor_col, age_col), drop = FALSE])
    colnames(age_unique) <- c("donor", "age")
    age_unique$donor <- as.character(age_unique$donor)
    age_unique$age <- as.numeric(age_unique$age)

    split_counts <- aggregate(mc_iter ~ donor, data = oof_all, FUN = length)
    colnames(split_counts)[2] <- "n_oof"

    resid_mean <- aggregate(resid ~ donor, data = oof_all, FUN = mean)
    resid_median <- aggregate(resid ~ donor, data = oof_all, FUN = stats::median)
    resid_sd <- aggregate(resid ~ donor, data = oof_all, FUN = stats::sd)

    colnames(resid_mean)[2] <- "resid_mean"
    colnames(resid_median)[2] <- "resid_median"
    colnames(resid_sd)[2] <- "resid_sd"

    pred_mean <- aggregate(pred ~ donor, data = oof_all, FUN = mean)
    colnames(pred_mean)[2] <- "pred_mean"

    donor_summary <- merge(age_unique, split_counts, by = "donor", all.x = TRUE)
    donor_summary <- merge(donor_summary, pred_mean, by = "donor", all.x = TRUE)
    donor_summary <- merge(donor_summary, resid_mean, by = "donor", all.x = TRUE)
    donor_summary <- merge(donor_summary, resid_median, by = "donor", all.x = TRUE)
    donor_summary <- merge(donor_summary, resid_sd, by = "donor", all.x = TRUE)

    avg_pred_metrics <- compute_age_metrics(donor_summary$pred_mean, donor_summary$age)

    list(
        donor_oof_predictions = oof_all,             # long form: one row per (donor, mc_iter) when held out
        donor_residual_summary = donor_summary,      # one row per donor: mean/median/sd residual + n_oof
        per_repeat_metrics = per_repeat_metrics,     # metrics per split
        avg_pred_metrics = avg_pred_metrics,         # metrics using pred_mean across repeats
        params = list(
            test_prop = test_prop,
            n_bins = n_bins,
            n_repeats = n_repeats,
            optimize_alpha = optimize_alpha,
            alpha_fixed = alpha_fixed,
            seed = seed
        )
    )
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


plot_model_performance <- function(cv_metrics, per_fold_metrics, test_metrics, cellType = NULL) {
    stopifnot(is.data.frame(cv_metrics), is.data.frame(per_fold_metrics), is.data.frame(test_metrics))

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
            title = paste0(
                "Performance",
                if (!is.null(cellType)) paste0(" (", cellType, ")")
            ),
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

    df_plot <- dplyr::semi_join(df, feature_summary, by = "feature")
    df_plot$feature <- factor(df_plot$feature, levels = rev(feature_summary$feature))

    median_lines <- dplyr::mutate(
        dplyr::filter(feature_summary, feature %in% df_plot$feature),
        feature = factor(feature, levels = rev(levels(df_plot$feature)))
    )

    ggplot2::ggplot(df_plot, ggplot2::aes(x = coef, y = feature, color = fold)) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
        ggplot2::geom_point(
            size = 2.2,
            alpha = 0.8,
            position = ggplot2::position_jitter(height = 0.15, width = 0)
        ) +
        ggplot2::geom_segment(
            data = median_lines,
            ggplot2::aes(
                x = med,
                xend = med,
                y = as.numeric(feature) - 0.4,
                yend = as.numeric(feature) + 0.4
            ),
            color = "black",
            linewidth = 0.7
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

get_output_basename <- function(pdf_file, default = "age_prediction_results") {

    if (is.null(pdf_file) || is.na(pdf_file) || !nzchar(pdf_file)) {
        return(default)
    }

    sub("\\.pdf$", "", basename(pdf_file))
}

extract_age_outputs <- function(r, cellType, region = NA_character_) {

    preds <- r$final_oof_predictions[, c("donor", "age", "pred")]
    preds <- cbind(cell_type = cellType, region = region, preds)

    coefs <- r$final_cv_model$final_model
    coefs <- cbind(cell_type = cellType, region = region, coefs)

    # Summary metrics (3 rows per cell type/region)
    mets <- rbind(
        cbind(set = "Final OOF (all donors)", r$final_oof_metrics),
        cbind(set = "Inner CV OOF (80% train)", r$cv_model$overall_metrics),
        cbind(set = "Outer holdout (20% donors)", r$test_set_metrics)
    )
    mets <- cbind(cell_type = cellType, region = region, mets)

    # Per-fold metrics (inner CV on 80% train)
    fold_mets <- r$cv_model$per_fold_metrics
    fold_mets <- cbind(cell_type = cellType, region = region, fold_mets)

    list(
        donor_predictions = preds,
        model_coefficients = coefs,
        model_metrics = mets,
        per_fold_metrics = fold_mets
    )
}

write_age_outputs_all <- function(outputs, result_dir, output_basename) {

    donor_predictions <- do.call(rbind, lapply(outputs, function(x) x$donor_predictions))
    model_coefficients <- do.call(rbind, lapply(outputs, function(x) x$model_coefficients))
    model_metrics <- do.call(rbind, lapply(outputs, function(x) x$model_metrics))
    per_fold_metrics <- do.call(rbind, lapply(outputs, function(x) x$per_fold_metrics))

    out_pred <- file.path(result_dir,
                          paste0(output_basename, "_donor_predictions.txt"))
    out_coef <- file.path(result_dir,
                          paste0(output_basename, "_model_coefficients.txt"))
    out_met  <- file.path(result_dir,
                          paste0(output_basename, "_model_metrics.txt"))
    out_fmet <- file.path(result_dir,
                          paste0(output_basename, "_model_per_fold_metrics.txt"))

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

    invisible(list(
        donor_predictions = out_pred,
        model_coefficients = out_coef,
        model_metrics = out_met,
        per_fold_metrics = out_fmet
    ))
}

#' Write donor-level age predictions for a cell type
#'
#' Writes donor-level predicted ages (and true ages) to a tab-delimited file for a
#' single cell type. This function uses the final full-dataset out-of-fold (OOF)
#' predictions stored in \code{r$final_oof_predictions}, which provide one
#' prediction per donor from a model that did not train on that donor.
#'
#' @param r Result list returned by \code{predict_age_celltype()}.
#' @param result_dir Output directory for result files.
#' @param cellType Cell type name (used in output file name and in a \code{cell_type} column).
#'
#' @return Invisibly returns the output file path.
#' @export
write_donor_age_predictions <- function(r, result_dir, cellType) {

    colsToKeep <- c("donor", "age", "pred")

    # Use FINAL full-dataset OOF predictions only
    df <- r$final_oof_predictions[, colsToKeep]

    # Add cell type for tracking
    df <- cbind(cell_type = cellType, df)

    out_file <- file.path(result_dir,
                          paste0("donor_age_predictions_", cellType, ".txt"))

    write.table(df, file = out_file,
                sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = TRUE)
}

#' Write age model coefficients for a cell type
#'
#' Writes the coefficient table for the all-data model trained on the full cohort.
#' The coefficients written by this function are intended for prediction on
#' external datasets and are stored in \code{r$final_cv_model$final_model}.
#'
#' The output includes the feature name, coefficient, per-feature training mean and
#' standard deviation used for z-scoring, and the model intercept.
#'
#' @param r Result list returned by \code{predict_age_celltype()}.
#' @param result_dir Output directory for result files.
#' @param cellType Cell type name (used in output file name and in a \code{cell_type} column).
#'
#' @return Invisibly returns the output file path.
#' @export
write_model_coefficients <- function(r, result_dir, cellType) {

    # Use coefficients from the ALL-DATA model (for external prediction)
    df <- cbind(cell_type = cellType, r$final_cv_model$final_model)

    out_file <- file.path(result_dir,
                          paste0("donor_age_model_coefficients_", cellType, ".txt"))

    write.table(df, file = out_file,
                sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = TRUE)
}

#' Write age model performance metrics for a cell type
#'
#' Writes summary and per-fold performance metrics for a cell type.
#'
#' The summary metrics file includes:
#' \itemize{
#'   \item \strong{Final OOF (all donors)}: metrics computed from full-dataset
#'         out-of-fold predictions (\code{r$final_oof_metrics}).
#'   \item \strong{Inner CV OOF (80\% train)}: metrics from inner cross-validation
#'         on the 80\% training split (\code{r$cv_model$overall_metrics}).
#'   \item \strong{Outer holdout (20\% donors)}: metrics computed on the held-out
#'         donor split (\code{r$test_set_metrics}).
#' }
#'
#' The per-fold metrics file contains fold-wise inner CV metrics from
#' \code{r$cv_model$per_fold_metrics}.
#'
#' @param r Result list returned by \code{predict_age_celltype()}.
#' @param result_dir Output directory for result files.
#' @param cellType Cell type name (used in output file name and in a \code{cell_type} column).
#'
#' @return Invisibly returns a character vector of length 2 with the summary and per-fold file paths.
#' @export
write_model_metrics <- function(r, result_dir, cellType) {

    out_file <- file.path(result_dir,
                          paste0("donor_age_model_metrics_", cellType, ".txt"))

    df <- rbind(
        cbind(set = "Final OOF (all donors)", r$final_oof_metrics),
        cbind(set = "Inner CV OOF (80% train)", r$cv_model$overall_metrics),
        cbind(set = "Outer holdout (20% donors)", r$test_set_metrics)
    )

    df <- cbind(cell_type = cellType, df)

    write.table(df, file = out_file,
                sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = TRUE)

    # Also write per-fold metrics (inner CV on 80% train)
    out_file2 <- file.path(result_dir,
                           paste0("donor_age_model_per_fold_metrics_", cellType, ".txt"))

    df2 <- r$cv_model$per_fold_metrics
    df2 <- cbind(cell_type = cellType, df2)

    write.table(df2, file = out_file2,
                sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = TRUE)
}

#########################
# MAKE TEST / TRAIN DATA
#########################

# Quantile bin helper
.quantile_bins <- function(x, n_bins = 5) {
    if (length(unique(na.omit(x))) <= 1L) return(factor(rep("bin1", length(x))))
    qs <- unique(stats::quantile(x, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
    cut(x, breaks = qs, include.lowest = TRUE, right = TRUE)
}

# Hamilton apportionment helper
.apportion_round <- function(n_per_bin, prop) {
    target <- n_per_bin * prop
    base   <- floor(target)
    rem    <- round(sum(target)) - sum(base)
    if (rem > 0) {
        frac <- target - base
        add  <- order(frac, decreasing = TRUE)[seq_len(min(rem, length(frac)))]
        base[add] <- base[add] + 1L
    }
    base
}

# 1) Donor-grouped, age-stratified holdout split
make_holdout_stratified <- function(dge, age_col = "age", donor_col = "donor",
                                    test_prop = 0.20, n_bins = 5, seed = 1) {
    stopifnot("DGEList" %in% class(dge), test_prop > 0, test_prop < 1)
    smp   <- dge$samples
    age   <- as.numeric(smp[[age_col]])
    names (age) <- as.character(smp[[donor_col]])

    #bin donors by age quantiles
    d_bins <- .quantile_bins(age, n_bins = n_bins)
    set.seed(seed)

    donors <- names(age)
    by_bin <- split(donors, d_bins)
    n_bin  <- vapply(by_bin, length, 1L)
    n_test_bin <- .apportion_round(n_bin, prop = test_prop)

    test_donors <- character(0)
    i <- 0L
    # samples from donor age bins
    for (b in names(by_bin)) {
        i <- i + 1L
        pool <- by_bin[[b]]
        if (length(pool) == 0) next
        k <- n_test_bin[i]
        if (k > 0) test_donors <- c(test_donors, sample(pool, k))
    }

    test_donors <- unique(test_donors)

    # return indices and convenient subsets
    dge_train<-dge[ , !dge$samples$donor %in% test_donors, keep.lib.sizes = TRUE]
    dge_test<-dge[ , dge$samples$donor %in% test_donors, keep.lib.sizes = TRUE]

    list(
        dge_train    = dge_train,
        dge_test     = dge_test
    )
}

# 2) Donor-grouped, age-stratified K-folds on a (possibly pre-subset) DGEList
make_cv_folds_stratified <- function(dge, age_col = "age", donor_col = "donor",
                                     K = 5, n_bins = 5, seed = 1) {
    stopifnot("DGEList" %in% class(dge), K >= 2)
    smp   <- dge$samples
    donor <- as.character(smp[[donor_col]])          # define donor
    age   <- as.numeric(smp[[age_col]])
    names(age) <- donor                               # donors == samples

    d_bins <- .quantile_bins(age, n_bins = n_bins)
    set.seed(seed)

    fold_by_donor <- integer(length(age)); names(fold_by_donor) <- names(age)
    for (b in levels(d_bins)) {
        pool <- names(age)[d_bins == b]
        if (length(pool) == 0) next
        pool <- sample(pool)
        fseq <- rep(1:K, length.out = length(pool))
        fold_by_donor[pool] <- fseq
    }

    fold_id <- as.integer(fold_by_donor[donor])       # map to sample order
    names(fold_id) <- rownames(smp)
    fold_id
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
collapse_by_donor <- function(dge, donor_col = "donor", keep_cols = character(0)) {
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

    # build samples data.frame with only keep_cols and validate uniqueness per donor
    if (length(keep_cols)) {
        missing <- setdiff(keep_cols, colnames(smp))
        if (length(missing)) stop("keep_cols not found in dge$samples: ", paste(missing, collapse = ", "))

        nD <- nlevels(f)
        out_list <- vector("list", length(keep_cols))
        names(out_list) <- keep_cols

        for (j in seq_along(keep_cols)) {
            colname <- keep_cols[j]
            colj <- smp[[colname]]
            # treat factors as character to avoid level coercion
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

    # assemble output
    dge_out <- dge
    dge_out$counts  <- C
    dge_out$samples <- samples_out
    dge_out
}

get_age_de_results<-function (cellType, age_de_results_dir, fdr_threshold=0.05) {
    if (is.null(age_de_results_dir))
        return (NULL)

    f=list.files(path=age_de_results_dir, pattern="age", full.names = TRUE)
    expectedFileName=paste(cellType, "age_DE_results.txt", sep="_")
    f=grep(expectedFileName, f, value=TRUE)
    if (length(f)!=1) {
        stop("Did not find exactly one DE results file for cell type ", cellType, " in directory ", age_de_results_dir)
    }
    de_results=utils::read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    de_results=de_results[de_results$adj.P.Val<=fdr_threshold,]
    return (de_results)
}

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

