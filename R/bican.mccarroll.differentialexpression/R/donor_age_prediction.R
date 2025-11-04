# library(edgeR)
# library(glmnet)
# library(ggplot2)
# library(logger)
# library (cowplot)
# library (dplyr)

# data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metacells"
# data_name="donor_rxn_DGEList"
#
# age_de_results_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age/cell_type"
# result_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/age_prediction"
# outPDFFile="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/age_prediction/age_prediction_results.pdf"
#
# #filtering to autosomal genes
# contig_yaml_file="/broad/mccarroll/software/metadata/individual_reference/GRCh38_ensembl_v43/GRCh38_ensembl_v43.contig_groups.yaml"
# reduced_gtf_file="/broad/mccarroll/software/metadata/individual_reference/GRCh38_ensembl_v43/GRCh38_ensembl_v43.reduced.gtf"

#filtering to gene functional types.  #TODO collapse this and autosomal processes.
# gtf_path="/broad/mccarroll/software/metadata/individual_reference/GRCh38_ensembl_v43/GRCh38_ensembl_v43.gtf"

# donor_col = "donor"
# age_col = "age"
# seed =12345; fdr_threshold=0.05; optimize_alpha=TRUE; alpha_fixed=0.5


#run Emi's data:
# data_dir="/broad/mccarroll/dropulation/analysis/cellarium_upload/SNAP200_freeze1/metacells"
# data_name="donor_rxn_DGEList"
#
# #still use the same set of overall features from BICAN.
# #I don't want to have to revisit this unless neccesary.
# age_de_results_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/sex_age/cell_type"
# result_dir="/broad/mccarroll/dropulation/analysis/age_predictions/SNAP200"
# outPDFFile="/broad/mccarroll/dropulation/analysis/age_predictions/SNAP200/age_prediction_results.pdf"


predict_age_by_celltype<-function (data_dir, data_name) {
    #validate the output directory exists
    if (!dir.exists(result_dir)) {
        stop("Result directory does not exist: ", result_dir)
    }

    #load the DGEList and prepare the data
    d=bican.mccarroll.differentialexpression::prepare_data_for_differential_expression(data_dir, data_name, randVars=c(), fixedVars=c())
    dge=d$dge

    #filter the genes to autosomal only
    dge=filter_dge_to_autosomes (dge, contig_yaml_file, reduced_gtf_file)

    #experimental: filter out lncRNAs and other non-gene symbols - this generally makes things a little worse
    # I guess those lncRNAs are doing something!
    #dge=filter_dge_to_gene_functions (dge, gene_functions=c("protein_coding"), gtf_path)

    cell_type_list=unique(dge$samples$cell_type)
    #cell_type_list="microglia" #hard coded for now.
    lineStr <- strrep("=", 80)

    if (!is.null(outPDFFile)) {
        pdf(outPDFFile)
    }

    #infer the columns you'll need to retain
    retained_features=c(donor_col, age_col)

    test_set_metrics_list=list()

    for (cellType in cell_type_list) {
        logger::log_info(lineStr)
        logger::log_info(paste("Learning donor age model from expression for cell type:", cellType))
        logger::log_info(lineStr)

        r=predict_age_celltype(cellType, dge, retained_features=retained_features, donor_col = donor_col, age_de_results_dir=age_de_results_dir, fdr_threshold=fdr_threshold, optimize_alpha=optimize_alpha, alpha_fixed=alpha_fixed, seed=seed)
        test_set_metrics_list[[cellType]]=r$test_set_metrics

        # QC plots
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

        print (p)
        print (coef1)
        if (!is.null(r$age_de_results))
            print (model_vs_de_plot)

        #write out the donor predictions, model coefficients, and metrics
        write_donor_age_predictions(r, result_dir, cellType)
        write_model_coefficients(r, result_dir, cellType)
        write_model_metrics(r, result_dir, cellType)
    }

    test_set_metrics_df <- do.call(rbind, test_set_metrics_list)
    mean_metrics <- apply (test_set_metrics_df, 2, mean)
    logger::log_info(
        "Mean metrics of held out test sets across all cell types: {paste(names(mean_metrics), sprintf('[%.3f]', mean_metrics), collapse = ', ')}"
    )

    if (!is.null(outPDFFile)) {
        dev.off()
    }

    # for each cell type

}

compare_age_residuals_to_donor_coefficients<-function () {
    #Load the DGEList, get the sample metadata
    #simplify to per-donor covariates
    #compare covariates to each cell type donor age residuals (prediction - chronological age)

    d=bican.mccarroll.differentialexpression::prepare_data_for_differential_expression(data_dir, data_name, randVars=c(), fixedVars=c())
    dge=d$dge

}

# The main function to predict age from expression for a given cell type
# Assumes all data has been materialized.
predict_age_celltype<-function (cellType, dge, retained_features=c("donor", "age"), donor_col = "donor", age_de_results_dir, fdr_threshold=0.05, optimize_alpha=FALSE, alpha_fixed=0.5, seed=1) {

    if (optimize_alpha) {
        logger::log_info("Optimizing alpha via grid search.")
    } else {
        logger::log_info(paste("Using fixed alpha =", alpha_fixed))
    }


    dge_cell <- dge[, dge$samples$cell_type == cellType, keep.lib.sizes = TRUE]
    #collapse across regions and muliple samples per donor
    dge_cell<-collapse_by_donor(dge_cell, donor_col = donor_col, keep_cols = retained_features)

    #filtering samples by library size
    r<- bican.mccarroll.differentialexpression::filter_by_libsize(dge_cell, threshold_sd = 1.96, bins = 50, strTitlePrefix = cellType)
    dge_cell<- r$dge

    #filter to the top 75% of highly expressed genes as a first pass.
    dge_cell<-bican.mccarroll.differentialexpression::filter_top_expressed_genes(dge_cell, gene_filter_frac = 0.75, verbose = TRUE)

    #filter to cpm cutoff of 1.
    r2=bican.mccarroll.differentialexpression::plot_logCPM_density_quantiles(dge_cell, cpm_cutoff = 1, logCPM_xlim = c(-5, 15), lower_quantile = 0.05, upper_quantile = 0.95, quantile_steps = 5)
    dge_cell=r2$filtered_dge

    #get age DE results (optional)
    age_de_results=get_age_de_results(cellType, age_de_results_dir, fdr_threshold=fdr_threshold)
    if (!is.null(age_de_results)) {
        genes_to_keep=intersect(rownames(dge_cell), rownames(age_de_results))
        logger::log_info(paste("Filtering to", length(genes_to_keep), "genes significant in age DE results at FDR <=", fdr_threshold))
        dge_cell=dge_cell[genes_to_keep,]
    }

    # start training the model by splitting away final train/test data.
    n_bins=5
    d=make_holdout_stratified(dge_cell, age_col = "age", donor_col = "donor", test_prop = 0.20, n_bins = n_bins, seed = seed)
    dge_train=d$dge_train
    dge_test=d$dge_test
    #validate train/test split are completely disjoint on donors
    if (length(intersect(dge_train$samples$donor, dge_test$samples$donor))>0) {
        stop("Train and test donors are not disjoint!")
    }

    logger::log_info(paste("Training data has", ncol(dge_train), "donors across", nrow(dge_train$counts), "genes"))
    logger::log_info(paste("Final holdout test data has", ncol(dge_test), "donors across", nrow(dge_test$counts), "genes"))

    k_fold_index=make_cv_folds_stratified(dge_train, age_col = "age", donor_col = "donor", K = 5, n_bins = n_bins, seed = seed)

    # Optionally optimize alpha via grid search
    if (optimize_alpha)
        r=train_enet_cv_optimize_alpha(dge_train, k_fold_index, age_col = "age", donor_col = "donor", seed = seed)
    else
        r=train_enet_cv(dge_train, k_fold_index, age_col = "age", donor_col = "donor", alpha_fixed = alpha_fixed, seed = seed)

    # the final model only needs the non-zero coefficients, but keep all for now.
    #r$final_model=r$final_model[which(r$final_model$coef != 0), ]
    #logger::log_info(paste("Final model has", nrow(r$final_model), "non-zero coefficients."))

    num_nonzero_coefs=sum(r$final_model$coef != 0)
    logger::log_info(paste("Final model has", num_nonzero_coefs, "non-zero coefficients."))

    #running predictions on 20% held out data.
    test_set_predictions=predict_age_from_dge(dge_test, r$final_model, prior.count = 1)

    metrics_test <- compute_age_metrics(test_set_predictions$pred, test_set_predictions$age)

    result=list(
        cell_type = cellType,
        dge_train = dge_train,
        dge_test = dge_test,
        cv_model = r,
        test_set_predictions = test_set_predictions,
        test_set_metrics = metrics_test,
        age_de_results=age_de_results
    )

    return (result)

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
# Returns: dataframe with donor, age, pred columns.
predict_age_from_dge <- function(dge_new, model_final, prior.count = 1) {
    stopifnot("DGEList" %in% class(dge_new))

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

    ggplot(df, aes(x = age, y = pred)) +
        geom_point(size = 2, alpha = 0.7, color = "steelblue") +
        geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
        geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.6, color = "red") +
        labs(
            title = titleStr,
            x = "Chronological Age",
            y = "Predicted Age"
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
    suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

    df <- dplyr::bind_rows(fold_models)
    df$fold <- factor(df$fold, levels = sort(unique(df$fold)))

    # Remove zero coefficients
    df <- df %>% filter(coef != 0)

    # Summarize feature-level medians
    feature_summary <- df %>%
        group_by(feature) %>%
        summarize(any_pos = any(coef > 0),
                  nnz = sum(coef != 0),
                  med = median(coef),
                  .groups = "drop") %>%
        filter((!positive_only | any_pos), nnz >= min_nonzero) %>%
        arrange(desc(med))

    if (!is.na(top_n) && is.finite(top_n)) feature_summary <- head(feature_summary, top_n)
    df_plot <- df %>% semi_join(feature_summary, by = "feature")
    df_plot$feature <- factor(df_plot$feature, levels = rev(feature_summary$feature))

    # Add median values to plotting dataframe
    median_lines <- feature_summary %>%
        filter(feature %in% df_plot$feature) %>%
        mutate(feature = factor(feature, levels = rev(feature)))

    ggplot(df_plot, aes(x = coef, y = feature, color = fold)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
        geom_point(size = 2.2, alpha = 0.8,
                   position = position_jitter(height = 0.15, width = 0)) +
        geom_segment(data = median_lines,
                     aes(x = med, xend = med,
                         y = as.numeric(feature) - 0.4,
                         yend = as.numeric(feature) + 0.4),
                     color = "black", linewidth = 0.7) +
        scale_color_brewer(palette = "Dark2") +
        labs(
            title = paste0("Elastic Net coefficients across folds",
                           if (!is.null(cellType)) paste0(" (", cellType, ")")),
            x = "Coefficient", y = NULL, color = "Fold"
        ) +
        theme_classic(base_size = 12) +
        theme(
            axis.text.y = element_text(size = 9),
            axis.ticks.y = element_blank(),
            axis.ticks.length.y = unit(0, "pt"),
            legend.position = "top",
            panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3)
        )

}

# Compare elastic net model coefficients with DE logFC (nonzero only)
#final_model=r$cv_model$final_model; age_de_results=r$age_de_results; cellType=cellType
plot_model_vs_de <- function(final_model, age_de_results, cellType = NULL) {
    suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tibble) })

    df <- final_model %>%
        filter(coef != 0) %>%
        select(feature, model_coef = coef) %>%
        inner_join(
            age_de_results %>%
                rownames_to_column("feature") %>%
                select(feature, logFC),
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

write_donor_age_predictions<-function (r, result_dir, cellType) {
    colsToKeep=c("donor", "age", "pred")
    test_set_predictions=r$test_set_predictions[,colsToKeep]
    train_set_predictions=r$cv_model$donor_age_predictions[,colsToKeep]
    #write out the predictions
    df=rbind(test_set_predictions, train_set_predictions)
    #add in the cell type for easier tracking later.
    df=cbind(cell_type=cellType, df)
    out_file=file.path(result_dir, paste0("donor_age_predictions_", cellType, ".txt"))
    write.table(df, file=out_file, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

write_model_coefficients<-function (r, result_dir, cellType) {
    df=cbind(cell_type=cellType, r$cv_model$final_model)
    out_file=file.path(result_dir, paste0("donor_age_model_coefficients_", cellType, ".txt"))
    write.table(df, file=out_file, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

write_model_metrics<-function (r, result_dir, cellType) {
    out_file=file.path(result_dir, paste0("donor_age_model_metrics_", cellType, ".txt"))
    df=rbind(
        cbind(set="Cross-validation", r$cv_model$overall_metrics),
        cbind(set="Held-out donors", r$test_set_metrics)
    )
    df=cbind(cell_type=cellType, df)
    write.table(df, file=out_file, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

    #also write per-fold metrics
    out_file2=file.path(result_dir, paste0("donor_age_model_per_fold_metrics_", cellType, ".txt"))
    df=r$cv_model$per_fold_metrics
    df=cbind(cell_type=cellType, df)
    write.table(r$cv_model$per_fold_metrics, file=out_file2, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

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
    if (is.null(gene_functions) | is.null(reduced_gtf_file)) {
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

