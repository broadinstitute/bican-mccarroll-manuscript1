# -------------------------------------------------------------------------
# Shared internal helpers for training elastic net age models
# -------------------------------------------------------------------------

#' Z-score (center/scale) a training design matrix
#'
#' Internal helper to z-score features using training-set mean and standard
#' deviation. Returns the scaled matrix along with the per-feature mean and
#' standard deviation needed to apply the same scaling to new data.
#'
#' @param M Numeric matrix of shape samples x features.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{X}: z-scored matrix (samples x features).
#'   \item \code{mu}: numeric vector of feature means.
#'   \item \code{s}: numeric vector of feature standard deviations (zeros set to 1).
#' }
#'
#' @keywords internal
.scale_train <- function(M) {
    mu <- colMeans(M)
    s <- apply(M, 2, stats::sd)
    s[s == 0] <- 1
    list(X = sweep(sweep(M, 2, mu, "-"), 2, s, "/"), mu = mu, s = s)
}

#' Apply training z-scoring parameters to a matrix
#'
#' Internal helper to center and scale a matrix using pre-computed training-set
#' parameters.
#'
#' @param M Numeric matrix of shape samples x features.
#' @param mu Numeric vector of feature means (training).
#' @param s Numeric vector of feature standard deviations (training).
#'
#' @return A numeric matrix of the same shape as \code{M}, z-scored using \code{mu}
#'   and \code{s}.
#'
#' @keywords internal
.scale_apply <- function(M, mu, s) {
    sweep(sweep(M, 2, mu, "-"), 2, s, "/")
}

#' Select (alpha, lambda) for elastic net by grid search over alpha
#'
#' Internal helper that performs a grid search over \code{alpha_grid} and uses
#' \code{glmnet::cv.glmnet()} to select \code{lambda.min} for each alpha. Returns
#' the alpha whose cross-validated loss is minimal, along with its lambda.
#'
#' @param X Numeric matrix (samples x features), typically already scaled.
#' @param y Numeric response vector (length = nrow(X)).
#' @param alpha_grid Numeric vector of candidate alpha values in [0, 1].
#' @param type.measure Character; passed to \code{glmnet::cv.glmnet()}.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{alpha}: selected alpha.
#'   \item \code{lambda_min}: selected lambda (cv lambda.min for selected alpha).
#'   \item \code{cvmin}: best (minimum) achieved cross-validation loss.
#' }
#'
#' @keywords internal
.select_alpha_lambda <- function(X, y, alpha_grid, type.measure = "mae", seed = NULL) {
    if (!is.null(seed)) set.seed(seed)

    best_cvmin <- Inf
    best_alpha <- NA_real_
    best_lambda <- NA_real_

    for (a in alpha_grid) {
        cv <- glmnet::cv.glmnet(
            X, y, family = "gaussian",
            alpha = a,
            standardize = FALSE,
            type.measure = type.measure
        )
        cvmin <- min(cv$cvm, na.rm = TRUE)
        if (is.finite(cvmin) && cvmin < best_cvmin) {
            best_cvmin <- cvmin
            best_alpha <- a
            best_lambda <- cv$lambda.min
        }
    }

    list(alpha = best_alpha, lambda_min = best_lambda, cvmin = best_cvmin)
}

#' Fit final elastic net model with fixed alpha on all provided data
#'
#' Internal helper to fit a final model on all samples in \code{X_all} given a
#' fixed alpha. Uses \code{glmnet::cv.glmnet()} to select \code{lambda.min} and
#' returns a model data.frame in the same format expected by
#' \code{predict_age_from_dge()}.
#'
#' @param X_all Numeric matrix (samples x genes), unscaled logCPM values.
#' @param y Numeric response vector (age), length = nrow(X_all).
#' @param alpha_fixed Numeric alpha in [0, 1].
#'
#' @return A list with:
#' \itemize{
#'   \item \code{final_model}: data.frame with columns \code{feature}, \code{coef},
#'     \code{mean_train}, \code{sd_train}, \code{intercept}.
#'   \item \code{final_alpha}: the alpha used (same as \code{alpha_fixed}).
#'   \item \code{final_lambda}: selected \code{lambda.min}.
#' }
#'
#' @keywords internal
.fit_final_model_fixed_alpha <- function(X_all, y, alpha_fixed) {

    Xs <- .scale_train(X_all)

    cv <- glmnet::cv.glmnet(
        Xs$X, y, family = "gaussian",
        alpha = alpha_fixed,
        standardize = FALSE,
        type.measure = "mae"
    )

    fit <- glmnet::glmnet(
        Xs$X, y, family = "gaussian",
        alpha = alpha_fixed,
        lambda = cv$lambda.min,
        standardize = FALSE
    )

    feats <- rownames(fit$beta)
    final_model <- data.frame(
        feature = feats,
        coef = as.numeric(fit$beta),
        mean_train = Xs$mu[feats],
        sd_train = Xs$s[feats],
        intercept = as.numeric(fit$a0),
        row.names = NULL,
        check.names = FALSE
    )

    list(final_model = final_model, final_alpha = alpha_fixed, final_lambda = cv$lambda.min)
}

#' Fit final elastic net model with alpha optimized by grid search
#'
#' Internal helper to fit a final model on all samples in \code{X_all}. Performs
#' a grid search over \code{alpha_grid}; for each alpha, selects \code{lambda.min}
#' by \code{glmnet::cv.glmnet()}, then chooses the alpha with lowest CV loss and
#' refits a final \code{glmnet} model at that (alpha, lambda).
#'
#' @param X_all Numeric matrix (samples x genes), unscaled logCPM values.
#' @param y Numeric response vector (age), length = nrow(X_all).
#' @param alpha_grid Numeric vector of candidate alpha values.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{final_model}: data.frame with columns \code{feature}, \code{coef},
#'     \code{mean_train}, \code{sd_train}, \code{intercept}.
#'   \item \code{final_alpha}: selected alpha.
#'   \item \code{final_lambda}: selected \code{lambda.min}.
#' }
#'
#' @keywords internal
.fit_final_model_opt_alpha <- function(X_all, y, alpha_grid, seed = NULL) {

    Xs <- .scale_train(X_all)

    sel <- .select_alpha_lambda(Xs$X, y, alpha_grid = alpha_grid, seed = seed)

    fit <- glmnet::glmnet(
        Xs$X, y, family = "gaussian",
        alpha = sel$alpha,
        lambda = sel$lambda_min,
        standardize = FALSE
    )

    feats <- rownames(fit$beta)
    final_model <- data.frame(
        feature = feats,
        coef = as.numeric(fit$beta),
        mean_train = Xs$mu[feats],
        sd_train = Xs$s[feats],
        intercept = as.numeric(fit$a0),
        row.names = NULL,
        check.names = FALSE
    )

    list(final_model = final_model, final_alpha = sel$alpha, final_lambda = sel$lambda_min)
}

#' Prepare donor age training design from a DGEList
#'
#' Internal helper to extract donor and age metadata from \code{dge_train$samples},
#' compute normalization factors, and build the log-CPM expression design matrix
#' used for modeling.
#'
#' @param dge_train A \code{DGEList}.
#' @param age_col Column name in \code{dge_train$samples} containing donor age.
#' @param donor_col Column name in \code{dge_train$samples} containing donor IDs.
#' @param prior.count Prior count used in \code{edgeR::cpm(log=TRUE)}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{dge_train}: normalized \code{DGEList} (norm factors computed).
#'   \item \code{smp}: \code{dge_train$samples}.
#'   \item \code{donor}: donor ID vector (character).
#'   \item \code{age}: age vector (numeric).
#'   \item \code{X_all}: numeric matrix (samples x genes) of log-CPM values.
#' }
#'
#' @keywords internal
.prepare_design_from_dge <- function(dge_train, age_col, donor_col, prior.count = 1) {

    stopifnot("DGEList" %in% class(dge_train))

    smp <- dge_train$samples
    donor <- as.character(smp[[donor_col]])
    age <- as.numeric(smp[[age_col]])

    if (is.null(dge_train$samples$norm.factors)) {
        dge_train <- edgeR::calcNormFactors(dge_train)
    }

    X_all <- t(edgeR::cpm(dge_train, log = TRUE, prior.count = prior.count))  # samples x genes

    list(dge_train = dge_train, smp = smp, donor = donor, age = age, X_all = X_all)
}

#' Train elastic net age model with K-fold cross-validation (fixed alpha)
#'
#' Trains an elastic net age prediction model using donor-level K-fold
#' cross-validation for out-of-fold (OOF) performance estimation. The final model
#' is always refit on all provided samples using \code{glmnet::cv.glmnet()} to
#' select \code{lambda.min}.
#'
#' Set \code{compute_oof = FALSE} to skip all fold-based OOF predictions and
#' per-fold bookkeeping; in this mode the function fits only the final model on
#' all provided data, which is useful for repeated Monte Carlo train/test splits.
#'
#' @param dge_train A \code{DGEList} with donor-level samples and ages stored in
#'   \code{dge_train$samples}.
#' @param k_fold_index Named integer vector mapping donor IDs to fold numbers
#'   1..K. Required when \code{compute_oof = TRUE}. Ignored when
#'   \code{compute_oof = FALSE}.
#' @param age_col Column name in \code{dge_train$samples} containing donor age.
#' @param donor_col Column name in \code{dge_train$samples} containing donor IDs.
#' @param alpha_fixed Elastic net alpha in [0, 1] (\code{0} ridge, \code{1} lasso).
#' @param compute_oof Logical; if \code{TRUE}, compute OOF predictions/metrics
#'   using \code{k_fold_index}. If \code{FALSE}, skip OOF computation and return
#'   only the refit final model and associated metadata.
#' @param seed Integer seed for reproducibility.
#' @param verbose Logical; if \code{TRUE}, log progress messages.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{per_fold_metrics}: data.frame of metrics per fold (or \code{NULL}
#'     if \code{compute_oof = FALSE}).
#'   \item \code{overall_metrics}: metrics on OOF predictions across all donors
#'     (or \code{NULL} if \code{compute_oof = FALSE}).
#'   \item \code{donor_age_predictions}: OOF prediction data.frame with columns
#'     \code{donor}, \code{age}, \code{pred}, \code{fold} (or \code{NULL} if
#'     \code{compute_oof = FALSE}).
#'   \item \code{fold_models}: list of per-fold model coefficient data.frames
#'     (or \code{NULL} if \code{compute_oof = FALSE}).
#'   \item \code{model_df}: combined long-form coefficient data.frame across folds
#'     (or \code{NULL} if \code{compute_oof = FALSE}).
#'   \item \code{final_model}: data.frame representing the final model fit on all
#'     provided data, including scaling parameters and intercept.
#'   \item \code{final_alpha}: alpha used for the final model (equals
#'     \code{alpha_fixed}).
#'   \item \code{final_lambda}: \code{lambda.min} selected by CV on all data.
#' }
#'
#' @export
train_enet_cv <- function(dge_train, k_fold_index = NULL,
                          age_col = "age", donor_col = "donor",
                          alpha_fixed = 0.5,
                          compute_oof = TRUE,
                          seed = 1, verbose=FALSE) {

    stopifnot("DGEList" %in% class(dge_train))
    set.seed(seed)

    prep <- .prepare_design_from_dge(dge_train, age_col = age_col, donor_col = donor_col, prior.count = 1)
    smp <- prep$smp
    donor <- prep$donor
    age <- prep$age
    X_all <- prep$X_all

    per_fold <- NULL
    overall <- NULL
    oof_df <- NULL
    model_df <- NULL
    fold_models <- NULL

    if (isTRUE(compute_oof)) {

        stopifnot(!is.null(k_fold_index))
        stopifnot(!is.null(names(k_fold_index)))

        fold_id <- as.integer(k_fold_index[donor])
        if (any(is.na(fold_id))) stop("Missing fold assignment for some donors")
        K <- max(fold_id, na.rm = TRUE)

        oof_pred <- rep(NA_real_, length(age))
        per_fold <- data.frame(
            fold = integer(0),
            n_train = integer(0),
            n_test = integer(0),
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

            Xtr_s <- .scale_train(X_all[tr, , drop = FALSE])
            Xte <- .scale_apply(X_all[te, , drop = FALSE], Xtr_s$mu, Xtr_s$s)
            ytr <- age[tr]
            yte <- age[te]

            cv <- glmnet::cv.glmnet(
                Xtr_s$X, ytr, family = "gaussian",
                alpha = alpha_fixed,
                standardize = FALSE,
                type.measure = "mae"
            )

            fit <- glmnet::glmnet(
                Xtr_s$X, ytr, family = "gaussian",
                alpha = alpha_fixed,
                lambda = cv$lambda.min,
                standardize = FALSE
            )

            pred <- as.numeric(predict(fit, Xte))
            oof_pred[te] <- pred

            m <- compute_age_metrics(pred, yte)
            per_fold[nrow(per_fold) + 1, ] <- cbind(
                fold = k,
                n_train = length(tr),
                n_test = length(te),
                lambda = cv$lambda.min,
                m
            )

            fold_predictions[[k]] <- data.frame(
                donor = donor[te], age = yte, pred = pred, fold = k,
                row.names = rownames(smp)[te],
                check.names = FALSE
            )

            feats_k <- rownames(fit$beta)
            fold_models[[k]] <- data.frame(
                feature = feats_k,
                coef = as.numeric(fit$beta),
                mean_train = Xtr_s$mu[feats_k],
                sd_train = Xtr_s$s[feats_k],
                intercept = as.numeric(fit$a0),
                fold = k,
                row.names = NULL,
                check.names = FALSE
            )
        }

        overall <- compute_age_metrics(oof_pred, age)
        oof_df <- do.call(rbind, fold_predictions)
        model_df <- do.call(rbind, fold_models)
    }

    # Always refit final model on all provided data
    if (verbose)
        logger::log_info("Refitting final model using all training data.")
    final_fit <- .fit_final_model_fixed_alpha(X_all, age, alpha_fixed = alpha_fixed)
    final_model <- final_fit$final_model
    final_model$cor <- compute_age_expression_correlation(X_all[, final_model$feature, drop = FALSE], age)
    if (verbose)
        logger::log_info("Final model lambda: {final_fit$final_lambda}")

    list(
        per_fold_metrics = per_fold,
        overall_metrics = overall,
        donor_age_predictions = oof_df,
        fold_models = fold_models,
        model_df = model_df,
        final_model = final_model,
        final_alpha = final_fit$final_alpha,
        final_lambda = final_fit$final_lambda
    )
}

#' Train elastic net age model with K-fold cross-validation (alpha grid search)
#'
#' Trains an elastic net age prediction model using donor-level K-fold
#' cross-validation for out-of-fold (OOF) performance estimation, while selecting
#' \code{alpha} by grid search. For each fold, \code{alpha} and \code{lambda} are
#' selected using training donors only, then a fold model is fit and evaluated on
#' held-out donors.
#'
#' The final model is always refit on all provided samples after selecting the
#' best \code{alpha} and \code{lambda} on the full dataset.
#'
#' Set \code{compute_oof = FALSE} to skip fold-based OOF prediction and return
#' only the final refit model (useful for repeated Monte Carlo splits).
#'
#' @param dge_train A \code{DGEList} with donor-level samples and ages stored in
#'   \code{dge_train$samples}.
#' @param k_fold_index Named integer vector mapping donor IDs to fold numbers
#'   1..K. Required when \code{compute_oof = TRUE}. Ignored when
#'   \code{compute_oof = FALSE}.
#' @param age_col Column name in \code{dge_train$samples} containing donor age.
#' @param donor_col Column name in \code{dge_train$samples} containing donor IDs.
#' @param alpha_grid Numeric vector of candidate alpha values.
#' @param compute_oof Logical; if \code{TRUE}, compute OOF predictions/metrics
#'   using \code{k_fold_index}. If \code{FALSE}, skip OOF computation.
#' @param seed Integer seed for reproducibility.
#' @param verbose Logical; if \code{TRUE}, log progress messages.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{per_fold_metrics}: data.frame of metrics per fold (or \code{NULL}
#'     if \code{compute_oof = FALSE}).
#'   \item \code{overall_metrics}: metrics on OOF predictions across all donors
#'     (or \code{NULL} if \code{compute_oof = FALSE}).
#'   \item \code{donor_age_predictions}: OOF prediction data.frame with columns
#'     \code{donor}, \code{age}, \code{pred}, \code{fold} (or \code{NULL} if
#'     \code{compute_oof = FALSE}).
#'   \item \code{fold_models}: list of per-fold model coefficient data.frames
#'     (or \code{NULL} if \code{compute_oof = FALSE}).
#'   \item \code{model_df}: combined long-form coefficient data.frame across folds
#'     (or \code{NULL} if \code{compute_oof = FALSE}).
#'   \item \code{final_model}: data.frame representing the final model fit on all
#'     provided data, including scaling parameters and intercept.
#'   \item \code{final_alpha}: selected alpha for the final model.
#'   \item \code{final_lambda}: selected \code{lambda.min} for the final model.
#' }
#'
#' @export
train_enet_cv_optimize_alpha <- function(dge_train, k_fold_index = NULL,
                                         age_col = "age", donor_col = "donor",
                                         alpha_grid = seq(0, 1, by = 0.1),
                                         compute_oof = TRUE,
                                         seed = 1, verbose=FALSE) {

    stopifnot("DGEList" %in% class(dge_train))
    set.seed(seed)

    prep <- .prepare_design_from_dge(dge_train, age_col = age_col, donor_col = donor_col, prior.count = 1)
    smp <- prep$smp
    donor <- prep$donor
    age <- prep$age
    X_all <- prep$X_all

    per_fold <- NULL
    overall <- NULL
    oof_df <- NULL
    model_df <- NULL
    fold_models <- NULL

    if (isTRUE(compute_oof)) {

        stopifnot(!is.null(k_fold_index))
        stopifnot(!is.null(names(k_fold_index)))

        fold_id <- as.integer(k_fold_index[donor])
        if (any(is.na(fold_id))) stop("Missing fold assignment for some donors")
        K <- max(fold_id, na.rm = TRUE)

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

            Xtr_s <- .scale_train(X_all[tr, , drop = FALSE])
            Xte <- .scale_apply(X_all[te, , drop = FALSE], Xtr_s$mu, Xtr_s$s)
            ytr <- age[tr]
            yte <- age[te]

            # Select alpha and lambda on the training donors only
            sel <- .select_alpha_lambda(Xtr_s$X, ytr, alpha_grid = alpha_grid, seed = seed + k)

            fit <- glmnet::glmnet(
                Xtr_s$X, ytr, family = "gaussian",
                alpha = sel$alpha,
                lambda = sel$lambda_min,
                standardize = FALSE
            )

            pred <- as.numeric(predict(fit, newx = Xte))
            oof_pred[te] <- pred

            m <- compute_age_metrics(pred, yte)
            per_fold[nrow(per_fold) + 1, ] <- data.frame(
                fold = k,
                n_train = length(tr),
                n_test = length(te),
                alpha = sel$alpha,
                lambda = sel$lambda_min,
                r = m[["r"]],
                median_abs_error = m[["median_abs_error"]],
                mean_abs_error = m[["mean_abs_error"]],
                check.names = FALSE
            )

            fold_predictions[[k]] <- data.frame(
                donor = donor[te], age = yte, pred = pred, fold = k,
                row.names = rownames(smp)[te],
                check.names = FALSE
            )

            feats_k <- rownames(fit$beta)
            fold_models[[k]] <- data.frame(
                feature = feats_k,
                coef = as.numeric(fit$beta),
                mean_train = Xtr_s$mu[feats_k],
                sd_train = Xtr_s$s[feats_k],
                intercept = as.numeric(fit$a0),
                fold = k,
                row.names = NULL,
                check.names = FALSE
            )
        }

        overall <- compute_age_metrics(oof_pred, age)
        oof_df <- do.call(rbind, fold_predictions)
        model_df <- do.call(rbind, fold_models)
    }

    # Always refit final model on all provided data (select alpha+lambda on full data)
    if (verbose)
        logger::log_info("Refitting final model using all training data.")
    final_fit <- .fit_final_model_opt_alpha(X_all, age, alpha_grid = alpha_grid, seed = seed)
    final_model <- final_fit$final_model
    final_model$cor <- compute_age_expression_correlation(X_all[, final_model$feature, drop = FALSE], age)
    if (verbose)
        logger::log_info("Final model alpha: {final_fit$final_alpha}, lambda: {final_fit$final_lambda}")

    list(
        per_fold_metrics = per_fold,
        overall_metrics = overall,
        donor_age_predictions = oof_df,
        fold_models = fold_models,
        model_df = model_df,
        final_model = final_model,
        final_alpha = final_fit$final_alpha,
        final_lambda = final_fit$final_lambda
    )
}

#' Compute donor residuals using repeated stratified 80/20 splits
#'
#' Runs repeated Monte Carlo cross-validation by repeatedly splitting donors into
#' stratified train/test sets, fitting an age model on the training donors, and
#' predicting ages for held-out donors. For each donor, out-of-fold (held-out)
#' predictions are accumulated across repeats, enabling per-donor residual
#' summaries and stability estimates with respect to training-set composition.
#'
#' This function is intended to produce donor-level residuals for downstream
#' biological interpretation (e.g., "older/younger than chronological age") while
#' reducing dependence on any single fold assignment.
#'
#' @param dge_cell A donor-level \code{DGEList} (one sample per donor), typically
#'   produced by collapsing multiple observations per donor.
#' @param age_col Column name in \code{dge_cell$samples} containing donor age.
#' @param donor_col Column name in \code{dge_cell$samples} containing donor IDs.
#' @param n_bins Number of age quantile bins used for stratified splitting.
#' @param test_prop Fraction of donors assigned to the test set in each repeat.
#' @param n_repeats Integer; number of repeated splits to perform.
#' @param optimize_alpha Logical; if \code{TRUE}, fit each split using
#'   \code{train_enet_cv_optimize_alpha(..., compute_oof = FALSE)}. If \code{FALSE},
#'   use \code{train_enet_cv(..., compute_oof = FALSE)} with \code{alpha_fixed}.
#' @param alpha_fixed Elastic net alpha used when \code{optimize_alpha = FALSE}.
#' @param alpha_grid Numeric vector of candidate alpha values used when
#'   \code{optimize_alpha = TRUE}.
#' @param seed Integer seed controlling the sequence of splits and model fitting.
#'   Each repeat uses \code{seed + mc_iter - 1}.
#' @param prior.count Prior count passed to \code{predict_age_from_dge()}.
#' @param verbose_every Integer; if positive, log progress every
#'   \code{verbose_every} repeats.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{donor_oof_predictions}: long-format data.frame with one row per
#'     held-out donor per repeat, containing \code{donor}, \code{age}, \code{pred},
#'     \code{resid} (pred - age), and \code{mc_iter}.
#'   \item \code{donor_residual_summary}: data.frame with one row per donor,
#'     containing \code{age}, \code{n_oof}, \code{pred_mean}, \code{resid_mean},
#'     \code{resid_median}, and \code{resid_sd}.
#'   \item \code{per_repeat_metrics}: data.frame of metrics for each repeat split.
#'   \item \code{avg_pred_metrics}: metrics computed using per-donor average
#'     predictions (\code{pred_mean}) across repeats.
#'   \item \code{params}: list of the key parameter values used.
#' }
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
                                                    alpha_grid = seq(0, 1, by = 0.1),
                                                    seed = 1,
                                                    prior.count = 1,
                                                    verbose_every = 10) {

    stopifnot("DGEList" %in% class(dge_cell))
    stopifnot(test_prop > 0, test_prop < 1)
    stopifnot(n_repeats >= 1)

    oof_list <- vector("list", n_repeats)
    metrics_list <- vector("list", n_repeats)

    for (mc_iter in seq_len(n_repeats)) {

        seed_i <- seed + mc_iter - 1L

        if (!is.null(verbose_every) && verbose_every > 0 && (mc_iter %% verbose_every == 0)) {
            logger::log_info(paste("Repeated split", mc_iter, "of", n_repeats))
        }

        h <- make_holdout_stratified(dge_cell,
                                     age_col = age_col,
                                     donor_col = donor_col,
                                     test_prop = test_prop,
                                     n_bins = n_bins,
                                     seed = seed_i)

        dge_train <- h$dge_train
        dge_test <- h$dge_test

        if (optimize_alpha) {
            fit_obj <- train_enet_cv_optimize_alpha(
                dge_train,
                k_fold_index = NULL,
                age_col = age_col,
                donor_col = donor_col,
                alpha_grid = alpha_grid,
                compute_oof = FALSE,
                seed = seed_i
            )
        } else {
            fit_obj <- train_enet_cv(
                dge_train,
                k_fold_index = NULL,
                age_col = age_col,
                donor_col = donor_col,
                alpha_fixed = alpha_fixed,
                compute_oof = FALSE,
                seed = seed_i
            )
        }

        pred_df <- predict_age_from_dge(dge_test, fit_obj$final_model, prior.count = prior.count)
        pred_df[["mc_iter"]] <- mc_iter
        pred_df[["resid"]] <- pred_df[["pred"]] - pred_df[["age"]]

        oof_list[[mc_iter]] <- pred_df[, c("donor", "age", "pred", "resid", "mc_iter")]

        m <- compute_age_metrics(pred_df$pred, pred_df$age)
        metrics_list[[mc_iter]] <- cbind(mc_iter = mc_iter, m)
    }

    oof_all <- do.call(rbind, oof_list)
    per_repeat_metrics <- do.call(rbind, metrics_list)

    age_unique <- unique(dge_cell$samples[, c(donor_col, age_col), drop = FALSE])
    colnames(age_unique) <- c("donor", "age")
    age_unique$donor <- as.character(age_unique$donor)
    age_unique$age <- as.numeric(age_unique$age)

    split_counts <- aggregate(mc_iter ~ donor, data = oof_all, FUN = length)
    colnames(split_counts)[2] <- "n_oof"

    pred_mean <- aggregate(pred ~ donor, data = oof_all, FUN = mean)
    colnames(pred_mean)[2] <- "pred_mean"

    resid_mean <- aggregate(resid ~ donor, data = oof_all, FUN = mean)
    colnames(resid_mean)[2] <- "resid_mean"

    resid_median <- aggregate(resid ~ donor, data = oof_all, FUN = stats::median)
    colnames(resid_median)[2] <- "resid_median"

    resid_sd <- aggregate(resid ~ donor, data = oof_all, FUN = stats::sd)
    colnames(resid_sd)[2] <- "resid_sd"

    donor_summary <- merge(age_unique, split_counts, by = "donor", all.x = TRUE)
    donor_summary <- merge(donor_summary, pred_mean, by = "donor", all.x = TRUE)
    donor_summary <- merge(donor_summary, resid_mean, by = "donor", all.x = TRUE)
    donor_summary <- merge(donor_summary, resid_median, by = "donor", all.x = TRUE)
    donor_summary <- merge(donor_summary, resid_sd, by = "donor", all.x = TRUE)

    avg_pred_metrics <- compute_age_metrics(donor_summary$pred_mean, donor_summary$age)

    list(
        donor_oof_predictions = oof_all,
        donor_residual_summary = donor_summary,
        per_repeat_metrics = per_repeat_metrics,
        avg_pred_metrics = avg_pred_metrics,
        params = list(
            test_prop = test_prop,
            n_bins = n_bins,
            n_repeats = n_repeats,
            optimize_alpha = optimize_alpha,
            alpha_fixed = alpha_fixed,
            alpha_grid = alpha_grid,
            seed = seed,
            prior.count = prior.count
        )
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
            sd_override <- apply(lcpm_override[common, , drop = FALSE], 1, stats::sd)
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


