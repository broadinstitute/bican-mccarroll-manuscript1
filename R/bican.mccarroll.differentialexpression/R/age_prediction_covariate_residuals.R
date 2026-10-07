# Residualize donor-level age predictions on donor/sample covariates (PMI, RQS,
# library size), to test whether coordinated deviations in predicted age across
# cell types and regions persist after accounting for technical covariates.

#' Load PMI, RQS and library size covariates from a saved DGEList
#'
#' Loads a DGEList with \code{loadDGEList()} (typically a PMI/RQS
#' complete-cases subset), which recomputes \code{lib.size} from the counts, and
#' returns one row per donor x cell type x region. Library size is summed over
#' the donor's samples (reactions) within each cell type x region, matching the
#' donor-level \code{num_umis} reported by the age prediction. PMI and RQS are
#' donor-level traits, so they must be constant across a donor's sample rows;
#' this is checked.
#'
#' @param data_dir Directory containing the DGEList files (see \code{loadDGEList()}).
#' @param data_prefix Filename prefix of the DGEList. Default
#'   \code{"age_pmi_rqs_complete_cases_DGEList"}.
#' @param pmi_col PMI column name in the samples table. Default \code{"pmi_hr"}.
#' @param rqs_col RQS column name in the samples table. Default \code{"RQS"}.
#'
#' @return A data.frame with columns \code{donor}, \code{cell_type},
#'   \code{region}, \code{pmi}, \code{rqs}, \code{lib_size}, restricted to
#'   rows with non-missing PMI and RQS.
#' @export
load_age_prediction_covariates <- function(data_dir,
                                           data_prefix = "age_pmi_rqs_complete_cases_DGEList",
                                           pmi_col = "pmi_hr",
                                           rqs_col = "RQS") {
  dge <- loadDGEList(dir = data_dir, prefix = data_prefix)
  s <- dge$samples

  needed <- c("donor", "cell_type", "region", pmi_col, rqs_col, "lib.size")
  missing_cols <- setdiff(needed, colnames(s))
  if (length(missing_cols) > 0) {
    stop("DGEList samples is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  s <- data.frame(
    donor = as.character(s$donor),
    cell_type = as.character(s$cell_type),
    region = as.character(s$region),
    pmi = as.numeric(s[[pmi_col]]),
    rqs = as.numeric(s[[rqs_col]]),
    lib_size = as.numeric(s$lib.size),
    stringsAsFactors = FALSE
  )
  s <- s[stats::complete.cases(s), , drop = FALSE]

  donor_traits <- unique(s[, c("donor", "pmi", "rqs")])
  if (anyDuplicated(donor_traits$donor)) {
    dups <- unique(donor_traits$donor[duplicated(donor_traits$donor)])
    stop(
      "PMI/RQS are not constant within donor for: ",
      paste(utils::head(dups, 5), collapse = ", "),
      call. = FALSE
    )
  }

  lib <- stats::aggregate(lib_size ~ donor + cell_type + region, data = s, FUN = sum)
  out <- merge(lib, donor_traits, by = "donor", all.x = TRUE, sort = FALSE)
  out <- out[, c("donor", "cell_type", "region", "pmi", "rqs", "lib_size")]
  rownames(out) <- NULL
  out
}

#' Remove PMI, RQS and library size effects from donor age predictions
#'
#' For every cell type x region, fits
#' \code{pred_mean ~ age + z_pmi + z_rqs + z_loglib} by ordinary least squares
#' across donors, where \code{z_*} are z-scores (within the cell type x region)
#' of PMI, RQS and \code{log10(lib_size)}, and \code{lib_size} is the donor's
#' library size for that cell type x region. Age is in the model so the covariate
#' effects are estimated over and above age, but only the covariate terms are
#' subtracted (\code{covariate_adjustment}, the sum of coefficient x z-score over
#' the three covariates). Because the z-scores are centered, the mean prediction
#' is unchanged and the adjusted values stay on the original age scale.
#'
#' The same adjustment is subtracted from \code{pred_mean}, \code{resid_mean},
#' \code{pred_mean_corrected} and \code{resid_mean_corrected}. The age-bias
#' (GAM) correction itself is not refit, so the adjusted columns can be passed
#' straight to the existing residual plotting functions, and differ from the
#' published ones only by the covariate adjustment.
#'
#' Only donors present in \code{covariates} are retained. Cell type x region
#' groups with fewer than \code{min_donors} such donors are dropped.
#'
#' @param donor_predictions data.frame of donor-level predictions, as read from
#'   \code{age_prediction_results_donor_predictions.txt}.
#' @param covariates data.frame from \code{load_age_prediction_covariates()}
#'   with columns \code{donor}, \code{cell_type}, \code{region}, \code{pmi},
#'   \code{rqs}, \code{lib_size}.
#' @param min_donors Minimum donors required to fit a group. Default 10.
#'
#' @return A list with:
#' \describe{
#'   \item{donor_predictions}{The input columns for the retained donors, with
#'     \code{pred_mean}, \code{resid_mean}, \code{pred_mean_corrected} and
#'     \code{resid_mean_corrected} adjusted, plus \code{covariate_adjustment}
#'     and \code{resid_mean_corrected_unadjusted} (the input
#'     \code{resid_mean_corrected} for the same donors). \code{resid_median} and
#'     \code{resid_sd} are dropped because they no longer correspond to the
#'     adjusted predictions.}
#'   \item{covariate_fits}{One row per cell type x region x covariate with the
#'     estimate (change in \code{pred_mean} per SD of the covariate, in the
#'     units of \code{pred_mean}), standard error, p-value, Benjamini-Hochberg adjusted p-value
#'     (\code{p_value_bh}, across all tests in the table) and number of donors.}
#' }
#' @export
residualize_age_predictions_on_covariates <- function(donor_predictions,
                                                      covariates,
                                                      min_donors = 10) {
  required <- c(
    "cell_type", "region", "donor", "age",
    "pred_mean", "resid_mean", "pred_mean_corrected", "resid_mean_corrected"
  )
  missing_cols <- setdiff(required, colnames(donor_predictions))
  if (length(missing_cols) > 0) {
    stop("donor_predictions is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  cov_cols <- c("donor", "cell_type", "region", "pmi", "rqs", "lib_size")
  if (!all(cov_cols %in% colnames(covariates))) {
    stop("covariates must have columns ", paste(cov_cols, collapse = ", "), call. = FALSE)
  }

  dp <- merge(
    donor_predictions, covariates[, cov_cols],
    by = c("donor", "cell_type", "region"), all = FALSE, sort = FALSE
  )
  dp <- dp[is.finite(dp$lib_size) & dp$lib_size > 0 & is.finite(dp$pred_mean), , drop = FALSE]

  group_key <- paste(dp$cell_type, dp$region, sep = "\r")
  groups <- split(seq_len(nrow(dp)), group_key)

  covariate_terms <- c("z_pmi", "z_rqs", "z_loglib")

  adj_list <- list()
  fit_list <- list()

  for (idx in groups) {
    d <- dp[idx, , drop = FALSE]
    ct <- d$cell_type[1]
    rg <- d$region[1]

    if (nrow(d) < min_donors) {
      logger::log_warn("Skipping {ct} {rg}: {nrow(d)} donors with complete covariates (< {min_donors})")
      next
    }

    d$z_pmi <- as.numeric(scale(d$pmi))
    d$z_rqs <- as.numeric(scale(d$rqs))
    d$z_loglib <- as.numeric(scale(log10(d$lib_size)))

    if (anyNA(d[, covariate_terms])) {
      logger::log_warn("Skipping {ct} {rg}: a covariate has zero variance")
      next
    }

    fit <- stats::lm(pred_mean ~ age + z_pmi + z_rqs + z_loglib, data = d)
    beta <- stats::coef(fit)[covariate_terms]
    adjustment <- as.numeric(as.matrix(d[, covariate_terms]) %*% beta)

    d$resid_mean_corrected_unadjusted <- d$resid_mean_corrected
    d$covariate_adjustment <- adjustment
    for (v in c("pred_mean", "resid_mean", "pred_mean_corrected", "resid_mean_corrected")) {
      d[[v]] <- d[[v]] - adjustment
    }

    adj_list[[length(adj_list) + 1]] <- d

    sm <- summary(fit)$coefficients[covariate_terms, , drop = FALSE]
    fit_list[[length(fit_list) + 1]] <- data.frame(
      cell_type = ct,
      region = rg,
      covariate = unname(c(z_pmi = "PMI", z_rqs = "RQS", z_loglib = "log10_lib_size")[covariate_terms]),
      estimate_per_sd = sm[, "Estimate"],
      std_error = sm[, "Std. Error"],
      p_value = sm[, "Pr(>|t|)"],
      n_donors = nrow(d),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }

  if (length(adj_list) == 0) {
    stop("No cell type x region group had enough donors with complete covariates", call. = FALSE)
  }

  out <- do.call(rbind, adj_list)
  out <- out[, setdiff(colnames(out), c("pmi", "rqs", "lib_size", covariate_terms, "resid_median", "resid_sd")), drop = FALSE]
  rownames(out) <- NULL

  covariate_fits <- do.call(rbind, fit_list)
  covariate_fits$p_value_bh <- stats::p.adjust(covariate_fits$p_value, method = "BH")

  list(
    donor_predictions = out,
    covariate_fits = covariate_fits
  )
}
