library (data.table)

#' Test in-group versus out-group eQTL slopes across cell types
#'
#' Fits one mixed model for each gene-variant pair in a gathered eQTL data set.
#' The model estimates cell-type-specific intercepts, a pooled out-group genotype
#' slope, and an in-group genotype-slope deviation while accounting for repeated
#' donors with a donor random intercept.
#'
#' @param gathered Either the list returned by `build_cluster_dataframe()` or a
#'   data.frame containing the gathered data. Supply either `gathered` or
#'   `gathered_file`, but not both.
#' @param gathered_file Optional path to a plain or gzipped TSV containing the
#'   gathered data.
#' @param covariate_mapping Optional covariate-mapping data.frame. This is only
#'   needed when `gathered` is a data.frame rather than the complete list returned
#'   by `build_cluster_dataframe()`.
#' @param covariate_mapping_file Optional path to a plain or gzipped TSV
#'   containing the covariate mapping.
#' @param force_factor_covariates Optional character vector of covariate columns
#'   that must be coerced to factors before fitting. Unknown columns are errors.
#' @param results_outfile Required writable path for the gene-level results TSV
#'   or TSV.GZ file. Existing files are overwritten.
#' @param comparators_outfile Required writable path for the long comparator TSV
#'   or TSV.GZ file. Existing files are overwritten.
#'
#' @param outPDF File path for an output PDF
#'
#' @return Invisibly returns a list with `results` and `comparators` data.frames.
#' @export
run_cell_type_specific_eqtl_tests <- function(
    gathered = NULL,
    gathered_file = NULL,
    covariate_mapping = NULL,
    covariate_mapping_file = NULL,
    force_factor_covariates = NULL,
    results_outfile,
    comparators_outfile=NULL,
    outPDF) {

    .require_testing_packages()
    .validate_output_path(results_outfile, "results_outfile")
    .validate_output_path(comparators_outfile, "comparators_outfile")

    logger::log_info("Fetching Data for analysis")

    inputs <- .resolve_testing_inputs(
        gathered = gathered,
        gathered_file = gathered_file,
        covariate_mapping = covariate_mapping,
        covariate_mapping_file = covariate_mapping_file
    )

    eqtl_data <- data.table::as.data.table(inputs$data)
    mapping <- data.table::as.data.table(inputs$covariate_mapping)

    prepared <- .prepare_testing_data(
        eqtl_data = eqtl_data,
        covariate_mapping = mapping,
        force_factor_covariates = force_factor_covariates
    )

    eqtl_data <- prepared$data
    covariates <- prepared$covariates

    pairs <- unique(eqtl_data[, .(gene, variant_id, cluster)])
    data.table::setorderv(pairs, c("gene", "variant_id"))

    result_rows <- vector("list", nrow(pairs))
    comparator_rows <- vector("list", nrow(pairs))

    logger::log_info("Running regressions")
    for (i in seq_len(nrow(pairs))) {
        pair <- pairs[i]
        pair_data <- eqtl_data[
            gene == pair$gene & variant_id == pair$variant_id
        ]

        analyzed <- .analyze_gene_variant(
            pair_data = pair_data,
            gene = pair$gene,
            variant_id = pair$variant_id,
            cluster = pair$cluster,
            covariates = covariates
        )

        result_rows[[i]] <- analyzed$result
        comparator_rows[[i]] <- analyzed$comparators
    }
    logger::log_info("Finished running regressions")

    results <- data.table::rbindlist(result_rows, use.names = TRUE, fill = TRUE)
    comparators <- data.table::rbindlist(
        comparator_rows,
        use.names = TRUE,
        fill = TRUE
    )

    results[, FDR := NA_real_]
    successful <- is.na(results$failure_reason) & !is.na(results$p_value)
    results[successful, FDR := stats::p.adjust(p_value, method = "BH")]

    data.table::setorderv(results, c("gene", "variant_id"))
    data.table::setorderv(
        comparators,
        c("gene", "variant_id", "cell_type", "region")
    )

    data.table::fwrite(
        results,
        results_outfile,
        sep = "\t",
        na = "NA",
        compress = "auto"
    )
    if (!is.null(comparators_outfile)) {
        data.table::fwrite(
            comparators,
            comparators_outfile,
            sep = "\t",
            na = "NA",
            compress = "auto"
        )
    }

    #plot results
    if (!is.null(outPDF)) {
        logger::log_info("Writing PDF Report")
        write_eqtl_summary_pdf(results, eqtl_data, outfile=outPDF, show_interaction=T)
    }

    logger::log_info("Finished Analysis")

    invisible(list(
        results = as.data.frame(results),
        comparators = as.data.frame(comparators)
    ))
}

.require_testing_packages <- function() {
    required <- c("data.table", "lme4", "lmerTest", "logger")
    missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]

    if (length(missing)) {
        stop(
            "Required package(s) not installed: ",
            paste(missing, collapse = ", "),
            call. = FALSE
        )
    }
}

.resolve_testing_inputs <- function(
    gathered,
    gathered_file,
    covariate_mapping,
    covariate_mapping_file) {

    if (!is.null(gathered) && !is.null(gathered_file)) {
        stop("Supply only one of gathered or gathered_file.", call. = FALSE)
    }
    if (is.null(gathered) && is.null(gathered_file)) {
        stop("Supply either gathered or gathered_file.", call. = FALSE)
    }

    data <- NULL
    mapping_from_list <- NULL

    if (!is.null(gathered_file)) {
        .validate_input_file(gathered_file, "gathered_file")
        data <- data.table::fread(
            gathered_file,
            sep = "\t",
            data.table = FALSE,
            showProgress = FALSE
        )
    } else if (is.list(gathered) && all(c("data", "covariate_mapping") %in% names(gathered))) {
        data <- gathered$data
        mapping_from_list <- gathered$covariate_mapping
    } else if (is.data.frame(gathered) || data.table::is.data.table(gathered)) {
        data <- gathered
    } else {
        stop(
            "gathered must be a data.frame or the list returned by build_cluster_dataframe().",
            call. = FALSE
        )
    }

    supplied_mapping_count <- sum(c(
        !is.null(mapping_from_list),
        !is.null(covariate_mapping),
        !is.null(covariate_mapping_file)
    ))

    if (supplied_mapping_count == 0L) {
        stop(
            "A covariate mapping must be supplied directly, by file, or inside gathered.",
            call. = FALSE
        )
    }
    if (supplied_mapping_count > 1L) {
        stop("Supply only one covariate mapping source.", call. = FALSE)
    }

    if (!is.null(covariate_mapping_file)) {
        .validate_input_file(covariate_mapping_file, "covariate_mapping_file")
        mapping <- data.table::fread(
            covariate_mapping_file,
            sep = "\t",
            data.table = FALSE,
            showProgress = FALSE
        )
    } else if (!is.null(covariate_mapping)) {
        mapping <- covariate_mapping
    } else {
        mapping <- mapping_from_list
    }

    if (!is.data.frame(data) && !data.table::is.data.table(data)) {
        stop("Gathered data must be a data.frame.", call. = FALSE)
    }
    if (!is.data.frame(mapping) && !data.table::is.data.table(mapping)) {
        stop("Covariate mapping must be a data.frame.", call. = FALSE)
    }

    list(data = data, covariate_mapping = mapping)
}

.prepare_testing_data <- function(
        eqtl_data,
        covariate_mapping,
        force_factor_covariates = NULL
) {
    covariates <- unique(covariate_mapping$output_column)

    factor_covariates <- unique(c(
        covariate_mapping$output_column[
            covariate_mapping$encoding == "factor"
        ],
        force_factor_covariates
    ))

    unknown_forced <- setdiff(force_factor_covariates, covariates)

    if (length(unknown_forced)) {
        stop(
            "Unknown forced factor covariates: ",
            paste(unknown_forced, collapse = ", "),
            call. = FALSE
        )
    }

    numeric_covariates <- setdiff(
        covariate_mapping$output_column[
            covariate_mapping$encoding == "numeric"
        ],
        factor_covariates
    )

    if (length(factor_covariates)) {
        eqtl_data[, (factor_covariates) := lapply(
            .SD,
            factor
        ), .SDcols = factor_covariates]
    }

    if (length(numeric_covariates)) {
        eqtl_data[, (numeric_covariates) := lapply(
            .SD,
            function(x) as.numeric(scale(x))
        ), .SDcols = numeric_covariates]
    }

    list(
        data = eqtl_data,
        covariates = covariates
    )
}

.analyze_gene_variant <- function(
    pair_data,
    gene,
    variant_id,
    cluster,
    covariates) {

    region_varies <- data.table::uniqueN(pair_data$region) > 1L
    complete_columns <- c(
        "expression_normalized",
        "genotype",
        "donor",
        "cell_type",
        "in_group",
        covariates
    )
    if (region_varies) {
        complete_columns <- c(complete_columns, "region")
    }

    complete_index <- stats::complete.cases(pair_data[, ..complete_columns])
    complete_data <- pair_data[complete_index]

    comparators <- .build_comparator_rows(
        pair_data = pair_data,
        complete_index = complete_index,
        gene = gene,
        variant_id = variant_id,
        cluster = cluster
    )

    n_in_rows <- nrow(complete_data[in_group == 1])
    n_out_rows <- nrow(complete_data[in_group == 0])

    if (n_in_rows == 0L) {
        logger::log_warn(
            "Skipping gene {gene} ({variant_id}): no complete-case in-group rows."
        )
        return(list(
            result = .failed_result(
                gene,
                variant_id,
                cluster,
                complete_data,
                "no_complete_case_in_group"
            ),
            comparators = comparators
        ))
    }

    if (n_out_rows == 0L) {
        return(list(
            result = .failed_result(
                gene,
                variant_id,
                cluster,
                complete_data,
                "no_complete_case_out_group"
            ),
            comparators = comparators
        ))
    }

    if (data.table::uniqueN(complete_data$genotype) < 2L) {
        return(list(
            result = .failed_result(
                gene,
                variant_id,
                cluster,
                complete_data,
                "no_genotype_variation"
            ),
            comparators = comparators
        ))
    }

    complete_data[, donor := droplevels(factor(donor))]
    complete_data[, cell_type := droplevels(factor(cell_type))]
    complete_data[, region := droplevels(factor(region))]

    formula <- .build_eqtl_formula(covariates, region_varies, complete_data)
    fit_result <- .fit_eqtl_model(formula, complete_data)

    if (!is.null(fit_result$failure_reason)) {
        return(list(
            result = .failed_result(
                gene,
                variant_id,
                cluster,
                complete_data,
                fit_result$failure_reason
            ),
            comparators = comparators
        ))
    }

    estimates <- .extract_eqtl_estimates(fit_result$fit)
    if (!is.null(estimates$failure_reason)) {
        return(list(
            result = .failed_result(
                gene,
                variant_id,
                cluster,
                complete_data,
                estimates$failure_reason
            ),
            comparators = comparators
        ))
    }

    result <- data.table::data.table(
        gene = gene,
        variant_id = variant_id,
        cluster = cluster,
        n_cell_types = data.table::uniqueN(complete_data$cell_type),
        n_in_group_cell_types = data.table::uniqueN(
            complete_data[in_group == 1, cell_type]
        ),
        n_comparator_cell_types = data.table::uniqueN(
            complete_data[in_group == 0, cell_type]
        ),
        n_unique_donors = data.table::uniqueN(complete_data$donor),
        pooled_out_group_slope = estimates$out_estimate,
        pooled_out_group_slope_se = estimates$out_se,
        pooled_out_group_slope_ci_low = estimates$out_ci_low,
        pooled_out_group_slope_ci_high = estimates$out_ci_high,
        in_group_slope = estimates$in_estimate,
        in_group_slope_se = estimates$in_se,
        in_group_slope_ci_low = estimates$in_ci_low,
        in_group_slope_ci_high = estimates$in_ci_high,
        in_group_vs_out_group_estimate = estimates$interaction_estimate,
        in_group_vs_out_group_se = estimates$interaction_se,
        in_group_vs_out_group_ci_low = estimates$interaction_ci_low,
        in_group_vs_out_group_ci_high = estimates$interaction_ci_high,
        t_stat = estimates$test_statistic,
        p_value = estimates$p_value,
        failure_reason = NA_character_
    )

    list(result = result, comparators = comparators)
}

.build_comparator_rows <- function(
    pair_data,
    complete_index,
    gene,
    variant_id,
    cluster) {

    x <- data.table::copy(pair_data)
    x[, complete_case := complete_index]

    x[, {
        n_complete <- data.table::uniqueN(donor[complete_case])
        all_expression_missing <- all(is.na(expression_normalized))

        list(
            gene = gene,
            variant_id = variant_id,
            cluster = cluster,
            in_group = as.integer(in_group[[1L]]),
            eligible = n_complete > 0L,
            reason = if (n_complete > 0L) {
                "eligible"
            } else if (all_expression_missing) {
                "all_expression_missing"
            } else {
                "no_complete_cases"
            },
            n_donors = n_complete
        )
    }, by = .(cell_type, region, cell_type_region)]
}

.build_eqtl_formula <- function(covariates, region_varies, complete_data) {
    fixed_terms <- c(
        "0 + cell_type",
        "genotype",
        "genotype:in_group",
        covariates
    )

    if (region_varies && data.table::uniqueN(complete_data$region) > 1L) {
        fixed_terms <- c(fixed_terms, "region")
    }

    stats::as.formula(
        paste(
            "expression_normalized ~",
            paste(c(fixed_terms, "(1 | donor)"), collapse = " + ")
        )
    )
}

.fit_eqtl_model <- function(formula, data) {
    warnings <- character()

    fit <- tryCatch(
        withCallingHandlers(
            suppressMessages(
                lmerTest::lmer(
                    formula,
                    data = data,
                    REML = FALSE
                )
            ),
            warning = function(w) {
                warnings <<- c(warnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            }
        ),
        error = function(e) e
    )

    if (inherits(fit, "error")) {
        return(list(
            fit = NULL,
            failure_reason = paste0(
                "model_error: ",
                conditionMessage(fit)
            ),
            singular_fit = NA,
            donor_variance = NA_real_,
            warnings = warnings
        ))
    }

    if (.has_nonconvergence(fit, warnings)) {
        return(list(
            fit = NULL,
            failure_reason = "non_convergence",
            singular_fit = NA,
            donor_variance = NA_real_,
            warnings = warnings
        ))
    }

    singular_fit <- lme4::isSingular(
        fit,
        tol = 1e-4
    )

    variance_components <- as.data.frame(
        lme4::VarCorr(fit)
    )

    donor_variance <- variance_components$vcov[
        variance_components$grp == "donor" &
            variance_components$var1 == "(Intercept)"
    ]

    donor_variance <- if (length(donor_variance)) {
        donor_variance[[1L]]
    } else {
        NA_real_
    }

    list(
        fit = fit,
        failure_reason = NULL,
        singular_fit = singular_fit,
        donor_variance = donor_variance,
        warnings = warnings
    )
}

.has_nonconvergence <- function(fit, captured_warnings) {
    opt_messages <- fit@optinfo$conv$lme4$messages
    opt_messages <- if (is.null(opt_messages)) character() else as.character(opt_messages)
    all_messages <- c(opt_messages, captured_warnings)

    if (!length(all_messages)) {
        return(FALSE)
    }

    pattern <- paste(
        c(
            "failed to converge",
            "unable to evaluate scaled gradient",
            "degenerate Hessian",
            "negative eigenvalue",
            "convergence code"
        ),
        collapse = "|"
    )

    any(grepl(pattern, all_messages, ignore.case = TRUE))
}

# .extract_eqtl_estimates <- function(fit) {
#     coefficient_names <- names(lme4::fixef(fit))
#     interaction_name <- intersect(
#         c("genotype:in_group", "in_group:genotype"),
#         coefficient_names
#     )
#
#     if (!("genotype" %in% coefficient_names)) {
#         return(list(failure_reason = "genotype_coefficient_not_estimable"))
#     }
#     if (length(interaction_name) != 1L) {
#         return(list(failure_reason = "interaction_coefficient_not_estimable"))
#     }
#
#     interaction_name <- interaction_name[[1L]]
#     beta <- lme4::fixef(fit)
#     covariance <- as.matrix(stats::vcov(fit))
#     coefficient_table <- coef(summary(fit))
#
#     if (!(interaction_name %in% rownames(coefficient_table))) {
#         return(list(failure_reason = "interaction_test_not_available"))
#     }
#
#     p_column <- grep("^Pr\\(", colnames(coefficient_table), value = TRUE)
#     statistic_column <- intersect(c("t value", "z value"), colnames(coefficient_table))
#
#     if (length(p_column) != 1L || length(statistic_column) != 1L) {
#         return(list(failure_reason = "lmerTest_statistics_not_available"))
#     }
#
#     out_estimate <- unname(beta[["genotype"]])
#     interaction_estimate <- unname(beta[[interaction_name]])
#     in_estimate <- out_estimate + interaction_estimate
#
#     out_variance <- covariance["genotype", "genotype"]
#     interaction_variance <- covariance[interaction_name, interaction_name]
#     in_variance <- out_variance + interaction_variance +
#         2 * covariance["genotype", interaction_name]
#
#     if (any(!is.finite(c(out_variance, interaction_variance, in_variance))) ||
#         any(c(out_variance, interaction_variance, in_variance) < 0)) {
#         return(list(failure_reason = "invalid_coefficient_variance"))
#     }
#
#     out_se <- sqrt(out_variance)
#     interaction_se <- sqrt(interaction_variance)
#     in_se <- sqrt(in_variance)
#     z <- stats::qnorm(0.975)
#
#     list(
#         failure_reason = NULL,
#         out_estimate = out_estimate,
#         out_se = out_se,
#         out_ci_low = out_estimate - z * out_se,
#         out_ci_high = out_estimate + z * out_se,
#         in_estimate = in_estimate,
#         in_se = in_se,
#         in_ci_low = in_estimate - z * in_se,
#         in_ci_high = in_estimate + z * in_se,
#         interaction_estimate = interaction_estimate,
#         interaction_se = interaction_se,
#         interaction_ci_low = interaction_estimate - z * interaction_se,
#         interaction_ci_high = interaction_estimate + z * interaction_se,
#         test_statistic = unname(
#             coefficient_table[interaction_name, statistic_column]
#         ),
#         p_value = unname(coefficient_table[interaction_name, p_column])
#     )
# }

.extract_eqtl_estimates <- function(fit) {
    beta <- lme4::fixef(fit)
    coefficient_names <- names(beta)

    interaction_name <- intersect(
        c("genotype:in_group", "in_group:genotype"),
        coefficient_names
    )

    if (!("genotype" %in% coefficient_names)) {
        return(list(
            failure_reason = "genotype_coefficient_not_estimable"
        ))
    }

    if (length(interaction_name) != 1L) {
        return(list(
            failure_reason = "interaction_coefficient_not_estimable"
        ))
    }

    interaction_name <- interaction_name[[1L]]

    out_contrast <- setNames(
        numeric(length(beta)),
        coefficient_names
    )
    out_contrast["genotype"] <- 1

    in_contrast <- out_contrast
    in_contrast[interaction_name] <- 1

    interaction_contrast <- setNames(
        numeric(length(beta)),
        coefficient_names
    )
    interaction_contrast[interaction_name] <- 1

    out_test <- lmerTest::contest1D(
        fit,
        out_contrast,
        confint = TRUE,
        level = 0.95
    )

    in_test <- lmerTest::contest1D(
        fit,
        in_contrast,
        confint = TRUE,
        level = 0.95
    )

    interaction_test <- lmerTest::contest1D(
        fit,
        interaction_contrast,
        confint = TRUE,
        level = 0.95
    )

    list(
        failure_reason = NULL,

        out_estimate = unname(out_test[["Estimate"]]),
        out_se = unname(out_test[["Std. Error"]]),
        out_ci_low = unname(out_test[["lower"]]),
        out_ci_high = unname(out_test[["upper"]]),

        in_estimate = unname(in_test[["Estimate"]]),
        in_se = unname(in_test[["Std. Error"]]),
        in_ci_low = unname(in_test[["lower"]]),
        in_ci_high = unname(in_test[["upper"]]),

        interaction_estimate =
            unname(interaction_test[["Estimate"]]),
        interaction_se =
            unname(interaction_test[["Std. Error"]]),
        interaction_ci_low =
            unname(interaction_test[["lower"]]),
        interaction_ci_high =
            unname(interaction_test[["upper"]]),

        t_stat = unname(interaction_test[["t value"]]),
        denominator_df = unname(interaction_test[["df"]]),
        p_value = unname(interaction_test[["Pr(>|t|)"]])
    )
}

.failed_result <- function(
    gene,
    variant_id,
    cluster,
    complete_data,
    failure_reason) {

    data.table::data.table(
        gene = gene,
        variant_id = variant_id,
        cluster = cluster,
        n_cell_types = data.table::uniqueN(complete_data$cell_type),
        n_in_group_cell_types = data.table::uniqueN(
            complete_data[in_group == 1, cell_type]
        ),
        n_comparator_cell_types = data.table::uniqueN(
            complete_data[in_group == 0, cell_type]
        ),
        n_unique_donors = data.table::uniqueN(complete_data$donor),
        pooled_out_group_slope = NA_real_,
        pooled_out_group_slope_se = NA_real_,
        pooled_out_group_slope_ci_low = NA_real_,
        pooled_out_group_slope_ci_high = NA_real_,
        in_group_slope = NA_real_,
        in_group_slope_se = NA_real_,
        in_group_slope_ci_low = NA_real_,
        in_group_slope_ci_high = NA_real_,
        in_group_vs_out_group_estimate = NA_real_,
        in_group_vs_out_group_se = NA_real_,
        in_group_vs_out_group_ci_low = NA_real_,
        in_group_vs_out_group_ci_high = NA_real_,
        test_statistic = NA_real_,
        p_value = NA_real_,
        failure_reason = failure_reason
    )
}

.validate_input_file <- function(path, argument_name) {
    if (length(path) != 1L || is.na(path) || !nzchar(path)) {
        stop(argument_name, " must be one non-empty path.", call. = FALSE)
    }
    if (!file.exists(path)) {
        stop(argument_name, " does not exist: ", path, call. = FALSE)
    }
}

.validate_output_path <- function(path, argument_name) {
    if (is.null(path))
        return (TRUE)

    if (missing(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
        stop(argument_name, " must be one non-empty path.", call. = FALSE)
    }

    directory <- dirname(path)
    if (!dir.exists(directory)) {
        stop(
            "Output directory does not exist for ", argument_name, ": ",
            directory,
            call. = FALSE
        )
    }

    probe <- tempfile(pattern = ".write_test_", tmpdir = directory)
    created <- file.create(probe)
    if (!created) {
        stop(
            "Output directory is not writable for ", argument_name, ": ",
            directory,
            call. = FALSE
        )
    }
    unlink(probe)
}

.require_columns <- function(x, required, label) {
    missing <- setdiff(required, names(x))
    if (length(missing)) {
        stop(
            label,
            " is missing required columns: ",
            paste(missing, collapse = ", "),
            call. = FALSE
        )
    }
}

# z-scores numeric covariates.
.prepare_covariates <- function(data, mapping) {
    required_mapping_columns <- c(
        "output_column",
        "encoding"
    )

    missing_mapping_columns <- setdiff(
        required_mapping_columns,
        names(mapping)
    )

    if (length(missing_mapping_columns)) {
        stop(
            "Covariate mapping is missing required columns: ",
            paste(missing_mapping_columns, collapse = ", "),
            call. = FALSE
        )
    }

    covariate_columns <- unique(mapping$output_column)

    missing_covariates <- setdiff(
        covariate_columns,
        names(data)
    )

    if (length(missing_covariates)) {
        stop(
            "Mapped covariates not found in data: ",
            paste(missing_covariates, collapse = ", "),
            call. = FALSE
        )
    }

    unsupported_encodings <- setdiff(
        unique(mapping$encoding),
        c("numeric", "factor")
    )

    if (length(unsupported_encodings)) {
        stop(
            "Unsupported covariate encodings: ",
            paste(unsupported_encodings, collapse = ", "),
            call. = FALSE
        )
    }

    factor_columns <- unique(
        mapping$output_column[mapping$encoding == "factor"]
    )

    numeric_columns <- unique(
        mapping$output_column[mapping$encoding == "numeric"]
    )

    overlap <- intersect(factor_columns, numeric_columns)

    if (length(overlap)) {
        stop(
            "Covariates have conflicting encodings: ",
            paste(overlap, collapse = ", "),
            call. = FALSE
        )
    }

    for (column in factor_columns) {
        data[[column]] <- factor(data[[column]])
    }

    for (column in numeric_columns) {
        x <- data[[column]]

        if (!is.numeric(x)) {
            stop(
                "Mapped numeric covariate '",
                column,
                "' is not numeric.",
                call. = FALSE
            )
        }

        x_sd <- stats::sd(x, na.rm = TRUE)

        if (!is.finite(x_sd) || x_sd == 0) {
            stop(
                "Mapped numeric covariate '",
                column,
                "' has zero or undefined standard deviation.",
                call. = FALSE
            )
        }

        data[[column]] <- as.numeric(scale(x))
    }

    list(
        data = data,
        covariates = covariate_columns,
        numeric_covariates = numeric_columns,
        factor_covariates = factor_columns
    )
}

##################
# VISUALIZATION
###################

plot_eqtl_interaction_result <- function(result_row) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.", call. = FALSE)
    }

    if (is.null(result_row) || nrow(result_row) != 1L) {
        stop("`result_row` must contain exactly one row.", call. = FALSE)
    }

    if ("failure_reason" %in% names(result_row)) {
        failure_value <- result_row[["failure_reason"]][[1L]]

        if (!is.na(failure_value) && nzchar(as.character(failure_value))) {
            return(NULL)
        }
    }

    required_columns <- c(
        "gene",
        "variant_id",
        "n_in_group_cell_types",
        "n_comparator_cell_types",
        "n_unique_donors",
        "pooled_out_group_slope",
        "pooled_out_group_slope_ci_low",
        "pooled_out_group_slope_ci_high",
        "in_group_slope",
        "in_group_slope_ci_low",
        "in_group_slope_ci_high",
        "in_group_vs_out_group_estimate",
        "p_value"
    )

    missing_columns <- setdiff(required_columns, names(result_row))

    if (length(missing_columns) > 0L) {
        stop(
            "Missing required columns in `result_row`: ",
            paste(missing_columns, collapse = ", "),
            call. = FALSE
        )
    }

    plot_df <- data.frame(
        group = c("Pooled out-group", "In-group"),
        estimate = c(
            as.numeric(result_row$pooled_out_group_slope[[1L]]),
            as.numeric(result_row$in_group_slope[[1L]])
        ),
        ci_low = c(
            as.numeric(result_row$pooled_out_group_slope_ci_low[[1L]]),
            as.numeric(result_row$in_group_slope_ci_low[[1L]])
        ),
        ci_high = c(
            as.numeric(result_row$pooled_out_group_slope_ci_high[[1L]]),
            as.numeric(result_row$in_group_slope_ci_high[[1L]])
        ),
        stringsAsFactors = FALSE
    )

    plot_df$group <- factor(
        plot_df$group,
        levels = rev(c("Pooled out-group", "In-group"))
    )

    finite_values <- c(
        0,
        plot_df$ci_low,
        plot_df$ci_high,
        plot_df$estimate
    )
    finite_values <- finite_values[is.finite(finite_values)]

    if (!length(finite_values)) {
        stop(
            "No finite estimates or confidence intervals are available to plot.",
            call. = FALSE
        )
    }

    x_limit <- max(abs(finite_values)) * 1.05

    if (!is.finite(x_limit) || x_limit <= 0) {
        x_limit <- 1
    }

    title_text <- paste0(
        as.character(result_row$gene[[1L]]),
        " - ",
        as.character(result_row$variant_id[[1L]])
    )

    subtitle_parts <- c(
        paste0(
            as.integer(result_row$n_in_group_cell_types[[1L]]),
            " in-group cell type",
            if (
                as.integer(result_row$n_in_group_cell_types[[1L]]) == 1L
            ) {
                ""
            } else {
                "s"
            }
        ),
        paste0(
            as.integer(result_row$n_comparator_cell_types[[1L]]),
            " out-group cell type",
            if (
                as.integer(result_row$n_comparator_cell_types[[1L]]) == 1L
            ) {
                ""
            } else {
                "s"
            }
        ),
        paste0(
            as.integer(result_row$n_unique_donors[[1L]]),
            " donors"
        )
    )

    if ("FDR" %in% names(result_row)) {
        fdr_value <- suppressWarnings(
            as.numeric(result_row$FDR[[1L]])
        )

        if (length(fdr_value) == 1L && is.finite(fdr_value)) {
            subtitle_parts <- c(
                paste0(
                    "FDR = ",
                    formatC(
                        fdr_value,
                        format = "e",
                        digits = 3
                    )
                ),
                subtitle_parts
            )
        }
    }

    subtitle_text <- paste(
        subtitle_parts,
        collapse = "; "
    )

    p_value <- suppressWarnings(
        as.numeric(result_row$p_value[[1L]])
    )

    p_text <- if (
        length(p_value) == 1L &&
        is.finite(p_value)
    ) {
        formatC(
            p_value,
            format = "e",
            digits = 3
        )
    } else {
        "NA"
    }

    caption_text <- paste0(
        "In-group - out-group = ",
        formatC(
            as.numeric(
                result_row$in_group_vs_out_group_estimate[[1L]]
            ),
            format = "f",
            digits = 3
        ),
        "; p = ",
        p_text
    )

    ggplot2::ggplot(
        plot_df,
        ggplot2::aes(
            x = estimate,
            y = group,
            color = group
        )
    ) +
        ggplot2::geom_vline(
            xintercept = 0,
            linetype = "dashed",
            linewidth = 0.5
        ) +
        ggplot2::geom_errorbar(
            ggplot2::aes(
                xmin = ci_low,
                xmax = ci_high
            ),
            orientation = "y",
            width = 0.15,
            linewidth = 0.8
        ) +
        ggplot2::geom_point(size = 3) +
        ggplot2::scale_color_manual(
            values = c(
                "Pooled out-group" = "blue",
                "In-group" = "darkgreen"
            ),
            guide = "none"
        ) +
        ggplot2::coord_cartesian(
            xlim = c(-x_limit, x_limit)
        ) +
        ggplot2::labs(
            title = title_text,
            subtitle = subtitle_text,
            caption = caption_text,
            x = "Genotype effect on normalized expression",
            y = NULL
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(
                face = "bold"
            ),
            panel.grid.minor = ggplot2::element_blank()
        )
}

plot_eqtl_interaction_result <- function(
        result_row,
        show_interaction = FALSE
) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.", call. = FALSE)
    }

    if (is.null(result_row) || nrow(result_row) != 1L) {
        stop(
            "`result_row` must contain exactly one row.",
            call. = FALSE
        )
    }

    if (
        length(show_interaction) != 1L ||
        is.na(show_interaction) ||
        !is.logical(show_interaction)
    ) {
        stop(
            "`show_interaction` must be TRUE or FALSE.",
            call. = FALSE
        )
    }

    if ("failure_reason" %in% names(result_row)) {
        failure_value <- result_row[["failure_reason"]][[1L]]

        if (
            !is.na(failure_value) &&
            nzchar(as.character(failure_value))
        ) {
            return(NULL)
        }
    }

    required_columns <- c(
        "gene",
        "variant_id",
        "n_in_group_cell_types",
        "n_comparator_cell_types",
        "n_unique_donors",
        "pooled_out_group_slope",
        "pooled_out_group_slope_ci_low",
        "pooled_out_group_slope_ci_high",
        "in_group_slope",
        "in_group_slope_ci_low",
        "in_group_slope_ci_high",
        "in_group_vs_out_group_estimate",
        "p_value"
    )

    if (show_interaction) {
        required_columns <- c(
            required_columns,
            "in_group_vs_out_group_ci_low",
            "in_group_vs_out_group_ci_high"
        )
    }

    missing_columns <- setdiff(
        required_columns,
        names(result_row)
    )

    if (length(missing_columns) > 0L) {
        stop(
            "Missing required columns in `result_row`: ",
            paste(missing_columns, collapse = ", "),
            call. = FALSE
        )
    }

    plot_df <- data.frame(
        group = c(
            "Pooled out-group",
            "In-group"
        ),
        estimate = c(
            as.numeric(
                result_row$pooled_out_group_slope[[1L]]
            ),
            as.numeric(
                result_row$in_group_slope[[1L]]
            )
        ),
        ci_low = c(
            as.numeric(
                result_row$pooled_out_group_slope_ci_low[[1L]]
            ),
            as.numeric(
                result_row$in_group_slope_ci_low[[1L]]
            )
        ),
        ci_high = c(
            as.numeric(
                result_row$pooled_out_group_slope_ci_high[[1L]]
            ),
            as.numeric(
                result_row$in_group_slope_ci_high[[1L]]
            )
        ),
        stringsAsFactors = FALSE
    )

    if (show_interaction) {
        interaction_row <- data.frame(
            group = "In-group minus out-group",
            estimate = as.numeric(
                result_row$in_group_vs_out_group_estimate[[1L]]
            ),
            ci_low = as.numeric(
                result_row$in_group_vs_out_group_ci_low[[1L]]
            ),
            ci_high = as.numeric(
                result_row$in_group_vs_out_group_ci_high[[1L]]
            ),
            stringsAsFactors = FALSE
        )

        plot_df <- rbind(
            plot_df,
            interaction_row
        )
    }

    group_order <- c(
        "Pooled out-group",
        "In-group"
    )

    if (show_interaction) {
        group_order <- c(
            group_order,
            "In-group minus out-group"
        )
    }

    plot_df$group <- factor(
        plot_df$group,
        levels = rev(group_order)
    )

    finite_values <- c(
        0,
        plot_df$ci_low,
        plot_df$ci_high,
        plot_df$estimate
    )

    finite_values <- finite_values[
        is.finite(finite_values)
    ]

    if (!length(finite_values)) {
        stop(
            paste0(
                "No finite estimates or confidence intervals ",
                "are available to plot."
            ),
            call. = FALSE
        )
    }

    x_limit <- max(
        abs(finite_values)
    ) * 1.05

    if (!is.finite(x_limit) || x_limit <= 0) {
        x_limit <- 1
    }

    title_text <- paste0(
        as.character(result_row$gene[[1L]]),
        " - ",
        as.character(result_row$variant_id[[1L]])
    )

    n_in_group <- as.integer(
        result_row$n_in_group_cell_types[[1L]]
    )

    n_out_group <- as.integer(
        result_row$n_comparator_cell_types[[1L]]
    )

    n_donors <- as.integer(
        result_row$n_unique_donors[[1L]]
    )

    subtitle_parts <- c(
        paste0(
            n_in_group,
            " in-group cell type",
            if (n_in_group == 1L) "" else "s"
        ),
        paste0(
            n_out_group,
            " out-group cell type",
            if (n_out_group == 1L) "" else "s"
        ),
        paste0(
            n_donors,
            " donors"
        )
    )

    if ("FDR" %in% names(result_row)) {
        fdr_value <- suppressWarnings(
            as.numeric(result_row$FDR[[1L]])
        )

        if (
            length(fdr_value) == 1L &&
            is.finite(fdr_value)
        ) {
            subtitle_parts <- c(
                paste0(
                    "FDR = ",
                    formatC(
                        fdr_value,
                        format = "e",
                        digits = 3
                    )
                ),
                subtitle_parts
            )
        }
    }

    subtitle_text <- paste(
        subtitle_parts,
        collapse = "; "
    )

    p_value <- suppressWarnings(
        as.numeric(result_row$p_value[[1L]])
    )

    p_text <- if (
        length(p_value) == 1L &&
        is.finite(p_value)
    ) {
        formatC(
            p_value,
            format = "e",
            digits = 3
        )
    } else {
        "NA"
    }

    interaction_estimate <- suppressWarnings(
        as.numeric(
            result_row$in_group_vs_out_group_estimate[[1L]]
        )
    )

    interaction_text <- if (
        length(interaction_estimate) == 1L &&
        is.finite(interaction_estimate)
    ) {
        formatC(
            interaction_estimate,
            format = "f",
            digits = 3
        )
    } else {
        "NA"
    }

    caption_text <- paste0(
        "In-group - out-group = ",
        interaction_text,
        "; p = ",
        p_text
    )

    color_values <- c(
        "Pooled out-group" = "blue",
        "In-group" = "darkgreen",
        "In-group minus out-group" = "purple4"
    )

    x_axis_label <- if (show_interaction) {
        "Effect estimate on normalized-expression scale"
    } else {
        "Genotype effect on normalized expression"
    }

    ggplot2::ggplot(
        plot_df,
        ggplot2::aes(
            x = estimate,
            y = group,
            color = group
        )
    ) +
        ggplot2::geom_vline(
            xintercept = 0,
            linetype = "dashed",
            linewidth = 0.5
        ) +
        ggplot2::geom_errorbar(
            ggplot2::aes(
                xmin = ci_low,
                xmax = ci_high
            ),
            orientation = "y",
            width = 0.15,
            linewidth = 0.8
        ) +
        ggplot2::geom_point(
            size = 3
        ) +
        ggplot2::scale_color_manual(
            values = color_values,
            guide = "none"
        ) +
        ggplot2::coord_cartesian(
            xlim = c(
                -x_limit,
                x_limit
            )
        ) +
        ggplot2::labs(
            title = title_text,
            subtitle = subtitle_text,
            caption = caption_text,
            x = x_axis_label,
            y = NULL
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(
                face = "bold"
            ),
            panel.grid.minor = ggplot2::element_blank()
        )
}

plot_gene_tpm_by_group <- function(
        geneName,
        eqtl_data,
        show_title = FALSE
) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.", call. = FALSE)
    }

    if (!requireNamespace("ggrepel", quietly = TRUE)) {
        stop("Package 'ggrepel' is required.", call. = FALSE)
    }

    required_columns <- c(
        "gene",
        "cell_type",
        "region",
        "cell_type_region",
        "in_group",
        "expression_tpm"
    )

    missing_columns <- setdiff(required_columns, names(eqtl_data))

    if (length(missing_columns)) {
        stop(
            "Missing required columns: ",
            paste(missing_columns, collapse = ", "),
            call. = FALSE
        )
    }

    eqtl_data <- data.table::as.data.table(eqtl_data)

    gene_data <- eqtl_data[
        gene == geneName & !is.na(expression_tpm),
        .SD,
        .SDcols = required_columns
    ]

    if (!nrow(gene_data)) {
        return(NULL)
    }

    gene_data[, expression_tpm := as.numeric(expression_tpm)]

    group_consistency <- gene_data[
        ,
        .(n_in_group_values = data.table::uniqueN(in_group, na.rm = TRUE)),
        by = cell_type_region
    ]

    if (any(group_consistency$n_in_group_values != 1L)) {
        stop(
            paste0(
                "At least one cell_type_region has inconsistent ",
                "in_group values."
            ),
            call. = FALSE
        )
    }

    summary_data <- gene_data[
        ,
        .(
            median_tpm = stats::median(
                expression_tpm,
                na.rm = TRUE
            )
        ),
        by = .(
            cell_type,
            region,
            cell_type_region,
            in_group
        )
    ]

    summary_data[, log10_median_tpm := log10(median_tpm + 1)]

    summary_data[
        ,
        group := ifelse(
            as.logical(in_group),
            "In-group",
            "Out-group"
        )
    ]

    summary_data[
        ,
        group := factor(
            group,
            levels = c("In-group", "Out-group")
        )
    ]

    regions <- unique(summary_data$region)
    regions <- regions[!is.na(regions)]

    if (length(regions) == 1L) {
        summary_data[, point_label := cell_type]
    } else {
        summary_data[, point_label := cell_type_region]
    }

    # Match the effect plot: out-group above in-group.
    summary_data[
        ,
        group_y := ifelse(
            group == "In-group",
            1,
            2
        )
    ]

    # Deterministic vertical offsets prevent points from overlapping while
    # keeping them within their in-group or out-group boxplot.
    summary_data[, point_y := group_y]

    for (group_name in levels(summary_data$group)) {
        index <- which(summary_data$group == group_name)

        if (!length(index)) {
            next
        }

        offsets <- if (length(index) == 1L) {
            0
        } else {
            seq(
                from = -0.35,
                to = 0.35,
                length.out = length(index)
            )
        }

        ordered_index <- index[
            order(summary_data$log10_median_tpm[index])
        ]

        summary_data$point_y[ordered_index] <-
            summary_data$group_y[ordered_index] + offsets
    }

    x_range <- range(
        summary_data$log10_median_tpm,
        finite = TRUE
    )

    x_width <- diff(x_range)

    if (!is.finite(x_width) || x_width == 0) {
        x_width <- max(abs(x_range), 1)
    }

    label_nudge <- 0.025 * x_width

    title_text <- if (show_title) {
        paste0(
            geneName,
            " - median TPM by cell type-region"
        )
    } else {
        NULL
    }

    ggplot2::ggplot() +
        ggplot2::geom_boxplot(
            data = summary_data,
            ggplot2::aes(
                x = log10_median_tpm,
                y = group_y,
                group = group,
                fill = group
            ),
            orientation = "y",
            width = 0.35,
            alpha = 0.25,
            outlier.shape = NA,
            linewidth = 0.7
        ) +
        ggplot2::geom_point(
            data = summary_data,
            ggplot2::aes(
                x = log10_median_tpm,
                y = point_y,
                color = cell_type_region
            ),
            size = 3
        ) +
        ggrepel::geom_text_repel(
            data = summary_data,
            ggplot2::aes(
                x = log10_median_tpm + label_nudge,
                y = point_y,
                label = point_label,
                color = cell_type_region
            ),
            direction = "y",
            hjust = 0,
            size = 3.5,
            max.overlaps = Inf,
            min.segment.length = 0,
            segment.alpha = 0.5,
            box.padding = 0.25,
            point.padding = 0.15,
            show.legend = FALSE
        ) +
        ggplot2::scale_fill_manual(
            values = c(
                "Out-group" = "blue",
                "In-group" = "darkgreen"
            ),
            guide = "none"
        ) +
        ggplot2::scale_y_continuous(
            breaks = c(1, 2),
            labels = c(
                "In-group",
                "Out-group"
            ),
            limits = c(0.65, 2.35)
        ) +
        ggplot2::scale_x_continuous(
            expand = ggplot2::expansion(
                mult = c(0.05, 0.25)
            )
        ) +
        ggplot2::labs(
            title = title_text,
            x = "Median log10(TPM + 1)",
            y = NULL
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(
                face = "bold"
            ),
            panel.grid.minor = ggplot2::element_blank(),
            legend.position = "none"
        )
}

plot_fdr_fraction_by_outgroup_count <- function(
        results,
        fdr_cutoff = 0.05
) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.", call. = FALSE)
    }

    if (!requireNamespace("scales", quietly = TRUE)) {
        stop("Package 'scales' is required.", call. = FALSE)
    }

    required_columns <- c(
        "gene",
        "n_comparator_cell_types",
        "FDR"
    )

    missing_columns <- setdiff(
        required_columns,
        names(results)
    )

    if (length(missing_columns)) {
        stop(
            "Missing required columns: ",
            paste(missing_columns, collapse = ", "),
            call. = FALSE
        )
    }

    significant_label <- paste0(
        "FDR < ",
        fdr_cutoff
    )

    nonsignificant_label <- paste0(
        "FDR >= ",
        fdr_cutoff
    )

    aggregate_significant_label <- paste0(
        "All testable: ",
        significant_label
    )

    aggregate_nonsignificant_label <- paste0(
        "All testable: ",
        nonsignificant_label
    )

    aggregate_x_label <- "All testable\n(>=1 out-group)"

    d <- as.data.frame(results)

    # Use one row per gene-variant pair when variant_id is present.
    id_columns <- intersect(
        c("gene", "variant_id"),
        names(d)
    )

    d <- d[
        !duplicated(d[, id_columns, drop = FALSE]),
        ,
        drop = FALSE
    ]

    d$n_comparator_cell_types <- as.integer(
        d$n_comparator_cell_types
    )

    d$FDR <- as.numeric(d$FDR)

    # Genes with no out-group expression cannot be tested.
    private_expression <- d[
        !is.na(d$n_comparator_cell_types) &
            d$n_comparator_cell_types == 0L,
        ,
        drop = FALSE
    ]

    # For categories with at least one out-group cell type, rows with
    # missing FDR are excluded from both numerator and denominator.
    tested <- d[
        !is.na(d$n_comparator_cell_types) &
            d$n_comparator_cell_types > 0L &
            !is.na(d$FDR),
        ,
        drop = FALSE
    ]

    tested$status <- ifelse(
        tested$FDR < fdr_cutoff,
        significant_label,
        nonsignificant_label
    )

    plot_data <- data.frame(
        n_comparator_cell_types = integer(),
        x_label = character(),
        status = character(),
        fill_group = character(),
        n_genes = integer(),
        total_genes = integer(),
        fraction = numeric(),
        stringsAsFactors = FALSE
    )

    if (nrow(tested)) {
        tested_counts <- stats::aggregate(
            gene ~ n_comparator_cell_types + status,
            data = tested,
            FUN = length
        )

        names(tested_counts)[
            names(tested_counts) == "gene"
        ] <- "n_genes"

        tested_totals <- stats::aggregate(
            n_genes ~ n_comparator_cell_types,
            data = tested_counts,
            FUN = sum
        )

        names(tested_totals)[
            names(tested_totals) == "n_genes"
        ] <- "total_genes"

        tested_counts <- merge(
            tested_counts,
            tested_totals,
            by = "n_comparator_cell_types",
            all.x = TRUE,
            sort = FALSE
        )

        tested_counts$fraction <- (
            tested_counts$n_genes /
                tested_counts$total_genes
        )

        tested_counts$x_label <- as.character(
            tested_counts$n_comparator_cell_types
        )

        tested_counts$fill_group <- tested_counts$status

        plot_data <- tested_counts[
            ,
            c(
                "n_comparator_cell_types",
                "x_label",
                "status",
                "fill_group",
                "n_genes",
                "total_genes",
                "fraction"
            )
        ]

        # Aggregate all testable genes across comparator-count categories.
        aggregate_counts <- stats::aggregate(
            gene ~ status,
            data = tested,
            FUN = length
        )

        names(aggregate_counts)[
            names(aggregate_counts) == "gene"
        ] <- "n_genes"

        aggregate_total <- sum(
            aggregate_counts$n_genes
        )

        aggregate_counts$n_comparator_cell_types <- NA_integer_
        aggregate_counts$x_label <- aggregate_x_label
        aggregate_counts$total_genes <- aggregate_total
        aggregate_counts$fraction <- (
            aggregate_counts$n_genes /
                aggregate_total
        )

        aggregate_counts$fill_group <- ifelse(
            aggregate_counts$status == significant_label,
            aggregate_significant_label,
            aggregate_nonsignificant_label
        )

        aggregate_counts <- aggregate_counts[
            ,
            c(
                "n_comparator_cell_types",
                "x_label",
                "status",
                "fill_group",
                "n_genes",
                "total_genes",
                "fraction"
            )
        ]

        plot_data <- rbind(
            plot_data,
            aggregate_counts
        )
    }

    if (nrow(private_expression)) {
        private_row <- data.frame(
            n_comparator_cell_types = 0L,
            x_label = "Private expression\n(0 out-groups)",
            status = "Not testable",
            fill_group = "Not testable",
            n_genes = nrow(private_expression),
            total_genes = nrow(private_expression),
            fraction = 1,
            stringsAsFactors = FALSE
        )

        plot_data <- rbind(
            private_row,
            plot_data
        )
    }

    if (!nrow(plot_data)) {
        return(NULL)
    }

    outgroup_values <- sort(
        unique(
            plot_data$n_comparator_cell_types[
                !is.na(plot_data$n_comparator_cell_types) &
                    plot_data$n_comparator_cell_types > 0L
            ]
        )
    )

    x_levels <- c(
        if (any(plot_data$n_comparator_cell_types == 0L, na.rm = TRUE)) {
            "Private expression\n(0 out-groups)"
        },
        as.character(outgroup_values),
        if (any(plot_data$x_label == aggregate_x_label)) {
            aggregate_x_label
        }
    )

    plot_data$x_label <- factor(
        plot_data$x_label,
        levels = x_levels
    )

    # ggplot stacks factor levels in reverse order. Significant categories
    # are listed after their nonsignificant counterparts so they appear at
    # the bottom of each bar.
    stack_levels <- c(
        "Not testable",
        nonsignificant_label,
        significant_label,
        aggregate_nonsignificant_label,
        aggregate_significant_label
    )

    plot_data$fill_group <- factor(
        plot_data$fill_group,
        levels = stack_levels
    )

    total_labels <- unique(
        plot_data[
            ,
            c(
                "x_label",
                "total_genes"
            )
        ]
    )

    total_labels$label <- paste0(
        "n = ",
        total_labels$total_genes
    )

    fill_values <- stats::setNames(
        c(
            "darkgreen",
            "grey75",
            "grey45",
            "grey75",
            "blue"
        ),
        c(
            significant_label,
            nonsignificant_label,
            "Not testable",
            aggregate_nonsignificant_label,
            aggregate_significant_label
        )
    )

    legend_order <- c(
        significant_label,
        nonsignificant_label,
        "Not testable",
        aggregate_significant_label,
        aggregate_nonsignificant_label
    )

    ggplot2::ggplot(
        plot_data,
        ggplot2::aes(
            x = x_label,
            y = fraction,
            fill = fill_group
        )
    ) +
        ggplot2::geom_col(
            width = 0.75
        ) +
        ggplot2::geom_text(
            ggplot2::aes(
                label = n_genes
            ),
            position = ggplot2::position_stack(
                vjust = 0.5
            ),
            size = 3.5
        ) +
        ggplot2::geom_text(
            data = total_labels,
            ggplot2::aes(
                x = x_label,
                y = 1.04,
                label = label
            ),
            inherit.aes = FALSE,
            size = 3.5
        ) +
        ggplot2::scale_y_continuous(
            labels = scales::label_percent(),
            breaks = seq(
                0,
                1,
                by = 0.2
            ),
            expand = ggplot2::expansion(
                mult = c(0, 0.08)
            )
        ) +
        ggplot2::scale_fill_manual(
            values = fill_values,
            breaks = legend_order,
            labels = legend_order,
            drop = FALSE
        ) +
        ggplot2::coord_cartesian(
            ylim = c(0, 1.09),
            clip = "off"
        ) +
        ggplot2::labs(
            x = "Out-group cell-type category",
            y = "Fraction of genes",
            fill = NULL
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major.x = ggplot2::element_blank(),
            legend.position = "top",
            axis.text.x = ggplot2::element_text(
                angle = 0,
                hjust = 0.5
            )
        )
}

plot_in_vs_out_group_slopes <- function(
        results,
        fdr_cutoff = 0.05,
        show_title = FALSE,
        label_significant = FALSE
) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.", call. = FALSE)
    }

    if (
        label_significant &&
        !requireNamespace("ggrepel", quietly = TRUE)
    ) {
        stop(
            "Package 'ggrepel' is required when label_significant = TRUE.",
            call. = FALSE
        )
    }

    required_columns <- c(
        "gene",
        "pooled_out_group_slope",
        "in_group_slope",
        "FDR",
        "n_comparator_cell_types"
    )

    missing_columns <- setdiff(
        required_columns,
        names(results)
    )

    if (length(missing_columns)) {
        stop(
            "Missing required columns: ",
            paste(missing_columns, collapse = ", "),
            call. = FALSE
        )
    }

    d <- as.data.frame(results)

    # Keep one row per gene-variant pair when variant_id is present.
    id_columns <- intersect(
        c("gene", "variant_id"),
        names(d)
    )

    if (length(id_columns)) {
        d <- d[
            !duplicated(d[, id_columns, drop = FALSE]),
            ,
            drop = FALSE
        ]
    }

    d$n_comparator_cell_types <- as.integer(
        d$n_comparator_cell_types
    )

    d$pooled_out_group_slope <- as.numeric(
        d$pooled_out_group_slope
    )

    d$in_group_slope <- as.numeric(
        d$in_group_slope
    )

    d$FDR <- as.numeric(d$FDR)

    # Plot genes that had at least one comparator and finite slope estimates.
    # Rows with missing FDR remain visible as "FDR unavailable".
    d <- d[
        !is.na(d$n_comparator_cell_types) &
            d$n_comparator_cell_types > 0L &
            is.finite(d$pooled_out_group_slope) &
            is.finite(d$in_group_slope),
        ,
        drop = FALSE
    ]

    if (!nrow(d)) {
        return(NULL)
    }

    significant_label <- paste0(
        "FDR < ",
        fdr_cutoff
    )

    nonsignificant_label <- paste0(
        "FDR >= ",
        fdr_cutoff
    )

    unavailable_label <- "FDR unavailable"

    d$significance <- ifelse(
        is.na(d$FDR),
        unavailable_label,
        ifelse(
            d$FDR < fdr_cutoff,
            significant_label,
            nonsignificant_label
        )
    )

    significance_levels <- c(
        significant_label,
        nonsignificant_label,
        unavailable_label
    )

    d$significance <- factor(
        d$significance,
        levels = significance_levels
    )

    if ("variant_id" %in% names(d)) {
        d$point_label <- paste0(
            d$gene,
            "\n",
            d$variant_id
        )
    } else {
        d$point_label <- d$gene
    }

    axis_values <- c(
        d$pooled_out_group_slope,
        d$in_group_slope,
        0
    )

    max_abs <- max(
        abs(axis_values),
        na.rm = TRUE
    )

    if (!is.finite(max_abs) || max_abs == 0) {
        max_abs <- 1
    }

    axis_limit <- max_abs * 1.05

    title_text <- if (show_title) {
        "In-group versus pooled out-group genotype effects"
    } else {
        NULL
    }

    color_values <- stats::setNames(
        c(
            "darkgreen",
            "grey60",
            "grey85"
        ),
        significance_levels
    )

    p <- ggplot2::ggplot(
        d,
        ggplot2::aes(
            x = pooled_out_group_slope,
            y = in_group_slope,
            color = significance
        )
    ) +
        ggplot2::geom_hline(
            yintercept = 0,
            linetype = "dashed",
            color = "grey50",
            linewidth = 0.6
        ) +
        ggplot2::geom_vline(
            xintercept = 0,
            linetype = "dashed",
            color = "grey50",
            linewidth = 0.6
        ) +
        ggplot2::geom_abline(
            intercept = 0,
            slope = 1,
            color = "black",
            linewidth = 0.7,
            alpha = 0.8
        ) +
        ggplot2::geom_point(
            size = 2.5,
            alpha = 0.85
        ) +
        ggplot2::scale_color_manual(
            values = color_values,
            breaks = significance_levels,
            drop = FALSE
        ) +
        ggplot2::coord_fixed(
            xlim = c(-axis_limit, axis_limit),
            ylim = c(-axis_limit, axis_limit)
        ) +
        ggplot2::labs(
            title = title_text,
            x = "Pooled out-group slope",
            y = "In-group slope",
            color = NULL
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            legend.position = "top",
            plot.title = ggplot2::element_text(
                face = "bold"
            )
        )

    if (label_significant) {
        label_data <- d[
            !is.na(d$FDR) &
                d$FDR < fdr_cutoff,
            ,
            drop = FALSE
        ]

        if (nrow(label_data)) {
            p <- p +
                ggrepel::geom_text_repel(
                    data = label_data,
                    ggplot2::aes(
                        label = point_label
                    ),
                    size = 3,
                    max.overlaps = Inf,
                    show.legend = FALSE
                )
        }
    }

    p
}

write_eqtl_summary_pdf <- function(
        results,
        eqtl_data,
        show_interaction=FALSE,
        outfile,
        width = 11,
        height = 8.5
) {
    if (!requireNamespace("cowplot", quietly = TRUE)) {
        stop("Package 'cowplot' is required.", call. = FALSE)
    }

    if (!dir.exists(dirname(outfile))) {
        stop(
            "Output directory does not exist: ",
            dirname(outfile),
            call. = FALSE
        )
    }

    results <- data.table::as.data.table(results)

    grDevices::pdf(
        file = outfile,
        width = width,
        height = height,
        onefile = TRUE
    )

    on.exit(
        grDevices::dev.off(),
        add = TRUE
    )

    # Summary page 1: fraction significant by number of out-group cell types.
    p_fdr <- plot_fdr_fraction_by_outgroup_count(
        results = results
    )

    if (!is.null(p_fdr)) {
        print(p_fdr)
    }

    # Summary page 2: in-group versus pooled out-group slopes.
    p_slopes <- plot_in_vs_out_group_slopes(
        results = results
    )

    if (!is.null(p_slopes)) {
        print(p_slopes)
    }

    # One page per successfully analyzed gene-variant pair.
    for (i in seq_len(nrow(results))) {
        result_row <- results[i]

        p_effect <- plot_eqtl_interaction_result(
            result_row = result_row,
            show_interaction=show_interaction
        )

        # Failed or otherwise unplottable model results return NULL.
        if (is.null(p_effect)) {
            next
        }

        p_tpm <- plot_gene_tpm_by_group(
            geneName = result_row$gene[[1L]],
            eqtl_data = eqtl_data,
            show_title = FALSE
        )

        if (is.null(p_tpm)) {
            print(p_effect)
            next
        }

        combined <- cowplot::plot_grid(
            p_effect,
            p_tpm,
            ncol = 1,
            align = "v",
            axis = "lr",
            rel_heights = c(1, 1)
        )

        print(combined)
    }

    invisible(outfile)
}
