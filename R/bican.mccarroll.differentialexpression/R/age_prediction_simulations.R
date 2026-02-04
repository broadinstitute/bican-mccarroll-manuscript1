# Framework:
# 1) Simulate data (separate functions)
# 2) Fit/plot (reusable analysis function that produces the same set of plots)

simulate_uniform_ages <- function(n = 200, age_min = 3, age_max = 10, seed = 1) {
    set.seed(seed)
    runif(n, min = age_min, max = age_max)
}

# Old simulation: "very bad" predictors (median + noise)
simulate_predictions_bad <- function(age, noise_sd = 1.0) {
    med_age <- median(age)
    pred1 <- med_age + rnorm(length(age), mean = 0, sd = noise_sd)
    pred2 <- med_age + rnorm(length(age), mean = 0, sd = noise_sd)
    list(pred1 = pred1, pred2 = pred2)
}

# New simulation: more realistic predictors:
# - predictions track age around the center with a slope (signal) and noise
# - tails get biased toward the center ("regression to the mean"), but less aggressively than a pure shrinkage model
#
# bias_strength controls how much the tails get pulled toward the median (0 = none, 1 = full pull within the tail region)
# tail_width controls what counts as "tail": e.g., 0.15 means outer 15% on each side
simulate_predictions_tail_bias <- function(age,
                                           slope = 0.7,
                                           noise_sd = 1.0,
                                           tail_width = 0.15,
                                           bias_strength = 0.4) {
    n <- length(age)
    med_age <- median(age)

    # baseline signal: centered linear relationship + noise
    pred1 <- med_age + slope * (age - med_age) + rnorm(n, mean = 0, sd = noise_sd)
    pred2 <- med_age + slope * (age - med_age) + rnorm(n, mean = 0, sd = noise_sd)

    # apply extra tail bias: pull predictions toward the median for donors in the extremes
    q_lo <- as.numeric(stats::quantile(age, probs = tail_width))
    q_hi <- as.numeric(stats::quantile(age, probs = 1 - tail_width))
    is_tail <- (age <= q_lo) | (age >= q_hi)

    # In the tails, move pred partway toward the median (less aggressive regression-to-mean than setting slope tiny)
    pred1[is_tail] <- pred1[is_tail] + bias_strength * (med_age - pred1[is_tail])
    pred2[is_tail] <- pred2[is_tail] + bias_strength * (med_age - pred2[is_tail])

    list(pred1 = pred1, pred2 = pred2)
}

# Reusable analysis: produces the same plot set as before + returns correlations and fitted curves.
analyze_predictions_with_gam_correction <- function(age,
                                                    pred1,
                                                    pred2,
                                                    gam_k = 5) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install ggplot2")
    if (!requireNamespace("mgcv", quietly = TRUE)) stop("Please install mgcv")

    # Raw residuals (pred - age)
    resid1 <- pred1 - age
    resid2 <- pred2 - age
    cor_raw <- cor(resid1, resid2)

    df1 <- data.frame(age = age, pred = pred1)
    df2 <- data.frame(age = age, pred = pred2)
    dfr <- data.frame(resid1 = resid1, resid2 = resid2)

    # Make R CMD CHECK Happy
    intercept <- slope <- NULL

    # Bad/realistic predictor plots with LM line in red (as before)
    p1 <- ggplot2::ggplot(df1, ggplot2::aes(x = age, y = pred)) +
        ggplot2::geom_point(alpha = 0.6) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        ggplot2::geom_smooth(method = "lm", se = FALSE, color = "red") +
        ggplot2::labs(
            title = "Predictor 1: age vs prediction",
            subtitle = "Red line: linear fit to data",
            x = "Chronological age",
            y = "Predicted age"
        ) +
        ggplot2::theme_classic(base_size = 16)

    p2 <- ggplot2::ggplot(df2, ggplot2::aes(x = age, y = pred)) +
        ggplot2::geom_point(alpha = 0.6) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        ggplot2::geom_smooth(method = "lm", se = FALSE, color = "red") +
        ggplot2::labs(
            title = "Predictor 2: age vs prediction",
            subtitle = "Red line: linear fit to data",
            x = "Chronological age",
            y = "Predicted age"
        ) +
        ggplot2::theme_classic(base_size = 16)

    p3 <- ggplot2::ggplot(dfr, ggplot2::aes(x = resid1, y = resid2)) +
        ggplot2::geom_point(alpha = 0.6) +
        ggplot2::geom_smooth(method = "lm", se = FALSE) +
        ggplot2::annotate(
            "text",
            x = -Inf, y = Inf,
            hjust = -0.1, vjust = 1.2,
            label = sprintf("cor = %.3f", cor_raw),
            size = 4.5
        ) +
        ggplot2::labs(
            title = "Raw residual correlation",
            x = "Residual 1 (pred1 - age)",
            y = "Residual 2 (pred2 - age)"
        ) +
        ggplot2::theme_classic(base_size = 16)

    # GAM correction: fit pred ~ s(age) and take distance from the curve
    df_gam <- data.frame(age = age, pred1 = pred1, pred2 = pred2)

    fit1 <- mgcv::gam(pred1 ~ s(age, k = gam_k, bs = "cr"), data = df_gam, method = "REML")
    fit2 <- mgcv::gam(pred2 ~ s(age, k = gam_k, bs = "cr"), data = df_gam, method = "REML")

    gam_fit1 <- stats::predict(fit1, newdata = data.frame(age = age))
    gam_fit2 <- stats::predict(fit2, newdata = data.frame(age = age))

    corr_resid1 <- pred1 - gam_fit1
    corr_resid2 <- pred2 - gam_fit2
    cor_corr <- cor(corr_resid1, corr_resid2)

    df1$gam_fit <- gam_fit1
    df2$gam_fit <- gam_fit2

    # MAKE R CMD CHECK Happy
    gam_fit <- pred <- NULL

    # Show GAM curve (red)
    p4 <- ggplot2::ggplot(df1, ggplot2::aes(x = age, y = pred)) +
        ggplot2::geom_point(alpha = 0.6) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        ggplot2::geom_line(ggplot2::aes(y = gam_fit), linewidth = 1, color = "red") +
        ggplot2::labs(
            title = "GAM fit through age vs predicted age (predictor 1)",
            x = "Chronological age",
            y = "Predicted age"
        ) +
        ggplot2::theme_classic(base_size = 16)

    p5 <- ggplot2::ggplot(df2, ggplot2::aes(x = age, y = pred)) +
        ggplot2::geom_point(alpha = 0.6) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        ggplot2::geom_line(ggplot2::aes(y = gam_fit), linewidth = 1, color = "red") +
        ggplot2::labs(
            title = "GAM fit through age vs predicted age (predictor 2)",
            x = "Chronological age",
            y = "Predicted age"
        ) +
        ggplot2::theme_classic(base_size = 16)

    df_corr <- data.frame(resid1 = corr_resid1, resid2 = corr_resid2)

    p6 <- ggplot2::ggplot(df_corr, ggplot2::aes(x = resid1, y = resid2)) +
        ggplot2::geom_point(alpha = 0.6) +
        ggplot2::geom_smooth(method = "lm", se = FALSE) +
        ggplot2::annotate(
            "text",
            x = -Inf, y = Inf,
            hjust = -0.1, vjust = 1.2,
            label = sprintf("cor = %.3f", cor_corr),
            size = 4.5
        ) +
        ggplot2::labs(
            title = "Residual correlation after GAM correction",
            x = "Corrected residual 1 (pred1 - gam_fit1)",
            y = "Corrected residual 2 (pred2 - gam_fit2)"
        ) +
        ggplot2::theme_classic(base_size = 16)

    print(p1); print(p2); print(p3)
    print(p4); print(p5); print(p6)

    invisible(list(
        cor_raw = cor_raw,
        cor_corrected = cor_corr,
        age = age,
        pred1 = pred1,
        pred2 = pred2,
        resid1 = resid1,
        resid2 = resid2,
        gam_fit1 = gam_fit1,
        gam_fit2 = gam_fit2,
        corr_resid1 = corr_resid1,
        corr_resid2 = corr_resid2,
        gam_model1 = fit1,
        gam_model2 = fit2
    ))
}

# Convenience wrappers

run_simulation_bad <- function(n = 200, age_min = 3, age_max = 10, noise_sd = 0.5, seed = 1, gam_k = 5) {
    age <- simulate_uniform_ages(n = n, age_min = age_min, age_max = age_max, seed = seed)
    preds <- simulate_predictions_bad(age = age, noise_sd = noise_sd)
    analyze_predictions_with_gam_correction(age = age, pred1 = preds$pred1, pred2 = preds$pred2, gam_k = gam_k)
}

run_simulation_tail_bias <- function(n = 200,
                                     age_min = 3,
                                     age_max = 10,
                                     slope = 0.7,
                                     noise_sd = 0.5,
                                     tail_width = 0.15,
                                     bias_strength = 0.4,
                                     seed = 1,
                                     gam_k = 5)
{
    age <- simulate_uniform_ages(n = n, age_min = age_min, age_max = age_max, seed = seed)
    preds <- simulate_predictions_tail_bias(age = age,
                                            slope = slope,
                                            noise_sd = noise_sd,
                                            tail_width = tail_width,
                                            bias_strength = bias_strength)
    analyze_predictions_with_gam_correction(age = age, pred1 = preds$pred1, pred2 = preds$pred2, gam_k = gam_k)
}

# Example runs:
# run_simulation_bad(seed = 1)
# run_simulation_tail_bias(seed = 1, slope = 0.8, noise_sd = 0.8, tail_width = 0.15, bias_strength = 0.4)
# run_simulation_tail_bias(seed = 1, slope = 0.8, noise_sd = 0.4, tail_width = 0.15, bias_strength = 0.3)
