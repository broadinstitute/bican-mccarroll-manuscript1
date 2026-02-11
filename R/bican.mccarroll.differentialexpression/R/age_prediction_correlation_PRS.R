# source ("R/age_prediction_correlation_cell_type_proportion.R")
# age_preds_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/age_prediction/age_prediction_region_alpha_0/age_prediction_results_alpha0_donor_predictions.txt"
# prs_file="/broad/mccarroll/giulio/gwas/bican_um1_wgs_gvs/score/bican_um1_wgs_gvs.scores.tsv"
# metadata_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/metacells/LEVEL_3/donor_rxn_DGEList_samples.tsv.gz"
# covar_cols = c("imputed_sex", "PC1", "PC2", "PC3", "PC4", "PC5", "apoe_score")
# min_donors=100
# cellTypeFile="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/metadata/cell_types_for_prs_test.txt"
# n_cores=4; #This can really slam memory usage if set too high.
# seed=1234
# outFile="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/age_prediction/age_prediction_region_alpha_0/age_prediction_PRS_correlation_results.txt"
# outPDF="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/age_prediction/age_prediction_region_alpha_0/age_prediction_PRS_correlation_results.pdf"

correlate_age_prediction_with_PRS<-function (age_preds_file, prs_file, metadata_file, cellTypeFile, covar_cols, min_donors=100, seed=1234) {
    set.seed(seed)

    data_list=prepare_data_for_PRS_regression(age_preds_file, prs_file, metadata_file, cellTypeFile)
    df=data_list$data
    prs_categories=data_list$prs_categories

    #quick test of the MSN result
    df2=df[df$cell_type=="MSN_D1" & df$region=="CaH", ]
    results_test=run_prs_assoc_all_celltype_region(df2, prs_categories=c("TS_2019"), outcome_col = "resid_mean_corrected", covar_cols = covar_cols, min_n = min_donors, n_permutation = 1e6, n_core = 1)

    results <- time_run_log(quote(
        run_prs_assoc_all_celltype_region(df, prs_categories=prs_categories, outcome_col = "resid_mean_corrected", covar_cols = covar_cols, min_n = min_donors, n_permutation = 1e5, n_core = n_cores)),
        label = paste("permutation run ", n_cores, " core")
    )

    #Add the FDR by PRS score
    results <- add_bh_by_prs(results, p_col="p_perm", out_col="p_fdr")

    write.table(results, file=outFile, sep="\t", quote=FALSE, row.names = FALSE)

    pdf(outPDF, width=11, height=11)
    #Do the permutations buy me anything? Not really!
    maxXY=(max(c(-log10(results$p), -log10(results$p_perm)), na.rm = TRUE))*1.1
    plot (-log10(results$p), -log10(results$p_perm), xlab="empiric p-value", ylab="permuted p", xlim=c(0,maxXY), ylim=c(0,maxXY), main="PRS correlation with age residuals: empirical vs permuted p-values")
    abline(0,1, col='red')


    #Look at TS_2019
    z=results[results$prs=="TS_2019",]
    z=z[order(z$p_fdr),]
    plot_fdr_heatmap(z, value_col="beta", use_neg_log10 = F)
    plot_fdr_heatmap(z, value_col="p_fdr", use_neg_log10 = T)
    plot_fdr_heatmap2(z, value_col = "beta", use_neg_log10 = FALSE, title = "TS_2019",
                      label_mode = "star", fdr_col = "p_fdr", fdr_thresh = 0.05, label_size = 8)
    plot_fdr_heatmap2(z, value_col = "beta", use_neg_log10 = FALSE, title = "TS_2019",
                      label_mode = "fdr", fdr_col = "p_fdr", fdr_thresh = 0.05, label_size=5)

    #Look at BIP2021
    z=results[results$prs=="BIP_2021",]
    plot_fdr_heatmap2(z, value_col = "beta", use_neg_log10 = FALSE, title = "BIP_2021",
                      label_mode = "fdr", fdr_col = "p_fdr", fdr_thresh = 0.05, label_size=5)


    dev.off()

    results_pass_fdr=results[results$p_fdr<0.05,]
    results[which.min(results$p_fdr_prs),]
    hist (-log10(results$p_fdr_prs), main="PRS correlation with residuals FDR p-values", xlab="-log10(FDR p-value)")
    abline(v=-log10(0.05), col="red", lty=2)
    hist (results$p_adj, main="PRS correlation with residuals FDR p-values", xlab="FDR p-value")
    abline(v=(0.05), col="red", lty=2)




    #predicted age z-score
    results=run_prs_assoc_all_celltype_region(df, prs_categories, outcome_col = "pred_mean_corrected_z", covar_cols = c(covar_cols, "age"))
    results$p_adj=p.adjust(results$p, method = "fdr")
    results_pass_fdr=results[results$p_adj<0.05,]

    #chronological age
    results_age=run_prs_assoc_all_celltype_region(df, prs_categories, outcome_col = "age", covar_cols = covar_cols)
    results_age$p_adj=p.adjust(results_age$p, method = "fdr")
    results_age[which.min(results_age$p_adj),]
    #the best FDR pvalues were 0.03-0.05. Let's be a little lenient and look at all results with FDR < 0.1
    #I did not write the 2nd sentence above, copilot filled that in.
    results_age_pass_fdr=results_age[results_age$p_adj<0.05,]
    hist (-log10(results_age$p_adj), main="PRS correlation with chronological age FDR p-values", xlab="-log10(FDR p-value)")

    sex_result=run_sex_assoc_all_celltype_region(df, outcome_col = "resid_mean_corrected",
                                                 sex_col = "imputed_sex",
                                                 covar_cols = c("PC1", "PC2", "PC3", "PC4", "PC5", "apoe_score", "age"),
                                                 min_n = 50)
    sex_result$p_adj=p.adjust(sex_result$p, method = "fdr")
    sex_result$p_adj_fmt <- formatC(sex_result$p_adj, format = "e", digits = 2)

    sex_result_pass_fdr=sex_result[sex_result$p_adj<0.05, ]

    p <- plot_sex_fdr_heatmap(sex_result, value_col = "p_adj", use_neg_log10 = TRUE)

    index=1
    plot_sex_outcome_celltype_region(df, cell_type = sex_result_pass_fdr[index,]$cell_type,
                                     region = sex_result_pass_fdr[index,]$region,
                                     outcome_col = "resid_mean_corrected", sex_col= "imputed_sex")

    plot_sex_fdr_heatmap(sex_result, value_col = "p_adj", use_neg_log10 = TRUE)
    plot_sex_fdr_heatmap(sex_result, value_col = "beta", use_neg_log10 = FALSE)


    #Using just age
    sex_result_chron_age=run_sex_assoc_all_celltype_region(df, outcome_col = "age",
                                                 sex_col = "imputed_sex",
                                                 covar_cols = c("PC1", "PC2", "PC3", "PC4", "PC5", "apoe_score"),
                                                 min_n = 50)
    sex_result_chron_age$p_adj=p.adjust(sex_result_chron_age$p, method = "fdr")

    sex_result_chron_age=sex_result_chron_age[sex_result_chron_age$p_adj<0.05, ]

    index=1
    plot_sex_outcome_celltype_region(df, cell_type = sex_result_pass_fdr[index,]$cell_type,
                                     region = sex_result_pass_fdr[index,]$region,
                                     outcome_col = "resid_mean_corrected", sex_col= "imputed_sex")

    #APOE - pretend it's a PRS and see how it correlates with age residuals

    results_apoe=run_prs_assoc_all_celltype_region(df, c("apoe_score"), outcome_col = "resid_mean_corrected", covar_cols = c("imputed_sex", "PC1", "PC2", "PC3", "PC4", "PC5"))
    results_apoe$p_adj=p.adjust(results_apoe$p, method = "fdr")
    results_pass_fdr=results_apoe[results_apoe$p_adj<0.05,]
    results_apoe[which.min(results_apoe$p_adj),]
    hist (results_apoe$p_adj, main="APOE correlation with residuals FDR p-values", xlab="FDR p-value", xlim=c(0,1))
    abline(v=0.05, col="red", lty=2)

    #Look at the TS scores in MSN_D1
    z=results[results$cell_type=="MSN_D1" & results$region=="CaH", ]
    z=z[order(z$p_adj),]
    hist (z$t, xlab="t-statistic", main="PRS in MSN_D1")
    abline (v=z[z$prs=="TS_2019",]$t, col='red')


}






##########################
# PRS association functions
###########################

run_prs_assoc_all_celltype_region <- function(df,
                                              prs_categories,
                                              outcome_col = "resid_mean_corrected",
                                              covar_cols = c("imputed_sex",
                                                             "PC1", "PC2", "PC3", "PC4", "PC5",
                                                             "apoe_score"),
                                              min_n = 50,
                                              n_permutation = 1e6,
                                              n_core=1L) {

    required_cols <- c("donor", "cell_type", "region", outcome_col, covar_cols)
    missing_required <- setdiff(required_cols, colnames(df))
    if (length(missing_required) > 0) {
        stop("Missing required columns in df: ",
             paste(missing_required, collapse = ", "))
    }

    pairs <- unique(df[, c("cell_type", "region"), drop = FALSE])
    pairs$cell_type <- as.character(pairs$cell_type)
    pairs$region <- as.character(pairs$region)

    res_list <- vector("list", nrow(pairs))
    k <- 0L

    for (i in seq_len(nrow(pairs))) {
        cell_type_i <- pairs$cell_type[i]
        region_i <- pairs$region[i]

        idx <- df$cell_type == cell_type_i & df$region == region_i
        n_rows_pair <- sum(idx)
        n_donors_pair <- length(unique(df$donor[idx]))

        logger::log_info(
            paste0(
                "Running PRS associations: cell_type=", cell_type_i,
                " region=", region_i,
                " n_prs=", length(prs_categories),
                " n_rows=", n_rows_pair,
                " n_donors=", n_donors_pair
            )
        )

        one <- run_prs_assoc_one_celltype_region(
            df = df,
            cell_type = cell_type_i,
            region = region_i,
            prs_categories = prs_categories,
            outcome_col = outcome_col,
            covar_cols = covar_cols,
            min_n = min_n,
            n_permutation=n_permutation,
            n_core=n_core
        )

        if (!is.null(one)) {
            k <- k + 1L
            res_list[[k]] <- one
        }
    }

    if (k == 0L) {
        return(NULL)
    }

    res <- do.call(rbind, res_list[seq_len(k)])
    rownames(res) <- NULL
    res
}




run_prs_assoc_one_celltype_region <- function(df,
                                              cell_type,
                                              region,
                                              prs_categories,
                                              outcome_col,
                                              covar_cols,
                                              min_n = 100,
                                              n_permutation = 1e6,
                                              n_core = 1L) {

    if (!is.numeric(n_core) || length(n_core) != 1L || is.na(n_core) || n_core < 1) {
        stop("n_core must be a single integer >= 1")
    }
    n_core <- as.integer(n_core)

    # Pre-filter once (cell_type/region/outcome/covars/complete cases)
    df_celltype_region <- prepare_prs_assoc_data(
        df = df,
        cell_type = cell_type,
        region = region,
        outcome_col = outcome_col,
        covar_cols = covar_cols,
        prs_name = prs_categories,
        min_n = min_n,
        verbose = TRUE
    )
    if (is.null(df_celltype_region)) {
        return(NULL)
    }

    worker <- function(prs_name) {
        run_prs_assoc_one_celltype_region_prs(
            df = df_celltype_region,
            cell_type = cell_type,
            region = region,
            prs_name = prs_name,
            outcome_col = outcome_col,
            covar_cols = covar_cols,
            min_n = min_n,
            n_permutation = n_permutation,
            verbose = FALSE
        )
    }

    if (n_core == 1L) {
        out_list <- lapply(prs_categories, worker)
    } else {
        out_list <- parallel::mclapply(
            prs_categories,
            worker,
            mc.cores = n_core
        )
    }

    out_list <- out_list[!vapply(out_list, is.null, logical(1))]
    if (length(out_list) == 0L) {
        return(NULL)
    }

    res <- do.call(rbind, out_list)
    rownames(res) <- NULL
    res
}


# cell_type="GABA_CGE_DFC"; region="DFC"; prs_name="BIP_2021"; outcome_col="resid_mean_corrected"; covar_cols = c("imputed_sex","PC1", "PC2", "PC3", "PC4", "PC5", "apoe_score")
# cell_type="GABA_CGE_DFC"; region="DFC"; prs_name="BIP_2021"; outcome_col="apoe_score"; covar_cols = c("imputed_sex","PC1", "PC2", "PC3", "PC4", "PC5")
# cell_type="MSN_D1"; region="CaH"; prs_name="TS_2019"; outcome_col="resid_mean_corrected"; covar_cols = c("imputed_sex","PC1", "PC2", "PC3", "PC4", "PC5", "apoe_score")

run_prs_assoc_one_celltype_region_prs <- function(df,
                                                  cell_type,
                                                  region,
                                                  prs_name,
                                                  outcome_col,
                                                  covar_cols,
                                                  min_n = 100,
                                                  n_permutation = 1e6,
                                                  verbose = FALSE) {

    # Further filter for this PRS (adds prs col + complete cases including PRS)
    d_cc <- prepare_prs_assoc_data(
        df = df,
        cell_type = cell_type,
        region = region,
        outcome_col = outcome_col,
        covar_cols = covar_cols,
        prs_name = prs_name,
        min_n = min_n,
        verbose = verbose
    )
    if (is.null(d_cc)) {
        return(NULL)
    }

    fml <- stats::reformulate(
        termlabels = c(prs_name, covar_cols),
        response = outcome_col
    )

    fit <- stats::lm(fml, data = d_cc)
    s <- summary(fit)
    coef_mat <- s$coefficients

    if (!prs_name %in% rownames(coef_mat)) {
        if (verbose) {
            logger::log_warn(
                paste0(
                    "Skipping PRS (term dropped): cell_type=", cell_type,
                    " region=", region,
                    " prs=", prs_name
                )
            )
        }
        return(NULL)
    }

    p_perm <- NA_real_

    if (!is.null(n_permutation) && is.numeric(n_permutation) && n_permutation > 0) {

        perm_fit <- permuco::aovperm(
            fml,
            data = d_cc,
            np = as.integer(n_permutation),
            method = "freedman_lane"
        )

        perm_tab <- summary(perm_fit)

        if (!is.data.frame(perm_tab) || !prs_name %in% rownames(perm_tab)) {
            stop("Could not extract permutation results for term: ", prs_name)
        }
        if (!"resampled P(>F)" %in% colnames(perm_tab)) {
            stop("Expected column not found in permuco summary: resampled P(>F)")
        }

        p_perm <- unname(perm_tab[prs_name, "resampled P(>F)"])
        #clean up the possibly large permutation object.
        rm (perm_fit)
        gc(verbose=FALSE)
    }

    data.frame(
        cell_type = as.character(cell_type),
        region = as.character(region),
        prs = prs_name,
        outcome = outcome_col,
        n = nrow(d_cc),
        beta = unname(coef_mat[prs_name, "Estimate"]),
        se = unname(coef_mat[prs_name, "Std. Error"]),
        t = unname(coef_mat[prs_name, "t value"]),
        p = unname(coef_mat[prs_name, "Pr(>|t|)"]),
        p_perm = p_perm,
        r2 = unname(s$r.squared),
        adj_r2 = unname(s$adj.r.squared),
        stringsAsFactors = FALSE
    )
}



prepare_prs_assoc_data <- function(df,
                                   cell_type,
                                   region,
                                   outcome_col,
                                   covar_cols,
                                   prs_name = NULL,  # NULL or character vector
                                   min_n = 100L,
                                   verbose = FALSE) {

    idx <- df$cell_type == cell_type & df$region == region
    if (!any(idx)) {
        if (verbose) {
            logger::log_warn(paste0("Skipping (no rows): cell_type=", cell_type,
                                    " region=", region))
        }
        return(NULL)
    }

    keep_cols <- unique(c(
        "donor", "cell_type", "region",
        outcome_col, covar_cols,
        prs_name
    ))

    d <- df[idx, keep_cols, drop = FALSE]

    if (!is.null(prs_name)) {
        missing_prs <- setdiff(prs_name, colnames(d))
        if (length(missing_prs) > 0) {
            if (verbose) {
                logger::log_warn(paste0("Skipping (missing PRS cols): cell_type=", cell_type,
                                        " region=", region,
                                        " prs_missing=", paste(missing_prs, collapse = ",")))
            }
            return(NULL)
        }
    }

    if ("imputed_sex" %in% colnames(d)) {
        d$imputed_sex <- factor(d$imputed_sex)
    }

    cc_cols <- unique(c(outcome_col, covar_cols))
    cc <- stats::complete.cases(d[, cc_cols, drop = FALSE])
    d_cc <- d[cc, , drop = FALSE]

    if (nrow(d_cc) < min_n) {
        if (verbose) {
            logger::log_warn(paste0("Skipping (too few complete cases): cell_type=", cell_type,
                                    " region=", region,
                                    " n_cc=", nrow(d_cc),
                                    " min_n=", min_n))
        }
        return(NULL)
    }

    d_cc
}




add_bh_by_prs <- function(results,
                          p_col = "p",
                          prs_col = "prs",
                          out_col = "p_adj_prs") {

    if (!is.data.frame(results)) {
        stop("results must be a data.frame")
    }
    if (!p_col %in% colnames(results)) {
        stop("p_col not found in results: ", p_col)
    }
    if (!prs_col %in% colnames(results)) {
        stop("prs_col not found in results: ", prs_col)
    }
    if (out_col %in% colnames(results)) {
        stop("out_col already exists in results: ", out_col)
    }

    p <- results[[p_col]]
    if (!is.numeric(p)) {
        stop("p_col must be numeric: ", p_col)
    }

    results[[out_col]] <- NA_real_

    # Apply BH within each PRS category; preserve NAs
    prs_vals <- results[[prs_col]]
    for (prs in unique(prs_vals)) {
        idx <- prs_vals == prs & !is.na(p)
        if (sum(idx) > 0) {
            results[[out_col]][idx] <- stats::p.adjust(p[idx], method = "BH")
        }
    }

    results
}





prepare_data_for_PRS_regression<-function (age_preds_file,prs_file, metadata_file, cellTypeFile=NULL) {
    age_preds=read.table(age_preds_file, header=TRUE, sep="\t", stringsAsFactors = FALSE)
    age_preds=age_preds[, c("cell_type", "region", "donor", "age", "pred_mean_corrected", "resid_mean_corrected"), drop=FALSE]

    if (!is.null(cellTypeFile)) {
        logger::log_info(paste0("Filtering to cell types in file: ", cellTypeFile))
        ct=read.table(cellTypeFile, header=FALSE, sep="\t", stringsAsFactors = FALSE)$V1
        age_preds=age_preds[age_preds$cell_type %in% ct, , drop=FALSE]

    }

    #add z-score scalings. - these should be by cell type / region
    age_preds=add_age_preds_zscores(age_preds)

    metadata=.parse_metadata(metadata_file)

    prs=parsePRSFile(prs_file)
    colnames(prs)[1]="donor"

    #get the list of PRS categories.
    prs_category_list=setdiff(colnames(prs), "donor")

    ## Merge age predictions with donor metadata
    df <- merge(
        age_preds,
        metadata,
        by = "donor",
        all.x = TRUE,
        sort = FALSE
    )

    ## Merge in PRS (many columns, donor-level)
    df <- merge(
        df,
        prs,
        by = "donor",
        all.x = TRUE,
        sort = FALSE
    )

    result=list(data=df, prs_categories=prs_category_list)
    return (result)
}

.parse_metadata<-function (metadata_file) {
    cols=c("donor", "PC1", "PC2", "PC3", "PC4", "PC5", "apoe_score", "imputed_sex")
    metadata=read.table(metadata_file, header=TRUE, sep="\t", stringsAsFactors = FALSE)
    md=unique(metadata[, cols, drop=FALSE])
    md=bican.mccarroll.differentialexpression::scale_PC_cols(md)
    md$imputed_sex=factor(md$imputed_sex)
    #apoe_score is intended to be linear so it should NOT be a factor!
    return (md)
}

#' Parse and simplify a PRS score table
#'
#' Read a tab-delimited PRS score file and return a simplified table containing
#' \code{sample_id} and one PRS score per trait. The input may contain multiple
#' PRS columns per trait corresponding to different ancestries and model types.
#'
#' Columns whose names contain \code{"graphPred"} are removed. For each remaining
#' trait, a single PRS column is selected by prioritizing EUR population PGSX
#' models first, followed by any PGSX model, and then EUR population \code{pgs}
#' or \code{blupx} models. If multiple columns match the same priority level, the
#' first matching column is used.
#'
#' Trait names are derived from column names by removing any \code{adj_} prefix
#' and truncating at the first occurrence of a \code{TRAIT_YYYY} pattern.
#'
#' @param inFile Path to a tab-delimited file with header. Must contain a
#'   \code{sample_id} column.
#'
#' @return A data.frame with \code{sample_id} and one selected PRS column per
#'   trait.
#'
#' @export
parsePRSFile <- function(inFile) {
    a <- read.table(inFile, header = TRUE, stringsAsFactors = FALSE, sep = "\t")

    all_cols <- colnames(a)

    # Step 1: Filter to PRS columns (starts with adj_) and exclude graphPred
    prs_cols <- all_cols[
            !grepl("graphPred", all_cols)
    ]

    # Step 2: Extract canonical trait name
    extract_trait <- function(colname) {
        stripped <- sub("^adj_", "", colname)
        # Remove ancestry & model suffix
        sub("([A-Za-z0-9]+_[0-9]{4}).*", "\\1", stripped)
    }

    # Map PRS columns to base traits
    trait_map <- setNames(sapply(prs_cols, extract_trait), prs_cols)

    # Step 3: Preference-based selection
    select_best_col <- function(trait) {
        cols <- names(trait_map[trait_map == trait])

        eur_pgsx <- grep("\\.EUR_.*pgsx", cols, value = TRUE)
        if (length(eur_pgsx) > 0) return(eur_pgsx[1])

        any_pgsx <- grep("pgsx", cols, value = TRUE)
        if (length(any_pgsx) > 0) return(any_pgsx[1])

        eur_alt <- grep("\\.EUR_.*(pgs|blupx)", cols, value = TRUE)
        if (length(eur_alt) > 0) return(eur_alt[1])

        return(cols[1])  # fallback
    }

    # Step 4: Apply selection
    unique_traits <- unique(trait_map)
    best_cols <- sapply(unique_traits, select_best_col, USE.NAMES = TRUE)

    # Step 5: Subset and rename
    simplified_names <- names(best_cols)
    a_selected <- a[, c(best_cols), drop = FALSE]
    colnames(a_selected) <- c(simplified_names)

    return(a_selected)
}

#############################
# SEX EFFECT functions
#############################

# run_sex_assoc_all_celltype_region <- function(df,
#                                               outcome_col = "resid_mean_corrected",
#                                               sex_col = "imputed_sex",
#                                               covar_cols = c("PC1", "PC2", "PC3", "PC4", "PC5", "apoe_score"),
#                                               min_n = 50) {
#
#     required_cols <- c("donor", "cell_type", "region", outcome_col, sex_col, covar_cols)
#     missing_required <- setdiff(required_cols, colnames(df))
#     if (length(missing_required) > 0) {
#         stop("Missing required columns in df: ",
#              paste(missing_required, collapse = ", "))
#     }
#
#     pairs <- unique(df[, c("cell_type", "region"), drop = FALSE])
#     pairs$cell_type <- as.character(pairs$cell_type)
#     pairs$region <- as.character(pairs$region)
#
#     res_list <- vector("list", nrow(pairs))
#     k <- 0L
#
#     for (i in seq_len(nrow(pairs))) {
#         cell_type_i <- pairs$cell_type[i]
#         region_i <- pairs$region[i]
#
#         idx <- df$cell_type == cell_type_i & df$region == region_i
#         n_rows_pair <- sum(idx)
#         n_donors_pair <- length(unique(df$donor[idx]))
#
#         logger::log_info(
#             paste0(
#                 "Running sex association: cell_type=", cell_type_i,
#                 " region=", region_i,
#                 " n_rows=", n_rows_pair,
#                 " n_donors=", n_donors_pair
#             )
#         )
#
#         one <- run_sex_assoc_one_celltype_region(
#             df = df,
#             cell_type = cell_type_i,
#             region = region_i,
#             outcome_col = outcome_col,
#             sex_col = sex_col,
#             covar_cols = covar_cols,
#             min_n = min_n
#         )
#
#         if (!is.null(one)) {
#             k <- k + 1L
#             res_list[[k]] <- one
#         }
#     }
#
#     if (k == 0L) {
#         return(NULL)
#     }
#
#     res <- do.call(rbind, res_list[seq_len(k)])
#     rownames(res) <- NULL
#     res
# }
#
#
# # cell_type="GABA_CGE_DFC"; region="DFC"; outcome_col="resid_mean_corrected"; sex_col="imputed_sex"; covar_cols = c("PC1", "PC2", "PC3", "PC4", "PC5", "apoe_score", "age")
# # cell_type="GABA_CGE_DFC"; region="DFC"; outcome_col="age"; sex_col="imputed_sex"; covar_cols = c("PC1", "PC2", "PC3", "PC4", "PC5", "apoe_score")
# run_sex_assoc_one_celltype_region <- function(df,
#                                               cell_type,
#                                               region,
#                                               outcome_col,
#                                               sex_col,
#                                               covar_cols,
#                                               min_n = 50) {
#
#     if (!outcome_col %in% colnames(df)) {
#         stop("outcome_col not found in df: ", outcome_col)
#     }
#     if (!sex_col %in% colnames(df)) {
#         stop("sex_col not found in df: ", sex_col)
#     }
#
#     missing_covars <- setdiff(covar_cols, colnames(df))
#     if (length(missing_covars) > 0) {
#         stop("Missing covariate columns in df: ",
#              paste(missing_covars, collapse = ", "))
#     }
#
#     idx <- df$cell_type == cell_type & df$region == region
#     if (!any(idx)) {
#         logger::log_warn(
#             paste0(
#                 "Skipping sex association (no rows): cell_type=", cell_type,
#                 " region=", region
#             )
#         )
#         return(NULL)
#     }
#
#     keep_cols <- unique(c("donor", "cell_type", "region", outcome_col, sex_col, covar_cols))
#     d <- df[idx, keep_cols, drop = FALSE]
#
#     cc <- stats::complete.cases(d[, c(outcome_col, sex_col, covar_cols), drop = FALSE])
#     n_cc <- sum(cc)
#     if (n_cc < min_n) {
#         logger::log_warn(
#             paste0(
#                 "Skipping sex association (too few complete cases): cell_type=", cell_type,
#                 " region=", region,
#                 " n_cc=", n_cc,
#                 " min_n=", min_n
#             )
#         )
#         return(NULL)
#     }
#
#     d_cc <- d[cc, , drop = FALSE]
#
#     d_cc[[sex_col]] <- factor(d_cc[[sex_col]])
#     if (nlevels(d_cc[[sex_col]]) < 2) {
#         logger::log_warn(
#             paste0(
#                 "Skipping sex association (only one sex level present): cell_type=", cell_type,
#                 " region=", region,
#                 " sex_levels=", paste(levels(d_cc[[sex_col]]), collapse = ",")
#             )
#         )
#         return(NULL)
#     }
#
#     fml <- stats::reformulate(
#         termlabels = c(sex_col, covar_cols),
#         response = outcome_col
#     )
#
#     fit <- stats::lm(fml, data = d_cc)
#     s <- summary(fit)
#     coef_mat <- s$coefficients
#
#     sex_terms <- grep(paste0("^", sex_col), rownames(coef_mat), value = TRUE)
#     if (length(sex_terms) < 1) {
#         logger::log_warn(
#             paste0(
#                 "Skipping sex association (sex term not in coefficients): cell_type=", cell_type,
#                 " region=", region
#             )
#         )
#         return(NULL)
#     }
#
#     sex_term <- sex_terms[1]
#
#     data.frame(
#         cell_type = as.character(cell_type),
#         region = as.character(region),
#         outcome = outcome_col,
#         predictor = sex_col,
#         term = sex_term,
#         sex_ref_level = levels(d_cc[[sex_col]])[1],
#         n = n_cc,
#         beta = unname(coef_mat[sex_term, "Estimate"]),
#         se = unname(coef_mat[sex_term, "Std. Error"]),
#         t = unname(coef_mat[sex_term, "t value"]),
#         p = unname(coef_mat[sex_term, "Pr(>|t|)"]),
#         r2 = unname(s$r.squared),
#         adj_r2 = unname(s$adj.r.squared),
#         stringsAsFactors = FALSE
#     )
# }
#
# plot_sex_outcome_celltype_region <- function(df,
#                                              cell_type,
#                                              region,
#                                              outcome_col = "resid_mean_corrected",
#                                              sex_col = "imputed_sex") {
#
#     ## R CMD CHECK safety
#     outcome <- sex <- NULL
#
#     if (!outcome_col %in% colnames(df)) {
#         stop("outcome_col not found in df: ", outcome_col)
#     }
#     if (!sex_col %in% colnames(df)) {
#         stop("sex_col not found in df: ", sex_col)
#     }
#
#     idx <- df$cell_type == cell_type & df$region == region
#     if (!any(idx)) {
#         stop("No rows for cell_type=", cell_type, " and region=", region)
#     }
#
#     d <- df[idx, c(outcome_col, sex_col), drop = FALSE]
#     colnames(d) <- c("outcome", "sex")
#
#     d$sex <- factor(d$sex)
#     d <- d[stats::complete.cases(d), , drop = FALSE]
#
#     if (nrow(d) == 0) {
#         stop("No complete cases after filtering")
#     }
#
#     ggplot2::ggplot(d, ggplot2::aes(x = sex, y = outcome)) +
#         ggplot2::geom_boxplot(outlier.shape = NA, linewidth = 0.7) +
#         ggplot2::geom_jitter(width = 0.15, alpha = 0.5, size = 1.2) +
#         ggplot2::labs(
#             x = "Sex",
#             y = outcome_col,
#             title = paste0(cell_type, " / ", region)
#         ) +
#         ggplot2::theme_bw()
# }
#
plot_fdr_heatmap <- function(result,
                                 value_col = "p_adj",
                                 use_neg_log10 = TRUE) {

    ## R CMD CHECK safety
    cell_type <- region <- value <- NULL

    required_cols <- c("cell_type", "region", value_col)
    missing_required <- setdiff(required_cols, colnames(result))
    if (length(missing_required) > 0) {
        stop("Missing required columns in sex_result: ",
             paste(missing_required, collapse = ", "))
    }

    d <- result[, c("cell_type", "region", value_col), drop = FALSE]
    colnames(d)[3] <- "value"

    d$cell_type <- as.character(d$cell_type)
    d$region <- as.character(d$region)

    # If there are duplicates (shouldn't be), keep the minimum adjusted p-value
    key <- paste(d$cell_type, d$region, sep = "\t")
    if (anyDuplicated(key)) {
        keep <- !duplicated(key)
        # Initialize with first occurrences
        d0 <- d[keep, , drop = FALSE]
        # Replace with minima across duplicates
        tapply_min <- tapply(d$value, key, min, na.rm = TRUE)
        d0$value <- unname(tapply_min[paste(d0$cell_type, d0$region, sep = "\t")])
        d <- d0
    }

    if (use_neg_log10) {
        d$value <- -log10(pmax(d$value, .Machine$double.xmin))
        fill_label <- paste0("-log10(", value_col, ")")
    } else {
        fill_label <- value_col
    }

    if (use_neg_log10) {
        p <- ggplot2::scale_fill_viridis_c(
            option = "magma",
            na.value = "grey90"
        )
    } else {
        max_abs <- max(abs(d$value), na.rm = TRUE)
        p <- ggplot2::scale_fill_gradient2(
            low = "blue",
            mid = "white",
            high = "red",
            midpoint = 0,
            limits = c(-max_abs, max_abs),
            na.value = "grey90"
        )
    }


    # Order factors so rows/cols are stable
    d$region <- factor(d$region, levels = sort(unique(d$region)))
    d$cell_type <- factor(d$cell_type, levels = sort(unique(d$cell_type)))

    ggplot2::ggplot(d, ggplot2::aes(x = cell_type, y = region, fill = value)) +
        ggplot2::geom_tile() +
        ggplot2::labs(x = "Cell type", y = "Region", fill = fill_label) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            panel.grid = ggplot2::element_blank()
        ) +
        p
}

plot_fdr_heatmap2 <- function(result,
                             value_col = "p_adj",
                             use_neg_log10 = TRUE,
                             title = NULL,
                             fdr_col = "p_fdr",
                             fdr_thresh = 0.05,
                             label_mode = c("none", "star", "fdr"),
                             fdr_digits = 2,
                             label_size = 3) {

    ## R CMD CHECK safety
    cell_type <- region <- value <- label <- NULL

    label_mode <- match.arg(label_mode)

    required_cols <- c("cell_type", "region", value_col)
    if (label_mode != "none") {
        required_cols <- unique(c(required_cols, fdr_col))
    }
    missing_required <- setdiff(required_cols, colnames(result))
    if (length(missing_required) > 0) {
        stop("Missing required columns in result: ",
             paste(missing_required, collapse = ", "))
    }

    d <- result[, required_cols, drop = FALSE]
    colnames(d)[colnames(d) == value_col] <- "value"

    d$cell_type <- as.character(d$cell_type)
    d$region <- as.character(d$region)

    # If there are duplicates (shouldn't be), keep the minimum adjusted p-value (or min fdr)
    key <- paste(d$cell_type, d$region, sep = "\t")
    if (anyDuplicated(key)) {
        keep <- !duplicated(key)
        d0 <- d[keep, , drop = FALSE]

        # value: if value_col is p-like we should take min; otherwise keep first
        # (retain your original behavior: min over duplicates)
        tapply_min_value <- tapply(d$value, key, min, na.rm = TRUE)
        d0$value <- unname(tapply_min_value[paste(d0$cell_type, d0$region, sep = "\t")])

        if (label_mode != "none") {
            tapply_min_fdr <- tapply(d[[fdr_col]], key, min, na.rm = TRUE)
            d0[[fdr_col]] <- unname(tapply_min_fdr[paste(d0$cell_type, d0$region, sep = "\t")])
        }

        d <- d0
    }

    if (use_neg_log10) {
        d$value <- -log10(pmax(d$value, .Machine$double.xmin))
        fill_label <- paste0("-log10(", value_col, ")")
    } else {
        fill_label <- value_col
    }

    # Label to overlay
    if (label_mode == "none") {
        d$label <- NA_character_
    } else if (label_mode == "star") {
        sig <- !is.na(d[[fdr_col]]) & d[[fdr_col]] < fdr_thresh
        d$label <- ifelse(sig, "*", "")
    } else if (label_mode == "fdr") {
        # print FDR only for significant cells; blank otherwise
        sig <- !is.na(d[[fdr_col]]) & d[[fdr_col]] < fdr_thresh
        d$label <- ""
        d$label[sig] <- format(d[[fdr_col]][sig],
                               scientific = TRUE,
                               digits = fdr_digits)
    }

    # Order factors so rows/cols are stable
    d$region <- factor(d$region, levels = sort(unique(d$region)))
    d$cell_type <- factor(d$cell_type, levels = sort(unique(d$cell_type)))

    # Scales
    if (use_neg_log10) {
        fill_scale <- ggplot2::scale_fill_viridis_c(
            option = "magma",
            na.value = "grey90"
        )
    } else {
        max_abs <- max(abs(d$value), na.rm = TRUE)
        fill_scale <- ggplot2::scale_fill_gradient2(
            low = "blue",
            mid = "white",
            high = "red",
            midpoint = 0,
            limits = c(-max_abs, max_abs),
            na.value = "grey90"
        )
    }

    g <- ggplot2::ggplot(d, ggplot2::aes(x = cell_type, y = region, fill = value)) +
        ggplot2::geom_tile() +
        ggplot2::labs(x = "Cell type", y = "Region", fill = fill_label, title = title) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            panel.grid = ggplot2::element_blank()
        ) +
        fill_scale

    if (label_mode != "none") {
        g <- g + ggplot2::geom_text(ggplot2::aes(label = label), size = label_size)
    }

    g
}


time_run_log <- function(expr, label = NULL) {
    t_start <- Sys.time()

    res <- eval(expr, envir = parent.frame())

    elapsed_sec <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

    msg <- if (is.null(label)) {
        paste0("Time taken: ", round(elapsed_sec, 2), " seconds")
    } else {
        paste0(label, " – time taken: ", round(elapsed_sec, 2), " seconds")
    }

    logger::log_info(msg)

    invisible(res)
}

#
#
