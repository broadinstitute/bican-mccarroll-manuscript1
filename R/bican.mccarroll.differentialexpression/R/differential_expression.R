
# library(bican.mccarroll.differentialexpression)
# library(variancePartition)
# library(Glimma)
#
# library(ggplot2)
# library(ggrepel)


# data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metacells"
# data_name="donor_rxn_DGEList"
# randVars=c("donor", "village")
# #note: "imputed_sex" moves to a fixed effect!
# fixedVars=c("age", "PC1", "PC2", "PC3", "PC4", "PC5", "pmi_hr", "pct_intronic", "frac_contamination", "imputed_sex", "toxicology_group", "single_cell_assay", "region", "biobank")
# contrast_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metadata/differential_expression_contrasts_all.txt"
# result_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/cell_type_results"
# cellTypeListFile=NULL

# Dropping PMI, HBCAC, toxicology from model.  Not all donors have PMI, PMI does not contribute very much to variance explained.
# Only running on sex and age to be more donor inclusive.
# data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metacells"
# data_name="donor_rxn_DGEList"
# randVars=c("donor", "village")
# fixedVars=c("age", "PC1", "PC2", "PC3", "PC4", "PC5", "pct_intronic", "frac_contamination", "imputed_sex", "single_cell_assay", "region", "biobank")
# contrast_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metadata/differential_expression_contrasts_sex_age.txt"
# result_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/cell_type_results_sex_age"
# outPDF=paste(result_dir, "volcano_plots.pdf", sep="/")
# cellTypeListFile=NULL

# region tests
# data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metacells"
# data_name="donor_rxn_DGEList"
# randVars=c("donor", "village")
# fixedVars=c("age", "PC1", "PC2", "PC3", "PC4", "PC5", "pct_intronic", "frac_contamination", "imputed_sex", "single_cell_assay", "region", "biobank")
# contrast_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metadata/differential_expression_contrasts_region.txt"
# result_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/region_comparisons"
# outPDF=paste(result_dir, "volcano_plots.pdf", sep="/")
# cellTypeListFile="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metadata/cell_types_for_region_analysis.txt"


# bican.mccarroll.differentialexpression::differential_expression(data_dir, data_name, randVars, fixedVars, contrast_file, cellTypeListFile, outPDF, result_dir)

#' Run differential expression analysis for each cell type in the DGEList.
#'
#' This uses a means model (~0 + fixedVars) and a random effects model for the specified random variables.
#' For each contrast group, the fixed effects are reordered so the contrast group is first, which
#' makes all levels of the contrast available for comparison.
#' @param data_dir Directory containing the DGEList data.
#' @param data_name Name of the DGEList data file (without extension).
#' @param randVars Vector of random effect variables.
#' @param fixedVars Vector of fixed effect variables.
#' @param contrast_file Path to the file containing contrast definitions.
#' @param cellTypeListFile A file containing an explicit list of cell types to test.  If NULL, all cell types in the DGEList will be tested.
#' @param outPDF Optional path to output PDF file for plots.
#' @param result_dir Directory to save the differential expression results.
#' @export
differential_expression <- function(data_dir, data_name, randVars, fixedVars, contrast_file, cellTypeListFile=NULL, outPDF=NULL, result_dir) {
    #load the DGEList and prepare the data
    d=bican.mccarroll.differentialexpression::prepare_data_for_differential_expression(data_dir, data_name, randVars, fixedVars)
    dge=d$dge; fixedVars=d$fixedVars; randVars=d$randVars

    dge=filter_dgelist_by_celltype_list(dge, cellTypeListFile)

    contrast_defs <- read.table(contrast_file, stringsAsFactors = FALSE, sep="\t", header=TRUE)

    # Variance Partition by cell type
    cell_type_list=unique(dge$samples$cell_type)
    #cellType="MSN_D1_matrix"
    line <- strrep("=", 80)

    plot_list= list()
    if (length(cell_type_list) > 0) {
        for (cellType in cell_type_list) {
            logger::log_info(line)
            logger::log_info(paste("Creating differential expression analysis for cell type:", cellType))
            logger::log_info(line)

            dge_cell <- dge[, dge$samples$cell_type == cellType, keep.lib.sizes = TRUE]
            #filtering samples by library size
            r<- filter_by_libsize(dge_cell, threshold_sd = 1.96, bins = 50, strTitlePrefix = cellType)
            dge_cell<- r$dge

            #filter to the top 75% of highly expressed genes as a first pass.
            dge_cell<-filter_top_expressed_genes(dge_cell, gene_filter_frac = 0.75, verbose = TRUE)
            #filter to cpm cutoff of 1.
            r2=plot_logCPM_density_quantiles(dge_cell, cpm_cutoff = 1, logCPM_xlim = c(-5, 15), lower_quantile = 0.05, upper_quantile = 0.95, quantile_steps = 5)
            dge_cell=r2$filtered_dge

            #run differential expression
            #this produces one list per contrast comparison.
            z<-differential_expression_one_cell_type(dge_cell, fixedVars, randVars, contrast_defs,
                                                     verbose = TRUE)
            #str(z)

            # flatten the results for summary and plotting
            # keep only data frames, keep ONLY inner names, preserve order
            z_flat <- do.call(c, lapply(unname(z), function(x) x[sapply(x, is.data.frame)]))

            #save the results
            for (contrast in names(z_flat)) {
                out=z_flat[[contrast]]
                n=paste(cellType, contrast, sep="_")
                outFile <- file.path(result_dir, paste0(n, "_DE_results.txt"))
                logger::log_info(paste("Saving results to:", outFile))
                write.table(out, file = outFile, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
            }

            #make a volcano plot for each contrast
            for (contrast in names(z_flat)) {
                n=paste(cellType, contrast, sep="_")
                df <- z_flat[[contrast]]
                if (nrow(df) > 0) {
                    p <- make_volcano(df, fdr_thresh = 0.05, lfc_thresh = 0,
                                      top_n_each = 10, title = paste(cellType, contrast))
                    plot_list[[n]] <- p
                }
            }
        }
    } else {
        logger::log_info("No cell types found in the DGEList samples.")
    }

    if (!is.null(outPDF)) {
        logger::log_info(paste("Saving all plots to PDF:", outPDF))
        grDevices::pdf(outPDF)
        pages=paginate_plots(plot_list, plots_per_page = 2)
        for (i in 1:length(pages)) {
            print(pages[[i]])
        }
        grDevices::dev.off()
    }

}



#Use ~0 + group for a "means" model: If you want a column for each group
#representing the mean expression for that group (and no intercept term),
#you can use this formula. This can be useful for certain types of comparisons.

#design <- model.matrix(~0 + group, data=data)

differential_expression_one_cell_type<-function (dge_cell, fixedVars, randVars, contrast_defs, verbose = TRUE, n_cores = parallel::detectCores() - 2) {
    #have to handle values that are not valid R column names.
    dge_cell$samples<- sanitize_levels(dge_cell$samples)
    #sanitize the contrast definitions to match the data
    #this broke with more complex contrasts so need to do it later.
    # contrast_defs<- sanitize_contrast_levels_old(contrast_defs)

    contrast_groups<-unique (contrast_defs$variable)
    topTables_all_list=list()
    for (contrast_group in contrast_groups) {
        logger::log_info(paste("Running differential expression for contrast group:", contrast_group))
        topTables_list<-differential_expression_one_cell_type_contrast_group(dge_cell, fixedVars, randVars, contrast_defs, contrast_group=contrast_group, verbose = TRUE, n_cores = n_cores)
        topTables_all_list[[contrast_group]]<-topTables_list
    }
    return(topTables_all_list)

}

differential_expression_one_cell_type_contrast_group<-function (dge_cell, fixedVars, randVars, contrast_defs,
                                                                contrast_group="toxicology_group", verbose = TRUE,
                                                                n_cores = parallel::detectCores() - 2) {

    #drop levels in a consistent way.
    #the sex encoding gets mangled somewhere internally in dream.
    dge_cell_this=dge_cell
    dge_cell_this$samples <- droplevels(dge_cell_this$samples)

    # Drop random effects if they have insufficient replication
    rv <- prune_random_effects_insufficient_replication(randVars, data=dge_cell_this$samples)

    # Ensure the contrast grouping variable is at the front of the fixed effects list
    fv=move_to_front(fixedVars, contrast_group)
    #drop fixed effects if they have insufficient replication
    fv <- drop_single_level_rand_effects(fv, metadata=dge_cell_this$samples, verbose = TRUE)

    rand_part <- paste0("(1|", rv, ")", collapse = " + ")
    fixed_part <- paste(fv, collapse = " + ")
    # the fixed effects formula for the design matrix.
    fixed_form <- stats::as.formula(paste(" ~ 0 +", fixed_part))
    #the fixed + random effects formula for the dream model.
    formula_str <- paste(fixed_part, rand_part, sep = " + ")
    full_form <- stats::as.formula(paste(" ~ 0 +", formula_str))

    design <- stats::model.matrix(fixed_form, data = dge_cell_this$samples)
    if (qr(design)$rank < ncol(design)) {
        stop("Design matrix is not full rank; consider dropping colinear variables.")
    }

    contrast_defs_this= contrast_defs[contrast_defs$variable == contrast_group, ]
    #sanitize_contrast_levels_old(contrast_defs_this)

    contrast_defs_this= sanitize_contrast_levels(contrast_defs_this, design, verbose=TRUE)

    contrast_matrix <- generate_contrasts_from_defs(contrast_defs_this, design)
    # variancePartition::plotContrasts(contrast_matrix)

    param <- BiocParallel::MulticoreParam(workers = n_cores)
    cell_type <- unique(dge_cell_this$samples$cell_type)

    ###################
    #After running voom, filter the extreme outlier genes and refit.
    ####################
    vobjDream <- variancePartition::voomWithDreamWeights(counts=dge_cell_this, formula=full_form, data=dge_cell_this$samples, BPPARAM = param)

    # Identify good genes
    genes_to_keep <- filter_high_weight_genes(vobjDream, dge_cell_this, quantile_threshold = 0.999)

    # Subset DGE and refit
    dge_cell_this <- dge_cell_this[genes_to_keep, ]
    vobjDream <- variancePartition::voomWithDreamWeights(dge_cell_this, full_form, data = dge_cell_this$samples, BPPARAM = param, plot=FALSE)

    #is this group a contrast group, or is it continuous (IE: age)
    #This defines L as the contrast matrix of there are levels being compared, or NULL for categories like age.
    has_contrasts_groups <- !all(is.na(contrast_defs_this$reference_level) & is.na(contrast_defs_this$comparison_level))
    L <- if (has_contrasts_groups) contrast_matrix else NULL

    fitmm <- capture_dream_warnings({
        variancePartition::dream(exprObj=vobjDream, formula=full_form, data=dge_cell_this$samples, BPPARAM = param, L=L)
    })

    fitmm <- variancePartition::eBayes(fitmm, trend = TRUE, robust = TRUE)
    #plotSA(fitmm, main=paste(cell_type, "mean-variance trend"))

    log_decide_tests_summary(fitmm, L=contrast_matrix, label = paste("DREAM DE summary for", contrast_group))

    topTables_list <- list()

    if (has_contrasts_groups) {
        for (contrast_name in colnames(contrast_matrix)) {
            topTables_list[[contrast_name]] <- limma::topTable(fitmm, coef = contrast_name, number = Inf)
        }
    } else {
        # If it's a continuous variable, we can just use the first coefficient
        topTables_list[[contrast_group]] <- limma::topTable(fitmm, coef = contrast_group, number = Inf)
    }

    return(topTables_list)
}


# Build limma contrasts from a data.frame that may include expressions like:
#   comparison_level = "(CaH + Pu)/2"
#   reference_level  = "NAC"
# Works for simple cases too (e.g., comparison=opioid, reference=control).
# For continuous tests, leave reference_level and comparison_level as NA.
generate_contrasts_from_defs <- function(contrast_defs, design_matrix) {
    stopifnot(is.data.frame(contrast_defs))
    stopifnot(is.matrix(design_matrix) || is.data.frame(design_matrix))

    design_cols <- colnames(design_matrix)

    # Translate a side ("expr") for a given factor "var" into a makeContrasts-friendly string
    translate_side <- function(expr, var, design_cols) {
        if (is.na(expr) || is.null(expr) || nchar(trimws(expr)) == 0) return("0")
        s <- gsub("\\s+", "", as.character(expr))

        var_pat  <- paste0("^", var)
        var_cols <- grep(var_pat, design_cols, value = TRUE)
        if (length(var_cols) == 0) stop("No design columns found for factor '", var, "'. Did you use '~ 0 + ", var, " + ...'?")
        levels_available <- sub(var_pat, "", var_cols)

        # First, rewrite numeric tokens to their sanitized level name (e.g. "1" -> "X1") iff that exists
        m <- gregexpr("[A-Za-z0-9_.-]+", s, perl = TRUE)
        toks <- regmatches(s, m)[[1]]
        if (length(toks)) {
            mapped <- vapply(toks, function(tok) {
                if (grepl("^[0-9]+(\\.[0-9]+)?$", tok)) {
                    sani <- make.names(tok)           # "1" -> "X1"
                    if (sani %in% levels_available) sani else tok
                } else tok
            }, character(1))
            regmatches(s, m)[[1]] <- mapped
        }

        # Now replace level names with full design column names (var + level)
        levels_available <- levels_available[order(nchar(levels_available), decreasing = TRUE)]
        for (lev in levels_available) {
            s <- gsub(paste0("(?<![A-Za-z0-9_.])", lev, "(?![A-Za-z0-9_.])"),
                      paste0(var, lev), s, perl = TRUE)
        }
        s
    }

    contrast_list <- list()

    for (i in seq_len(nrow(contrast_defs))) {
        row <- contrast_defs[i, ]
        cname <- as.character(row$contrast_name)
        var   <- as.character(row$variable)
        ref   <- row$reference_level
        comp  <- row$comparison_level

        # Continuous covariate: both NA
        if ((is.na(ref) || length(ref) == 0) && (is.na(comp) || length(comp) == 0)) {
            # expect a column exactly named <var> in the design (continuous term)
            if (!(var %in% design_cols)) {
                stop("Continuous term '", var, "' not found in design columns.")
            }
            contrast_list[[cname]] <- var
            next
        }

        # Factor/composite case: translate both sides
        comp_str <- translate_side(comp, var, design_cols)
        ref_str  <- translate_side(ref,  var, design_cols)

        # Build final contrast expression
        # Handle zero on either side to keep expressions tidy
        if (identical(ref_str, "0")) {
            contrast_expr <- comp_str
        } else if (identical(comp_str, "0")) {
            contrast_expr <- paste0("0 - (", ref_str, ")")
        } else {
            contrast_expr <- paste0("(", comp_str, ") - (", ref_str, ")")
        }

        contrast_list[[cname]] <- contrast_expr
    }

    # Construct the contrast matrix
    contrast_matrix <- do.call(
        limma::makeContrasts,
        args = c(contrast_list, list(levels = design_matrix))
    )

    contrast_matrix
}



sanitize_levels <- function(df, exclude = character()) {
    for (col in setdiff(names(df), exclude)) {
        if (is.factor(df[[col]])) {
            levels(df[[col]]) <- make.names(levels(df[[col]]))
        } else if (is.character(df[[col]])) {
            df[[col]] <- make.names(df[[col]])
        }
    }
    df
}



# Simple, design-aware sanitizer:
# - Maps tokens to actual factor levels present in the design (~ 0 + variable).
# - If any unknown token appears in an expression, that side becomes NA.
# - After sanitizing both sides: drop rows where exactly one side is NA.
# - Keep rows where both sides are NA (continuous) or both valid (factor contrast).
sanitize_contrast_levels <- function(contrast_defs, design_matrix, verbose = TRUE) {
    stopifnot(is.data.frame(contrast_defs))
    stopifnot(is.matrix(design_matrix) || is.data.frame(design_matrix))

    design_cols <- colnames(design_matrix)

    sanitize_expr_for_var <- function(expr, var) {
        # NA means "no expression" (allowed for continuous rows or 0-side)
        if (length(expr) == 0 || is.na(expr)) return(NA_character_)
        s <- as.character(expr)

        # pull suffixes for this var from design (means coding: ~ 0 + var)
        var_pat  <- paste0("^", var)
        var_cols <- grep(var_pat, design_cols, value = TRUE)
        if (!length(var_cols)) return(NA_character_)  # no columns → treat as invalid for this var
        lvl_suffixes <- sub(var_pat, "", var_cols)

        # tokenize: contiguous level-like chunks; leave operators/parens
        m <- gregexpr("[A-Za-z0-9_.-]+", s, perl = TRUE)
        toks <- regmatches(s, m)[[1]]
        if (!length(toks)) return(NA_character_)  # nothing meaningful → invalid

        had_unknown <- FALSE
        mapped <- vapply(toks, function(tok) {
            # numeric literal? could be a level like "1"/"2" *or* a true number (e.g., "/2")
            if (grepl("^[0-9]+(\\.[0-9]+)?$", tok)) {
                # prefer direct suffix match, else try make.names("1")->"X1"
                if (tok %in% lvl_suffixes) return(tok)
                tok_sani <- make.names(tok)
                if (tok_sani %in% lvl_suffixes) return(tok_sani)
                return(tok)  # true numeric constant
            }
            # text token already a suffix?
            if (tok %in% lvl_suffixes) return(tok)
            # try sanitized text
            tok_sani <- make.names(tok)
            if (tok_sani %in% lvl_suffixes) return(tok_sani)

            had_unknown <<- TRUE
            tok  # keep to reinsert (we’ll null out the whole side below)
        }, character(1))

        if (had_unknown) return(NA_character_)  # this side invalid → caller will drop the row

        # put mapped tokens back (still suffixes, not full var+suffix)
        regmatches(s, m)[[1]] <- mapped
        s
    }

    out <- contrast_defs
    # force character columns and sanitize per-row
    out$reference_level  <- NA_character_
    out$comparison_level <- NA_character_
    for (i in seq_len(nrow(out))) {
        var <- as.character(contrast_defs$variable[i])
        out$reference_level[i]  <- sanitize_expr_for_var(contrast_defs$reference_level[i],  var)
        out$comparison_level[i] <- sanitize_expr_for_var(contrast_defs$comparison_level[i], var)
    }

    # Decide which rows to keep:
    # - keep if both sides NA (continuous)
    # - keep if both sides non-NA (valid factor contrast)
    # - drop if exactly one side is NA
    both_na     <- is.na(out$reference_level) & is.na(out$comparison_level)
    both_non_na <- !is.na(out$reference_level) & !is.na(out$comparison_level)
    keep <- both_na | both_non_na

    dropped <- out$contrast_name[!keep]
    if (verbose && length(dropped)) {
        message(sprintf("Dropping %d contrast(s) due to missing/unknown levels: %s",
                        length(dropped), paste(dropped, collapse = ", ")))
    }

    out[keep, , drop = FALSE]
}


move_to_front <- function(vec, name_to_move) {
    if (!(name_to_move %in% vec)) {
        stop(sprintf("'%s' not found in vector", name_to_move))
    }
    c(name_to_move, setdiff(vec, name_to_move))
}

log_decide_tests_summary <- function(fit, L, label = "DE summary") {
    contrast_names <- colnames(L)
    out <- summary(limma::decideTests(fit)[, contrast_names, drop = FALSE])
    logger::log_info("{label}:\n{paste(capture.output(print(out)), collapse = '\n')}")
}

# Filter genes with extreme voom weights
filter_high_weight_genes <- function(vobj, dge, quantile_threshold = 0.999, verbose = TRUE) {
    stopifnot(!is.null(vobj), !is.null(dge))

    weights <- vobj$weights
    gene_names <- rownames(dge$counts)

    if (nrow(weights) != length(gene_names)) {
        stop("Number of genes in weights does not match gene names from DGEList.")
    }

    max_weights <- apply(weights, 1, max, na.rm = TRUE)
    threshold <- stats::quantile(max_weights, quantile_threshold, na.rm = TRUE)

    keep <- max_weights < threshold
    n_dropped <- sum(!keep)
    n_total <- length(keep)

    if (verbose) {
        message(sprintf(
            "Filtering %d of %d genes (%.2f%%) with extreme weights (quantile threshold %.3f -> %.2e)",
            n_dropped, n_total, 100 * n_dropped / n_total, quantile_threshold, threshold
        ))
    }

    return(gene_names[keep])
}

capture_dream_warnings <- function(expr) {
    warning_msgs <- character()

    result <- withCallingHandlers(
        expr,
        warning = function(w) {
            msg <- conditionMessage(w)
            if (grepl("Model failed to converge with .*negative eigenvalue", msg)) {
                warning_msgs <<- c(warning_msgs, msg)
                invokeRestart("muffleWarning")  # Suppress the warning
            }
        }
    )

    failed_count <- length(warning_msgs)
    if (failed_count > 0) {
        message(sprintf("%d genes failed to converge", failed_count))
    }

    return(result)
}

summarize_top_tables_for_celltype <- function(topTables_list,
                                              lfc_threshold = 0.5,
                                              pval_threshold = 0.05) {
    # Initialize result holder
    summary_list <- list()

    for (contrast_group in names(topTables_list)) {
        group_results <- topTables_list[[contrast_group]]
        for (contrast_name in names(group_results)) {
            df <- group_results[[contrast_name]]

            # Define status per gene
            status <- rep("NotSig", nrow(df))
            status[df$logFC >= lfc_threshold & df$adj.P.Val < pval_threshold] <- "Up"
            status[df$logFC <= -lfc_threshold & df$adj.P.Val < pval_threshold] <- "Down"

            # Tabulate and align
            tab <- table(factor(status, levels = c("Down", "NotSig", "Up")))
            summary_list[[contrast_name]] <- tab
        }
    }

    # Combine to a matrix
    summary_matrix <- do.call(cbind, summary_list)
    return(summary_matrix)
}



make_volcano <- function(df,
                             fdr_thresh = 0.05,
                             lfc_thresh = 0,
                             top_n_each = 10,
                             title = NULL,
                             point_alpha = 0.6,
                             add_counts_inset = TRUE) {

    if (!all(c("logFC", "P.Value", "adj.P.Val") %in% names(df)))
        stop("df must contain: logFC, P.Value, adj.P.Val")
    if (is.null(rownames(df)))
        stop("Row names must contain gene symbols for labeling.")

    d <- df
    d$neglog10FDR <- -log10(pmax(d$adj.P.Val, .Machine$double.xmin))
    d$sig <- d$adj.P.Val <= fdr_thresh & abs(d$logFC) >= lfc_thresh

    d$dir <- ifelse(d$sig & d$logFC < 0, "down",
                    ifelse(d$sig & d$logFC > 0, "up", "ns"))
    # Fix display order to match your legend preference
    d$dir <- factor(d$dir, levels = c("down", "ns", "up"))

    # labels: top |logFC| among FDR-significant hits, split by direction
    sig_up   <- d[d$dir == "up",   , drop = FALSE]
    sig_down <- d[d$dir == "down", , drop = FALSE]
    ord_up   <- if (nrow(sig_up))   order(-abs(sig_up$logFC), sig_up$adj.P.Val)   else integer(0)
    ord_down <- if (nrow(sig_down)) order(-abs(sig_down$logFC), sig_down$adj.P.Val) else integer(0)
    lab_up   <- if (length(ord_up))   utils::head(sig_up[ord_up,   , drop = FALSE], top_n_each) else d[0,]
    lab_down <- if (length(ord_down)) utils::head(sig_down[ord_down, , drop = FALSE], top_n_each) else d[0,]
    labs <- rbind(lab_up, lab_down)
    if (nrow(labs)) labs$label <- rownames(labs)

    # counts (ordered down/ns/up)
    n_down <- sum(d$dir == "down", na.rm = TRUE)
    n_ns   <- sum(d$dir == "ns",   na.rm = TRUE)
    n_up   <- sum(d$dir == "up",   na.rm = TRUE)

    col_map <- c("down" = "steelblue3", "ns" = "gray70", "up" = "firebrick2")

    legend_labels <- c(
        down = sprintf("down (%s)", format(n_down, big.mark = ",")),
        ns   = sprintf("ns (%s)",   format(n_ns,   big.mark = ",")),
        up   = sprintf("up (%s)",   format(n_up,   big.mark = ","))
    )

    # Make R CMD HAPPY
    logFC<- neglog10FDR <- label <- NULL

    p <- ggplot2::ggplot(d, ggplot2::aes(x = logFC, y = neglog10FDR, color = dir)) +
        ggplot2::geom_point(alpha = point_alpha, size = 1.2) +
        ggplot2::scale_color_manual(
            values = c("down" = "steelblue3", "ns" = "gray70", "up" = "firebrick2"),
            breaks = c("down", "ns", "up"),
            labels = legend_labels,
            drop   = FALSE
        ) +
        ggplot2::geom_hline(yintercept = -log10(fdr_thresh), linetype = "dashed", linewidth = 0.4) +
        ggplot2::labs(
            title = title,
            x = "log2 fold change",
            y = expression(-log[10]("FDR")),
            color = "Direction"
        ) +
        ggplot2::theme_classic(base_size = 12)

    # Symmetric x-axis: ±(1.1 * max |logFC|), also respects lfc_thresh
    xmax <- max(abs(d$logFC), lfc_thresh, na.rm = TRUE)
    if (is.finite(xmax) && xmax > 0) {
        xlim <- c(-1, 1) * (1.1 * xmax)
        p <- p + ggplot2::scale_x_continuous(limits = xlim,
                                             expand = ggplot2::expansion(mult = 0))
    }


    if (lfc_thresh > 0) {
        p <- p + ggplot2::geom_vline(xintercept = c(-lfc_thresh, lfc_thresh),
                                     linetype = "dotted", linewidth = 0.4)
    }

    if (nrow(labs) > 0) {
        p <- p + ggrepel::geom_text_repel(
            data = labs,
            ggplot2::aes(x = logFC, y = neglog10FDR, label = label),
            size = 3,
            color = "black",
            min.segment.length = 0,
            max.overlaps = 10000,
            box.padding = 0.3,
            point.padding = 0.2,
            show.legend = FALSE  # keep legend showing points
        )
    }


    return (p)
}

paginate_plots <- function(plots, plots_per_page = 2) {
    if (!length(plots)) return(invisible(list()))
    ppp <- as.integer(plots_per_page)
    if (is.na(ppp) || ppp < 1) stop("plots_per_page must be a positive integer.")

    # pad with blanks so length is a multiple of ppp
    n_pad <- (ppp - (length(plots) %% ppp)) %% ppp
    if (n_pad > 0) plots <- c(plots, rep(list(cowplot::ggdraw()), n_pad))

    idx <- split(seq_along(plots), ceiling(seq_along(plots) / ppp))

    pages <- lapply(idx, function(ii) {
        cowplot::plot_grid(
            plotlist = plots[ii],
            ncol = 1, nrow = ppp,
            align = "v",
            rel_heights = rep(1, ppp)
        )
    })

    return (pages)
}

# test<-function () {
#     df1=read.table("/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/cell_type_results_sex_age/astrocyte_age_DE_results.txt", header=T, stringsAsFactors=F, sep="\t")
#     df2=read.table("/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/differential_expression/cell_type_results_sex_age/astrocyte_female_vs_male_DE_results.txt", header=T, stringsAsFactors=F, sep="\t")
#
#     l=list(df1=df1, df2=df2)
#
#     plot_list= list()
#
#     p <- make_volcano(df1, fdr_thresh = 0.05, lfc_thresh = 0,
#                       top_n_each = 10, title = paste(cellType, contrast))
#     plot_list[[1]] <- p
#
#     p <- make_volcano(df2, fdr_thresh = 0.05, lfc_thresh = 0,
#                       top_n_each = 10, title = paste(cellType, contrast))
#     plot_list[[2]] <- p
#
#     outPDF="/downloads/test.pdf"
#     if (!is.null(outPDF)) {
#         logger::log_info(paste("Saving all plots to PDF:", outPDF))
#         grDevices::pdf(outPDF)
#         pages=paginate_plots(plot_list, plots_per_page = 2)
#         for (i in 1:length(pages)) {
#             print(pages[[i]])
#         }
#         grDevices::dev.off()
#     }
#
# }

#quickly regenerate the PDF from the files in the result directory.
generate_pdf_from_files<-function (result_dir, outPDF) {
    files=list.files(result_dir, pattern="_DE_results.txt", full.names = TRUE, recursive = TRUE)
    plot_list= list()
    i=1
    for (f in files) {
        df=read.table(f, header=T, stringsAsFactors=F, sep="\t")
        n=gsub("_DE_results.txt", "", basename(f))
        p <- make_volcano(df, fdr_thresh = 0.05, lfc_thresh = 0,
                          top_n_each = 10, title = n)
        plot_list[[i]] <- p
        i=i+1
    }

    if (!is.null(outPDF)) {
        logger::log_info(paste("Saving all plots to PDF:", outPDF))
        grDevices::pdf(outPDF)
        pages=paginate_plots(plot_list, plots_per_page = 2)
        for (i in 1:length(pages)) {
            print(pages[[i]])
        }
        grDevices::dev.off()
    }
}

filter_dgelist_by_celltype_list<-function (dge, cellTypeListFile=NULL) {
    if (is.null(cellTypeListFile)) {
        return(dge)
    }

    size_prefilter <- dim(dge)[2]

    cell_type_list=read.table(cellTypeListFile, stringsAsFactors = FALSE, sep="\t", header=FALSE)$V1
    idx=dge$samples$cell_type %in% cell_type_list
    if (any(is.na(idx))) {
        stop("Some cell types in the list are not present in the DGEList samples.")
    }
    dge_filtered <- dge[, idx, keep.lib.sizes = TRUE]
    size_postfilter <- dim(dge_filtered)[2]
    logger::log_info(paste("Filtered DGEList from", size_prefilter, "to", size_postfilter, "metacells based on cell type list."))
    dge_filtered$samples$cell_type <- factor(dge_filtered$samples$cell_type, levels = cell_type_list)
    return(dge_filtered)
}
