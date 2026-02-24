#' # styler::tidyverse_style(scope = I(c("spaces", "indention", "tokens")))
#'
#' # libraries -------------------------------------------------------------------------------------------------------
#'
#' library(data.table)
#' library(bican.mccarroll.differentialexpression)
#' library(TRADEtools)
#' library(pheatmap)
#' library(viridisLite)
#' library(sandwich)
#' library(lmtest)
#' library(Matrix)
#' library(matrixStats)
#' library(cluster)
#' library(fgsea)
#'
#'
#'
#' # define cell type ordering and groups ----------------------------------------------------------------------------
#'
#' ## paths ----
#'
#' # ct_file = "C:/Users/sburger/Documents/manuscripts/bican/cell_types_use.txt"
#' ct_file <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/sburger_tmp/cell_types_use.txt"
#'
#'
#' ## run ----
#'
#' cell_types_use <- scan(ct_file, what = character())
#'
#' #TODO: Not used
#' # dfc_cell_types <- c(
#' #   "glutamatergic_L23IT",
#' #   "glutamatergic_L4IT",
#' #   "glutamatergic_L5IT",
#' #   "glutamatergic_L6IT",
#' #   "GABA_PVALB",
#' #   "GABA_SST",
#' #   "GABA_VIP",
#' #   "GABA_LAMP5"
#' # )
#'
#' # striatal_cell_types <- c(
#' #   "MSN_D1_matrix",
#' #   "MSN_D1_striosome",
#' #   "MSN_D2_matrix",
#' #   "MSN_D2_striosome",
#' #   "GABA_TAC3-PLPP4",
#' #   "GABA_PTHLH-PVALB"
#' # )
#'
#' non_neuron_types <- c("astrocyte", "OPC", "oligodendrocyte", "microglia")
#' neuron_types <- cell_types_use[!(cell_types_use %in% non_neuron_types)]
#'
#' region_order <- c("CaH", "Pu", "NAC", "ic", "DFC")
#'
#'
#'
#' # load DE results -------------------------------------------------------------------------------------------------
#'
#' ## functions ----
#'
#' #'
#' format_de_results <- function(df, gene_to_chr) {
#'
#'   dt <- as.data.table(df)
#'   colnames(dt) <- c("gene", "cell_type", "region", "test", "log_fc", "ave_expr", "t", "p_value", "adj_p_val", "b", "z_std")
#'
#'   dt <- merge(gene_to_chr, dt, by = "gene", all.y = TRUE) # add chr
#'   dt[, log_fc_se := log_fc / t] # add log_fc_se column
#'
#'   dt
#'
#' }
#'
#'
#' ## paths ----
#'
#' # de_results_dir = "C:/Users/sburger/Documents/r_io/server_copy/differential_expression/results"
#' de_results_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results"
#'
#' de_dir <- paste0(de_results_dir, "/LEVEL_3/sex_age/cell_type")
#' de_region_subset_dir <- paste0(de_results_dir, "/LEVEL_3/sex_age/cell_type_subset_region")
#' de_region_interaction_dir <- paste0(de_results_dir, "/LEVEL_3/sex_age/cell_type_region_interaction_absolute_effects")
#'
#' # gene_to_chr_file <- "C:/Users/sburger/Documents/manuscripts/bican/gene_to_chromosome.txt"
#' gene_to_chr_file <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/sburger_tmp/gene_to_chromosome.txt"
#'
#'
#' ## run ----
#'
#' gene_to_chr <- fread(gene_to_chr_file)
#' gene_to_chr <- gene_to_chr[chr %in% c(1:22, "X", "Y", "M")]
#'
#' # regions combined
#' de_age <- format_de_results(bican.mccarroll.differentialexpression::parse_de_inputs(de_dir, "age", ct_file), gene_to_chr)
#' de_sex <- format_de_results(bican.mccarroll.differentialexpression::parse_de_inputs(de_dir, "female_vs_male", ct_file), gene_to_chr)
#'
#' # region subset
#' de_rs_age <- format_de_results(parse_de_inputs(de_region_subset_dir, "age", ct_file), gene_to_chr)
#' de_rs_sex <- format_de_results(parse_de_inputs(de_region_subset_dir, "female_vs_male", ct_file), gene_to_chr)
#'
#' # region interaction
#' de_ri_age <- format_de_results(parse_de_inputs(de_region_interaction_dir, "age", ct_file), gene_to_chr)
#'
#'
#' # volcano plots ---------------------------------------------------------------------------------------------------
#'
#' ## functions ----
#'
#' #'
#' de_volcano_plot <- function(de_dt, cell_type_use, region_use, fdr_cutoff = 0.05, abs_log_fc_cutoff = log2(1.05)) {
#'
#'   if (is.na(region_use)) {
#'     dt <- de_dt[cell_type == cell_type_use & is.na(region)]
#'   } else {
#'     dt <- de_dt[cell_type == cell_type_use & region == region_use]
#'   }
#'
#'   range <- dt[, max(abs(log_fc))]
#'   range <- c(-range, range)
#'
#'   par(pty = "s", xpd = FALSE)
#'   dt[, plot(log_fc, -log10(adj_p_val), pch = 20, xlim = range, col = "lightgrey",
#'             xlab = "Effect size, log2", ylab = "Adjusted p-value, -log10", main = "")]
#'   dt[adj_p_val < fdr_cutoff, points(log_fc, -log10(adj_p_val), pch = 20, xlim = range, col = "cornflowerblue",
#'                                     xlab = "Effect size, log2", ylab = "Adjusted p-value, -log10", main = "")]
#'   title(main = cell_type_use, adj = 0)
#'   abline(h = -log10(fdr_cutoff), lty = 2)
#'   abline(v = c(-abs_log_fc_cutoff, abs_log_fc_cutoff), lty = 2)
#'
#'   up <- dt[log_fc > abs_log_fc_cutoff & adj_p_val < fdr_cutoff, .N]
#'   down <- dt[log_fc < -abs_log_fc_cutoff & adj_p_val < fdr_cutoff, .N]
#'   p <- binom.test(up, up + down, 0.5)$p.value
#'
#'   par(xpd = NA)
#'   if (p < 0.05) {
#'     p <- formatC(p, format = "e", digits = 1)
#'     legend("topright", inset = c(0, -0.16),
#'            c(paste("proportion ↑ =", round(up / (up + down), 2)), paste("p-value =", p)), bty = "n")
#'   }
#'
#' }
#'
#'
#' ## run ----
#'
#' for (i in cell_types_use) {
#'
#'   de_volcano_plot(de_age, cell_type_use = i, region_use = NA)
#'
#' }
#'
#'
#'
#' # scatter plots ---------------------------------------------------------------------------------------------------
#'
#' ## functions ----
#'
#' #'
#' de_scatter_plot <- function(
#'     de_dt,
#'     cell_type_a,
#'     cell_type_b,
#'     region_a = NA,
#'     region_b = NA,
#'     fdr_cutoff = 0.05,
#'     add_fit = TRUE
#' ) {
#'
#'   if (is.na(region_a)) {
#'     x <- de_dt[cell_type == cell_type_a & is.na(region)]
#'   } else {
#'     x <- de_dt[cell_type == cell_type_a & region == region_a]
#'   }
#'
#'   if (is.na(region_b)) {
#'     y <- de_dt[cell_type == cell_type_b & is.na(region)]
#'   } else {
#'     y <- de_dt[cell_type == cell_type_b & region == region_b]
#'   }
#'
#'   name_a <- paste0(toupper(substr(cell_type_a, 1, 1)), substr(cell_type_a, 2, nchar(cell_type_a)))
#'   name_b <- paste0(toupper(substr(cell_type_b, 1, 1)), substr(cell_type_b, 2, nchar(cell_type_b)))
#'
#'   m <- merge(x, y, by = c("chr", "gene"))
#'
#'   range <- m[adj_p_val.x < 0.05 | adj_p_val.y < 0.05, max(abs(c(log_fc.x, log_fc.y)))]
#'   range <- c(-range, range)
#'
#'   par(pty = "s", mar = c(6, 6, 4, 2))
#'
#'   m[, plot(log_fc.x, log_fc.y, pch = 20, col = "lightgrey", xlim = range, ylim = range,
#'            xlab = c("Effect size, log2", paste(name_a, region_a, sep = ", ")),
#'            ylab = c(paste(name_b, region_b, sep = ", "), "Effect size, log2"))]
#'   m[(adj_p_val.x < 0.05 | adj_p_val.y < 0.05), points(log_fc.x, log_fc.y, pch = 20, col = "cornflowerblue")]
#'
#'   abline(h = 0, v = 0, lty = 2)
#'   abline(0, 1, lty = 2)
#'
#'   c <- m[adj_p_val.x < 0.05 | adj_p_val.y < 0.05, cor.test(log_fc.x, log_fc.y, method = "spearman")]
#'   # p = c$p.value # these tests are always highly significant
#'   # p = formatC(p, format = "e", digits = 1)
#'   rho_sqrd <- round(c$estimate^2, 2)
#'   legend("topleft", legend = bquote(rho^2 * " = " * .(rho_sqrd)), bty = "n")
#'
#'   if (rho_sqrd > 0.5 & add_fit) {
#'     fit <- lm(log_fc.y ~ log_fc.x, data = m[adj_p_val.x < 0.05 & adj_p_val.y < 0.05])
#'     ct <- coeftest(fit, vcov = vcovHC(fit, type = "HC1"))
#'     b <- ct["log_fc.x", "Estimate"]
#'     b_se <- ct["log_fc.x", "Std. Error"]
#'     b_ci <- b + c(-1, 1) * 1.96 * b_se
#'     fit_color <- "tomato"
#'     if (b_ci[1] >= 1 & b_ci[1] <= 1) fit_color <- "lightgrey"
#'     abline(fit, lty = 2, col = fit_color)
#'     legend("bottomright", legend = bquote(beta * " = " * .(round(b_ci[1], 2)) * " - " * .(round(b_ci[2], 2))), bty = "n")
#'   }
#'
#' }
#'
#'
#' ## run ----
#'
#' pdf("/downloads/tmp/age_DE_scatterplots.pdf")
#' if (TRUE) {
#'
#'   #TODO: there are 3 region combinations, this is "within cell type"
#'   de_scatter_plot(de_ri_age, "MSN_D1_matrix", "MSN_D1_matrix", "CaH", "Pu")
#'   de_scatter_plot(de_ri_age, "MSN_D1_matrix", "MSN_D1_matrix", "CaH", "NAC")
#'   de_scatter_plot(de_ri_age, "MSN_D1_matrix", "MSN_D1_matrix", "Pu", "NAC")
#'   de_scatter_plot(de_ri_age, "MSN_D2_matrix", "MSN_D2_matrix", "CaH", "Pu")
#'   de_scatter_plot(de_ri_age, "MSN_D2_matrix", "MSN_D2_matrix", "CaH", "NAC")
#'   de_scatter_plot(de_ri_age, "MSN_D2_matrix", "MSN_D2_matrix", "Pu", "NAC")
#'   de_scatter_plot(de_ri_age, "MSN_D1_striosome", "MSN_D1_striosome", "CaH", "Pu")
#'   de_scatter_plot(de_ri_age, "MSN_D1_striosome", "MSN_D1_striosome", "CaH", "NAC")
#'   de_scatter_plot(de_ri_age, "MSN_D1_striosome", "MSN_D1_striosome", "Pu", "NAC")
#'   de_scatter_plot(de_ri_age, "MSN_D2_striosome", "MSN_D2_striosome", "CaH", "Pu")
#'   de_scatter_plot(de_ri_age, "MSN_D2_striosome", "MSN_D2_striosome", "CaH", "NAC")
#'   de_scatter_plot(de_ri_age, "MSN_D2_striosome", "MSN_D2_striosome", "Pu", "NAC")
#'
#'   de_scatter_plot(de_ri_age, "astrocyte", "astrocyte", "CaH", "Pu")
#'   de_scatter_plot(de_ri_age, "astrocyte", "astrocyte", "CaH", "NAC")
#'   de_scatter_plot(de_ri_age, "astrocyte", "astrocyte", "Pu", "NAC")
#'   de_scatter_plot(de_ri_age, "OPC", "OPC", "CaH", "Pu")
#'   de_scatter_plot(de_ri_age, "OPC", "OPC", "CaH", "NAC")
#'   de_scatter_plot(de_ri_age, "OPC", "OPC", "Pu", "NAC")
#'   de_scatter_plot(de_ri_age, "oligodendrocyte", "oligodendrocyte", "CaH", "Pu")
#'   de_scatter_plot(de_ri_age, "oligodendrocyte", "oligodendrocyte", "CaH", "NAC")
#'   de_scatter_plot(de_ri_age, "oligodendrocyte", "oligodendrocyte", "Pu", "NAC")
#'   de_scatter_plot(de_ri_age, "microglia", "microglia", "CaH", "Pu")
#'   de_scatter_plot(de_ri_age, "microglia", "microglia", "CaH", "NAC")
#'   de_scatter_plot(de_ri_age, "microglia", "microglia", "Pu", "NAC")
#'
#'   de_scatter_plot(de_ri_age, "astrocyte", "astrocyte", "CaH", "DFC")
#'   de_scatter_plot(de_ri_age, "OPC", "OPC", "CaH", "DFC")
#'   de_scatter_plot(de_ri_age, "oligodendrocyte", "oligodendrocyte", "CaH", "DFC")
#'   de_scatter_plot(de_ri_age, "microglia", "microglia", "CaH", "DFC")
#'
#'   de_scatter_plot(de_age, "MSN_D1_matrix", "MSN_D2_matrix")
#'   de_scatter_plot(de_age, "MSN_D1_striosome", "MSN_D2_striosome")
#'   de_scatter_plot(de_age, "MSN_D1_matrix", "MSN_D1_striosome")
#'   de_scatter_plot(de_age, "MSN_D2_matrix", "MSN_D2_striosome")
#'
#'   de_scatter_plot(de_age, "astrocyte", "OPC")
#'   de_scatter_plot(de_age, "astrocyte", "oligodendrocyte")
#'   de_scatter_plot(de_age, "astrocyte", "microglia")
#'   de_scatter_plot(de_age, "OPC", "oligodendrocyte")
#'
#'   de_scatter_plot(de_age, "MSN_D1_matrix", "GABA_PTHLH-PVALB")
#'   de_scatter_plot(de_ri_age, "MSN_D1_matrix", "glutamatergic_L23IT", "CaH", "DFC")
#'
#' }
#' dev.off()
#'
#'
#'
#' # correlation heatmap ---------------------------------------------------------------------------------------------
#'
#' ## functions ----
#'
#' #'
#' create_cor_mat <- function(de_dt, cell_types_use, regions_use, non_neuron_types, fdr_cutoff = 0.05) {
#'
#'   de_dt <- copy(de_dt)
#'
#'   de_dt <- de_dt[cell_type %in% cell_types_use]
#'   de_dt <- de_dt[region %in% regions_use]
#'
#'   de_dt <- de_dt[region != "ic" | cell_type %in% non_neuron_types]
#'
#'   de_dt[, cell_type := factor(cell_type, levels = cell_types_use)]
#'   de_dt[, region := factor(region, levels = regions_use)]
#'   de_dt[, cr := paste(cell_type, region, sep = "__")]
#'
#'   setorderv(de_dt, c("cell_type", "region"))
#'
#'   n <- length(unique(de_dt$cr))
#'   out_mat <- matrix(NA, nrow = n, ncol = n, dimnames = list(unique(de_dt$cr), unique(de_dt$cr)))
#'   for (i in unique(de_dt$cr)) {
#'     message(i)
#'     for (j in unique(de_dt$cr)) {
#'       if (i == j) { out_mat[i, j] <- 1
#'       next }
#'
#'       a <- de_dt[cr == i]
#'       b <- de_dt[cr == j]
#'       m <- merge(a, b, by = "gene")
#'
#'       c <- m[adj_p_val.x < fdr_cutoff | adj_p_val.y < fdr_cutoff, cor.test(log_fc.x, log_fc.y, method = "spearman")]
#'       r <- sign(c$estimate) * c$estimate^2
#'
#'       out_mat[i, j] <- r
#'     }
#'   }
#'
#'   out_mat
#'
#' }
#'
#'
#' ## run ----
#'
#' cor_mat_main <- create_cor_mat(de_ri_age, cell_types_use, regions_use = c("CaH", "DFC"), non_neuron_types)
#' cor_mat_supp <- create_cor_mat(de_ri_age, cell_types_use, regions_use = region_order, non_neuron_types)
#'
#' pheatmap(
#'   cor_mat_main,
#'   breaks = seq(-1, 1, length.out = 101),
#'   color = colorRampPalette(c("steelblue", "white", "darkorange"))(100),
#'   clustering_method = "complete"
#' )
#'
#' pheatmap(
#'   cor_mat_supp,
#'   breaks = seq(-1, 1, length.out = 101),
#'   color = colorRampPalette(c("steelblue", "white", "darkorange"))(100),
#'   clustering_method = "complete"
#' )
#'
#'
#'
#' # load cell metadata ----------------------------------------------------------------------------------------------
#'
#' ## paths ----
#'
#' # cell_metadata_file = "C:/Users/sburger/Documents/r_io/server_copy/metadata/CAP_freeze_3/CAP_cell_metadata.annotated.txt.gz"
#' cell_metadata_file <- "/broad/bican_um1_mccarroll/RNAseq/analysis/cellarium_upload/CAP_freeze_3/CAP_cell_metadata.annotated.txt.gz"
#'
#'
#' ## run ----
#'
#' cell_metadata <- fread(cell_metadata_file)
#'
#' donor_ages_dt <- unique(cell_metadata[, .(donor_external_id, age)])
#' donor_ages <- as.numeric(donor_ages_dt$age)
#' names(donor_ages) <- donor_ages_dt$donor_external_id
#' donor_ages <- sort(donor_ages)
#'
#'
#'
#' # load and summarize metacells ------------------------------------------------------------------------------------
#'
#' ## functions ----
#'
#' #'
#' load_metacells <- function(path, regions_use = c("CaH", "Pu", "NAC", "ic", "DFC")) {
#'
#'   dt <- fread(path)
#'   names(dt)[1] <- "gene"
#'   mat <- as.matrix(dt, rownames = "gene")
#'   mat <- as(mat, "dgCMatrix")
#'
#'   col_metadata <- as.data.table(do.call(rbind, strsplit(colnames(mat), "__")))
#'   colnames(col_metadata) <- c("donor", "village", "cell_type", "region", "single_cell_assay")
#'
#'   print(cell_types_use)
#'   col_mask <- col_metadata[, cell_type %in% cell_types_use & region %in% regions_use]
#'   col_metadata <- col_metadata[cell_type %in% cell_types_use & region %in% regions_use]
#'   mat <- mat[, col_mask]
#'
#'   col_metadata[, group := paste(donor, cell_type, region, sep = "__")]
#'   f <- col_metadata[, factor(group)]
#'
#'   # design matrix
#'   dm <- sparseMatrix(i = seq_along(f), j = as.integer(f), x = 1, dims = c(length(f), nlevels(f)), dimnames = list(NULL, levels(f)))
#'
#'   # aggregate columns
#'   mat <- mat %*% dm
#'
#'   # normalize
#'   cs <- colSums(mat)
#'   d <- Diagonal(x = 1e6 / cs)
#'   mat_tpm <- mat %*% d
#'   colnames(mat_tpm) <- colnames(mat)
#'
#'   col_metadata <- as.data.table(do.call(rbind, strsplit(colnames(mat_tpm), "__")))
#'   colnames(col_metadata) <- c("donor", "cell_type", "region")
#'
#'   list(metacells = mat_tpm, col_metadata = col_metadata)
#'
#' }
#'
#' #'
#' summarize_metacells <- function(metacells, col_metadata, donor_ages) {
#'
#'   # avoid modifying caller's data.table by reference
#'   col_metadata <- data.table::copy(col_metadata)
#'
#'   # map donor to age
#'   col_metadata[, age := unname(donor_ages[donor])]
#'
#'   col_metadata[, age_bin := cut(
#'     age,
#'     breaks = c(-Inf, 39, 49, 59, 69, 79, 89, Inf),
#'     labels = c("30", "40", "50", "60", "70", "80", "90"),
#'     right = TRUE
#'   )]
#'
#'   # group columns by cell_type, region
#'   col_metadata[, cr := paste(cell_type, region, sep = "__")]
#'   groups <- split(seq_len(nrow(col_metadata)), col_metadata$cr)
#'
#'   genes <- rownames(metacells)
#'   age_bins <- c("30", "40", "50", "60", "70", "80", "90")
#'
#'   res_list <- vector("list", length(groups))
#'   names(res_list) <- names(groups)
#'
#'   i <- 0
#'   for (gn in names(groups)) {
#'     message(gn)
#'     i <- i + 1
#'     cols <- groups[[gn]]
#'
#'     sub_sp <- metacells[, cols, drop = FALSE]
#'     sub_dense <- as.matrix(sub_sp)
#'
#'     # existing summary across all donors/metacells in this (cell_type, region) group
#'     st <- row_stats_block_fast(sub_dense)
#'
#'     cr_parts <- data.table::tstrsplit(gn, "__", fixed = TRUE)
#'     ct <- cr_parts[[1]][1]
#'     rg <- cr_parts[[2]][1]
#'
#'     dt_out <- data.table::data.table(
#'       gene = genes,
#'       cell_type = ct,
#'       region = rg,
#'       median = as.numeric(st$median),
#'       mad = as.numeric(st$mad),
#'       q_10 = as.numeric(st$q_10),
#'       q_25 = as.numeric(st$q_25),
#'       q_75 = as.numeric(st$q_75),
#'       q_90 = as.numeric(st$q_90),
#'       n_donors = as.integer(st$n) # keep your current meaning
#'     )
#'
#'     # add age-binned medians
#'     bins_for_cols <- col_metadata$age_bin[cols]
#'     for (b in age_bins) {
#'       idx <- which(bins_for_cols == b)
#'       colname <- paste0("median_", b)
#'
#'       if (length(idx) == 0L) {
#'         dt_out[[colname]] <- NA_real_
#'       } else if (length(idx) == 1L) {
#'         # median of 1 column is the column itself
#'         dt_out[[colname]] <- as.numeric(sub_dense[, idx])
#'       } else {
#'         dt_out[[colname]] <- matrixStats::rowMedians(sub_dense[, idx, drop = FALSE])
#'       }
#'     }
#'
#'     res_list[[i]] <- dt_out
#'   }
#'
#'   out <- data.table::rbindlist(res_list, use.names = TRUE, fill = TRUE)
#'   data.table::setorder(out, gene, cell_type, region)
#'   out
#'
#' }
#'
#' #'
#' row_stats_block_fast <- function(mat) {
#'
#'   stopifnot(requireNamespace("matrixStats", quietly = TRUE))
#'   med <- matrixStats::rowMedians(mat, na.rm = TRUE)
#'   mad <- matrixStats::rowMads(mat, na.rm = TRUE, constant = 1.4826)
#'   q <- matrixStats::rowQuantiles(mat, probs = c(0.10, 0.25, 0.75, 0.90), na.rm = TRUE, type = 7)
#'   n <- apply(mat, 1, function(x) sum(!is.na(x)))
#'   list(
#'     median = med,
#'     mad = mad,
#'     q_10 = q[, 1],
#'     q_25 = q[, 2],
#'     q_75 = q[, 3],
#'     q_90 = q[, 4],
#'     n = n
#'   )
#'
#' }
#'
#' #'
#' split_metacells_by_cell_type_region <- function(metacells, col_metadata, donor_ages) {
#'
#'   # don't modify caller by reference
#'   col_metadata <- data.table::copy(col_metadata)
#'
#'   # map donor to age
#'   col_metadata[, age := unname(donor_ages[donor])]
#'
#'   # group key
#'   col_metadata[, cr := paste(cell_type, region, sep = "__")]
#'   groups <- split(seq_len(nrow(col_metadata)), col_metadata$cr)
#'
#'   out <- vector("list", length(groups))
#'   names(out) <- names(groups)
#'
#'   i <- 0
#'   for (gn in names(groups)) {
#'     i <- i + 1
#'     cols <- groups[[gn]]
#'
#'     sub_sp <- metacells[, cols, drop = FALSE]
#'
#'     donors <- col_metadata$donor[cols]
#'
#'     # ages for sorting
#'     ages <- unname(donor_ages[donors])
#'     ord <- order(ages, donors)
#'     sub_sp <- sub_sp[, ord, drop = FALSE]
#'     donors <- donors[ord]
#'
#'     colnames(sub_sp) <- donors
#'
#'     out[[i]] <- sub_sp
#'   }
#'
#'   out
#'
#' }
#'
#'
#' ## paths ----
#'
#' # metacells_file = "C:/Users/sburger/Documents/r_io/server_copy/metacells/LEVEL_3/donor_rxn_DGEList_counts.tsv.gz"
#' metacells_file <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/metacells/LEVEL_3/donor_rxn_DGEList_counts.tsv.gz"
#'
#'
#' ## run ----
#'
#' # load metacells
#' tmp <- load_metacells(metacells_file)
#' metacells <- tmp$metacells
#' col_metadata <- tmp$col_metadata
#'
#' # summarize gene expression using metacells; includes median TPM per decade
#' metacell_summary <- summarize_metacells(metacells, col_metadata, donor_ages)
#'
#' # split metacells out by cell-type-region
#' metacell_cr_list <- split_metacells_by_cell_type_region(metacells, col_metadata, donor_ages)
#'
#'
#'
#' # prepare DE matrices ---------------------------------------------------------------------------------------------
#'
#' ## functions ----
#'
#' #'
#' prepare_de_matrices <- function(de_dt,
#'                                 metacell_summary,
#'                                 cell_types_use,
#'                                 fdr_cutoff = 0.01,
#'                                 abs_lfc_cutoff = log2(1.05),
#'                                 min_tpm = 10,
#'                                 regions_use = c("CaH", "DFC")) {
#'
#'   de_dt <- data.table::copy(de_dt)
#'
#'   de_dt <- de_dt[cell_type %in% cell_types_use]
#'
#'   de_dt <- merge(de_dt, metacell_summary[region == "CaH", .(gene, cell_type, median_tpm_ca = median)], by = c("gene", "cell_type"), all.x = TRUE)
#'   de_dt <- merge(de_dt, metacell_summary[region == "DFC", .(gene, cell_type, median_tpm_dfc = median)], by = c("gene", "cell_type"), all.x = TRUE)
#'
#'   de_dt <- de_dt[median_tpm_ca > min_tpm | median_tpm_dfc > min_tpm]
#'
#'   lfc_mat <- as.matrix(dcast(de_dt, gene ~ cell_type, value.var = "log_fc"), rownames = "gene")[, cell_types_use]
#'   fdr_mat <- as.matrix(dcast(de_dt, gene ~ cell_type, value.var = "adj_p_val"), rownames = "gene")[, cell_types_use]
#'
#'   n_sig <- rowSums((fdr_mat < fdr_cutoff) & (abs(lfc_mat) > abs_lfc_cutoff), na.rm = TRUE)
#'   table(n_sig > 0)
#'   lfc_mat <- lfc_mat[n_sig > 0, ]
#'   fdr_mat <- fdr_mat[n_sig > 0, ]
#'
#'   # setting missing values to 0 to avoids NA issues
#'   lfc_mat[is.na(lfc_mat)] <- 0
#'   fdr_mat[is.na(fdr_mat)] <- 1
#'
#'   lfc_mat_z <- scale(lfc_mat, scale = TRUE, center = TRUE)
#'
#'   list(lfc_mat = lfc_mat, fdr_mat = fdr_mat, lfc_mat_z = lfc_mat_z)
#'
#' }
#'
#' #'
#' prepare_region_lfc_matrix <- function(de_dt, genes_use, cell_types_use, regions_use) {
#'
#'   de_dt <- data.table::copy(de_dt)
#'
#'   de_dt <- de_dt[cell_type %in% cell_types_use & region %in% regions_use]
#'
#'   lfc_mat <- as.matrix(dcast(de_dt, gene ~ cell_type + region, value.var = "log_fc", sep = "__"), rownames = "gene")
#'   lfc_mat <- lfc_mat[genes_use, ]
#'   tmp <- as.data.table(do.call(rbind, strsplit(colnames(lfc_mat), "__")))
#'   tmp[, V1 := factor(V1, levels = cell_types_use)]
#'   tmp[, V2 := factor(V2, levels = regions_use)]
#'   col_order <- order(tmp$V1, tmp$V2)
#'   lfc_mat <- lfc_mat[, col_order]
#'
#'   lfc_mat
#'
#' }
#'
#'
#' ## run ----
#'
#' de_age_mat_list <- prepare_de_matrices(de_age, metacell_summary, cell_types_use)
#'
#' de_ri_age_flc_mat <- prepare_region_lfc_matrix(
#'   de_dt=de_ri_age,
#'   genes_use = rownames(de_age_mat_list$lfc_mat),
#'   cell_types_use = cell_types_use,
#'   regions_use = c("CaH", "DFC")
#' )
#'
#'
#'
#' # k-means clustering heatmap --------------------------------------------------------------------------------------
#'
#' ## functions ----
#'
#' #'
#' plot_k_means_silhouette <- function(mat) {
#'
#'   d <- dist(mat)
#'
#'   ks <- 10:30
#'   avg_sil <- sapply(ks, function(k) {
#'
#'     set.seed(42)
#'     km <- stats::kmeans(mat, centers = k, nstart = 200, iter.max = 20)
#'     sil <- cluster::silhouette(km$cluster, d)
#'     mean(sil[, "sil_width"])
#'
#'   })
#'
#'   plot(
#'     ks, avg_sil, type = "b",
#'     xlab = "Number of clusters",
#'     ylab = "Average silhouette width"
#'   )
#'
#' }
#'
#' #'
#' plot_k_means_heatmap <- function(k_means_mat, lfc_mat, scaling_factor) {
#'
#'   set.seed(42)
#'   k <- stats::kmeans(k_means_mat, centers = 19, nstart = 200, iter.max = 20)
#'
#'   k_use <- factor(k$cluster, levels = c(2, 10, 6, 3, 5, 14, 9, 13, 4, 1, 15, 19, 18, 7, 17, 8, 11, 12))
#'   gn <- names(k_use)
#'   k_use <- as.numeric(k_use)
#'   names(k_use) <- gn
#'   k_use <- k_use[!is.na(k_use)]
#'   gene_order <- names(k_use)[order(k_use)]
#'   boundaries <- cumsum(table(k_use))
#'
#'   pheatmap(t(lfc_mat[gene_order, ] * scaling_factor),
#'            cluster_cols = FALSE,
#'            cluster_rows = FALSE,
#'            breaks = seq(from = -1, to = 1, length.out = 101),
#'            colorRampPalette(c("steelblue", "white", "darkorange"))(100),
#'            show_colnames = FALSE,
#'            gaps_col = boundaries
#'   )
#'
#'   k_use
#'
#' }
#'
#'
#' ## run ----
#'
#' plot_k_means_silhouette(de_age_mat_list$lfc_mat_z)
#'
#' gene_clusters <- plot_k_means_heatmap(de_age_mat_list$lfc_mat_z, de_age_mat_list$lfc_mat, scaling_factor = 5)
#' plot_k_means_heatmap(de_age_mat_list$lfc_mat_z, de_ri_age_flc_mat, scaling_factor = 5)
#'
#'
#'
#' # write lightweight DE outputs ------------------------------------------------------------------------------------
#'
#' ## functions ----
#'
#' #'
#' write_de_lite <- function(
#'     de_dt,
#'     metacell_summary,
#'     donor_ages,
#'     cell_type_use,
#'     region_use,
#'     out_name,
#'     out_dir,
#'     n_top = 200,
#'     fdr_thresh = 0.05
#' ) {
#'
#'   de_dt <- copy(de_dt)
#'   metacell_summary <- copy(metacell_summary)
#'
#'   # subset DE to requested cell type
#'   de_dt <- de_dt[cell_type == cell_type_use]
#'
#'   # medians for requested cell type + region
#'   medians <- metacell_summary[
#'     cell_type == cell_type_use & region == region_use,
#'     .(gene, median_30, median_40, median_50, median_60, median_70, median_80, median_90)
#'   ]
#'
#'   # rename medians columns
#'   setnames(medians, gsub("^median$", "median_tpm", names(medians)))
#'   idx_dec <- which(!(names(medians) %in% c("gene")))
#'   setnames(medians, idx_dec, paste0(names(medians)[idx_dec], "s"))
#'
#'   # build up table
#'   up_genes <- de_dt[adj_p_val < fdr_thresh & log_fc > 0, .(gene, log_fc, log_fc_se, t_stat = t, fdr = adj_p_val)]
#'   up_genes <- merge(up_genes, medians, by = "gene", all.x = TRUE)
#'   setorderv(up_genes, "t_stat", -1)
#'   up_genes <- up_genes[seq_len(min(n_top, nrow(up_genes)))]
#'
#'   # build down table
#'   down_genes <- de_dt[adj_p_val < fdr_thresh & log_fc < 0, .(gene, log_fc, log_fc_se, t_stat = t, fdr = adj_p_val)]
#'   down_genes <- merge(down_genes, medians, by = "gene", all.x = TRUE)
#'   setorderv(down_genes, "t_stat", 1)
#'   down_genes <- down_genes[seq_len(min(n_top, nrow(down_genes)))]
#'
#'   # capture summary messages as text
#'   summary_lines <- character(0)
#'
#'   add_line <- function(...) summary_lines <<- c(summary_lines, paste0(...))
#'
#'   add_line("Input summary")
#'   add_line("cell type: ", cell_type_use)
#'   add_line("brain region: ", region_use)
#'   n_donors <- metacell_summary[cell_type == cell_type_use & region == region_use, n_donors][1]
#'   add_line("N donors: ", n_donors)
#'   add_line("donor age range: ", min(donor_ages), " - ", max(donor_ages))
#'   add_line("donors median age: ", median(donor_ages))
#'   str <- convert_ages_to_decade_string(donor_ages)
#'   add_line("donors per decade bin: ", str)
#'   add_line("")
#'   add_line("Results summary")
#'   add_line("genes tested: ", de_dt[, .N])
#'   n_sig <- de_dt[adj_p_val < fdr_thresh, .N]
#'   n_up <- de_dt[adj_p_val < fdr_thresh & log_fc > 0, .N]
#'   n_dn <- de_dt[adj_p_val < fdr_thresh & log_fc < 0, .N]
#'   add_line("FDR < ", fdr_thresh, " total: ", n_sig)
#'   add_line("FDR < ", fdr_thresh, " up: ", n_up)
#'   add_line("FDR < ", fdr_thresh, " down: ", n_dn)
#'   pct_down <- (n_dn / n_sig) * 100
#'   add_line("percent down FDR < ", fdr_thresh, ": ", ifelse(is.na(pct_down), "NA", round(pct_down, 3)))
#'
#'   # out paths
#'   summary_path <- file.path(out_dir, paste0(out_name, "_summary.txt"))
#'   up_path <- file.path(out_dir, paste0(out_name, "_up_genes.txt"))
#'   down_path <- file.path(out_dir, paste0(out_name, "_down_genes.txt"))
#'
#'   summary_text <- paste(
#'     paste(summary_lines, collapse = "\n"),
#'     "",
#'     "Additional files",
#'     basename(up_path),
#'     basename(down_path),
#'     sep = "\n"
#'   )
#'   writeLines(summary_text, con = summary_path)
#'
#'   fwrite(up_genes, file = up_path, sep = "\t", quote = FALSE, na = "NA")
#'   fwrite(down_genes, file = down_path, sep = "\t", quote = FALSE, na = "NA")
#'
#'   invisible(list(
#'     summary_path = summary_path,
#'     up_path = up_path,
#'     down_path = down_path,
#'     up_genes = up_genes,
#'     down_genes = down_genes,
#'     summary_lines = summary_lines
#'   ))
#'
#' }
#'
#' #'
#' convert_ages_to_decade_string <- function(donor_ages) {
#'
#'   dd <- substr(donor_ages, 1, 1)
#'   dd[dd == "2"] <- "3"
#'   dd <- as.numeric(dd)
#'   dd <- dd * 10
#'   tbl <- table(dd)
#'   str <- paste0(names(tbl), "s N = ", tbl, collapse = ", ")
#'
#'   str
#'
#' }
#'
#'
#' ## run ----
#'
#' # write summaries and gene tables for selected cell types
#' de_lite_dir <- "/downloads/tmp/de_lite_old"
#' if (TRUE) {
#'
#'   write_de_lite(
#'     de_dt = de_age,
#'     metacell_summary = metacell_summary,
#'     donor_ages = donor_ages,
#'     cell_type_use = "MSN_D1_matrix",
#'     region_use = "CaH",
#'     out_name = "MSN_D1_matrix",
#'     out_dir = de_lite_dir,
#'     n_top = 200,
#'     fdr_thresh = 0.05
#'   )
#'
#'   write_de_lite(
#'     de_dt = de_age,
#'     metacell_summary = metacell_summary,
#'     donor_ages = donor_ages,
#'     cell_type_use = "glutamatergic_L23IT",
#'     region_use = "DFC",
#'     out_name = "glutamatergic_L23IT",
#'     out_dir = de_lite_dir,
#'     n_top = 200,
#'     fdr_thresh = 0.05
#'   )
#'
#'   write_de_lite(
#'     de_dt = de_age,
#'     metacell_summary = metacell_summary,
#'     donor_ages = donor_ages,
#'     cell_type_use = "GABA_PTHLH-PVALB",
#'     region_use = "CaH",
#'     out_name = "GABA_PTHLH-PVALB",
#'     out_dir = de_lite_dir,
#'     n_top = 200,
#'     fdr_thresh = 0.05
#'   )
#'
#'   write_de_lite(
#'     de_dt = de_age,
#'     metacell_summary = metacell_summary,
#'     donor_ages = donor_ages,
#'     cell_type_use = "astrocyte",
#'     region_use = "CaH",
#'     out_name = "astrocyte",
#'     out_dir = de_lite_dir,
#'     n_top = 200,
#'     fdr_thresh = 0.05
#'   )
#'
#'   write_de_lite(
#'     de_dt = de_age,
#'     metacell_summary = metacell_summary,
#'     donor_ages = donor_ages,
#'     cell_type_use = "OPC",
#'     region_use = "CaH",
#'     out_name = "OPC",
#'     out_dir = de_lite_dir,
#'     n_top = 200,
#'     fdr_thresh = 0.05
#'   )
#'
#'   write_de_lite(
#'     de_dt = de_age,
#'     metacell_summary = metacell_summary,
#'     donor_ages = donor_ages,
#'     cell_type_use = "oligodendrocyte",
#'     region_use = "CaH",
#'     out_name = "oligodendrocyte",
#'     out_dir = de_lite_dir,
#'     n_top = 200,
#'     fdr_thresh = 0.05
#'   )
#'
#'   write_de_lite(
#'     de_dt = de_age,
#'     metacell_summary = metacell_summary,
#'     donor_ages = donor_ages,
#'     cell_type_use = "microglia",
#'     region_use = "CaH",
#'     out_name = "microglia",
#'     out_dir = de_lite_dir,
#'     n_top = 200,
#'     fdr_thresh = 0.05
#'   )
#'
#' }
#'
#'
#'
#' # run GSEA --------------------------------------------------------------------------------------------------------
#'
#' ## functions ----
#'
#' #'
#' run_gsea <- function(de_dt, fgsea_cell_types, gmt_files, seed=42) {
#'   set.seed(seed)
#'   gsea_results <- list()
#'   for (i in fgsea_cell_types) {
#'
#'     message(i)
#'
#'     ranks <- de_dt[cell_type == i, t] # using t-statistic to rank genes
#'     names(ranks) <- de_dt[cell_type == i, gene]
#'     ranks <- sort(ranks)
#'
#'     for (j in gmt_files) {
#'
#'       message(j)
#'
#'       pathways <- gmtPathways(j) # EDIT: would be faster to load all of these once; but fgsea is rate limiting
#'
#'       new_results <- fgsea(pathways = pathways,
#'                            stats = ranks,
#'                            minSize = 15,
#'                            maxSize = 500)
#'
#'       # colnames(new_results) = tolower(colnames(new_results)) # leaving these alone for now; would be nice if they were snake_case like everything else
#'
#'       new_results$cell_type <- i
#'       new_results$gmt <- basename(j)
#'
#'       setorderv(new_results, "pval")
#'
#'       gsea_results[[length(gsea_results) + 1]] <- new_results
#'
#'     }
#'
#'   }
#'
#'   gsea_results <- rbindlist(gsea_results)
#'
#'   gsea_results[, id := 1:nrow(gsea_results)]
#'
#' }
#'
#' #'
#' write_gsea_lite <- function(gsea_results, out_dir) {
#'
#'   for (i in sort(unique(gsea_results$cell_type))) {
#'
#'     message(i)
#'
#'     out_pos <- gsea_results[cell_type == i & grepl("c5.go", gmt, fixed = TRUE) & NES > 0 & padj < 0.05]
#'     out_pos <- out_pos[, .(pathway, p_adj = padj, nes = NES, size, gmt)]
#'     write.table(out_pos, paste0(out_dir, "/", i, "_pos_gsea.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
#'
#'     out_neg <- gsea_results[cell_type == i & grepl("c5.go", gmt, fixed = TRUE) & NES < 0 & padj < 0.05]
#'     out_neg <- out_neg[, .(pathway, p_adj = padj, nes = NES, size, gmt)]
#'     write.table(out_neg, paste0(out_dir, "/", i, "_neg_gsea.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
#'
#'   }
#'
#' }
#'
#'
#' ## run ----
#'
#' fgsea_cell_types <- c(
#'   "MSN_D1_matrix",
#'   "glutamatergic_L23IT",
#'   "GABA_PTHLH-PVALB",
#'   "astrocyte",
#'   "OPC",
#'   "oligodendrocyte",
#'   "microglia"
#' )
#'
#' # gmt_files <- c("C:/Users/sburger/Documents/gene_sets/msigdb_2025/c5.go.bp.v2025.1.Hs.symbols.gmt",
#' #                "C:/Users/sburger/Documents/gene_sets/msigdb_2025/c5.go.cc.v2025.1.Hs.symbols.gmt",
#' #                "C:/Users/sburger/Documents/gene_sets/msigdb_2025/c5.go.mf.v2025.1.Hs.symbols.gmt",
#' #                list.files("C:/Users/sburger/Documents/gene_sets/enrichr", pattern = ".txt", full.names = TRUE))
#' # abridged set we'll likey be using; I did upload the enrichr GMT files to the same parent dir
#' gmt_files <- sort(list.files("/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/sburger_tmp/gmt_files/C5", full=T))
#'
#' # EDIT: I just let this go overnight, but it would make more sense run in parallel
#' gsea_results <- run_gsea(de_age, fgsea_cell_types, gmt_files)
#' # saveRDS(gsea_results, "bican_age_de_gsea_results.rds")
#' # gsea_results = readRDS("bican_age_de_gsea_results.rds")
#'
#' # gsea_lite_dir = "C:/Users/sburger/Documents/manuscripts/bican/gsea_lite"
#' gsea_lite_dir <- "/downloads/tmp/de_gsea_old"
#' write_gsea_lite(gsea_results, gsea_lite_dir)
#'
#'
#'
#' # plot donor-level GEX --------------------------------------------------------------------------------------------
#'
#' ## functions ----
#'
#' #'
#' plot_donor_gex_age_heatmap <- function(metacells, gs, donor_ages, gs_gaps = NULL, cluster_gs = FALSE, transpose = FALSE) {
#'
#'   if (!all(gs %in% rownames(metacells))) { message("WARNING: missing ", gs[!(gs %in% rownames(metacells))]) }
#'   gs <- gs[gs %in% rownames(metacells)]
#'   metacells <- metacells[gs, ]
#'
#'   p_10 <- apply(metacells, 1, quantile, probs = 0.10, na.rm = TRUE)
#'   p_90 <- apply(metacells, 1, quantile, probs = 0.90, na.rm = TRUE)
#'   range <- p_90 - p_10
#'   metacells <- sweep(metacells, 1, p_10, "-")
#'   metacells <- sweep(metacells, 1, range, "/")
#'   metacells[metacells < 0] <- 0
#'   metacells[metacells > 1] <- 1
#'
#'   idx <- which(diff(donor_ages[colnames(metacells)] %/% 10) != 0)
#'
#'   if (cluster_gs) {
#'     dist_mat <- as.dist(1 - cor(t(as.matrix(metacells)), method = "spearman"))
#'     cluster_arg <- hclust(dist_mat, method = "average")
#'   } else {cluster_arg <- FALSE}
#'
#'   if (transpose) {
#'     pheatmap(t(metacells), cluster_rows = FALSE, cluster_cols = cluster_arg, color = viridis(100), show_rownames = FALSE,
#'              gaps_col = gs_gaps, gaps_row = idx)
#'   } else {
#'     pheatmap(metacells, cluster_rows = cluster_arg, cluster_cols = FALSE, color = viridis(100), show_colnames = FALSE,
#'              gaps_col = idx, gaps_row = gs_gaps)
#'   }
#'
#' }
#'
#' #'
#' plot_donor_gex_age_scatterplot <- function(exp_vector, donor_ages, main = "") {
#'
#'   plot(donor_ages[names(exp_vector)], exp_vector, pch = 20, xlab = "Age", ylab = "Expression, TPM", main = main)
#'   # abline(v = 60, lty = 2, col = "lightgrey")
#'   # abline(h = 10, lty = 2, col = "lightgrey")
#'
#' }
#'
#'
#' ## run ----
#'
#' # plot a multi-gene heatmap
#' gs <- c(
#'   "HOMER1", "RIMS1", "NLGN4X", "IL1RAPL1", "SLITRK2", "LRRN1", "CSMD2",
#'   "PPP1R1B", "PPP1R1A", "ADCY1", "ADCY3", "PDE10A", "PDE4B", "PDE7A", "PRKCB",
#'   "KCNQ3", "KCNAB1", "RELN", "DAB1",
#'   "PTGDR", "PTGER2", "IFNLR1", "TLR5", "SERPING1", "IL17RB",
#'   "FKBP5", "BLVRA",
#'   "ATG9B", "SLC17A5",
#'   "TGFBR3", "ACVR2B", "NOTCH2", "CHRD",
#'   "WDR19", "BBS1", "KIZ", "TTC29", "DNAI7"
#' )
#' gs_gaps <- c(7, 15, 19, 25, 27, 29, 33)
#'
#' dev.off()
#'
#' plot_donor_gex_age_heatmap(
#'   metacell_cr_list$MSN_D1_matrix__CaH,
#'   gs = gs,
#'   donor_ages = donor_ages,
#'   gs_gaps = gs_gaps
#' )
#'
#' # plot a single gene scatterplot
#' plot_donor_gex_age_scatterplot(metacell_cr_list$MSN_D1_matrix__CaH["LINC01735", ], donor_ages, main = "")
#'
#' # plot a metagene scatterplot
#' microglia_priming <- c("MS4A6A", "MS4A4A", "MS4A4E", "CD163", "GPNMB", "IL15", "AIM2", "NOD2", "TRIM14", "TNFAIP8", "DOCK5", "ANXA2", "CEACAM1", "PLXNA1", "PLXNC1", "GAS7", "ESR1", "IFIH1", "DDX60L", "HERC5", "PARP9", "DTX3L", "HLA-C", "ERAP2", "BTN3A2")
#' plot_donor_gex_age_scatterplot(
#'   colSums(metacell_cr_list$microglia__CaH[microglia_priming, ]),
#'   donor_ages,
#'   main = ""
#' )
