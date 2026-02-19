#' # styler::tidyverse_style(scope = I(c("spaces", "indention", "tokens")))
#'
#' # environment -----------------------------------------------------------------------------------------------------
#'
#' library(data.table)
#' library(Matrix)
#' library(TRADEtools)
#' library(pheatmap)
#' library(viridisLite)
#' library(fgsea)
#' library(bican.mccarroll.differentialexpression)
#'
#' if (FALSE) {
#'
#'   remotes::install_github(
#'     "broadinstitute/bican-mccarroll-manuscript1",
#'     subdir = "R/bican.mccarroll.differentialexpression",
#'     dependencies = FALSE
#'   )
#'
#'   ls(getNamespace("bican.mccarroll.differentialexpression"))
#'
#' }
#'
#'
#'
#'
#'
#' # functions -------------------------------------------------------------------------------------------------------
#'
#' #'
#' format_de_results = function(df, gene_to_chr) {
#'
#'   dt = as.data.table(df)
#'   colnames(dt) = c("gene", "cell_type", "region", "test", "log_fc", "ave_expr", "t", "p_value", "adj_p_val", "b", "z_std")
#'
#'   dt <- merge(gene_to_chr, dt, by = "gene", all.y = TRUE) # add chr
#'   dt[, log_fc_se := log_fc / t] # add log_fc_se column
#'
#'   dt
#'
#' }
#'
#' #'
#' load_metacells = function(path, regions_use = c("CaH", "DFC")) {
#'
#'   #TODO: what is this doing?  I'm guessing it's an error.
#'   #path = metacells_path
#'
#'   dt <- fread(path)
#'   names(dt)[1] = "gene"
#'   mat = as.matrix(dt, rownames = "gene")
#'   mat <- as(mat, "dgCMatrix")
#'
#'   col_metadata = as.data.table(do.call(rbind, strsplit(colnames(mat), "__")))
#'   colnames(col_metadata) = c("donor", "village", "cell_type", "region", "single_cell_assay")
#'
#'   print(cell_types_use)
#'   col_mask = col_metadata[, cell_type %in% cell_types_use & region %in% regions_use]
#'   col_metadata = col_metadata[cell_type %in% cell_types_use & region %in% regions_use]
#'   mat = mat[, col_mask]
#'
#'   col_metadata[, group := paste(donor, cell_type, region, sep = "__")]
#'   f <- col_metadata[, factor(group)]
#'
#'   # design matrix
#'   dm = sparseMatrix(i = seq_along(f), j = as.integer(f), x = 1, dims = c(length(f), nlevels(f)), dimnames = list(NULL, levels(f)))
#'
#'   # aggregate columns
#'   mat <- mat %*% dm
#'
#'   # normalize
#'   cs <- colSums(mat)
#'   d <- Diagonal(x = 1e6 / cs)
#'   mat_tpm <- mat %*% d
#'   colnames(mat_tpm) = colnames(mat)
#'
#'   col_metadata = as.data.table(do.call(rbind, strsplit(colnames(mat_tpm), "__")))
#'   colnames(col_metadata) = c("donor", "cell_type", "region")
#'
#'   list(metacells = mat_tpm, col_metadata = col_metadata)
#'
#' }
#'
#' #TODO: never called?
#' summarize_metacells <- function(metacells, col_metadata) {
#'
#'   # group columns by cell_type, region
#'   col_metadata[, cr := paste(cell_type, region, sep = "__")]
#'   groups <- split(seq_len(nrow(col_metadata)), col_metadata$cr)
#'
#'   genes <- rownames(metacells)
#'
#'   res_list <- vector("list", length(groups))
#'   names(res_list) <- names(groups)
#'
#'   # iterate groups; convert only the group's sub-matrix to dense; manageable if group sizes are reasonable
#'   i <- 0
#'   for (gn in names(groups)) {
#'     message(gn)
#'
#'     i <- i + 1
#'     cols <- groups[[gn]]
#'
#'     sub_sp <- metacells[, cols, drop = FALSE]
#'     sub_dense <- as.matrix(sub_sp)
#'     st <- row_stats_block_fast(sub_dense)
#'
#'     cr_parts <- data.table::tstrsplit(gn, "__", fixed = TRUE)
#'     ct <- cr_parts[[1]][1]
#'     rg <- cr_parts[[2]][1]
#'
#'     dt_out <- data.table::data.table(gene = genes,
#'                                      cell_type = ct,
#'                                      region = rg,
#'                                      median = as.numeric(st$median),
#'                                      mad = as.numeric(st$mad),
#'                                      q_10 = as.numeric(st$q_10),
#'                                      q_25 = as.numeric(st$q_25),
#'                                      q_75 = as.numeric(st$q_75),
#'                                      q_90 = as.numeric(st$q_90),
#'                                      n_donors = as.integer(st$n))
#'
#'     res_list[[i]] <- dt_out
#'   }
#'
#'   out <- data.table::rbindlist(res_list, use.names = TRUE)
#'   data.table::setorder(out, gene, cell_type, region)
#'
#'   out
#'
#' }
#'
#' #TODO: never called?
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
#' run_trade = function(de_dt) {
#'
#'   combinations <- unique(de_dt[, .(test, cell_type, region)])
#'   nrow(combinations)
#'
#'   trade_results <- data.table()
#'   for (i in 1:nrow(combinations)) {
#'     test_use <- unlist(combinations[i, "test"])
#'     cell_type_use <- unlist(combinations[i, "cell_type"])
#'     region_use <- unlist(combinations[i, "region"])
#'
#'     message(test_use, " ", cell_type_use, " ", region_use)
#'
#'     if(is.na(region_use)) {
#'       tmp <- de_dt[test == test_use & cell_type == cell_type_use & is.na(region)]
#'     } else {
#'       tmp <- de_dt[test == test_use & cell_type == cell_type_use & region == region_use]
#'     }
#'
#'     if(nrow(tmp) == 0) { message("no DE results; skipping"); next }
#'
#'     tmp_a <- tmp[chr %in% 1:22]
#'     out_a <- TRADE(
#'       mode = "univariate",
#'       results1 = tmp_a,
#'       results2 = NULL,
#'       annot_table = NULL,
#'       log2FoldChange = "log_fc",
#'       lfcSE = "log_fc_se",
#'       pvalue = "p_value",
#'       model_significant = TRUE, # EDIT: explore the metrics this adds
#'       genes_exclude = NULL,
#'       estimate_sampling_covariance = FALSE,
#'       covariance_matrix_set = "combined",
#'       component_varexplained_threshold = 0,
#'       weight_nocorr = 1,
#'       n_sample = NULL,
#'       verbose = FALSE
#'     )
#'
#'     tmp_x <- tmp[chr == "X"]
#'     out_x <- TRADE(
#'       mode = "univariate",
#'       results1 = tmp_x,
#'       log2FoldChange = "log_fc",
#'       lfcSE = "log_fc_se",
#'       pvalue = "p_value"
#'     )
#'
#'     new_row <- data.table(
#'       test = test_use, cell_type = cell_type_use, region = region_use,
#'       trade_twi = out_a$distribution_summary$transcriptome_wide_impact,
#'       trade_degs = out_a$distribution_summary$Me,
#'       trade_twi_x = out_x$distribution_summary$transcriptome_wide_impact,
#'       trade_degs_x = out_x$distribution_summary$Me
#'     )
#'
#'     trade_results <- rbind(trade_results, new_row)
#'   }
#'
#'   trade_results
#'
#' }
#'
#'
#'
#' # paths -----------------------------------------------------------------------------------------------------------
#'
#' # de results
#' # de_dir = "C:/Users/sburger/Documents/r_io/server_copy/differential_expression/results/LEVEL_3/sex_age/cell_type"
#' # de_region_subset_dir = "C:/Users/sburger/Documents/r_io/server_copy/differential_expression/results/LEVEL_3/sex_age/cell_type_subset_region"
#' # de_region_interaction_dir = "C:/Users/sburger/Documents/r_io/server_copy/differential_expression/results/LEVEL_3/sex_age/cell_type_region_interaction_absolute_effects"
#'
#' de_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results/LEVEL_3/sex_age/cell_type"
#' de_region_subset_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results/LEVEL_3/sex_age/cell_type_subset_region"
#' de_region_interaction_dir = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results/LEVEL_3/sex_age/cell_type_region_interaction_absolute_effects"
#'
#' # metacells
#' # metacells_file = "C:/Users/sburger/Documents/r_io/server_copy/metacells/LEVEL_3/donor_rxn_DGEList_counts.tsv.gz"
#' #TODO this isn't referenced anywhere?
#' metacells_file = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/metacells/LEVEL_3/donor_rxn_DGEList_counts.tsv.gz"
#'
#' # gene to chromosome table
#' #gene_to_chr <- fread("C:/Users/sburger/Documents/manuscripts/bican/gene_to_chromosome.txt")
#' gene_to_chr <- fread("/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/metadata/gene_to_chromosome.txt")
#' gene_to_chr <- gene_to_chr[chr %in% c(1:22, "X", "Y", "M")]
#'
#' # cell metadata
#' cell_metadata = fread("/broad/bican_um1_mccarroll/RNAseq/analysis/cellarium_upload/CAP_freeze_3/CAP_cell_metadata.annotated.txt.gz")
#'
#' # cell types
#' #ct_file = "C:/Users/sburger/Documents/manuscripts/bican/cell_types_use.txt"
#' ct_file = "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/metadata/cell_types_for_de_filtering_plot.txt"
#' cell_types_use = scan(ct_file, what = character())
#'
#'
#' # define and order cell types, regions ----------------------------------------------------------------------------
#'
#' print(cell_types_use)
#'
#' #TODO: NEVER REFERENCED
#' dfc_cell_types <- c(
#'   "glutamatergic_L23IT",
#'   "glutamatergic_L4IT",
#'   "glutamatergic_L5IT",
#'   "glutamatergic_L6IT",
#'   "GABA_PVALB",
#'   "GABA_SST",
#'   "GABA_VIP",
#'   "GABA_LAMP5"
#' )
#'
#' #TODO: NEVER REFERENCED.
#' striatal_cell_types <- c(
#'   "MSN_D1_matrix",
#'   "MSN_D1_striosome",
#'   "MSN_D2_matrix",
#'   "MSN_D2_striosome",
#'   "GABA_TAC3-PLPP4",
#'   "GABA_PTHLH-PVALB"
#' )
#'
#' non_neuron_types <- c("astrocyte", "OPC", "oligodendrocyte", "microglia")
#' neuron_types <- cell_types_use[!(cell_types_use %in% non_neuron_types)]
#' neuron_types
#'
#' region_order <- c("CaH", "Pu", "NAC", "ic", "DFC")
#'
#'
#'
#' # load DE results -------------------------------------------------------------------------------------------------
#'
#' # regions combined
#' de_age = format_de_results(bican.mccarroll.differentialexpression::parse_de_inputs(de_dir, "age", ct_file), gene_to_chr)
#' de_sex = format_de_results(bican.mccarroll.differentialexpression::parse_de_inputs(de_dir, "female_vs_male", ct_file), gene_to_chr)
#'
#' # region subset
#' de_rs_age = format_de_results(parse_de_inputs(de_region_subset_dir, "age", ct_file), gene_to_chr)
#' de_rs_sex = format_de_results(parse_de_inputs(de_region_subset_dir, "female_vs_male", ct_file), gene_to_chr)
#'
#' # region interaction
#' de_ri_age = format_de_results(parse_de_inputs(de_region_interaction_dir, "age", ct_file), gene_to_chr)
#'
#'
#'
#' # TRADE -----------------------------------------------------------------------------------------------------------
#'
#' ## run and save ----
#'
#' # slow
#' trade_age = run_trade(de_age)
#' trade_sex = run_trade(de_sex)
#' trade_rs_age = run_trade(de_rs_age)
#' trade_ri_age = run_trade(de_ri_age)
#' trade_list = list(age = trade_age, sex = trade_sex, region_subset_age = trade_rs_age, region_interaction_age = trade_ri_age)
#' # saveRDS(trade_list, "C:/Users/sburger/Documents/manuscripts/bican/trade_list.rds")
#'
#'
#' z1=de_rs_age[de_rs_age$cell_type == "MSN_D1_matrix" & de_rs_age$region == "Pu"]
#' z2=de_rs_age[de_rs_age$cell_type == "astrocyte" & de_rs_age$region == "Pu"]
#' z3=de_rs_age[de_rs_age$cell_type == "astrocyte" & de_rs_age$region == "NAC"]
#' # run_trade(z1)
#' # run_trade(z2)
#' # run_trade(z3)
#'
#' ## load; if previously saved ---
#'
#' # trade_list = readRDS("C:/Users/sburger/Documents/manuscripts/bican/trade_list.rds")
#'
#'
#' ## plots ----
#'
#' #'
#' trade_heatmap = function(trade_results, cell_types_use, region_order) {
#'
#'   trade_results_wide <- dcast(trade_results, cell_type ~ region, value.var = "trade_twi")
#'   trade_results_wide <- trade_results_wide[cell_type %in% cell_types_use]
#'   trade_results_wide[, cell_type := factor(cell_type, levels = cell_types_use)]
#'   setorderv(trade_results_wide, "cell_type")
#'   trade_results_wide[!(cell_type %in% non_neuron_types), ic := NA]
#'   trade_results_wide <- as.matrix(trade_results_wide, rownames = "cell_type")
#'   trade_results_wide <- trade_results_wide[, region_order]
#'   pheatmap(trade_results_wide, cluster_rows = FALSE, cluster_cols = FALSE, na_col = "white", color = mako(100))
#'
#' }
#'
#' trade_heatmap(trade_list$region_subset_age, cell_types_use, region_order)
#' trade_heatmap(trade_list$region_interaction_age, cell_types_use, region_order)
#'
#' #'
#' trade_barplot = function(trade_results, cell_types_use) {
#'
#'   tmp = as.matrix(trade_results[, .(cell_type, trade_twi)], rownames = "cell_type")
#'   tmp = tmp[cell_types_use[length(cell_types_use):1], ]
#'   par(mar = c(7, 10, 2, 2), mgp = c(5, 1, 0))
#'   barplot(tmp, col = "black", las = 2, horiz = TRUE, xlab = "TRADE transcriptome-wide impact (TWI)")
#'
#' }
#'
#' trade_barplot(trade_list$age, cell_types_use)
#' # trade_barplot(trade_list$sex, cell_types_use)
