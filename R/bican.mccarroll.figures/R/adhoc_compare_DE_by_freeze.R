# old_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results/LEVEL_3/sex_age/cell_type"
# new_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/results/LEVEL_6/sex_age/cell_type"
#
# old_cell_types <- c(
#     "MSN_D1_matrix",
#     "MSN_D1_striosome",
#     "MSN_D2_matrix",
#     "MSN_D2_striosome",
#     "glutamatergic_L23IT",
#     "glutamatergic_L4IT",
#     "glutamatergic_L5IT",
#     "glutamatergic_L6IT",
#     "GABA_TAC3-PLPP4",
#     "GABA_PTHLH-PVALB",
#     "GABA_PVALB",
#     "GABA_SST",
#     "GABA_VIP",
#     "GABA_LAMP5",
#     "astrocyte",
#     "OPC",
#     "oligodendrocyte",
#     "microglia"
# )
#
# new_cell_types <- c(
#     "SPN_D1_matrix",
#     "SPN_D1_striosome",
#     "SPN_D2_matrix",
#     "SPN_D2_striosome",
#     "glut_L23_IT",
#     "glut_L4_IT",
#     "glut_L5_IT",
#     "glut_L6_IT",
#     "GABA_TAC3-PLPP4",
#     "GABA_PTHLH-PVALB",
#     "GABA_PVALB",
#     "GABA_SST",
#     "GABA_VIP",
#     "GABA_LAMP5",
#     "astrocyte",
#     "OPC",
#     "oligodendrocyte",
#     "microglia"
# )
#
# stopifnot(length(old_cell_types) == length(new_cell_types))
#
# cell_type_map <- data.frame(
#     old = old_cell_types,
#     new = new_cell_types,
#     stringsAsFactors = FALSE
# )
#
# compare_one <- function(old_cell_type, new_cell_type) {
#     old_file <- file.path(old_dir, paste0(old_cell_type, "__age_DE_results.txt"))
#     new_file <- file.path(new_dir, paste0(new_cell_type, "__age_DE_results.txt"))
#
#     message("Comparing: ", old_cell_type, " vs ", new_cell_type)
#
#     o <- read.table(
#         old_file,
#         header = TRUE,
#         stringsAsFactors = FALSE,
#         sep = "\t"
#     )
#
#     n <- read.table(
#         new_file,
#         header = TRUE,
#         stringsAsFactors = FALSE,
#         sep = "\t"
#     )
#
#     missing_from_new <- setdiff(rownames(o), rownames(n))
#     missing_from_old <- setdiff(rownames(n), rownames(o))
#
#     shared_genes <- intersect(rownames(o), rownames(n))
#
#     if (length(shared_genes) == 0) {
#         stop(old_cell_type, " vs ", new_cell_type, ": no shared genes")
#     }
#
#     message("  Genes in CAP 3 only:   ", length(missing_from_new))
#     message("  Genes in CAP 3.1 only: ", length(missing_from_old))
#     message("  Shared genes:          ", length(shared_genes))
#
#     if (length(missing_from_new) > 0) {
#         message("  Example CAP 3 only genes: ", paste(head(missing_from_new, 10), collapse = ", "))
#     }
#
#     if (length(missing_from_old) > 0) {
#         message("  Example CAP 3.1 only genes: ", paste(head(missing_from_old, 10), collapse = ", "))
#     }
#
#     o <- o[shared_genes, ]
#     n <- n[shared_genes, ]
#
#     plot(
#         o$t,
#         n$t,
#         xlab = "CAP 3 t-stat",
#         ylab = "CAP 3.1 t-stat",
#         main = paste(old_cell_type, "vs", new_cell_type, "t-stat")
#     )
#     abline(0, 1, col = "red")
#
#     plot(
#         o$logFC,
#         n$logFC,
#         xlab = "CAP 3 logFC",
#         ylab = "CAP 3.1 logFC",
#         main = paste(old_cell_type, "vs", new_cell_type, "logFC")
#     )
#     abline(0, 1, col = "red")
#
#     data.frame(
#         old_cell_type = old_cell_type,
#         new_cell_type = new_cell_type,
#         n_genes = nrow(o),
#         cor_t = cor(o$t, n$t, use = "complete.obs"),
#         cor_logFC = cor(o$logFC, n$logFC, use = "complete.obs"),
#         max_abs_diff_t = max(abs(o$t - n$t), na.rm = TRUE),
#         max_abs_diff_logFC = max(abs(o$logFC - n$logFC), na.rm = TRUE),
#         mean_abs_diff_t = mean(abs(o$t - n$t), na.rm = TRUE),
#         mean_abs_diff_logFC = mean(abs(o$logFC - n$logFC), na.rm = TRUE),
#         stringsAsFactors = FALSE
#     )
# }
#
# runMeFreezeComparison<-function () {
#     pdf("/downloads/CAP3_vs_CAP3.1_age_DE_comparison.pdf", width = 10, height = 5)
#
#     par(mfrow = c(1, 2))
#
#     summary_list <- lapply(seq_len(nrow(cell_type_map)), function(i) {
#         compare_one(cell_type_map$old[i], cell_type_map$new[i])
#     })
#
#     dev.off()
#
#     summary_df <- do.call(rbind, summary_list)
#
#     write.table(
#         summary_df,
#         file = "/downloads/CAP3_vs_CAP3.1_age_DE_comparison_summary.tsv",
#         sep = "\t",
#         quote = FALSE,
#         row.names = FALSE
#     )
#
#     summary_df
# }
#
