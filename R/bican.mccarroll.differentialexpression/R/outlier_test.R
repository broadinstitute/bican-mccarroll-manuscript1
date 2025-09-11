#BiocManager::install("dreamlet")
#BiocManager::install("scran")
#BiocManager::install("DEFormats")
# Dependencies
# library(SingleCellExperiment)
# library(S4Vectors)
# library(SummarizedExperiment)
# library(dreamlet)
# library (DEFormats)
#
# data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metacells"
# data_name="donor_rxn_DGEList"


# testMe<-function (data_dir, data_name) {
#     dge=bican.mccarroll.differentialexpression::loadDGEList(data_dir, prefix = data_name)
#     dge$samples$cell_type_region<-factor(paste(dge$samples$cell_type, dge$samples$region, sep="_"))
#     dge$samples$n_cells=dge$samples$num_nuclei
#
#     factor_vars <- c("imputed_sex", "biobank", "single_cell_assay")
#     numeric_vars <- c("age","PC1","PC2","PC3","PC4","PC5")
#     covariate_cols <- c(numeric_vars, factor_vars)
#
#     form <- if (length(covariate_cols)) stats::reformulate(covariate_cols) else ~ 1
#
#     #for a given cell type and region, create the pseudo-bulk SCE object
#     #ctr=unique(dge$samples$cell_type_region)[1]
#     for (ctr in unique(dge$samples$cell_type_region)) {
#         message(ctr)
#
#         # Create a SingleCellExperiment object for this cell type and region
#         sce <- make_sce_for_ctr(dge, ctr,
#                                 id_col = "sample_name",
#                                 cluster_col = "cell_type_region",
#                                 ncells_col = "num_nuclei")
#
#
#         dge_sub=dge[,dge$samples$cell_type_region==ctr]
#         counts=list(dge_sub$counts)
#         names(counts)=ctr
#         geneDF=DataFrame(rownames (dge_sub$counts))
#         colData=DataFrame(dge_sub$samples)
#         sce<-SingleCellExperiment(counts, colData=colData)
#         cellCounts(sce)
#
#         res.proc <- processAssays(sce, formula = form, assays=ctr)
#
#     }
#
#
#     factor_vars <- c("imputed_sex", "biobank", "single_cell_assay")
#     numeric_vars <- c("age","PC1","PC2","PC3","PC4","PC5")
#
#     covariate_cols <- c(numeric_vars, factor_vars)
#
#     pb<- dge_to_pb_sce(dge,
#                        id_col = "sample_name",
#                        cluster_col = "cell_type_region",
#                        covariate_cols = covariate_cols,
#                        assay_names = NULL,
#                        ncells_col="num_nuclei")
#
#
#     # (A) Coerce types explicitly (adjust lists as needed)
#     factor_vars <- c("imputed_sex", "biobank", "single_cell_assay")
#     numeric_vars <- c("age","PC1","PC2","PC3","PC4","PC5")
#
#     for (nm in intersect(factor_vars, colnames(SummarizedExperiment::colData(pb)))) {
#         SummarizedExperiment::colData(pb)[[nm]] <- as.factor(SummarizedExperiment::colData(pb)[[nm]])
#     }
#     for (nm in intersect(numeric_vars, colnames(SummarizedExperiment::colData(pb)))) {
#         SummarizedExperiment::colData(pb)[[nm]] <- suppressWarnings(as.numeric(SummarizedExperiment::colData(pb)[[nm]]))
#     }
#
#     # (B) Start from your desired set; keep only those that exist
#     covariate_cols <- intersect(covariate_cols, colnames(SummarizedExperiment::colData(pb)))
#
#     # (C) Drop columns that are unusable *globally* (all NA or zero variance)
#     cd <- as.data.frame(SummarizedExperiment::colData(pb))
#     usable_global <- vapply(covariate_cols, function(v) {
#         x <- cd[[v]]
#         x <- x[!is.na(x)]
#         if (length(x) < 2) return(FALSE)
#         if (is.factor(x)) return(length(unique(x)) > 1)
#         # numeric/character
#         ux <- unique(x)
#         length(ux) > 1
#     }, logical(1))
#     covariate_cols <- covariate_cols[usable_global]
#
#     # (D) Optional: tighten further by checking *per assay* variability among samples with nonzero counts
#     # If this eliminates everything, we’ll fall back to ~1 automatically.
#     assays <- SummarizedExperiment::assayNames(pb)
#     cols_with_signal <- function(v) {
#         for (a in assays) {
#             # samples used for this assay are all columns of pb
#             # but we only want ones with nonzero library size in this assay
#             M <- SummarizedExperiment::assay(pb, a)
#             libsz <- colSums(M)
#             keep <- libsz > 0
#             if (!any(keep)) next
#             x <- cd[[v]][keep]
#             x <- x[!is.na(x)]
#             if (length(x) >= 2) {
#                 if (is.factor(cd[[v]])) {
#                     if (length(unique(x)) > 1) return(TRUE)
#                 } else {
#                     if (length(unique(x)) > 1) return(TRUE)
#                 }
#             }
#         }
#         FALSE
#     }
#     if (length(covariate_cols)) {
#         covariate_cols <- covariate_cols[vapply(covariate_cols, cols_with_signal, logical(1))]
#     }
#
#     # (E) Build a safe formula; fallback to ~1 if nothing remains
#     form <- if (length(covariate_cols)) stats::reformulate(covariate_cols) else ~ 1
#
#     # (F) (Optional) Quick per-assay diagnostics so you can see what’s happening
#     diag <- lapply(assays, function(a) {
#         M <- SummarizedExperiment::assay(pb, a)
#         libsz <- colSums(M)
#         keep <- libsz > 0
#         present <- vapply(covariate_cols, function(v) {
#             x <- cd[[v]][keep]
#             nonNA <- sum(!is.na(x))
#             varies <- length(unique(x[!is.na(x)])) > 1
#             paste0("n=", nonNA, "; varies=", varies)
#         }, character(1))
#         data.frame(assay=a, t(present), check.names = FALSE)
#     })
#     diag_df <- do.call(rbind, diag)
#     print(diag_df)  # inspect which covariates actually vary within each assay
#
#     # (G) Now run dreamlet preprocessing
#     res.proc <- processAssays(pb, formula = form)
#
#     # And then outliers
#     out_res <- outlierByAssay(res.proc, assays)
#
#
# }
#
#
#
#
# test<-function () {
#     #BiocManager::install("muscat")
#     library(muscat)
#     library(SingleCellExperiment)
#     library (dreamlet)
#
#     data(example_sce)
#     pb <- aggregateToPseudoBulk(example_sce,
#                                 assay = "counts",
#                                 cluster_id = "cluster_id",
#                                 sample_id = "sample_id",
#                                 verbose = FALSE
#     )
#     res.proc <- processAssays(pb, ~group_id)
#     outlierByAssay( res.proc, c("B cells", "CD14+ Monocytes"))
#
# }
#
#
#
# make_sce_for_ctr <- function(dge, ctr,
#                              id_col = "sample_name",
#                              cluster_col = "cell_type_region",
#                              ncells_col = "n_cells") {
#     stopifnot(inherits(dge, "DGEList"))
#     idx <- dge$samples[[cluster_col]] == ctr
#     stopifnot(any(idx))
#
#     # Subset to this cell_type_region
#     M <- dge$counts[, idx, drop = FALSE]
#     samp_sub <- dge$samples[idx, , drop = FALSE]
#
#     # Aggregate to one column per sample (if multiple columns exist)
#     lab <- as.character(samp_sub[[id_col]])
#     M_by_sample <- t(rowsum(t(M), group = lab, reorder = FALSE))  # genes × unique(sample)
#
#     # Build colData: one row per sample, using the first row per sample
#     uS <- colnames(M_by_sample)
#     split_list <- split.data.frame(samp_sub, samp_sub[[id_col]], drop = TRUE)
#     one_per <- lapply(split_list, function(df) df[1, , drop = FALSE])
#     cd <- do.call(rbind, one_per)
#     cd <- cd[uS, , drop = FALSE]                # order to match columns
#     rownames(cd) <- uS
#
#     # Build n_cells vector per sample for this assay (sum within sample)
#     n_per_sample <- tapply(as.numeric(samp_sub[[ncells_col]]), lab, sum, default = 0)
#     n_per_sample <- n_per_sample[uS]            # align order
#     ncells_list <- lapply(uS, function(s) { v <- as.numeric(n_per_sample[s]); names(v) <- as.character(ctr); v })
#     names(ncells_list) <- uS
#
#     # SingleCellExperiment with one assay named by ctr
#     sce <- SingleCellExperiment(
#         assays  = setNames(list(M_by_sample), as.character(ctr)),
#         colData = S4Vectors::DataFrame(cd)
#     )
#     int_colData(sce) <- S4Vectors::DataFrame(n_cells = I(ncells_list))
#     sce
# }
