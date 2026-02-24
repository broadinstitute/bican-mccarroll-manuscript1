## ============================================================
## Task: Run GSEA (age) and write lightweight outputs
## ============================================================

## -----------------------
## Parameters
## -----------------------

ct_file <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/sburger_tmp/cell_types_use.txt"
gene_to_chr_file <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/sburger_tmp/gene_to_chromosome.txt"

de_results_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results"
de_dir <- base::paste0(de_results_dir, "/LEVEL_3/sex_age/cell_type")

test <- "age"

fgsea_cell_types <- c(
    "MSN_D1_matrix",
    "glutamatergic_L23IT",
    "GABA_PTHLH-PVALB",
    "astrocyte",
    "OPC",
    "oligodendrocyte",
    "microglia"
)

gmt_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/sburger_tmp/gmt_files/C5"
gmt_files <- sort(list.files(gmt_dir, full.names = TRUE))

minSize <- 15
maxSize <- 500

gsea_lite_dir <- "/downloads/tmp/de_gsea_new"
padj_thresh <- 0.05
gmt_pattern <- "c5.go"

## Optional caching
#save_rds <- TRUE
#gsea_rds <- base::file.path(gsea_lite_dir, "bican_age_de_gsea_results.rds")

## -----------------------
## Execution
## -----------------------

if (!base::dir.exists(gsea_lite_dir)) {
    dir.create(gsea_lite_dir, recursive = TRUE)
}

gene_to_chr <- read_gene_to_chr(gene_to_chr_file)
de_age <- read_de_results(de_dir, test, ct_file, gene_to_chr)

gsea_results <- run_gsea(
    de_dt = de_age,
    fgsea_cell_types = fgsea_cell_types,
    gmt_files = gmt_files
)

write_gsea_lite(
    gsea_results = gsea_results,
    out_dir = gsea_lite_dir,
    padj_thresh = padj_thresh,
    gmt_pattern = gmt_pattern
)
