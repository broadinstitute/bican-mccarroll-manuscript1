# Set up the metadata needed to compare the Kana et al. 2026 (bioRxiv) CaH
# CPS external dataset (converted by kana_2026_biorxiv.R in this directory)
# against BICAN's CaH region, following the same manifest/cell-type-list
# conventions used for the other external_comparison_* datasets.
#
# Unlike the other external comparisons (Gabitto 2024, PMID_39227716,
# PMID_40903571, SNAP200), Kana uses its own striatal cell type nomenclature
# rather than BICAN/SEA-AD names, so a small hand-curated name map is the
# single source of truth for which six cell types are compared and how they
# are paired. Everything else in this script is derived from that map.

deRoot <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression"

kanaVoomDir  <- file.path(deRoot, "external_comparison_kana_2026_biorxiv/voom-like")
kanaMetaDir  <- file.path(deRoot, "external_comparison_kana_2026_biorxiv/metadata")
bicanCaHDir  <- file.path(deRoot, "results/LEVEL_6/sex_age/cell_type_region_interaction_absolute_effects")
sharedMetaDir <- file.path(deRoot, "metadata")

# Kana (striatal) cell type -> BICAN cell type, for the six cell types this
# comparison is restricted to.
cell_type_map <- data.frame(
    kana = c(
        "STRd_D1_Matrix_MSN",
        "STRd_D1_Striosome_MSN",
        "STRd_D2_Matrix_MSN",
        "STRd_D2_Striosome_MSN",
        "STR_FS_PTHLH-PVALB_GABA",
        "STR_TAC3-PLPP4_GABA"
    ),
    bican = c(
        "SPN_D1_matrix",
        "SPN_D1_striosome",
        "SPN_D2_matrix",
        "SPN_D2_striosome",
        "GABA_PTHLH-PVALB",
        "GABA_TAC3-PLPP4"
    ),
    stringsAsFactors = FALSE
)

dir.create(kanaMetaDir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. Cell type name map
# ---------------------------------------------------------------------------

map_file <- file.path(kanaMetaDir, "kana_2026_bican_cell_type_map.tsv")
utils::write.table(cell_type_map, map_file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# ---------------------------------------------------------------------------
# 2. Sign-test overlap manifest (BICAN CaH age vs. Kana CaH ad_cps)
# ---------------------------------------------------------------------------

bican_file <- function(bican_ct) file.path(bicanCaHDir, paste0(bican_ct, "__CaH__age_DE_results.txt"))
kana_file  <- function(kana_ct)  file.path(kanaVoomDir, paste0(kana_ct, "__CaH__ad_cps_DE_results.txt"))

manifest <- data.frame(
    cell_type = cell_type_map$bican,
    bican = vapply(cell_type_map$bican, bican_file, character(1)),
    KANA_2026 = vapply(cell_type_map$kana, kana_file, character(1)),
    stringsAsFactors = FALSE
)

missing_files <- c(manifest$bican[!file.exists(manifest$bican)], manifest$KANA_2026[!file.exists(manifest$KANA_2026)])
if (length(missing_files) > 0L) {
    stop(
        "kana_2026_comparison_setup: missing DE file(s) referenced by the cell type map: ",
        paste(missing_files, collapse = ", "),
        call. = FALSE
    )
}

manifest_file <- file.path(kanaMetaDir, "kana_2026_bican_ad_cps_de_overlap_manifest.tsv")
utils::write.table(manifest, manifest_file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# ---------------------------------------------------------------------------
# 3. Cell type list files (shared differential_expression/metadata/, matching
#    where cell_types_for_sea_ad_mtg_mmc_plots.txt lives)
# ---------------------------------------------------------------------------

# The 6 Kana names (restricted, Kana-only figures).
kana_ct_file <- file.path(sharedMetaDir, "cell_types_for_kana_2026_plots.txt")
writeLines(cell_type_map$kana, kana_ct_file)

# All 66 Kana names present in voom-like/ (all-cell-types Kana-only figures).
kana_all_files <- list.files(kanaVoomDir, pattern = "__CaH__ad_cps_DE_results\\.txt$")
kana_all_cts <- sort(unique(sub("__CaH__ad_cps_DE_results\\.txt$", "", kana_all_files)))
if (length(kana_all_cts) == 0L) {
    stop("kana_2026_comparison_setup: no Kana voom-like DE files found in '", kanaVoomDir, "'.", call. = FALSE)
}
kana_all_ct_file <- file.path(sharedMetaDir, "cell_types_for_kana_2026_all_plots.txt")
writeLines(kana_all_cts, kana_all_ct_file)

# The 6 BICAN names (post-rename filtering/ordering for the BICAN-vs-Kana
# correlation heatmap and overlap enrichment figures).
bican_ct_file <- file.path(sharedMetaDir, "cell_types_for_kana_2026_bican_plots.txt")
writeLines(cell_type_map$bican, bican_ct_file)

logger::log_info(
    "Kana 2026 comparison setup complete: cell type map ({nrow(cell_type_map)} rows), ",
    "overlap manifest ({nrow(manifest)} rows), and 3 cell type list files written."
)
