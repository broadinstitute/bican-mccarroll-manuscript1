# Produce a file with start_distance (variant pos - TSS, in bp) for each
# gene-variant pair in the index_snp_slope_matrix.
#
# Reads start_distance directly from the tensorQTL index files, then
# subsets to the ~9,899 index SNPs. No need to re-run get_egene_union_pairs.
#
# Usage:
#   Rscript run_get_index_snp_start_distance.R <eqtl_dir> <region_cell_type_path> <index_snp_matrix_path> <output_path>

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
eqtl_dir              <- args[1]
region_cell_type_path <- args[2]
index_snp_matrix_path <- args[3]
output_path           <- args[4]

# ---- 1. Read index SNP list ----
index_dt <- fread(index_snp_matrix_path, select = c("phenotype_id", "variant_id"))
message(sprintf("Index SNPs: %d gene-variant pairs", nrow(index_dt)))

# ---- 2. Read start_distance from tensorQTL index files ----
region_cell_type_dt <- fread(region_cell_type_path)

dist_list <- vector("list", nrow(region_cell_type_dt))
for (i in seq_len(nrow(region_cell_type_dt))) {
    ct  <- region_cell_type_dt$cell_type[i]
    reg <- region_cell_type_dt$region[i]
    subdir <- paste0(ct, "__", reg)
    eqtl_file <- file.path(eqtl_dir, subdir, paste0(subdir, ".cis_qtl.txt.gz"))

    message(sprintf("Reading: %s", subdir))
    dt <- fread(eqtl_file, select = c("phenotype_id", "variant_id", "start_distance"))
    dist_list[[i]] <- dt
}

all_dist <- rbindlist(dist_list)
all_dist <- unique(all_dist, by = c("phenotype_id", "variant_id"))
message(sprintf("Unique gene-variant pairs with start_distance: %d", nrow(all_dist)))

# ---- 3. Join: index SNPs + start_distance ----
result <- merge(index_dt, all_dist, by = c("phenotype_id", "variant_id"), all.x = TRUE)

n_matched <- sum(!is.na(result$start_distance))
n_missing <- sum(is.na(result$start_distance))
message(sprintf("Matched: %d / %d (missing: %d)", n_matched, nrow(result), n_missing))

# ---- 4. Write output ----
fwrite(result, output_path, sep = "\t")
message(sprintf("Written to: %s", output_path))
