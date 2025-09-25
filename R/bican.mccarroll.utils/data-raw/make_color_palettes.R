# Create cell type palettes

# Load cell type to HBMA cell type mapping
cell_type2hbma_fn <- "data-raw/cell_type2hbma.txt"
cell_type2hbma_df <- read.table(cell_type2hbma_fn, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Load HBMA cell type to color mapping
hbma_cluser_annotation_fn <- "data-raw/hbma_cluster_annotation_term.csv"
hbma_cluser_annotation_df <- read.table(hbma_cluser_annotation_fn, header = TRUE, sep = ",", stringsAsFactors = FALSE, fill = TRUE, comment.char = "")

# Merge
cell_type2hbma_df <- merge(cell_type2hbma_df, hbma_cluser_annotation_df[,c("label", "color_hex_triplet")], all.x = TRUE, by.x = "hbma_label", by.y = "label")
cell_type2hbma_df <- cell_type2hbma_df[,c("cell_type", "color_hex_triplet")]

if (any(is.na(cell_type2hbma_df$color_hex_triplet))) {
  stop("ERROR: 1 or more cell types do not have a HBMA analog.")
}

# Create ggplot color and fill scales
hbma_cell_type_palette <- cell_type2hbma_df$color_hex_triplet
names(hbma_cell_type_palette) <- cell_type2hbma_df$cell_type

# Save
usethis::use_data(hbma_cell_type_palette, overwrite = TRUE)
