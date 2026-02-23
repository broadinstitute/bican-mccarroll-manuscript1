


ctp_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/cell_type_proportions/LEVEL_1"
ctp_annotations <- data.table::fread(
  file.path(ctp_dir, "donor_region.annotation.cell_type_proportions.txt"),
)



# 
cell_types <- unique(ctp_annotations$annotation)
interneurons <- c(grep("GABA", cell_types, value = TRUE),
                  "striatal_cholinergic")

projection_neurons <- c(grep("MSN", cell_types, value = TRUE),
                        grep("IT", cell_types, value = TRUE),
                        "cortical_glutamatergic_L5ET",
                        "cortical_glutamatergic_L56NP",
                        "cortical_glutamatergic_L6",
                        "VTR-HTH_glut")
                        

setdiff(cell_types, c(interneurons, projection_neurons))
                        


regroup_cell_types <- function(long_ctp_df, cell_type_column, groupings){
  
}
print("testing!")
