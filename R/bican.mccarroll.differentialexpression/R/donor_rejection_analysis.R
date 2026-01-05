
data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metacells"
data_name="donor_rxn_DGEList"
randVars=c("donor", "imputed_sex", "biobank", "single_cell_assay", "region", "hbcac_status", "toxicology_group")
fixedVars=c("age", "PC1", "PC2", "PC3", "PC4", "PC5", "pmi_hr", "pct_intronic", "frac_contamination")

additionalDonorMetadata=c("/broad/bican_um1_mccarroll/RNAseq/analysis/cellarium_upload/CAP_freeze_2/donor_metadata_rejection.txt")

classify_donor_rejections<-function () {
    logger::log_info(paste("Loading DGEList from:", data_dir, "with prefix:", data_name))
    dge=bican.mccarroll.differentialexpression::loadDGEList(data_dir, prefix = data_name)

    #merge in additional donor metadata
    if (!is.null(additionalDonorMetadata)) {
        logger::log_info(sprintf("Merging additional donor metadata from: %s", additionalDonorMetadata))

        donor_metadata <- utils::read.table(additionalDonorMetadata, header = TRUE, sep = "\t")

        if ("donor" %in% colnames(donor_metadata)) {
            m <- match(dge$samples$donor, donor_metadata$donor)
            for (cn in setdiff(names(donor_metadata), "donor")) {
                dge$samples[[cn]] <- donor_metadata[[cn]][m]
            }
        } else {
            logger::log_warn("No 'donor' column found in additional metadata; skipping merge.")
        }
    }
}
