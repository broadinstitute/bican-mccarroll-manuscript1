# Are donor predicted ages or their residuals correlated across cell types?
# If variance in prediction is random, then the residuals should not correlate.
# If there is a biological signal, then the residuals should correlate.

model_file_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/age_prediction"

compare_age_residuals<-function () {
    model_predictions=load_model_predictions(model_file_dir)
    res_mat <- compute_residual_matrix(results, cellType = cellType)

}


load_model_predictions<-function (model_file_dir) {
    model_files=list.files(model_file_dir, pattern="donor_age_predictions*", full.names = TRUE)

    #x=model_files[1]
    parseOne<-function (x) {
        a=read.table(x, header=TRUE, sep="\t", stringsAsFactors = FALSE)
        return (a)
    }

    all_models=lapply(model_files, parseOne)
    logger::log_info(paste("Loaded model predictions for "), length(all_models), " cell types")
    all_models=do.call(rbind, all_models)
    return(all_models)
}
