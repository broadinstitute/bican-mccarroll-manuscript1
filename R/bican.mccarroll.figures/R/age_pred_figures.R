

# Generation of raw data for all plots for age prediction.
# These results will be cached.
# This returns a list of dataframes that can be used for plotting.
run_age_prediction<-function (data_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/metacells/LEVEL_3",
                              data_name="donor_rxn_DGEList",
                              age_de_results_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results/LEVEL_3/sex_age/cell_type_region_interaction_absolute_effects",
                              result_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/age_prediction/age_prediction_region_alpha_0",
                              contig_yaml_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/metadata/GRCh38_ensembl_v43.contig_groups.yaml",
                              reduced_gtf_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/metadata/GRCh38_ensembl_v43.reduced.gtf.gz",
                              n_cores=12,
                              data_cache_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/figure_repository/data_cache"
                              ) {

    cache_dir <- file.path(data_cache_dir, "age_prediction_cache")
    if (!dir.exists(cache_dir)) {
        dir.create(cache_dir, recursive = TRUE)
    }

    if (dir.exists(cache_file)) {
        logger::log_info("Using cached data from {cache_file}")
        #read in all the files and return
    } else { # process the data as usual, write to the cache
        #clustering_min_genes and num_clusters don't do anything when outDir is NULL
        logger::log_info("No cached data from {cache_dir} regenerating data from sources.  This can take a while")
        bican.mccarroll.differentialexpression::predict_age_by_celltype_region(
            data_dir = data_dir,
            data_name = data_name,
            age_de_results_dir = age_de_results_dir,
            outPDFFile=NULL,
            result_dir = result_dir,
            contig_yaml_file = contig_yaml_file,
            reduced_gtf_file = reduced_gtf_file,
            optimize_alpha = FALSE,
            alpha_fixed = 0,
            n_cores = n_cores)
    }




}
