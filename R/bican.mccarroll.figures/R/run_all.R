#' This runs all of the plotting functions in the package with the same configured paths.
#' This is used internally to generate all figures, but is not a public method
#' because it hard-codes the options.
run_all <- function() {
  options(
    bican.mccarroll.figures.data_root_dir =
      "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis",
    bican.mccarroll.figures.out_dir =
      "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/figure_repository",
    bican.mccarroll.figures.cache_dir =
      "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/figure_repository/data_cache"
  )

  bican.mccarroll.figures::plot_sample_covariate_correlations()
  bican.mccarroll.figures::plot_eqtl_filtering_trajectories()
  bican.mccarroll.figures::plot_eqtl_filtering_examples()
  bican.mccarroll.figures::plot_de_filtering_trajectories()
  bican.mccarroll.figures::plot_de_filtering_examples()
  bican.mccarroll.figures::age_prediction_mean_residual_correlation_plots()
  bican.mccarroll.figures::age_prediction_error_plots()
  bican.mccarroll.figures::age_prediction_residual_corr_and_jaccard_heatmaps_region()
  bican.mccarroll.figures::age_prediction_residual_corr_and_jaccard_heatmaps_region(region = "DFC")
  bican.mccarroll.figures::age_prediction_residual_corr_and_jaccard_heatmaps_cell_type()
  bican.mccarroll.figures::age_prediction_corrected_residual_pairwise_scatter_region()
  bican.mccarroll.figures::age_prediction_uncorrected_residual_pairwise_scatter_region()
  bican.mccarroll.figures::age_prediction_examples()
  bican.mccarroll.figures::plot_trade_analysis_bican()
  bican.mccarroll.figures::plot_trade_analysis_bican_sea_ad_mtg_mmc()
  bican.mccarroll.figures::plot_trade_analysis_PMID_39402379()
  bican.mccarroll.figures::plot_kmeans_age()
  bican.mccarroll.figures::plot_ancestry_pca()

  bican.mccarroll.figures::plot_donor_gex_vs_age()
  bican.mccarroll.figures::plot_de_cor_heatmaps_age()
  bican.mccarroll.figures::plot_de_cor_heatmap_bican_sea_ad_vs_pmid_39402379()

  bican.mccarroll.figures::de_sex_age_scatter_plots()
  bican.mccarroll.figures::plot_de_volcano()
  bican.mccarroll.figures::plot_sex_de_by_chromosome_bican()
  bican.mccarroll.figures::plot_de_summary_barplot_bican_sea_ad_mtg_mmc()
  bican.mccarroll.figures::plot_de_summary_barplot_pmid_39402379()

  bican.mccarroll.figures::plot_de_sign_test_snap200_age_bican()
  bican.mccarroll.figures::plot_de_sign_test_pmid_40903571_bican()
  bican.mccarroll.figures::plot_de_sign_test_pmid_39227716_bican()
  bican.mccarroll.figures::plot_de_sign_test_pmid_39402379_bican()
  bican.mccarroll.figures::plot_de_sign_test_snap200_sex_bican()

  bican.mccarroll.figures::run_private_eqtl_tests_bican()
  bican.mccarroll.figures::summarize_private_eqtl_tests_bican()
  bican.mccarroll.figures::run_private_eqtl_tests_downsampled_bican()
  bican.mccarroll.figures::summarize_private_eqtl_tests_downsampled_bican()
}
