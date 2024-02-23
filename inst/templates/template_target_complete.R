library(targets)
library(tarchetypes)
library(moiraine)

## add any package you need to use in the pipeline here
tar_option_set(
  packages = c(
    "moiraine",
    "MOFA2",
    "mixOmics",
    "readr",
    "tibble",
    "tidyr",
    "dplyr",
    "ggplot2",
    "patchwork"
  )
)

list(

  # Data import ----------------------------------------------------------------

  ## Data import using a target factory
  import_dataset_csv_factory(
    files = c(
      system.file("extdata/genomics_dataset.csv", package = "moiraine"),
      system.file("extdata/transcriptomics_dataset.csv", package = "moiraine"),
      system.file("extdata/metabolomics_dataset.csv", package = "moiraine")
    ),
    col_ids = c("marker", "gene_id", "sample_id"),
    features_as_rowss = c(TRUE, TRUE, FALSE),
    target_name_suffixes = c("geno", "transcripto", "metabo")
  ),

  ## Genomics features metadata file
  tar_target(
    fmetadata_file_geno,
    system.file("extdata/genomics_features_info.csv", package = "moiraine"),
    format = "file"
  ),

  ## Genomics features metadata import
  tar_target(
    fmetadata_geno,
    import_fmetadata_csv(
      fmetadata_file_geno,
      col_id = "marker",
      col_types = c("chromosome" = "c")
    )
  ),


  ## Metabolomics features metadata import
  import_fmetadata_csv_factory(
    files = c(
      system.file("extdata/metabolomics_features_info.csv", package = "moiraine")
    ),
    col_ids = c("feature_id"),
    target_name_suffixes = c("metabo")
  ),

  ## Transcriptomics features metadata import
  import_fmetadata_gff_factory(
    files = system.file("extdata/bos_taurus_gene_model.gff3", package = "moiraine"),
    feature_types = "genes",
    add_fieldss = c("Name", "description"),
    target_name_suffixes = "transcripto"
  ),

  ## Samples metadata import
  import_smetadata_csv_factory(
    files = system.file("extdata/samples_info.csv", package = "moiraine"),
    col_ids = "animal_id",
    target_name_suffixes = "all"
  ),

  ## Creating omics sets for each dataset
  create_omics_set_factory(
    datasets = c(data_geno, data_transcripto, data_metabo),
    omics_types = c("genomics", "transcriptomics", "metabolomics"),
    features_metadatas = c(fmetadata_geno, fmetadata_transcripto, fmetadata_metabo),
    samples_metadatas = c(smetadata_all, smetadata_all, smetadata_all)
  ),

  ## Creating the MultiDataSet object
  tar_target(
    mo_set,
    create_multiomics_set(
      list(set_geno,
           set_transcripto,
           set_metabo)
    )
  ),

  # Inspecting datasets --------------------------------------------------------

  ## Creating a density plot for each dataset
  tar_target(
    density_plots,
    plot_density_data(
      mo_set,
      combined = FALSE,
      scales = "free"
    )
  ),

  ## Plotting the relationship between features mean and standard deviation
  ## for each dataset
  tar_target(
    mean_sd_plots,
    plot_meansd_data(mo_set)
  ),

  ## Assessing missing values
  tar_target(
    n_missing_values,
    check_missing_values(mo_set)
  ),

  # Data pre-processing --------------------------------------------------------
  ## Applying transformations to the datasets
  transformation_datasets_factory(
    mo_set_de,
    c("rnaseq" = "vst-deseq2",
      "metabolome" = "logx"),
    log_bases = 2,
    pre_log_functions = zero_to_half_min(),
    standardize = FALSE,
    transformed_data_name = "mo_set_transformed"
  ),

  ## Density plot for each transformed dataset
  tar_target(
    density_plots_transformed,
    plot_density_data(
      mo_set_transformed,
      combined = FALSE,
      scales = "free"
    )
  ),

  ## Plotting the mean-SD trend for transformed each dataset
  tar_target(
    mean_sd_plots_transformed,
    plot_meansd_data(mo_set_transformed)
  ),

  ## Summary table of the transformations applied
  tar_target(
    transformation_summary,
    get_table_transformations(transformations_runs_list)
  ),

  ## Running a PCA on each dataset
  pca_complete_data_factory(
    mo_set_transformed,
    complete_data_name = "mo_set_complete"
  ),

  ## PCA screeplots
  tar_target(
    pca_screeplots,
    plot_screeplot_pca(pca_runs_list)
  ),

  ## PCA sample plots
  tar_target(
    pca_sample_plots,
    plot_samples_coordinates_pca(
      pca_runs_list,
      datasets = "snps",
      pcs = 1:3,
      mo_data = mo_set_de,
      colour_upper = "geno_comp_cluster",
      shape_upper = "status",
      colour_lower = "feedlot"
    )
  ),

  # Data pre-filtering ---------------------------------------------------------
  ## Unsupervised feature selection based on MAD score
  feature_preselection_mad_factory(
    mo_set_complete,
    to_keep_ns = c("snps" = 1000, "rnaseq" = 1000),
    with_ties = TRUE,
    filtered_set_target_name = "mo_presel_unsupervised"
  ),

  ## Diagnostic plot for MAD-based feature selection
  tar_target(
    preselection_mad_plot,
    plot_feature_preselection_mad(individual_mad_values)
  ),

  ## Supervised feature selection based on bruising groups
  feature_preselection_splsda_factory(
    mo_set_complete,
    group = "status",
    to_keep_ns = c("snps" = 1000, "rnaseq" = 1000),
    filtered_set_target_name = "mo_presel_supervised"
  ),

  ## Diagnostic plot for sPLS-DA based feature selection
  tar_target(
    preselection_splsda_plot,
    plot_feature_preselection_splsda(individual_splsda_perf)
  ),

  # sPLS integration -----------------------------------------------------------
  ## Creating sPLS input
  tar_target(
    spls_input,
    get_input_spls(
      mo_presel_supervised,
      mode = "canonical",
      datasets = c("rnaseq", "metabolome")
    )
  ),

  ## Initial PLS run with no feature selection and large number of components
  tar_target(
    spls_novarsel,
    spls_run(
      spls_input,
      ncomp = 4
    )
  ),

  ## Cross-validation for number of components
  tar_target(
    spls_perf_res,
    mixOmics::perf(
      spls_novarsel,
      validation = "Mfold",
      folds = 10,
      nrepeat = 10,
      cpus = 3
    )
  ),

  ## Plotting cross-validation results (for number of components)
  ## Can try criterion = 'Q2.total', 'cor.tpred', 'cor.upred', 'RSS.tpred',
  ## 'RSS.upred' (but avoid 'RSS' and 'PRESS')
  tar_target(
    spls_perf_plot,
    plot(spls_perf_res, criterion = "Q2.total")
  ),

  ## Selected value for ncomp
  tar_target(
    spls_optim_ncomp,
    spls_get_optim_ncomp(spls_perf_res, min_ncomp = 2)
  ),

  ## Cross-validation for number of features to retain
  tar_target(
    spls_tune_res,
    spls_tune(
      spls_input,
      ncomp = spls_optim_ncomp,
      keepX = seq(10, 100, 10),
      keepY = seq(10, 100, 10),
      validation = "Mfold",
      folds = 10,
      nrepeat = 5,
      measure = "cor",
      cpus = 3
    )
  ),

  ## Plotting cross-validation results (for number of features)
  tar_target(
    spls_tune_plot,
    spls_plot_tune(spls_tune_res)
  ),

  ## Final sPLS run
  tar_target(
    spls_final_run,
    spls_run(
      spls_input,
      ncomp = spls_optim_ncomp,
      keepX = spls_tune_res$choice.keepX,
      keepY = spls_tune_res$choice.keepY
    )
  ),

  ## Generating standardised output
  tar_target(
    spls_output,
    get_output(spls_final_run)
  ),

  ## Percentage of variance explained
  tar_target(
    spls_plot_variance_explained,
    plot_variance_explained(spls_output)
  ),

  ## Samples score matrix plot
  tar_target(
    spls_samples_scores_plot,
    plot_samples_score(
      spls_output,
      mo_data = mo_set_complete,
      colour_upper = "status",
      scale_colour_upper = palette_status,
      shape_upper = "gender",
      colour_lower = "feedlot",
      scale_colour_lower = palette_feedlot
    ) +
      theme(legend.box = "vertical")
  ),

  ## Distribution of features weight
  tar_target(
    spls_features_weight_distribution,
    plot_features_weight_distr(spls_output)
  ),

  ## Plot of top contributing features
  tar_target(
    spls_top_features_plot,
    plot_top_features(
      spls_output,
      mo_data = mo_set_complete,
      label_cols = list(
        "rnaseq" = "Name",
        "metabolome" = "name"
      )
    )
  ),

  ## Table of selected features
  tar_target(
    spls_selected_features,
    get_selected_features(spls_output)
  ),

  # sO2PLS integration ---------------------------------------------------------

  ## Creating sO2PLS input
  tar_target(
    omicspls_input,
    get_input_omicspls(
      mo_presel_supervised,
      datasets = c("rnaseq", "metabolome")
    )
  ),

  ## Adjusted cross-validation for number of components
  tar_target(
    so2pls_cv_adj,
    so2pls_crossval_o2m_adjR2(
      omicspls_input,
      a = 1:5,
      ax = seq(0, 10, by = 2),
      ay = seq(0, 10, by = 2),
      nr_folds = 10,
      nr_cores = 6,
      seed = 127
    )
  ),
  tar_target(
    so2pls_cv_adj_res,
    so2pls_get_optim_ncomp_adj(so2pls_cv_adj)
  ),

  ## Plotting adjusted cross-validation results
  tar_target(
    so2pls_cv_adj_plot,
    so2pls_plot_cv_adj(so2pls_cv_adj)
  ),

  ## Standard cross-validation for number of components
  tar_target(
    so2pls_cv,
    so2pls_crossval_o2m(
      omicspls_input,
      so2pls_cv_adj,
      nr_folds = 10,
      nr_cores = 6,
      seed = 356
    )
  ),
  tar_target(
    so2pls_cv_res,
    so2pls_get_optim_ncomp(so2pls_cv)
  ),

  ## Plotting standard cross-validation results
  tar_target(
    so2pls_cv_plot,
    so2pls_plot_cv(so2pls_cv)
  ),

  ## Cross-validation for sparsity parameters
  tar_target(
    so2pls_cv_sparsity,
    so2pls_crossval_sparsity(
      omicspls_input,
      n = so2pls_cv_res["n"],
      nx = so2pls_cv_res["nx"],
      ny = so2pls_cv_res["ny"],
      nr_folds = 10,
      keepx_seq = c(seq(5, 30, 5), seq(40, 100, 10)),
      keepy_seq = c(seq(5, 40, 5))
    )
  ),
  tar_target(
    so2pls_cv_sparsity_res,
    so2pls_get_optim_keep(so2pls_cv_sparsity)
  ),

  ## Plotting the results of the cross-validation for the number of features
  ## to retain from each dataset for the different joint components
  tar_target(
    so2pls_cv_sparsity_plot,
    so2pls_plot_cv_sparsity(so2pls_cv_sparsity)
  ),

  ## Extracting sparsity results in table format
  tar_target(
    so2pls_cv_sparsity_table,
    so2pls_print_cv_sparsity(so2pls_cv_sparsity_res)
  ),

  ## Final sO2PLS run
  tar_target(
    so2pls_final_run,
    so2pls_o2m(
      omicspls_input,
      so2pls_cv_res,
      so2pls_cv_sparsity_res
    )
  ),

  ## Summary plot of percentage of variance explained
  tar_target(
    so2pls_summary_plot,
    so2pls_plot_summary(so2pls_final_run)
  ),

  ## Screeplot
  tar_target(
    so2pls_screeplot,
    so2pls_screeplot(so2pls_final_run)
  ),

  ## Comparison of samples score for joint components
  tar_target(
    so2pls_joint_components_comparison_plot,
    so2pls_compare_samples_joint_components(
      so2pls_final_run,
      mo_data = mo_set_de,
      colour_by = "status",
      shape_by = "feedlot"
    )
  ),

  ## Coefficient plot for joint components
  tar_target(
    so2pls_joint_components_coefficients_plot,
    so2pls_plot_joint_components_coefficients(so2pls_final_run)
  ),

  ## Joint component samples score plot
  tar_target(
    so2pls_joint_components_samples_score_plot,
    so2pls_plot_samples_joint_components(
      so2pls_final_run,
      mo_data = mo_set_de,
      colour_upper = "status",
      scale_colour_upper = scale_colour_brewer(palette = "Paired"),
      shape_upper = "feedlot"
    ) +
      theme(legend.box = "vertical")
  ),

  ## Specific components samples score plot
  tar_target(
    so2pls_specific_components_samples_score_plot,
    so2pls_plot_samples_specific_components(
      so2pls_final_run,
      mo_data = mo_set_de,
      colour_upper = "feedlot",
      scale_colour_upper = scale_colour_brewer(palette = "Paired"),
      colour_lower = "rnaseq_batch",
      shape_upper = "gender"
    ) |>
      map(\(x) x + theme(legend.box = "vertical"))
  ),

  ## Generating standardised output
  tar_target(
    so2pls_output,
    get_output(so2pls_final_run)
  ),

  ## Percentage of variance explained
  tar_target(
    so2pls_plot_variance_explained,
    plot_variance_explained(so2pls_output, ncol = 1) +
      theme(axis.text.x = element_text(size = 9, angle = 30, hjust = 1))
  ),

  ## Samples score matrix plot
  tar_target(
    so2pls_samples_scores_plot,
    plot_samples_score(
      so2pls_output,
      latent_dimensions = "joint component 1",
      mo_data = mo_set_complete,
      colour_upper = "status",
      scale_colour_upper = palette_status,
      shape_upper = "gender"
    )
  ),

  ## Distribution of features weight
  tar_target(
    so2pls_features_weight_distribution,
    plot_features_weight_distr(so2pls_output) +
      plot_layout(ncol = 2)
  ),

  ## Plot of top contributing features
  tar_target(
    so2pls_top_features_plot,
    plot_top_features(
      so2pls_output,
      mo_data = mo_set_complete,
      label_cols = list(
        "rnaseq" = "Name",
        "metabolome" = "name"
      )
    )
  ),

  ## Table of selected features
  tar_target(
    so2pls_selected_features,
    get_selected_features(
      so2pls_output,
      latent_dimensions = "joint component 1"
    )
  ),

  # MOFA integration -----------------------------------------------------------

  ## Creating MOFA input
  tar_target(
    mofa_input,
    get_input_mofa(
      mo_presel_supervised,
      options_list = list(
        data_options = list(scale_views = TRUE),
        model_options = list(likelihoods = c(
          "snps" = "poisson",
          "rnaseq" = "gaussian",
          "metabolome" = "gaussian")
        ),
        training_options = list(seed = 43)
      ),
      only_common_samples = FALSE
    )
  ),

  ## Overview plot of the samples in each dataset
  tar_target(
    mofa_input_plot,
    plot_data_overview(mofa_input)
  ),

  ## Training MOFA model
  tar_target(
    mofa_trained,
    run_mofa(
      mofa_input,
      save_data = TRUE,
      use_basilisk = TRUE
    )
  ),

  ## Formatting MOFA output
  tar_target(
    mofa_output,
    get_output(mofa_trained)
  ),

  ## Plots of variance explained
  tar_target(
    mofa_var_explained_plot,
    plot_variance_explained(
      mofa_trained,
      x = "view",  ## datasets on the x-axis
      y = "factor" ## factors on the y-axis
    )
  ),
  tar_target(
    mofa_total_var_explained_plot,
    plot_variance_explained(
      mofa_trained,
      x = "view",
      y = "factor",
      plot_total = TRUE
    )[[2]]
  ),

  ## Plot of factors correlation with covariates
  tar_target(
    mofa_factors_covariates_cor_plot,
    mofa_plot_cor_covariates(mofa_trained)
  ),

  ## Generating standardised output
  tar_target(
    mofa_output,
    get_output(mofa_trained)
  ),

  ## Percentage of variance explained
  tar_target(
    mofa_plot_variance_explained,
    plot_variance_explained(mofa_output, ncol = 1) +
      theme(axis.text.x = element_text(size = 9, angle = 30, hjust = 1))
  ),

  ## Samples score matrix plot
  tar_target(
    mofa_samples_scores_plot,
    plot_samples_score(
      mofa_output,
      latent_dimensions = paste("Factor", 1:4),
      mo_data = mo_set_complete,
      colour_upper = "status",
      scale_colour_upper = palette_status,
      shape_upper = "gender",
      colour_lower = "geno_comp_cluster",
      scale_colour_lower = palette_geno_comp
    ) +
      theme(legend.box = "vertical")
  ),

  ## Distribution of features weight
  tar_target(
    mofa_features_weight_distribution,
    plot_features_weight_distr(
      mofa_output,
      latent_dimensions = paste("Factor", 1:4)
    ) +
      plot_layout(ncol = 1)
  ),

  ## Plot of top contributing features
  tar_target(
    mofa_top_features_plot,
    plot_top_features(
      mofa_output,
      mo_data = mo_set_complete,
      label_cols = list(
        "rnaseq" = "Name",
        "metabolome" = "name"
      ),
      latent_dimensions = paste("Factor", 1:2)
    )
  ),

  ## Table of top contributing features
  tar_target(
    mofa_top_features,
    get_top_features(
      mofa_output,
      min_importance = 0.8,
      mo_data = mo_set_complete
    )
  ),

  # DIABLO integration ---------------------------------------------------------

  ## Creating the DIABLO input
  tar_target(
    diablo_input,
    get_input_mixomics_supervised(
      mo_presel_supervised,
      group = "status"
    )
  ),

  ## Running sPLS on each dataset to construct the design matrix
  diablo_pairwise_pls_factory(diablo_input),

  ## Initial DIABLO run with no feature selection and large number of components
  tar_target(
    diablo_novarsel,
    diablo_run(
      diablo_input,
      diablo_design_matrix,
      ncomp = 7
    )
  ),

  ## Cross-validation for number of components
  tar_target(
    diablo_perf_res,
    mixOmics::perf(
      diablo_novarsel,
      validation = "Mfold",
      folds = 10,
      nrepeat = 10,
      cpus = 3
    )
  ),

  ## Plotting cross-validation results (for number of components)
  tar_target(
    diablo_perf_plot,
    diablo_plot_perf(diablo_perf_res)
  ),

  ## Selected value for ncomp
  tar_target(
    diablo_optim_ncomp,
    diablo_get_optim_ncomp(diablo_perf_res)
  ),

  ## Cross-validation for number of features to retain
  tar_target(
    diablo_tune_res,
    diablo_tune(
      diablo_input,
      diablo_design_matrix,
      ncomp = diablo_optim_ncomp,
      validation = "Mfold",
      folds = 10,
      nrepeat = 5,
      dist = "centroids.dist",
      cpus = 3
    )
  ),

  ## Plotting cross-validation results (for number of features)
  tar_target(
    diablo_tune_plot,
    diablo_plot_tune(diablo_tune_res)
  ),

  ## Final DIABLO run
  tar_target(
    diablo_final_run,
    diablo_run(
      diablo_input,
      diablo_design_matrix,
      ncomp = diablo_optim_ncomp,
      keepX = diablo_tune_res$choice.keepX
    )
  ),

  ## Generating standardised output
  tar_target(
    diablo_output,
    get_output(diablo_final_run)
  ),

  ## Percentage of variance explained
  tar_target(
    diablo_plot_variance_explained,
    plot_variance_explained(diablo_output, ncol = 2)
  ),

  ## Samples score matrix plot
  tar_target(
    diablo_samples_scores_plot,
    plot_samples_score(
      diablo_output,
      mo_data = mo_set_complete,
      colour_upper = "status",
      scale_colour_upper = palette_status,
      shape_upper = "gender",
      colour_lower = "rnaseq_batch",
      scale_colour_lower = palette_rnaseq_batch
    ) +
      theme(legend.box = "vertical")
  ),

  ## Distribution of features weight
  tar_target(
    diablo_features_weight_distribution,
    plot_features_weight_distr(
      diablo_output
    ) +
      plot_layout(ncol = 1)
  ),

  ## Plot of top contributing features
  tar_target(
    diablo_top_features_plot,
    plot_top_features(
      diablo_output,
      mo_data = mo_set_complete,
      label_cols = list(
        "rnaseq" = "Name",
        "metabolome" = "name"
      ),
      latent_dimensions = paste("Component", 1:2)
    )
  ),

  ## Table of top contributing features
  tar_target(
    diablo_selected_features,
    get_selected_features(diablo_output)
  ),

  # Results evaluation ---------------------------------------------------------

  ## Evaluating DIABLO selected features against single-omics results
  tar_target(
    diablo_selected_vs_singleomics_table,
    evaluate_feature_selection_table(
      diablo_output,
      mo_data = mo_set_complete,
      col_names = list(
        "snps" = "qtl_type",
        "rnaseq" = "de_signif",
        "metabolome" = "de_signif"
      )
    )
  ),

  ## Plotting DIABLO features weight against single-omics results
  tar_target(
    diablo_features_weight_vs_singleomics_plot,
    plot_features_weight_covariate(
      diablo_output,
      mo_data = mo_set_complete,
      covariate = list(
        "snps" = "qtl_type",
        "rnaseq" = "de_status",
        "metabolome" = "de_status"
      ),
      remove_null_weight = TRUE
    )
  ),

  ## Genes GO annotation file
  tar_target(
    rnaseq_go_terms_file,
    system.file(
      "extdata/transcriptomics_go_annotation.csv",
      package = "moiraine"
    ),
    format = "file"
  ),

  ## Genes GO annotation data-frame
  tar_target(
    rnaseq_go_df,
    read_csv(rnaseq_go_terms_file) |>
      filter(go_domain == "Biological process")
  ),

  ## GO term sets
  tar_target(
    go_sets,
    make_feature_sets_from_df(
      rnaseq_go_df,
      col_id = "gene_id",
      col_set = "go_id"
    )
  ),

  ## Filtering GO term sets against measured features
  tar_target(
    go_sets_filtered,
    reduce_feature_sets_data(go_sets, mo_set_complete)
  ),

  ## Checking genes GO term sets against datasets
  tar_target(
    go_sets_check,
    check_feature_sets(
      go_sets_filtered,
      mo_set_complete,
      datasets = "rnaseq"
    )
  ),

  ## Table of information about GO terms
  tar_target(
    go_sets_info,
    rnaseq_go_df |>
      dplyr::select(go_id, go_name) |>
      dplyr::distinct()
  ),

  ## MOFA latent components enrichment analysis
  tar_target(
    mofa_enrichment_results,
    evaluate_method_enrichment(
      mofa_output,
      go_sets_filtered,
      datasets = "rnaseq",
      latent_dimensions = paste("Factor", 1:3),
      use_abs = TRUE,
      min_set_size = 10,
      add_missing_features = TRUE,
      mo_data = mo_set_complete,
      sets_info_df = go_sets_info,
      col_set = "go_id"
    )
  ),

  ## Plotting features weight for GO term 'GO:0031424'
  tar_target(
    mofa_enrichment_go0031424_plot,
    plot_features_weight_set(
      mofa_output,
      go_sets_filtered[["GO:0031424"]],
      set_name = "GO:0031424 (keratinization)",
      features_metric = "importance",
      datasets = "rnaseq",
      latent_dimensions = paste("Factor", c(1, 3)),
      point_alpha = 0.2
    )
  ),

  ## Assessing DIABLO samples clustering
  tar_target(
    diablo_silhouette,
    compute_samples_silhouette(
      diablo_output,
      mo_set_complete,
      "status"
    )
  ),

  # Results comparison ---------------------------------------------------------
  ## Creating a list of integration methods output objects
  tar_target(
    output_list,
    list(spls_output, so2pls_output, mofa_output, diablo_output)
  ),

  ## Heatmap for comparison of integration methods output
  tar_target(
    comparison_methods_heatmap_plot,
    comparison_heatmap_corr(output_list)
  ),

  ## Correlation plot for comparison of sPLS and sO2PLS results
  tar_target(
    mofa_so2pls_correlation_plot,
    comparison_plot_correlation(output_list[2:1])
  ),

  ## Comparison of samples score for MOFA factor 1 and DIABLO component 1
  tar_target(
    mofa_so2pls_samples_score_comparison_plot,
    plot_samples_score_pair(
      output_list[3:4],
      list("MOFA" = "Factor 1", "DIABLO" = "Component 1"),
      mo_data = mo_set_complete,
      colour_by = "status"
    ) +
      scale_colour_brewer(palette = "Set1")
  ),

  ## Comparison of features weight for MOFA factor 1 and sO2PLS joint component 1
  tar_target(
    mofa_so2pls_features_weight_comparison_plot,
    plot_features_weight_pair(
      output_list[3:4],
      list("MOFA" = "Factor 1", "DIABLO" = "Component 1"),
      mo_data = mo_set_complete,
      label_cols = list(
        "rnaseq" = "Name",
        "metabolome" = "name"
      )
    )
  ),

  ## Table of features' consensus importance for MOFA factor 1 and sO2PLS joint
  ## component 1
  tar_target(
    mofa_so2pls_features_weight_consensus_importance,
    compute_consensus_importance(
      output_list[3:4],
      list("MOFA" = "Factor 1", "DIABLO" = "Component 1")
    ) |>
      left_join(
        get_features_labels(
          mo_set,
          list("rnaseq" = "Name",
               "metabolome" = "name")
        ),
        by = c("dataset", "feature_id")
      )
  )
)
