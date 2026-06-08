# Package index

## Project management

Getting started with a new data integration project

- [`create_moiraine_pipeline()`](https://plant-food-research-open.github.io/moiraine/reference/create_moiraine_pipeline.md)
  : Creates a target script file from template

## Data import

Importing data into R

- [`import_dataset_csv()`](https://plant-food-research-open.github.io/moiraine/reference/import_dataset_csv.md)
  : Import a dataset from a csv file
- [`import_fmetadata_csv()`](https://plant-food-research-open.github.io/moiraine/reference/import_fmetadata_csv.md)
  : Import feature metadata from a csv file
- [`import_fmetadata_gff()`](https://plant-food-research-open.github.io/moiraine/reference/import_fmetadata_gff.md)
  : Import features metadata from a GFF/GTF file
- [`import_smetadata_csv()`](https://plant-food-research-open.github.io/moiraine/reference/import_smetadata_csv.md)
  : Import samples metadata from a csv file
- [`import_dataset_csv_factory()`](https://plant-food-research-open.github.io/moiraine/reference/import_dataset_csv_factory.md)
  : Target factory for csv datasets import
- [`import_fmetadata_csv_factory()`](https://plant-food-research-open.github.io/moiraine/reference/import_fmetadata_csv_factory.md)
  : Target factory for csv features metadata import
- [`import_fmetadata_gff_factory()`](https://plant-food-research-open.github.io/moiraine/reference/import_fmetadata_gff_factory.md)
  : Target factory for GFF/GTF features metadata import
- [`import_smetadata_csv_factory()`](https://plant-food-research-open.github.io/moiraine/reference/import_smetadata_csv_factory.md)
  : Target factory for csv samples metadata import

## Omics and multi-omics sets creation

Creating omics sets

- [`create_omics_set()`](https://plant-food-research-open.github.io/moiraine/reference/create_omics_set.md)
  : Create a Biobase set object to store omics data
- [`create_omics_set_factory()`](https://plant-food-research-open.github.io/moiraine/reference/create_omics_set_factory.md)
  : Target factory for omics sets creation
- [`create_multiomics_set()`](https://plant-food-research-open.github.io/moiraine/reference/create_multiomics_set.md)
  : Create a MultiDataSet object to store multi-omics data
- [`add_omics_set()`](https://plant-food-research-open.github.io/moiraine/reference/add_omics_set.md)
  : Adds an omics set to a MultiDataSet object
- [`MetabolomeSet`](https://plant-food-research-open.github.io/moiraine/reference/MetabolomeSet-class.md)
  [`MetabolomeSet-class`](https://plant-food-research-open.github.io/moiraine/reference/MetabolomeSet-class.md)
  : Class to contain objects describing high-throughput metabolomics
  assays.
- [`PhenotypeSet`](https://plant-food-research-open.github.io/moiraine/reference/PhenotypeSet-class.md)
  [`PhenotypeSet-class`](https://plant-food-research-open.github.io/moiraine/reference/PhenotypeSet-class.md)
  : Class to contain objects describing phenotypic assays.
- [`add_metabo()`](https://plant-food-research-open.github.io/moiraine/reference/add_metabo-methods.md)
  : Method to add a MetabolomeSet to a MultiDataSet object.
- [`add_pheno()`](https://plant-food-research-open.github.io/moiraine/reference/add_pheno-methods.md)
  : Method to add a PhenotypeSet to a MultiDataSet object.

## Multi-omics sets (MultiDataSet objects)

### Querying

Querying the MultiDataSet object

- [`n_features()`](https://plant-food-research-open.github.io/moiraine/reference/n_features.md)
  : Number of features in each dataset of MultiDataSet object
- [`n_samples()`](https://plant-food-research-open.github.io/moiraine/reference/n_samples.md)
  : Number of samples in each dataset of MultiDataSet object
- [`get_features()`](https://plant-food-research-open.github.io/moiraine/reference/get_features.md)
  : Get feature IDs from MultiDataSet
- [`get_samples()`](https://plant-food-research-open.github.io/moiraine/reference/get_samples.md)
  : Get sample IDs from MultiDataSet
- [`get_datasets()`](https://plant-food-research-open.github.io/moiraine/reference/get_datasets.md)
  : Get multi-omics measurement datasets
- [`get_dataset_matrix()`](https://plant-food-research-open.github.io/moiraine/reference/get_dataset_matrix.md)
  : Get multi-omics dataset as matrix
- [`get_features_metadata()`](https://plant-food-research-open.github.io/moiraine/reference/get_features_metadata.md)
  : Get features metadata dataframes from MultiDataSet
- [`get_samples_metadata()`](https://plant-food-research-open.github.io/moiraine/reference/get_samples_metadata.md)
  : Get samples metadata dataframes from MultiDataSet
- [`get_samples_metadata_combined()`](https://plant-food-research-open.github.io/moiraine/reference/get_samples_metadata_combined.md)
  : Get combined samples metadata data-frame from MultiDataSet
- [`check_missing_values()`](https://plant-food-research-open.github.io/moiraine/reference/check_missing_values.md)
  : Check for missing values in MultiDataSet
- [`get_features_labels()`](https://plant-food-research-open.github.io/moiraine/reference/get_features_labels.md)
  : Get feature labels
- [`join_features_metadata()`](https://plant-food-research-open.github.io/moiraine/reference/join_features_metadata.md)
  : Join feature metadata to table
- [`join_samples_metadata()`](https://plant-food-research-open.github.io/moiraine/reference/join_samples_metadata.md)
  : Join samples metadata to table

### Plotting

Plotting properties of omics datasets

- [`plot_samples_upset()`](https://plant-food-research-open.github.io/moiraine/reference/plot_samples_upset.md)
  : Upset plot of samples
- [`plot_density_data()`](https://plant-food-research-open.github.io/moiraine/reference/plot_density_data.md)
  : Per-dataset density plot for MultiDataSet object
- [`plot_meansd_data()`](https://plant-food-research-open.github.io/moiraine/reference/plot_meansd_data.md)
  : Per-dataset mean-sd trend plots for MultiDataSet object
- [`plot_data_covariate()`](https://plant-food-research-open.github.io/moiraine/reference/plot_data_covariate.md)
  : Plots omics data vs sample covariate
- [`plot_data_heatmap()`](https://plant-food-research-open.github.io/moiraine/reference/plot_data_heatmap.md)
  : Plots omics data as heatmap

### Modifying

Modifying or subsetting datasets or metadata

- [`replace_dataset()`](https://plant-food-research-open.github.io/moiraine/reference/replace_dataset.md)
  : Replace matrix dataset within a MultiDataSet object
- [`round_dataset()`](https://plant-food-research-open.github.io/moiraine/reference/round_dataset.md)
  : Round values in omics dataset from MultiDataSet object
- [`add_features_metadata()`](https://plant-food-research-open.github.io/moiraine/reference/add_features_metadata.md)
  : Adding data-frame to features metadata
- [`add_samples_metadata()`](https://plant-food-research-open.github.io/moiraine/reference/add_samples_metadata.md)
  : Adding data-frame to samples metadata
- [`subset_features()`](https://plant-food-research-open.github.io/moiraine/reference/subset_features.md)
  : Subset a MultiDataSet object by feature

## PCA

PCA and missing values imputation on multi-omics set

- [`run_pca()`](https://plant-food-research-open.github.io/moiraine/reference/run_pca.md)
  : Run PCA on MultiDataSet
- [`run_pca_matrix()`](https://plant-food-research-open.github.io/moiraine/reference/run_pca_matrix.md)
  : Run PCA on matrix
- [`get_complete_data()`](https://plant-food-research-open.github.io/moiraine/reference/get_complete_data.md)
  : Get MultiDataSet object with imputed values
- [`pca_complete_data_factory()`](https://plant-food-research-open.github.io/moiraine/reference/pca_complete_data_factory.md)
  : Target factory for PCA run and missing values imputation on each
  omics dataset
- [`plot_screeplot_pca()`](https://plant-food-research-open.github.io/moiraine/reference/plot_screeplot_pca.md)
  : Screeplots for single-omics PCA
- [`plot_samples_coordinates_pca()`](https://plant-food-research-open.github.io/moiraine/reference/plot_samples_coordinates_pca.md)
  : Samples score plots for single-omics PCA
- [`get_pca_arguments()`](https://plant-food-research-open.github.io/moiraine/reference/get_pca_arguments.md)
  : Extract arguments used in PCA run

## Transformation

Omics datasets transformation and normalisation

- [`transform_bestNormalise_auto()`](https://plant-food-research-open.github.io/moiraine/reference/transform_bestNormalise_auto.md)
  : Applies the bestNormalize function to rows of a matrix
- [`transform_bestNormalise_manual()`](https://plant-food-research-open.github.io/moiraine/reference/transform_bestNormalise_manual.md)
  : Applies a normalisation method from bestNormalize to rows of a
  matrix
- [`transform_dataset()`](https://plant-food-research-open.github.io/moiraine/reference/transform_dataset.md)
  : Applies a transformation to a dataset from a MultiDataSet object
- [`transform_logx()`](https://plant-food-research-open.github.io/moiraine/reference/transform_logx.md)
  : Applies a log-x transformation to matrix
- [`transform_vsn()`](https://plant-food-research-open.github.io/moiraine/reference/transform_vsn.md)
  : Applies Variance Stabilising Normalisation (vsn) to matrix
- [`transform_vst()`](https://plant-food-research-open.github.io/moiraine/reference/transform_vst.md)
  : Applies Variance Stabilising Transformation (DESeq2) to matrix
- [`transformation_datasets_factory()`](https://plant-food-research-open.github.io/moiraine/reference/transformation_datasets_factory.md)
  : Target factory for datasets transformation
- [`get_transformed_data()`](https://plant-food-research-open.github.io/moiraine/reference/get_transformed_data.md)
  : Get MultiDataSet with transformed data
- [`get_table_transformations()`](https://plant-food-research-open.github.io/moiraine/reference/get_table_transformations.md)
  : Get table with transformation applied to each dataset
- [`zero_to_half_min()`](https://plant-food-research-open.github.io/moiraine/reference/zero_to_half_min.md)
  : Replace zeros with half-min in matrix

## Prefiltering

Features prefiltering for multi-omics set

### Unsupervised (Median Absolute Variation)

- [`select_features_mad()`](https://plant-food-research-open.github.io/moiraine/reference/select_features_mad.md)
  : Select features based on Median Absolute Deviation from MultiDataSet
- [`select_features_cov()`](https://plant-food-research-open.github.io/moiraine/reference/select_features_cov.md)
  : Select features based on Coefficient of Variation from MultiDataSet
- [`select_features_mad_matrix()`](https://plant-food-research-open.github.io/moiraine/reference/select_features_mad_matrix.md)
  : Select features based on Median Absolute Deviation from matrix
- [`select_features_cov_matrix()`](https://plant-food-research-open.github.io/moiraine/reference/select_features_cov_matrix.md)
  : Select features based on Coefficient of Variation from matrix
- [`get_filtered_dataset_variability()`](https://plant-food-research-open.github.io/moiraine/reference/get_filtered_dataset_variability.md)
  : Get filtered MultiDataSet object based on variability measure
- [`feature_preselection_mad_factory()`](https://plant-food-research-open.github.io/moiraine/reference/feature_preselection_mad_factory.md)
  : Target factory for feature preselection based on Median Absolute
  Deviation
- [`feature_preselection_cov_factory()`](https://plant-food-research-open.github.io/moiraine/reference/feature_preselection_cov_factory.md)
  : Target factory for feature preselection based on Coefficient of
  Variation
- [`plot_feature_preselection_mad()`](https://plant-food-research-open.github.io/moiraine/reference/plot_feature_preselection_mad.md)
  : Diagnostics plots for MAD-based feature preselection
- [`plot_feature_preselection_cov()`](https://plant-food-research-open.github.io/moiraine/reference/plot_feature_preselection_cov.md)
  : Diagnostics plots for COV-based feature preselection

### Supervised (sPLS-DA)

- [`get_input_splsda()`](https://plant-food-research-open.github.io/moiraine/reference/get_input_splsda.md)
  : Generate sPLS-DA input data (for mixomics)
- [`perf_splsda()`](https://plant-food-research-open.github.io/moiraine/reference/perf_splsda.md)
  : Assess optimal number of components for sPLS-DA on omics dataset
  from MultiDataSet object
- [`run_splsda()`](https://plant-food-research-open.github.io/moiraine/reference/run_splsda.md)
  : Performs sPLS-DA on omics dataset from MultiDataSet object
- [`get_filtered_dataset_splsda()`](https://plant-food-research-open.github.io/moiraine/reference/get_filtered_dataset_splsda.md)
  : Get filtered MultiDataSet object based on sPLS-DA runs
- [`feature_preselection_splsda_factory()`](https://plant-food-research-open.github.io/moiraine/reference/feature_preselection_splsda_factory.md)
  : Target factory for feature preselection based on sPLS-DA
- [`plot_feature_preselection_splsda()`](https://plant-food-research-open.github.io/moiraine/reference/plot_feature_preselection_splsda.md)
  : Diagnostics plots for sPLS-DA-based feature preselection

## Supervised integration

Integration of datasets aiming to discriminate samples based on an
outcome of interest

### DIABLO

Supervised integration with the DIABLO method from mixOmics

- [`get_input_mixomics_supervised()`](https://plant-food-research-open.github.io/moiraine/reference/get_input_mixomics_supervised.md)
  : Generate mixomics input data for supervised methods
- [`run_pairwise_pls()`](https://plant-food-research-open.github.io/moiraine/reference/run_pairwise_pls.md)
  : Pairwise PLS datasets comparison
- [`diablo_generate_design_matrix()`](https://plant-food-research-open.github.io/moiraine/reference/diablo_generate_design_matrix.md)
  : Generate DIABLO design matrix
- [`diablo_get_optim_ncomp()`](https://plant-food-research-open.github.io/moiraine/reference/diablo_get_optim_ncomp.md)
  : Get the optimal ncomp value
- [`diablo_get_pairwise_pls_corr()`](https://plant-food-research-open.github.io/moiraine/reference/diablo_get_pairwise_pls_corr.md)
  : Get pairwise correlations from PLS run
- [`diablo_get_params()`](https://plant-food-research-open.github.io/moiraine/reference/diablo_get_params.md)
  : Get parameters from DIABLO run
- [`diablo_get_wa_coord()`](https://plant-food-research-open.github.io/moiraine/reference/diablo_get_wa_coord.md)
  : Get weighted average coordinates
- [`diablo_pairwise_pls_factory()`](https://plant-food-research-open.github.io/moiraine/reference/diablo_pairwise_pls_factory.md)
  : Target factory for pairwise PLS and design matrix estimation for
  DIABLO run
- [`diablo_plot()`](https://plant-food-research-open.github.io/moiraine/reference/diablo_plot.md)
  : Plots DIABLO output
- [`diablo_plot_circos()`](https://plant-food-research-open.github.io/moiraine/reference/diablo_plot_circos.md)
  : Plots DIABLO circos plot
- [`diablo_plot_perf()`](https://plant-food-research-open.github.io/moiraine/reference/diablo_plot_perf.md)
  : Plots DIABLO perf results
- [`diablo_plot_tune()`](https://plant-food-research-open.github.io/moiraine/reference/diablo_plot_tune.md)
  : Plots DIABLO tune results
- [`diablo_plot_var()`](https://plant-food-research-open.github.io/moiraine/reference/diablo_plot_var.md)
  : Plots DIABLO features correlation circle
- [`diablo_predefined_design_matrix()`](https://plant-food-research-open.github.io/moiraine/reference/diablo_predefined_design_matrix.md)
  : Generate a design matrix for DIABLO
- [`diablo_run()`](https://plant-food-research-open.github.io/moiraine/reference/diablo_run.md)
  : Runs DIABLO algorithm
- [`diablo_table_optim_keepX()`](https://plant-food-research-open.github.io/moiraine/reference/diablo_table_optim_keepX.md)
  : Formatted table with optimal keepX values
- [`diablo_tune()`](https://plant-food-research-open.github.io/moiraine/reference/diablo_tune.md)
  : Tunes keepX arg for DIABLO

### MEFISTO

Supervised integration with the MEFISTO method from MOFA2 (for
time-series or spatially resolved data)

- [`get_input_mefisto()`](https://plant-food-research-open.github.io/moiraine/reference/get_input_mefisto.md)
  : Generate MEFISTO input data

## Unsupervised integration

Integration of datasets aiming at assessing variation common to the
datasets

### MOFA

Unsupervised integration with the MOFA method from MOFA2

- [`get_input_mofa2()`](https://plant-food-research-open.github.io/moiraine/reference/get_input_mofa2.md)
  : Generate input data for MOFA2 package
- [`get_input_mofa()`](https://plant-food-research-open.github.io/moiraine/reference/get_input_mofa.md)
  : Generate MOFA input data
- [`mofa_get_weights()`](https://plant-food-research-open.github.io/moiraine/reference/mofa_get_weights.md)
  : Get feature weights from MOFA object
- [`mofa_plot_cor_covariates()`](https://plant-food-research-open.github.io/moiraine/reference/mofa_plot_cor_covariates.md)
  : Plots the correlation between factors and covariates for MOFA

### sO2PLS

Unsupervised integration of 2 datasets with the sO2PLS method from
omicsPLS

- [`get_input_omicspls()`](https://plant-food-research-open.github.io/moiraine/reference/get_input_omicspls.md)
  : Generate omicsPLS input data
- [`so2pls_compare_samples_joint_components()`](https://plant-food-research-open.github.io/moiraine/reference/so2pls_compare_samples_joint_components.md)
  : Compares sO2PLS samples joint component scores between the two
  datasets
- [`so2pls_crossval_o2m()`](https://plant-food-research-open.github.io/moiraine/reference/so2pls_crossval_o2m.md)
  : Wrapper for OmicsPLS::crossval_o2m function
- [`so2pls_crossval_o2m_adjR2()`](https://plant-food-research-open.github.io/moiraine/reference/so2pls_crossval_o2m_adjR2.md)
  : Wrapper for OmicsPLS::crossval_o2m_adjR2 function
- [`so2pls_crossval_sparsity()`](https://plant-food-research-open.github.io/moiraine/reference/so2pls_crossval_sparsity.md)
  : Perform cross-validation to find the optimal number of
  features/groups to keep for each joint component for sO2PLS
- [`so2pls_get_components()`](https://plant-food-research-open.github.io/moiraine/reference/so2pls_get_components.md)
  : Get list of latent components from sO2PLS results
- [`so2pls_get_optim_keep()`](https://plant-food-research-open.github.io/moiraine/reference/so2pls_get_optim_keep.md)
  : Extract optimal number of features to keep from cross-validation
  results for sO2PLS
- [`so2pls_get_optim_ncomp()`](https://plant-food-research-open.github.io/moiraine/reference/so2pls_get_optim_ncomp.md)
  : Extract optimal number of components from cross-validation results
  for sO2PLS
- [`so2pls_get_optim_ncomp_adj()`](https://plant-food-research-open.github.io/moiraine/reference/so2pls_get_optim_ncomp_adj.md)
  : Extract optimal number of components from adjusted cross-validation
  results for sO2PLS
- [`so2pls_get_variance_explained()`](https://plant-food-research-open.github.io/moiraine/reference/so2pls_get_variance_explained.md)
  : Percentage of variance explained for sO2PLS
- [`so2pls_get_wa_coord()`](https://plant-food-research-open.github.io/moiraine/reference/so2pls_get_wa_coord.md)
  : Computes average sample coordinates for sO2PLS joint components
- [`so2pls_o2m()`](https://plant-food-research-open.github.io/moiraine/reference/so2pls_o2m.md)
  : Wrapper for OmicsPLS::o2m function
- [`so2pls_plot_cv()`](https://plant-food-research-open.github.io/moiraine/reference/so2pls_plot_cv.md)
  : Plots cross-validation results for sO2PLS
- [`so2pls_plot_cv_adj()`](https://plant-food-research-open.github.io/moiraine/reference/so2pls_plot_cv_adj.md)
  : Plot adjusted cross-validation results for sO2PLS
- [`so2pls_plot_cv_sparsity()`](https://plant-food-research-open.github.io/moiraine/reference/so2pls_plot_cv_sparsity.md)
  : Plot sparsity cross-validation results for sO2PLS
- [`so2pls_plot_joint_components_coefficients()`](https://plant-food-research-open.github.io/moiraine/reference/so2pls_plot_joint_components_coefficients.md)
  : Plots sO2PLS contributions between datasets joint components
- [`so2pls_plot_samples_joint_components()`](https://plant-food-research-open.github.io/moiraine/reference/so2pls_plot_samples_joint_components.md)
  : Plots sO2PLS joint components samples scores
- [`so2pls_plot_samples_specific_components()`](https://plant-food-research-open.github.io/moiraine/reference/so2pls_plot_samples_specific_components.md)
  : Plots sO2PLS specific components samples scores
- [`so2pls_plot_summary()`](https://plant-food-research-open.github.io/moiraine/reference/so2pls_plot_summary.md)
  : Plot summary of sO2PLS run
- [`so2pls_print_cv_adj()`](https://plant-food-research-open.github.io/moiraine/reference/so2pls_print_cv_adj.md)
  : Print adjusted cross-validation results for sO2PLS
- [`so2pls_print_cv_sparsity()`](https://plant-food-research-open.github.io/moiraine/reference/so2pls_print_cv_sparsity.md)
  : Print sparsity cross-validation results for sO2PLS
- [`so2pls_screeplot()`](https://plant-food-research-open.github.io/moiraine/reference/so2pls_screeplot.md)
  : Screeplot sO2PLS run

### sPLS

Unsupervised integration of 2 datasets with the sPLS method from
mixOmics

- [`get_input_spls()`](https://plant-food-research-open.github.io/moiraine/reference/get_input_spls.md)
  : Generate sPLS input data (for mixomics)
- [`get_input_mixomics_unsupervised()`](https://plant-food-research-open.github.io/moiraine/reference/get_input_mixomics_unsupervised.md)
  : Generate mixomics input data for unsupervised methods
- [`spls_get_optim_ncomp()`](https://plant-food-research-open.github.io/moiraine/reference/spls_get_optim_ncomp.md)
  : Select the optimal ncomp from sPLS cross-validation results
- [`spls_get_params()`](https://plant-food-research-open.github.io/moiraine/reference/spls_get_params.md)
  : Get parameters from sPLS run
- [`spls_get_wa_coord()`](https://plant-food-research-open.github.io/moiraine/reference/spls_get_wa_coord.md)
  : Computes average sample coordinates for sPLS components
- [`spls_plot_tune()`](https://plant-food-research-open.github.io/moiraine/reference/spls_plot_tune.md)
  : Displays results of sPLS tuning
- [`spls_plot_var()`](https://plant-food-research-open.github.io/moiraine/reference/spls_plot_var.md)
  : Plots sPLS features correlation circle
- [`spls_run()`](https://plant-food-research-open.github.io/moiraine/reference/spls_run.md)
  : Run sPLS algorithm
- [`spls_tune()`](https://plant-food-research-open.github.io/moiraine/reference/spls_tune.md)
  : Performs cross-validation for mixomics sPLS to select optimal keepX
  and keepY

## Standardised method output

Functions to get and query the results of an integration method as a
standardised R object

- [`get_output()`](https://plant-food-research-open.github.io/moiraine/reference/get_output.md)
  [`get_output_pca()`](https://plant-food-research-open.github.io/moiraine/reference/get_output.md)
  [`get_output_splsda()`](https://plant-food-research-open.github.io/moiraine/reference/get_output.md)
  [`get_output_spls()`](https://plant-food-research-open.github.io/moiraine/reference/get_output.md)
  [`get_output_diablo()`](https://plant-food-research-open.github.io/moiraine/reference/get_output.md)
  [`get_output_mofa2()`](https://plant-food-research-open.github.io/moiraine/reference/get_output.md)
  [`get_output_so2pls()`](https://plant-food-research-open.github.io/moiraine/reference/get_output.md)
  : Extract output of integration method in standard format
- [`get_latent_dimensions()`](https://plant-food-research-open.github.io/moiraine/reference/get_latent_dimensions.md)
  : Get latent dimensions levels from dimension reduction output
- [`get_top_features()`](https://plant-food-research-open.github.io/moiraine/reference/get_top_features.md)
  : Extract top features
- [`get_selected_features()`](https://plant-food-research-open.github.io/moiraine/reference/get_selected_features.md)
  : Extract selected features

### Plotting functions

- [`plot_variance_explained()`](https://plant-food-research-open.github.io/moiraine/reference/plot_variance_explained.md)
  : Plot of variance explained
- [`plot_samples_score()`](https://plant-food-research-open.github.io/moiraine/reference/plot_samples_score.md)
  : Plots sample scores as scatterplot matrix
- [`plot_samples_score_pair()`](https://plant-food-research-open.github.io/moiraine/reference/plot_samples_score_pair.md)
  : Plots sample scores as a scatterplot
- [`plot_samples_score_covariate()`](https://plant-food-research-open.github.io/moiraine/reference/plot_samples_score_covariate.md)
  : Plots sample scores against covariate
- [`plot_features_weight_distr()`](https://plant-food-research-open.github.io/moiraine/reference/plot_features_weight_distr.md)
  : Plots features weight distribution
- [`plot_features_weight_pair()`](https://plant-food-research-open.github.io/moiraine/reference/plot_features_weight_pair.md)
  : Plots features weight as a scatterplot
- [`plot_features_weight_covariate()`](https://plant-food-research-open.github.io/moiraine/reference/plot_features_weight_covariate.md)
  : Plots features weight against covariate
- [`plot_top_features()`](https://plant-food-research-open.github.io/moiraine/reference/plot_top_features.md)
  : Plots top features importance

## Methods evaluation

Evaluating the results of an integration method against prior knowledge

- [`evaluate_feature_selection_table()`](https://plant-food-research-open.github.io/moiraine/reference/evaluate_feature_selection_table.md)
  : Evaluate feature selection against features label
- [`make_feature_sets_from_df()`](https://plant-food-research-open.github.io/moiraine/reference/make_feature_sets_from_df.md)
  : Makes list of feature sets from data-frame
- [`make_feature_sets_from_fm()`](https://plant-food-research-open.github.io/moiraine/reference/make_feature_sets_from_fm.md)
  : Makes list of feature sets from features metadata
- [`check_feature_sets()`](https://plant-food-research-open.github.io/moiraine/reference/check_feature_sets.md)
  : Checks features assignment to sets
- [`reduce_feature_sets_data()`](https://plant-food-research-open.github.io/moiraine/reference/reduce_feature_sets_data.md)
  : Reduce feature sets to match multi-omics dataset
- [`evaluate_method_enrichment()`](https://plant-food-research-open.github.io/moiraine/reference/evaluate_method_enrichment.md)
  : Enrichment analysis for integration results
- [`plot_features_weight_set()`](https://plant-food-research-open.github.io/moiraine/reference/plot_features_weight_set.md)
  : Plots features weight in/not in a set
- [`compute_samples_silhouette()`](https://plant-food-research-open.github.io/moiraine/reference/compute_samples_silhouette.md)
  : Computes samples silhouette score from method output

## Methods comparison

Comparison of the results from several integration methods

- [`get_samples_score_correlation()`](https://plant-food-research-open.github.io/moiraine/reference/get_samples_score_correlation.md)
  : Get samples score correlation
- [`get_features_weight_correlation()`](https://plant-food-research-open.github.io/moiraine/reference/get_features_weight_correlation.md)
  : Get features weight correlation
- [`comparison_heatmap_corr()`](https://plant-food-research-open.github.io/moiraine/reference/comparison_heatmap_corr.md)
  : Heatmap of correlation between latent dimensions
- [`comparison_plot_correlation()`](https://plant-food-research-open.github.io/moiraine/reference/comparison_plot_correlation.md)
  : Correlation plot between latent components
- [`compute_consensus_importance()`](https://plant-food-research-open.github.io/moiraine/reference/compute_consensus_importance.md)
  : Computes consensus feature importance
- [`consensus_importance_metric()`](https://plant-food-research-open.github.io/moiraine/reference/consensus_importance_metric.md)
  : Calculate features importance score
- [`show_consensus_metrics()`](https://plant-food-research-open.github.io/moiraine/reference/show_consensus_metrics.md)
  : Illustrates importance consensus metrics

## Running times assessment

Assessment and comparison of the pipeline running time

- [`get_targets_running_time()`](https://plant-food-research-open.github.io/moiraine/reference/get_targets_running_time.md)
  [`plot_running_time()`](https://plant-food-research-open.github.io/moiraine/reference/get_targets_running_time.md)
  : Extract/plot running time of functions used for different
  integration methods
- [`get_method_functions()`](https://plant-food-research-open.github.io/moiraine/reference/get_method_functions.md)
  : Get functions used for each integration method
- [`aggr_patterns_fct()`](https://plant-food-research-open.github.io/moiraine/reference/aggr_patterns_fct.md)
  : Aggregate regexp patterns for a search

## Other (utils)

Other miscellaneous functions

### Helper MultiDataSet object

- [`check_is_multidataset()`](https://plant-food-research-open.github.io/moiraine/reference/check_is_multidataset.md)
  : Checks whether object is MultiDataSet
- [`check_input_multidataset()`](https://plant-food-research-open.github.io/moiraine/reference/check_input_multidataset.md)
  : Check a MultiDataSet input

### Helper plotting functions

- [`plot_correlation_matrix()`](https://plant-food-research-open.github.io/moiraine/reference/plot_correlation_matrix.md)
  : Plot correlation matrix
- [`plot_correlation_matrix_full()`](https://plant-food-research-open.github.io/moiraine/reference/plot_correlation_matrix_full.md)
  : Plots a full correlation matrix (corrplot-style)
- [`ggpairs_custom()`](https://plant-food-research-open.github.io/moiraine/reference/ggpairs_custom.md)
  : ggpairs plot with custom colours
- [`plot_x_wrapper()`](https://plant-food-research-open.github.io/moiraine/reference/plot_x_wrapper.md)
  : Wrapper to create plot
- [`plot_x_continuous()`](https://plant-food-research-open.github.io/moiraine/reference/plot_x_continuous.md)
  : Scatter plot function
- [`plot_x_discrete()`](https://plant-food-research-open.github.io/moiraine/reference/plot_x_discrete.md)
  : Violin plot function function

### Misc.

- [`hclust_matrix_rows()`](https://plant-food-research-open.github.io/moiraine/reference/hclust_matrix_rows.md)
  : Hierarchical clustering of matrix rows
- [`options_list_as_tibble()`](https://plant-food-research-open.github.io/moiraine/reference/options_list_as_tibble.md)
  : Returns options list as a tibble
- [`is_equal_or_null()`](https://plant-food-research-open.github.io/moiraine/reference/is_equal_or_null.md)
  : Check null or equality
