# moiraine 1.0.1

- New functions `get_targets_running_time()` and `plot_running_time()` to view the running time of each target associated with an integration method. 

- In `comparison_heatmap_corr()`, the user can now choose the legend position through the `legend_position` parameter.

- Package here has been removed from dependencies (not needed).

- `where()` function now imported from tidyselect instead of dplyr (as it required a newer version of dplyr).

- Fixed typo in samples metadata file, samples with no value for "rnaseq_batch" variable now have `NA` rather than `"BNA"` values. 

- `perf_splsda()`, `run_splsda()`, `diablo_tune()`, `spls_tune()` and `so2pls_crossval_sparsity()` now have a `seed` argument (hopefully self-explanatory :)). Accordingly, `feature_preselection_splsda_factory` now has arguments `seed_perf` and `seed_run` to pass on seeds to `perf_splsda()` and `run_splsda()`. 

- `create_multiomis_set()` now returns an error if some feature IDs are used across different omics sets. This is to prevent errors further down the line when visualising or subsetting the multi-omics data.

- Fixed a bug which triggered an error when applying `get_output()` to DIABLO results with only one latent component.
