# Changelog

## moiraine 1.0.1

- Fixed bug in
  [`so2pls_get_optim_ncomp()`](https://plant-food-research-open.github.io/moiraine/reference/so2pls_get_optim_ncomp.md)
  which was breaking in presence of NAs.

- New functions
  [`get_targets_running_time()`](https://plant-food-research-open.github.io/moiraine/reference/get_targets_running_time.md)
  and
  [`plot_running_time()`](https://plant-food-research-open.github.io/moiraine/reference/get_targets_running_time.md)
  to view the running time of each target associated with an integration
  method.

- In
  [`comparison_heatmap_corr()`](https://plant-food-research-open.github.io/moiraine/reference/comparison_heatmap_corr.md),
  the user can now choose the legend position through the
  `legend_position` parameter.

- Package here has been removed from dependencies (not needed).

- `where()` function now imported from tidyselect instead of dplyr (as
  it required a newer version of dplyr).

- Fixed typo in samples metadata file, samples with no value for
  “rnaseq_batch” variable now have `NA` rather than `"BNA"` values.

- [`perf_splsda()`](https://plant-food-research-open.github.io/moiraine/reference/perf_splsda.md),
  [`run_splsda()`](https://plant-food-research-open.github.io/moiraine/reference/run_splsda.md),
  [`diablo_tune()`](https://plant-food-research-open.github.io/moiraine/reference/diablo_tune.md),
  [`spls_tune()`](https://plant-food-research-open.github.io/moiraine/reference/spls_tune.md)
  and
  [`so2pls_crossval_sparsity()`](https://plant-food-research-open.github.io/moiraine/reference/so2pls_crossval_sparsity.md)
  now have a `seed` argument (hopefully self-explanatory :)).
  Accordingly, `feature_preselection_splsda_factory` now has arguments
  `seed_perf` and `seed_run` to pass on seeds to
  [`perf_splsda()`](https://plant-food-research-open.github.io/moiraine/reference/perf_splsda.md)
  and
  [`run_splsda()`](https://plant-food-research-open.github.io/moiraine/reference/run_splsda.md).

- `create_multiomis_set()` now returns an error if some feature IDs are
  used across different omics sets. This is to prevent errors further
  down the line when visualising or subsetting the multi-omics data.

- Fixed a bug which triggered an error when applying
  [`get_output()`](https://plant-food-research-open.github.io/moiraine/reference/get_output.md)
  to DIABLO results with only one latent component.
