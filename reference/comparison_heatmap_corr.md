# Heatmap of correlation between latent dimensions

Constructs a heatmap displaying the correlation between the latent
dimensions constructed with several integration methods. The lower
triangle of the heatmap displays the correlation between features
weight, while the upper triangle shows the correlation between samples
score. Each triangle of the matrix is reordered separately to show the
most highly correlated dimensions next to each other.

## Usage

``` r
comparison_heatmap_corr(
  output_list,
  latent_dimensions = NULL,
  include_missing_features = FALSE,
  legend_ncol = length(output_list),
  legend_position = "bottom",
  legend_title = "correlation"
)
```

## Arguments

- output_list:

  List of integration methods output each generated via the
  [`get_output()`](https://plant-food-research-open.github.io/moiraine/reference/get_output.md)
  function. If named, the names will be used to annotate the plot. See
  details.

- latent_dimensions:

  Named list, where each element is a character vector giving the latent
  dimensions to retain in the corresponding element of `output_list`.
  Names must match those of `output_list`. Can be used to filter latent
  dimensions only in certain elements from output_list (see examples).
  If `NULL` (default value), all latent dimensions will be used.

- include_missing_features:

  Logical, see
  [`get_features_weight_correlation()`](https://plant-food-research-open.github.io/moiraine/reference/get_features_weight_correlation.md)
  for details. Default value is `FALSE`.

- legend_ncol:

  Integer, number of columns in the legend. Default value is set to the
  length of `output_list`.

- legend_position:

  Character, position of the legend. Should be one of `"bottom"`
  (default), `"top"`, `"left"` or `"right`.

- legend_title:

  Character, name to give to the heatmap colour legend. Intended to use
  to shorten the legend title if it goes out of frame.

## Value

a
[`ComplexHeatmap::Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html).

## Details

If `output_list` is unnamed, the different elements in the list will be
differentiated by the name of the method used to produce them (e.g.
DIABLO, sO2PLS, etc). In order to compare different results from a same
integration method (e.g. DIABLO applied to the full vs pre-filtered
data), it is possible to assign names to the elements of `output_list`
(see examples). These names will be used in place of the method name in
the plot to identify where the latent dimensions come from.

## Examples

``` r
if (FALSE) {
## Comparing the output from DIABLO, sO2PLS and MOFA

res <- list(
  get_output_diablo(diablo_res), ## diablo_res: output from diablo_run()
  get_output_so2pls(so2pls_res), ## so2pls_res: output from so2pls_o2m()
  get_output_mofa2(mofa_res) ## mofa_res: output from run_mofa
)

comparison_heatmap_corr(res)

## Selecting only some factors from a MOFA run for the comparison
## (for the other methods, all latent dimensions will be retained)
comparison_heatmap_corr(
  res,
  latent_dimensions = list(
    "MOFA" = paste0("Factor ", 1:3)
  )
)

## Comparing two different results from a same integration method -
## diablo_run_full and diablo_run_prefiltered would both be output
## from the diablo_run() function.

res <- list(
  "DIABLO full" = get_output_diablo(diablo_run_full),
  "DIABLO prefiltered" = get_output_diablo(diablo_run_prefiltered)
)

comparison_heatmap_corr(res)
}
```
