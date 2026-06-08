# Computes consensus feature importance

Computes the consensus feature importance from features weight obtained
with different integration methods (considering features importance for
one latent component per integration method), or for different latent
dimensions constructed by an integration method.

## Usage

``` r
compute_consensus_importance(
  output_list,
  latent_dimensions,
  metric = "geometric",
  include_missing_features = FALSE
)
```

## Arguments

- output_list:

  List of integration methods output, each generated via the
  [`get_output()`](https://plant-food-research-open.github.io/moiraine/reference/get_output.md)
  function, or a single integration method output (from
  [`get_output()`](https://plant-food-research-open.github.io/moiraine/reference/get_output.md)).

- latent_dimensions:

  Named list (if `output_list` is a list), where each element is a
  character giving the latent dimension to retain in the corresponding
  element of `output_list` (1 value). If `output_list` is a single
  output object, needs to be instead a character vector giving the
  latent dimensions to retain.

- metric:

  Character, one of the metrics to use to compute the consensus score.
  Can be one of `'min'`, `'max'`, `'average'`, `'product'`, `'l2'` (for
  L2-norm), `'geometric'` (for geometric mean) or `'harmonic'` (for
  harmonic mean). Default value is `'geometric'`. Names must match those
  of `output_list`.

- include_missing_features:

  Logical, whether features missing in some of the output should be
  included in the calculation (see Details). Default value is `FALSE`.

## Value

A tibble giving the consensus importance of each feature.

## Details

If `include_missing_features` is `FALSE` (default behaviour), and some
features are present in the output of one integration method but not the
other (e.g. because a different pre-filtering was applied to the input
data of the two methods), these features will be ignored. This does not
mean that features that were selected by one method but not the other
are discarded; in that case the feature will be assigned a weight of 0
for the method that did not select it. This is the recommended
behaviour, should only be changed in specific scenarios (e.g. to check
whether using all features in a dataset vs doing a variance-based
preselection affect which features are deemed most important). If
`include_missing_features` is `TRUE`, missing features will be assigned
a weight of 0. Note that the geometric and harmonic means only work for
strictly positive values. Therefore, all importance scores of 0 are
replaced with an offset when computing these metrics. The offset is
calculated per dataset, and corresponds to the minimum non-null
importance score observed across all features in the dataset (and across
all latent dimensions), divided by 2. The calculation of the offset is
done before removing missing features (if
`include_missing_features = FALSE`) so that results are consistent
between the two options for `include_missing_features`.
