# Get features weight correlation

Constructs the correlation matrix between the features weight of the
latent dimensions obtained with different integration methods.

## Usage

``` r
get_features_weight_correlation(output_list, include_missing_features = FALSE)
```

## Arguments

- output_list:

  List of integration methods output, each generated via the
  [`get_output()`](https://plant-food-research-open.github.io/moiraine/reference/get_output.md)
  function. If named, the names will be added at the beginning of each
  latent dimension' label. If unnamed, the name of the integration
  method will be used instead.

- include_missing_features:

  Logical, whether features missing in some of the output should be
  included in the calculation (see Details). Default value is `FALSE`.

## Value

A correlation matrix.

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
a weight of 0.
