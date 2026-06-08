# Get samples score correlation

Constructs the correlation matrix between the samples score of the
latent dimensions obtained with different integration methods.

## Usage

``` r
get_samples_score_correlation(output_list)
```

## Arguments

- output_list:

  List of integration methods output, each generated via the
  [`get_output()`](https://plant-food-research-open.github.io/moiraine/reference/get_output.md)
  function. If named, the names will be added at the beginning of each
  latent dimension' label. If unnamed, the name of the integration
  method will be used instead.

## Value

A correlation matrix.
