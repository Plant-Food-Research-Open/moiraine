# Calculate features importance score

Given a vector of features importance (between 1 and 0), returns a
consensus importance score.

## Usage

``` r
consensus_importance_metric(x, metric, offset = 0)
```

## Arguments

- x:

  Numeric vector of importance values (between 0 and 1).

- metric:

  Character, one of the metrics to use to compute the consensus score.
  Can be one of `'min'`, `'max'`, `'average'`, `'product'`, `'l2'` (for
  L2-norm), `'geometric'` (for geometric mean) or `'harmonic'` (for
  harmonic mean).

- offset:

  Numeric (strictly positive), used to replace zero values to compute
  the geometric and harmonic mean. If `0` (default value), the zero
  values will be ignored when calculating the geometric and harmonic
  mean. Accepting vector of values to facilitate use whithin
  [`dplyr::summarise()`](https://dplyr.tidyverse.org/reference/summarise.html).

## Value

A numeric value, the importance consensus score.
