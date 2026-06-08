# Get the optimal ncomp value

Selects the optimal `comp` value (number of components to compute) from
a DIABLO cross-validation run, for a given error measure and distance
metric.

## Usage

``` r
diablo_get_optim_ncomp(
  perf_res,
  measure = "Overall.BER",
  distance = "centroids.dist",
  min_ncomp = 1
)
```

## Arguments

- perf_res:

  The cross-validation results, computed with
  [`perf`](https://rdrr.io/pkg/mixOmics/man/perf.html).

- measure:

  The error measure for which to obtain the optimal value; possible
  values are `'Overall.ER'` and `'Overall.BER'`. Default value is
  `'Overall.BER'`.

- distance:

  The distance metric for which to obtain the optimal value; possible
  values are `'max.dist'`, `'centroids.dist'` and `'mahalanobis.dist'`.
  Default value is `'centroids.dist'`.

- min_ncomp:

  Integer, the minimum `ncomp` value to be returned. Default value is 1,
  i.e. this argument does not play a role in selecting the `comp` value.
  Can be useful if we want at least 2 latent components for final plots.

## Value

An integer, the optimal value of `ncomp` to use for the DIABLO run.
