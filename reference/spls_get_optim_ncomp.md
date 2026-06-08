# Select the optimal ncomp from sPLS cross-validation results

Select the optimal number of components to compute in a sPLS run from
the cross-validation results obtained with
[`mixOmics::perf()`](https://rdrr.io/pkg/mixOmics/man/perf.html) on a
PLS or sPLS result, using the the mean `Q2.total` values. Note that this
function is experimental, the corresponding diagnostic plots should be
considered when selecting the optimal number of components to use.

## Usage

``` r
spls_get_optim_ncomp(spls_perf, thr = 0.0975, min_ncomp = 1)
```

## Arguments

- spls_perf:

  List, the result of
  [`mixOmics::perf()`](https://rdrr.io/pkg/mixOmics/man/perf.html).

- thr:

  Numeric, the threshold to be used for the `Q2` values. Default value
  is 0.0975.

- min_ncomp:

  Integer, the minimum `ncomp` value to be returned. Default value is 1,
  i.e. this argument does not play a role in selecting the `comp` value.
  Can be useful if we want at least 2 latent components for final plots.

## Value

An integer, the optimal number of components to use for the sPLS run.

## Details

The selection is made as follows:

- If all `Q2` values are below the threshold specified with `thr`, the
  number of components yielding the highest `Q2` value is selected.

- If all `Q2` values are above the threshold, the number of components
  yielding the lowest `Q2` value is selected.

- If the `Q2` values are increasing, the number of components `n` is
  selected such that `n+1` is the smallest number of components with a
  `Q2` value above the threshold.

- If the `Q2` values are decreasing, the number of components `n` is
  selected such that `n+1` is the smallest number of components with a
  `Q2` value below the threshold.
