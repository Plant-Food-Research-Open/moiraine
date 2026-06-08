# Formatted table with optimal keepX values

Produces a nicely formatted table with the optimal number of features to
select from each dataset for a DIABLO run according to the results of a
cross-validation analysis. Used mostly for producing a nice report.

## Usage

``` r
diablo_table_optim_keepX(tune_res)
```

## Arguments

- tune_res:

  The cross-validation results, computed with
  [`diablo_tune`](https://plant-food-research-open.github.io/moiraine/reference/diablo_tune.md).

## Value

A tibble with a Dataset column, one column for each latent component and
a Total column giving the number of features to retain from the
corresponding dataset across all latent components.
