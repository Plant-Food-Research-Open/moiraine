# Displays results of sPLS tuning

Displays the results of cross-validation to tune the number of
components to retain in each dataset for a sPLS run. Similar to
[`mixOmics::plot.tune.spls()`](https://rdrr.io/pkg/mixOmics/man/plot.tune.html).

## Usage

``` r
spls_plot_tune(spls_tune_res)
```

## Arguments

- spls_tune_res:

  The result of
  [`spls_tune()`](https://plant-food-research-open.github.io/moiraine/reference/spls_tune.md)
  or
  [`mixOmics::tune.spls()`](https://rdrr.io/pkg/mixOmics/man/tune.spls.html).

## Value

A ggplot.

## Details

The plot displays the correlation or RSS between the latent components
obtained with the corresponding values of `keepX` (x-axis) and `keepY`
(y-axis) and the latent components of the full model (i.e. that retains
all features). The colour of the points shows the mean correlation/RSS
across the cross-validation folds, and the size of the points' shadow
(in gray) represents the coefficient of variation (COV) of the
correlation/RSS, i.e. standard error divided by mean. The point
corresponding with the optimal value for `keepX` and `keepY` is
indicating with a red border.
