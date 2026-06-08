# Plots DIABLO output

Displays the samples coordinates for a given latent component across the
datasets. This is a copy of the
[`plotDiablo`](https://rdrr.io/pkg/mixOmics/man/plotDiablo.html)
function, the only difference is that the margins have been increased to
accomodate for a title.

## Usage

``` r
diablo_plot(
  diablo_res,
  ncomp = 1,
  legend = TRUE,
  legend.ncol,
  col.per.group = NULL
)
```

## Arguments

- diablo_res:

  The output from
  [`block.splsda`](https://rdrr.io/pkg/mixOmics/man/block.splsda.html)
  or
  [`diablo_run`](https://plant-food-research-open.github.io/moiraine/reference/diablo_run.md).

- ncomp:

  Integer, the latent component to plot.

- legend:

  Logical, should the legend be added to the plot? Default value is
  `TRUE`.

- legend.ncol:

  Integer, the number of columns for the legend. If none specified, will
  be calculated as `min(5, nlevels(diablo_res$Y))`.

- col.per.group:

  Named character vector, provides the colours to use for each
  phenotypic group. Names must match the levels of `diablo_res$Y`.

## Value

None.
