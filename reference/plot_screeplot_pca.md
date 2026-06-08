# Screeplots for single-omics PCA

Produces a scree plot (percentage of variance explained by each
principal component) for the PCA run on each omics dataset, using
ggplot2.

## Usage

``` r
plot_screeplot_pca(pca_result, cumulative = FALSE, datasets = NULL)
```

## Arguments

- pca_result:

  List of PCA runs result on each of the datasets, each computed with
  the
  [`run_pca()`](https://plant-food-research-open.github.io/moiraine/reference/run_pca.md)
  function.

- cumulative:

  Logical, should the cumulative variance be plotted? Default is
  `FALSE`.

- datasets:

  Optional, character vector with the names of the datasets to plot. If
  `NULL` (default value), all datasets will be plotted.

## Value

A ggplot2 plot.
