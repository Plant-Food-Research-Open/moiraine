# Samples score plots for single-omics PCA

Produces a pairwise samples score plot for the PCA run on each omics
dataset, using
[`GGally::ggpairs()`](https://ggobi.github.io/ggally/reference/ggpairs.html).

## Usage

``` r
plot_samples_coordinates_pca(pca_result, pcs = NULL, datasets = NULL, ...)
```

## Arguments

- pca_result:

  List of PCA results on each of the datasets, each computed with the
  [`run_pca()`](https://plant-food-research-open.github.io/moiraine/reference/run_pca.md)
  function.

- pcs:

  Integer vector or named list of integer vectors, the principal
  components to display for each dataset. If integer vector (e.g.
  `1:5`), will be used for all datasets. Alternatively, a different set
  of PCs can be specified through a named list (e.g.
  `list('snps' = 1:4, 'rnaseq' = 1:5)`). The length of the list must
  match the number of datasets to be displayed, and the names must match
  the dataset names. Default value is `NULL`, i.e. all principal
  components will be plotted for each dataset.

- datasets:

  Optional, character vector of datasets for which the plots should be
  created.

- ...:

  Other arguments passed to
  [`plot_samples_score()`](https://plant-food-research-open.github.io/moiraine/reference/plot_samples_score.md).

## Value

A list of ggmatrix plots (or a single ggmatrix plot if `pcs` has length
1).

## Examples

``` r
if (FALSE) { # \dontrun{
## Default: plotting all PCs for all datasets
plot_samples_coordinates_pca(pca_result)

## Plotting only the first 3 PCs for each dataset
plot_samples_coordinates_pca(
  pca_result,
  pcs = 1:3
)

## Plotting the first 3 PCs for the genomics dataset, 4 PCs for the
## transcriptomics dataset, 5 PCs for the metabolomics dataset
plot_samples_coordinates_pca(
  pca_result,
  pcs = list(
    "snps" = 1:3,
    "rnaseq" = 1:4,
    "metabolome" = 1:5
  )
)

## Plotting the first 3 PCs for the genomics and transcriptomics datasets
plot_samples_coordinates_pca(
  pca_result,
  pcs = 1:3,
  datasets = c("snps", "rnaseq")
)

# Plotting the first 3 PCs for the genomics dataset and 4 PCs for the
## transcriptomics dataset (no plot for the metabolomics dataset)
plot_samples_coordinates_pca(
  pca_result,
  pcs = list(
    "snps" = 1:3,
    "rnaseq" = 1:4
  ),
  datasets = c("snps", "rnaseq")
)
} # }
```
