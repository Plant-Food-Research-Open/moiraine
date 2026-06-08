# Extract arguments used in PCA run

Extracts the list of arguments used for each PCA run from a list of PCA
results, and formats them into a tibble.

## Usage

``` r
get_pca_arguments(pca_result)
```

## Arguments

- pca_result:

  The result of a PCA run on each of the datasets, computed with the
  [`pcaMethods::pca()`](https://rdrr.io/pkg/pcaMethods/man/pca.html)
  function.

## Value

A tibble with the following columns: `"Omics dataset"`,
`"PCA method used"`, `"Number of Principal Components computed"`,
`"Scaling applied"` and `"Dataset centered"`.
