# Generate DIABLO design matrix

Generates a design matrix for the DIABLO algorithm, based on the
correlation between datasets inferred from pairwise PLS runs.

## Usage

``` r
diablo_generate_design_matrix(
  cormat,
  threshold = 0.8,
  low_val = 0.1,
  high_val = 1,
  y_val = 1
)
```

## Arguments

- cormat:

  The correlation matrix between datasets, obtained from
  `diablo_get_pairwise_pls_corr`.

- threshold:

  Numeric, correlation value above which datasets are considered as
  highly correlated (see Details). Default value is 0.8.

- low_val:

  Numeric, value in the design matrix for datasets that are not highly
  correlated. Default value is 0.1.

- high_val:

  Numeric, value in the design matrix for datasets that are highly
  correlated. Default value is 1.

- y_val:

  Numeric, value in the design matrix between datasets and the outcome
  (Y). Default value is 1.

## Value

A numeric matrix, to be used as design matrix when running
[`block.plsda`](https://rdrr.io/pkg/mixOmics/man/block.plsda.html), with
one row per dataset and one column per dataset.

## Details

Use a threshold to detect pairs of datasets that are highly correlated.
For those pairs of datasets, the corresponding cell in the design matrix
is set to `high_val`. For pairs of datasets with a correlation below the
threshold, the corresponding cell in the design matrix is set to
`low_val`.
