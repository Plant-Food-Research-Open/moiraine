# Get feature weights from MOFA object

Extracts the feature weights from a trained MOFA or MEFISTO model (from
the MOFA2 package).

## Usage

``` r
mofa_get_weights(
  object,
  views = "all",
  factors = "all",
  abs = FALSE,
  scale = "none",
  as.data.frame = TRUE
)
```

## Arguments

- object:

  A trained MOFA object.

- views:

  Character or integer vector, name or index of the views (i.e.
  datasets) for which the feature weights should be extracted. Default
  value is "all", i.e. all datasets are considered.

- factors:

  Character or integer vector, name or index of the factors for which
  the feature weights should be extracted. Default value is "all", i.e.
  all factors are considered.

- abs:

  Logical, should the absolute value of the weights be returned? Default
  value is `FALSE`.

- scale:

  Character, the type of scaling that should be performed on the feature
  weights. Possible values are `'none'`, `'by_view'`, `'by_factor'` or
  `'overall'` (see Details). Default value is `'none'`.

- as.data.frame:

  Logical, whether the function should return a long data-frame instead
  of a list of matrices. Default value is `TRUE`.

## Value

By default, returns a tibble with columns `view`, `feature`, `factor`,
`value`. Alternatively, if `as.data.frame = FALSE`, returns a list of
matrices, one per view, with features as rows and factors as columns.

## Details

Scaling options:

- `scale = 'none'`: no scaling performed;

- `scale = 'by_view'`: weights are divided by the maximum absolute
  weight in the corresponding view/dataset;

- `scale = 'by_factor'`: weights are divided by the maximum absolute
  weight in the corresponding factor;

- `scale = 'overall'`: weights are divided by the maximum absolute
  weight across all views/datasets and factors considered.
