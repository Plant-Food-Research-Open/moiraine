# Performs cross-validation for mixomics sPLS to select optimal keepX and keepY

Peforms cross-validation to assess the optimal number of features to
retain in each dataset for a sPLS run (implemented in the `mixOmics`
package).

## Usage

``` r
spls_tune(
  spls_input,
  keepX = NULL,
  keepY = NULL,
  cpus = NULL,
  seed = NULL,
  ...
)
```

## Arguments

- spls_input:

  A `mixOmics` input object created with
  [`get_input_mixomics_unsupervised()`](https://plant-food-research-open.github.io/moiraine/reference/get_input_mixomics_unsupervised.md).

- keepX:

  Numeric vector, values for the number of features to retain from
  dataset X to test. Default value is `NULL` (default sequence of values
  will be used, see details).

- keepY:

  Numeric vector, values for the number of features to retain from
  dataset Y to test. Default value is `NULL` (default sequence of values
  will be used, see details).

- cpus:

  Integer, the number of CPUs to use when running the code in parallel.
  For advanced users, see the `BPPARAM` argument of
  [`mixOmics::tune.spls()`](https://rdrr.io/pkg/mixOmics/man/tune.spls.html).

- seed:

  Integer, seed to use. Default is `NULL`, i.e. no seed is set inside
  the function.

- ...:

  Further arguments passed to
  [`mixOmics::tune.spls()`](https://rdrr.io/pkg/mixOmics/man/tune.spls.html).

## Value

A list (see
[`mixOmics::tune.spls()`](https://rdrr.io/pkg/mixOmics/man/tune.spls.html)).

## Details

If no value is provided for `keepX` or `keepY`, the sequence
`seq(5, 30, 5)` will be used, truncated to only retain values that are
inferior or equal to the number of features in the X or Y dataset.
