# Run sPLS algorithm

Runs the sPLS algorithm
([`mixOmics::spls()`](https://rdrr.io/pkg/mixOmics/man/spls.html)) from
the mixOmics package.

## Usage

``` r
spls_run(spls_input, ...)
```

## Arguments

- spls_input:

  A mixOmics input object created with
  [`get_input_spls()`](https://plant-food-research-open.github.io/moiraine/reference/get_input_spls.md).

- ...:

  Arguments passed to
  [`mixOmics::spls()`](https://rdrr.io/pkg/mixOmics/man/spls.html).

## Value

An object of class `mixo.spls` (if `keepX` and/or `keepY` arguments were
provided) or `mix.pls` (if they were not), see
[`mixOmics::spls()`](https://rdrr.io/pkg/mixOmics/man/spls.html) and
[`mixOmics::pls()`](https://rdrr.io/pkg/mixOmics/man/pls.html).
