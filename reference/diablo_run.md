# Runs DIABLO algorithm

Runs the DIABLO algorithm
([`block.splsda`](https://rdrr.io/pkg/mixOmics/man/block.splsda.html))
from the `mixOmics` package.

## Usage

``` r
diablo_run(mixomics_data, design_matrix, ...)
```

## Arguments

- mixomics_data:

  A `mixOmics` input object created with
  [`get_input_mixomics_supervised`](https://plant-food-research-open.github.io/moiraine/reference/get_input_mixomics_supervised.md).

- design_matrix:

  Either numeric matrix created through
  [`diablo_generate_design_matrix`](https://plant-food-research-open.github.io/moiraine/reference/diablo_generate_design_matrix.md),
  or character (accepted values are `'null'`, `'weighted_full'`,
  `'full'`). See Details.

- ...:

  Arguments to be passed to the
  [`block.splsda`](https://rdrr.io/pkg/mixOmics/man/block.splsda.html)
  function.

## Value

An object of class `block.splsda` (if `keepX` argument was provided) or
`block.splsda` (if it was not), see
[`mixOmics::block.splsda()`](https://rdrr.io/pkg/mixOmics/man/block.splsda.html)
and
[`mixOmics::block.plsda()`](https://rdrr.io/pkg/mixOmics/man/block.plsda.html).

## Details

The `design_matrix` argument can either be a custom design matrix (for
example as constructed via the
[`diablo_generate_design_matrix`](https://plant-food-research-open.github.io/moiraine/reference/diablo_generate_design_matrix.md)
function); or a character indicating the type of design matrix to
generate. Possible values include:

- `'null'`: Off-diagonal elements of the design matrix are set to 0;

- `'weighted_full'`: Off-diagonal elements of the design matrix are set
  to 0.1;

- `'full'`: Off-diagonal elements of the design matrix are set to 1.
