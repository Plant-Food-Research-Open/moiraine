# Tunes keepX arg for DIABLO

Performs cross-validation to estimate the optimal number of features to
retain from each dataset for a DIABLO run.

## Usage

``` r
diablo_tune(
  mixomics_data,
  design_matrix,
  keepX_list = NULL,
  cpus = NULL,
  seed = NULL,
  ...
)
```

## Arguments

- mixomics_data:

  A `mixOmics` input object created with
  [`get_input_mixomics_supervised()`](https://plant-food-research-open.github.io/moiraine/reference/get_input_mixomics_supervised.md).

- design_matrix:

  Either numeric matrix created through `diablo_generate_design_matrix`,
  or character (accepted values are `'null'`, `'weighted_full'`,
  `'full'`). See Details.

- keepX_list:

  Named list, gives for each omics dataset in the mixOmics input (i.e.
  excluding the response Y) a vector of values to test (i.e. number of
  features to return from this dataset). If `NULL` (default), a standard
  grid will be applied for each dataset and latent component, testing
  values: `seq(5, 30, 5)`.

- cpus:

  Integer, the number of CPUs to use when running the code in parallel.
  For advanced users, see the `BPPARAM` argument of
  [`mixOmics::tune.block.splsda()`](https://rdrr.io/pkg/mixOmics/man/tune.block.splsda.html).

- seed:

  Integer, seed to use. Default is `NULL`, i.e. no seed is set inside
  the function.

- ...:

  Arguments to be passed to the
  [`mixOmics::tune.block.splsda()`](https://rdrr.io/pkg/mixOmics/man/tune.block.splsda.html)
  function.

## Value

A list, see
[`mixOmics::tune.block.splsda()`](https://rdrr.io/pkg/mixOmics/man/tune.block.splsda.html).

## Details

The
``` design_matrix`` argument can either be a custom design matrix (for example as constructed via the  ```diablo_generate_design_matrix\`
function); or a character indicating the type of design matrix to
generate. Possible values include:

- `'null'`: Off-diagonal elements of the design matrix are set to 0;

- `'weighted_full'`: Off-diagonal elements of the design matrix are set
  to 0.1;

- `'full'`: Off-diagonal elements of the design matrix are set to 1.
