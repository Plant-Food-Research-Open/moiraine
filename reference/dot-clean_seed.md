# Clean seed for sPLS-DA preselection factory

Checks seed input arguments for the
[`feature_preselection_splsda_factory()`](https://plant-food-research-open.github.io/moiraine/reference/feature_preselection_splsda_factory.md)
function.

## Usage

``` r
.clean_seed(x, ds)
```

## Arguments

- x:

  Integer vector of seeds to use.

- ds:

  Character vector of dataset names.

## Value

`x` with names if it didn't have them or an error.
