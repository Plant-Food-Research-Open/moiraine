# Check that variable names corresponds to columns in samples metadata

Checks whether a variable name corresponds to a column in the samples
metadata of the corresponding dataset. If one value is provided, will be
used for all datasets.

## Usage

``` r
.check_input_var_smetadata_common(x, mo_data)
```

## Arguments

- x:

  Character, name of the column from the samples metadata.

- mo_data:

  A `MultiDataSet` object containing samples information for the
  datasets. Should be checked with
  [`check_input_multidataset()`](https://plant-food-research-open.github.io/moiraine/reference/check_input_multidataset.md).

## Value

Nothing. Will throw an error if need be.
