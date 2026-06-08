# Generate omicsPLS input data

Creates an object that can be used as input for the omicsPLS package. It
contains the omics datasets restricted to common samples. Each dataset
is feature-centred.

## Usage

``` r
get_input_omicspls(mo_data, datasets = names(mo_data), scale_data = FALSE)
```

## Arguments

- mo_data:

  A
  [`MultiDataSet-class`](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

- datasets:

  Character vector of length 2, the names of the datasets from `mo_data`
  to include in the analysis.

- scale_data:

  Boolean, should the datasets be scaled? Default value is `FALSE`.

## Value

A list, in which each element corresponds to one omics dataset, with
samples as rows and features as columns.
