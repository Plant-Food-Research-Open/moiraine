# Generate mixomics input data for unsupervised methods

Creates an object that can be used as input for the MixOmics package. It
contains the omics datasets restricted to common samples.

## Usage

``` r
get_input_mixomics_unsupervised(mo_data, datasets = names(mo_data))
```

## Arguments

- mo_data:

  A
  [`MultiDataSet-class`](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

- datasets:

  Character vector, the names of the datasets from `mo_data` to include
  in the analysis.

## Value

A list, in which each element corresponds to one omics dataset, with
samples as rows and features as columns.
