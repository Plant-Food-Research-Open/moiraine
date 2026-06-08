# Generate mixomics input data for supervised methods

Creates an object that can be used as input for the MixOmics package. It
contains the omics datasets restricted to common samples (with no
missing group information) and an outcome group for each sample.

## Usage

``` r
get_input_mixomics_supervised(mo_data, group, datasets = names(mo_data))
```

## Arguments

- mo_data:

  A
  [`MultiDataSet-class`](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

- group:

  Character, the column name in the samples metadata data-frames to use
  as samples group (use
  [`get_samples_metadata`](https://plant-food-research-open.github.io/moiraine/reference/get_samples_metadata.md)
  to view the samples information data-frame for each dataset). This
  column should be either of type factor, character or integer.

- datasets:

  Character vector, the names of the datasets from `mo_data` to include
  in the analysis.

## Value

A list, in which each element corresponds to one omics dataset, with
samples as rows and features as columns. The `Y` element is a factor
vector with the outcome group of each sample.

## Details

Samples with missing values in the `group` column of the sample metadata
will be removed from the dataset.
