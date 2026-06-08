# Adding data-frame to samples metadata

Adds information from a data-frame to the samples metadata of a
MultiDataSet object.

## Usage

``` r
add_samples_metadata(mo_data, df, datasets = NULL)
```

## Arguments

- mo_data:

  A
  [MultiDataSet::MultiDataSet](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

- df:

  A tibble or data-frame of samples information, with at least column
  `id` (giving the sample IDs), and one row per sample ID.

- datasets:

  Character vector, name of the datasets to which the samples
  information should be added. If `NULL` (default value), the
  information will be added to the samples metadata of all datasets.

## Value

A MultiDataSet object, with info from `df` added to its corresponding
samples metadata.
