# Join samples metadata to table

Adds samples metadata information to a table containing sample IDs.

## Usage

``` r
join_samples_metadata(df, mo_data, datasets = NULL)
```

## Arguments

- df:

  Data-frame or tibble with a column `id` containing sample IDs.

- mo_data:

  A
  [MultiDataSet::MultiDataSet](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

- datasets:

  Character vector, name(s) of datasets from which the samples metadata
  should be extracted. If `NULL` (default value), information from all
  datasets will be used.

## Value

The `df` table with additional columns containing information about the
samples from the samples metadata table.
