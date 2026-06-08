# Join feature metadata to table

Adds features metadata information to a table containing feature IDs.

## Usage

``` r
join_features_metadata(df, mo_data)
```

## Arguments

- df:

  Data-frame of tibble with a column `feature_id` containing feature
  IDs.

- mo_data:

  A
  [MultiDataSet::MultiDataSet](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

## Value

The `df` table with additional columns containing information about the
features from the features metadata table.
