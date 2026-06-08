# Adds features label to data-frame

Adds the features label to a data-frame for plotting. Can be extracted
from the features metadata of a `MultiDataSet` object; otherwise use the
feature IDs as label. If some labels are missing, feature IDs will be
used instead.

## Usage

``` r
.add_features_labels_toplot(toplot, label_cols, mo_data, truncate = NULL)
```

## Arguments

- toplot:

  The data-frame to which the labels should be added.

- label_cols:

  Character or named list of character, giving for each dataset the name
  of the column in the corresponding features metadata to use as label.
  If one value, will be used for all datasets. If list, the names must
  correspond to the names of the datasets in `mo_data`. If a dataset is
  missing from the list or no value is provided, feature IDs will be
  used as labels. Alternatively, use `feature_id` to get the feature IDs
  as labels.

- mo_data:

  A `MultiDataSet` object. Only used if `label_cols` is not `NULL`.

- truncate:

  Integer, width to which the labels should be truncated (to avoid
  issues with very long labels in plots). If `NULL` (default value), no
  truncation will be performed.

## Value

the `toplot` data-frame with an additional column `label`.
