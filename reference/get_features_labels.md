# Get feature labels

Extracts the feature labels from a MultiDataSet object given the name of
the column from the feature metadata of each dataset containing the
feature labels for this dataset.

## Usage

``` r
get_features_labels(mo_data, label_cols = "feature_id", truncate = NULL)
```

## Arguments

- mo_data:

  A
  [MultiDataSet::MultiDataSet](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

- label_cols:

  Character or named list of character, giving for each dataset the name
  of the column in the corresponding features metadata to use as label.
  If one value, will be used for all datasets. If list, the names must
  correspond to the names of the datasets in `mo_data`. Default value is
  `feature_id`, i.e. by default the ID of the features will be used as
  label.

- truncate:

  Integer, width to which the labels should be truncated (to avoid
  issues with very long labels in plots). If `NULL` (default value), no
  truncation will be performed.

## Value

A tibble with columns `dataset`, `feature_id` and `label`.

## Examples

``` r
if (FALSE) { # \dontrun{
## This works if each dataset in mo_data has in their features metadata table
## a column called `name` that contains the feature labels.
get_features_labels(mo_data, label_cols = "name")

## If instead we want to use a different column for each dataset:
get_features_labels(
  mo_data,
  label_cols = list(
    "snps" = "feature_id",
    "rnaseq" = "gene_name",
    "metabolome" = "comp_formula"
  )
)

## If we want to use the feature IDs as labels for the genomics dataset,
## we can simply remove it from the list (this is equivalent to the example
## above):
get_features_labels(
  mo_data,
  label_cols = list(
    "rnaseq" = "gene_name",
    "metabolome" = "comp_formula"
  )
)
} # }
```
