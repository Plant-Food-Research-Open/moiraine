# Makes list of feature sets from data-frame

Creates a list of feature sets from an annotation data-frame.

## Usage

``` r
make_feature_sets_from_df(annotation_df, col_id, col_set)
```

## Arguments

- annotation_df:

  A data-frame of feature annotation in long format, with at least a
  column of feature ID and a column giving the set to which the feature
  belongs. If a feature belongs to more than one set, there should be a
  row for each of these sets.

- col_id:

  Character, name of the column in `annotation_df` data-frame that
  contains the features ID.

- col_set:

  Character, name of the column `annotation_df` data-frame that contains
  the sets ID.

## Value

a named list, where each element corresponds to a set, and contains a
vector of features ID of all features belonging to that set.
