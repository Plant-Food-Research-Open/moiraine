# Plots omics data as heatmap

For a given set of features, plots their value across the samples as a
heatmap.

## Usage

``` r
plot_data_heatmap(
  mo_data,
  features,
  center = FALSE,
  scale = FALSE,
  samples = NULL,
  only_common_samples = FALSE,
  samples_info = NULL,
  features_info = NULL,
  colours_list = NULL,
  label_cols = NULL,
  truncate = NULL,
  legend_title_size = 10,
  legend_text_size = 10,
  legend_ncol = 1,
  ...
)
```

## Arguments

- mo_data:

  A
  [MultiDataSet::MultiDataSet](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

- features:

  Character vector, the ID of features to show in the plot.

- center:

  Logical, whether the data should be centered (feature-wise). Default
  value is `FALSE`.

- scale:

  Logical, whether the data should be scaled (feature-wise). Default
  value is `FALSE`.

- samples:

  Character vector, the ID of samples to include in the plot. If `NULL`
  (default), all samples will be used.

- only_common_samples:

  Logical, whether only samples that are present in all datasets should
  be plotted. Default value is `FALSE`.

- samples_info:

  Character vector, column names from the samples metadata tables of the
  datasets to be represented in the plot as samples annotation.

- features_info:

  Named list of character vectors, where each element corresponds to a
  dataset, and gives the column names from the features metadata of the
  dataset to be represented in the plot as features annotation. The
  names of the list must correspond to dataset names in the `mo_data`
  object.

- colours_list:

  Named list, where each element gives the colour palette to use in the
  samples or features annotation. Names must match values in
  `samples_info` vector and elements of `features_info` list. For
  continuous palettes, must use
  [`circlize::colorRamp2()`](https://rdrr.io/pkg/circlize/man/colorRamp2.html)
  function (see the [ComplexHeatmap reference
  book](https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html#simple-annotation)).

- label_cols:

  Character or named list of character, giving for each dataset the name
  of the column in the corresponding features metadata to use as label.
  If one value, will be used for all datasets. If list, the names must
  correspond to the names of the datasets in `mo_data`. If a dataset is
  missing from the list or no value is provided, feature IDs will be
  used as labels. Alternatively, use `feature_id` to get the feature IDs
  as labels.

- truncate:

  Integer, width to which the labels should be truncated (to avoid
  issues with very long labels in plots). If `NULL` (default value), no
  truncation will be performed.

- legend_title_size:

  Integer, size in points of legend title.

- legend_text_size:

  Integer, size in points of legend elements text.

- legend_ncol:

  Integer, number of columns in the legend. Default value is `1`.

- ...:

  Additional arguments passed to the
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
  function.

## Value

A
[ComplexHeatmap::Heatmap](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap-class.html)
object.

## Examples

``` r
if (FALSE) { # \dontrun{
## Selecting at random 3 features from each dataset
random_features <- get_features(mo_set) |>
  map(sample, size = 3, replace = FALSE) |>
  unlist() |>
  unname()

plot_data_heatmap(
  mo_set,
  random_features,
  center = TRUE,
  scale = TRUE,
  show_column_names = FALSE,
  only_common_samples = TRUE,
  samples_info = c("status", "day_on_feed"),
  features_info = c("chromosome"),
  colours_list = list(
    "status" = c("Control" = "gold", "BRD" = "navyblue"),
    "day_on_feed" = colorRamp2(c(5, 65), c("white", "pink3"))
  ),
  label_cols = list(
     "rnaseq" = "Name",
    "metabolome" = "name"
  )
 )
} # }
```
