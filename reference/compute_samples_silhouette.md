# Computes samples silhouette score from method output

Calculates the samples silhouette width from the results of a dimension
reduction method, according to some samples grouping from the samples
metadata of a `MultiDataSet` object.

## Usage

``` r
compute_samples_silhouette(
  method_output,
  mo_data,
  group_col,
  latent_dimensions = NULL,
  distance_metric = c("euclidean", "maximum", "manhattan", "canberra", "binary",
    "minkowski")
)
```

## Arguments

- method_output:

  method_output Integration method output generated via the
  [`get_output()`](https://plant-food-research-open.github.io/moiraine/reference/get_output.md)
  function.

- mo_data:

  A
  [`MultiDataSet-class`](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

- group_col:

  Character, name of the column in one of the samples metadata table
  from `mo_data` containing the samples grouping to be used.

- latent_dimensions:

  Character vector, the latent dimensions to use for computing distance
  between samples. If `NULL` (default value), all latent dimensions will
  be used.

- distance_metric:

  Character, name of the metric to use when computing the distance
  between samples from their coordinates in the latent dimensions. This
  will be passed to the
  [`stats::dist()`](https://rdrr.io/r/stats/dist.html) function. Options
  include: `"euclidean"` (default value), `"maximum"`, `"manhattan"`,
  `"canberra"`, `"binary"` or `"minkowski"`.

## Value

A list with two elements:

- `samples_silhouette`: a tibble giving for each sample (`sample_id`
  column) the group to which it belongs (`group` column), its closest
  (other) group in the space spanned by the latent dimensions
  (`neighbour_group`), and its silhouette width (`silhouette_width`
  column).

- `groups_average_silhouette`: a tibble giving for each samples group
  (`group` column) its average silhouette width (`group_average_width`
  column).

## Details

The samples silhouette width and groups average width are calculated
using the
[`cluster::silhouette()`](https://rdrr.io/pkg/cluster/man/silhouette.html)
function.
