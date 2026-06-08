# Enrichment analysis for integration results

Performs an enrichment analysis for each latent dimension in an
integration result, based on user-defined feature sets. The enrichment
analysis is done with the
[`gage::gage()`](https://rdrr.io/pkg/gage/man/gage.html) function from
the [`gage`](https://bioconductor.org/packages/gage/) package, using the
features' signed importance score as features metric.

## Usage

``` r
evaluate_method_enrichment(
  method_output,
  feature_sets,
  datasets = NULL,
  latent_dimensions = NULL,
  use_abs = TRUE,
  rank_test = FALSE,
  min_set_size = 5,
  add_missing_features = FALSE,
  mo_data = NULL,
  sets_info_df = NULL,
  col_set = NULL
)
```

## Arguments

- method_output:

  Integration method output generated via the
  [`get_output()`](https://plant-food-research-open.github.io/moiraine/reference/get_output.md)
  function.

- feature_sets:

  Named list, where each element corresponds to a feature set, and
  contains a vector of features ID of all features belonging to that
  set.

- datasets:

  Character vector, the names of the datasets to consider in the
  enrichment analysis. If `NULL` (default value), features from all
  datasets will be included in the analysis.

- latent_dimensions:

  Character vector, the latent dimensions for which an enrichment
  analysis should be performed. If `NULL` (default value), all latent
  dimensions will be analysed.

- use_abs:

  Logical, whether to use the absolute value of the features metric to
  perform the enrichment. If `TRUE` (default value), it allows to
  highlight feature sets in which the features have high
  weight/importance score, both positive and negative. If `FALSE`, it
  will instead highlight feature sets in which the weights all have the
  same sign (coordinated change).

- rank_test:

  Logical, whether a non-parametric Wilcoxon Mann-Whitney test should be
  used instead of the default two-sample t-test (i.e. based on features
  rank rather than their metric). Default value is `FALSE`.

- min_set_size:

  Integer, the minimum number of features in a set required in order to
  compute an enrichment score for the set. Default value is 5.

- add_missing_features:

  Logical, whether features that are in a multi-omics dataset (provided
  through the `mo_data` argument) but don't have a weight in the
  integration results (e.g. because they were not selected in the
  pre-processing step) should be added in the results. If `TRUE`
  (default value), they will be added with an importance score of 0.

- mo_data:

  A
  [`MultiDataSet-class`](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object. If `add_missing_features` is true, all features in the
  multi-omics dataset with no weight in the integration method result
  will be added with an importance score of 0.

- sets_info_df:

  Data-frame, information about the feature sets that will be added to
  the enrichment results. If `NULL` (default value), no information will
  be added to the results.

- col_set:

  Character, name of the column in `sets_info_df` containing the set
  IDs. Should match the names of the `feature_sets` list.

## Value

a tibble of enrichment results.

## Details

When `add_missing_features` is `TRUE` (default behaviour) and a
MultiDataSet object is passed through the `mo_data` argument, features
present in the multi-omics dataset but absent in the integration
method's results will be added to the method's result with a weight of
0. This make sure that if, from a set of 30 features, 25 of these
features were removed during the feature pre-selection stage, the
enrichment considers that these 25 features were not given high weights
by the method. Otherwise, if `add_missing_features` is `FALSE`, these 25
features will be ignored, and so the enrichment analysis may find that
one latent dimension is enriched for this particular set, even though
there only are 5 features out of 30 from the set that contribute to the
latent dimension. Also note that multiple-testing correction is applied
at the latent dimension level, and there is no correction across the
latent dimensions.

When setting `use_abs` to `FALSE`, for each latent dimension, their
enrichment for the features test is tested twice: once for enrichment in
features with positive weight/importance, and once for features with
negative weight/importance score. This will be indicated in the
`direction` column of the resulting tibble.

Note that we built this function using the gage vignette on [RNA-Seq
Data Pathway and Gene-set Analysis
Workflow](https://bioconductor.org/packages/release/bioc/vignettes/gage/inst/doc/RNA-seqWorkflow.pdf),
section 7.1.
