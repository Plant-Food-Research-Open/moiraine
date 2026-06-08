# Assess optimal number of components for sPLS-DA on omics dataset from MultiDataSet object

Performs cross-validation for a PLS-DA run (implemented in the
`mixOmics` package) on an omics dataset from a `MultiDataSet` object.
This allows to estimate the optimal number of latent components to
construct. This is intended for feature preselection in the omics
dataset (see examples below).

## Usage

``` r
perf_splsda(
  splsda_input,
  ncomp_max = 5,
  validation = "Mfold",
  folds = 5,
  nrepeat = 50,
  measure = "BER",
  distance = "centroids.dist",
  cpus = 1,
  progressBar = TRUE,
  seed = NULL
)
```

## Arguments

- splsda_input:

  Input for the sPLS-DA functions from mixOmics, created with
  [`get_input_splsda()`](https://plant-food-research-open.github.io/moiraine/reference/get_input_splsda.md).

- ncomp_max:

  Integer, the maximum number of latent components to test when
  estimating the number of latent components to use. Default value is
  `5`.

- validation:

  Character, which cross-validation method to use, can be one of
  `"Mfold"` or `"loo"` (see
  [`mixOmics::perf()`](https://rdrr.io/pkg/mixOmics/man/perf.html)).
  Default value is `"Mfold"`.

- folds:

  Integer, number of folds to use in the M-fold cross-validation (see
  [`mixOmics::perf()`](https://rdrr.io/pkg/mixOmics/man/perf.html)).
  Default value is 5.

- nrepeat:

  Integer, number of times the cross-validation is repeated (see
  [`mixOmics::perf()`](https://rdrr.io/pkg/mixOmics/man/perf.html)).

- measure:

  Performance measure used to select the optimal value of `ncomp`, can
  be one of `"BER"` or `"overall"` (see
  [`mixOmics::perf()`](https://rdrr.io/pkg/mixOmics/man/perf.html)).
  Default value is `"BER"`.

- distance:

  Distance metric used to select the optimal value of `ncomp`, can be
  one of `"max.dist"`, `"centroids.dist"` or `"mahalanobis.dist"` (see
  [`mixOmics::perf()`](https://rdrr.io/pkg/mixOmics/man/perf.html)).
  Default value is `"centroids.dist"`.

- cpus:

  Integer, number of cpus to use.

- progressBar:

  Logical, whether to display a progress bar during the optimisation of
  `ncomp`. Default value is `TRUE`.

- seed:

  Integer, seed to use. Default is `NULL`, i.e. no seed is set inside
  the function.

## Value

A list as per the output of the
[`mixOmics::perf()`](https://rdrr.io/pkg/mixOmics/man/perf.html)
function, with the following additional elements:

- `dataset_name`: the name of the dataset analysed;

- `group`: column name in the samples information data-frame used as
  samples group;

- `optim_ncomp`: the optimal number of latent components as per the
  `measure` and `distance` specified;

- `optim_measure`: the measure used to select the optimal number of
  latent components;

- `optim_distance`: the distance metric used to select the optimal
  number of latent components. In addition, the name of the dataset
  analysed and the column name in the samples information data-frame
  used as samples group as stored as attributes `dataset_name` and
  `group`, respectively.

## Details

This function uses the
[`mixOmics::plsda()`](https://rdrr.io/pkg/mixOmics/man/plsda.html) and
[`mixOmics::perf()`](https://rdrr.io/pkg/mixOmics/man/perf.html)
function from the `mixOmics` package.
