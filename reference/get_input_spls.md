# Generate sPLS input data (for mixomics)

Creates an object that can be used as input for the (s)PLS functions
from the mixOmics package. It contains the two omics datasets selected,
restricted to common samples.

## Usage

``` r
get_input_spls(mo_data, mode, datasets = names(mo_data), multilevel = NULL)
```

## Arguments

- mo_data:

  A
  [MultiDataSet::MultiDataSet](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

- mode:

  Character, the mode of PLS to use of the analysis (see [sPLS
  documentation](http://mixomics.org/methods/spls/)). Should be one of
  `'regression'`, `'canonical'`, `'invariant'` or `'classic'`.

- datasets:

  Character vector of length 2, the names of the datasets from `mo_data`
  to include in the analysis.

- multilevel:

  Character vector of length 1 or 3 to be used as information about
  repeated measurements. See Details. Default value is `NULL`, i.e. the
  multilevel option will not be used.

## Value

A list, in which each element corresponds to one omics dataset, with
samples as rows and features as columns. The mode to use for the
analysis is stored in the `mode` attribute of the returned object.

## Details

`multilevel` argument: enables the multilevel option (see [mixOmics
site](http://mixomics.org/methods/multilevel/)) to deal with repeated
measurements.
[`mixOmics::spls()`](https://rdrr.io/pkg/mixOmics/man/spls.html) enables
one- and two-factor decomposition. For one-factor decomposition,
`multilevel` argument should be the name of the column in the samples
metadata that gives the ID of the observation units (e.g. the ID of the
subjects that were measured several times). The resulting design matrix
(stored in the `multilevel` argument of the returned object) will be a
data-frame with one column which gives the ID (as integer) of the
observation units corresponding to each sample in the omics datasets.
For two-factor decomposition, `multilevel` should be of length 3. The
first value, similarly to the one-factor decomposition, should be the
name of the column in the samples metadata that gives the ID of the
observation units (e.g. the ID of the subjects that were measured
several times). The second and third values should be the name of the
columns in the samples metadata that give the two factors considered.
The resulting design matrix (stored in the `multilevel` argument of the
returned object) will be a data-frame with three columns: the first
column gives the ID (as integer) of the observation units corresponding
to each sample in the omics datasets; the second and third columns give
the levels of the two factors.
