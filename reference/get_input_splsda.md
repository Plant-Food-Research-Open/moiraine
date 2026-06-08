# Generate sPLS-DA input data (for mixomics)

Creates an object that can be used as input for the (s)PLS-DA functions
from the mixOmics package. It contains the omics dataset as well as the
samples group membership in a list.

## Usage

``` r
get_input_splsda(mo_data, dataset_name, group, multilevel = NULL)
```

## Arguments

- mo_data:

  A
  [`MultiDataSet-class`](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

- dataset_name:

  Character, name of the dataset from `mo_data` to analyse.

- group:

  Character, the column name in the samples information data-frame to
  use as samples group (use
  [`get_samples_metadata`](https://plant-food-research-open.github.io/moiraine/reference/get_samples_metadata.md)
  to view the samples information data-frame for a omics dataset).

- multilevel:

  Character vector of length 1 or 3 to be used as information about
  repeated measurements. See Details. Default value is `NULL` (no
  repeated measurements).

## Value

A list, in which the first element corresponds to the omics dataset,
with samples as rows and features as columns, and the second element
(named `'Y'`) is a named factor vector, giving for each sample its
group. The name of the dataset to be analysed is stored in the
`dataset_name` attribute of the returned object.

## Details

`multilevel` argument: enables the multilevel option (see [mixOmics
site](http://mixomics.org/methods/multilevel/)) to deal with repeated
measurements.
[`mixOmics::splsda()`](https://rdrr.io/pkg/mixOmics/man/splsda.html)
enables one- and two-factor decomposition. For one-factor decomposition,
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
