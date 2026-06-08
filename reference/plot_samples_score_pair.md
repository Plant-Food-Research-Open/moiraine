# Plots sample scores as a scatterplot

Plots the samples score from a pair of latent dimensions from a
dimension reduction analysis as a scatterplot.

## Usage

``` r
plot_samples_score_pair(
  method_output,
  latent_dimensions,
  mo_data = NULL,
  colour_by = NULL,
  shape_by = NULL
)
```

## Arguments

- method_output:

  A single Integration method output generated via the
  [`get_output()`](https://plant-food-research-open.github.io/moiraine/reference/get_output.md)
  function, or a list of two integration method outputs.

- latent_dimensions:

  Character vector of length 2 (if `method_output` is a single output
  object), or named list of length 2 (if `method_output` is a list of
  two output objects). In the first case, gives the name of the two
  latent dimensions that should be represented. In the second case, the
  names of the list should correspond to the names of the methods, and
  the values should each be a character giving for the corresponding
  method the name of the latent dimension to display.

- mo_data:

  A `MultiDataSet` object (will be used to extract samples information).

- colour_by:

  Character, name of column in one of the samples metadata tables from
  `mo_data` to use for colouring samples in the plot. Default value is
  `NULL`.

- shape_by:

  Character, name of column in one of the samples metadata tables from
  `mo_data` to use a shape for the samples in the plot. Default value is
  `NULL`.

## Value

A ggplot.

## Examples

``` r
if (FALSE) { # \dontrun{
## Let diablo_res be the output from a DIABLO analysis, and mofa_res the
## output from a MOFA analysis. Let mo_set be the corresponding MultiDataSet object.
output_diablo <- get_output_diablo(diablo_res)
output_mofa <- get_output_mofa2(mofa_res)

## Scatterplot of the first two DIABLO components
plot_samples_score_pair(output_diablo, c("Component 1", "Component 2"))

## Adding samples information to the plot - here 'Time' and 'Treatment' should
## two columns in the samples metadata of one of the datasets in mo_set
plot_samples_score_pair(
  output_diablo,
  c("Component 1", "Component 2"),
  mo_data <- mo_set,
  colour_by = "Time",
  shape_by = "Treatment"
)

## Comparing the first MOFA factor to the first DIABLO component
plot_samples_score_pair(
  list(output_diablo, output_mofa),
  list("DIABLO" = "Component 1", "MOFA" = "Factor 1"),
  mo_data <- mo_set,
  colour_by = "Time",
  shape_by = "Treatment"
)

## Giving custom names to the methods
plot_samples_score_pair(
  list("DIABLO prefiltered" = output_diablo, "MOFA full" = output_mofa),
  list("DIABLO prefiltered" = "Component 1", "MOFA full" = "Factor 1")
)
} # }
```
