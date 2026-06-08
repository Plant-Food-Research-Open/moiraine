# Extract/plot running time of functions used for different integration methods

In the plot, the coloured blocks in each bar each represent one of the
steps (targets) that correspond to an integration method. The text above
each bar shows the function that took the longest to run for each
integration method.

## Usage

``` r
get_targets_running_time(
  target_patterns = c(sPLS = "^spls_", sO2PLS = "^so2pls_", MOFA = "^mofa_", DIABLO =
    "^diablo_"),
  patterns_to_methods = c("sPLS", "sO2PLS", "MOFA", "DIABLO"),
  target_exclude_patterns = NULL,
  method_functions = get_method_functions()
)

plot_running_time(
  target_patterns = c(sPLS = "^spls_", sO2PLS = "^so2pls_", MOFA = "^mofa_", DIABLO =
    "^diablo_"),
  patterns_to_methods = c("sPLS", "sO2PLS", "MOFA", "DIABLO"),
  target_exclude_patterns = NULL,
  method_functions = get_method_functions()
)
```

## Arguments

- target_patterns:

  Named character vector, regex patterns used to extract targets by
  their names.

- patterns_to_methods:

  Character vector of same length as `target_patterns`, gives the
  integration method to which each pattern from `target_patterns`
  corresponds.

- target_exclude_patterns:

  Character vector, regex used to exclude some targets from the
  analysis.

- method_functions:

  Named list of character vector, giving for each integration method
  (i.e. values in `patterns_to_methods`) the name of the functions used
  in the pipeline to be considered when computing the running time.
  Names should match values in `patterns_to_methods`.

## Value

Either a tibble of running time for each function used as part of each
integration method (`get_targets_running_time`), or a ggplot
(`plot_running_time`).

## Examples

``` r
if (FALSE) { # \dontrun{
## By default uses all targets whose name starts with "spls_" as part of the
## "sPLS" method, etc
get_targets_running_time()
plot_running_time()

## If we ran two versions of MOFA, e.g. one with supervised input and the
## other with unsupervised input
get_targets_running_time(target_exclude_patterns = "^mofa_unsupervised")
plot_running_time(target_exclude_patterns = "^mofa_unsupervised")

## If instead we want to compare the two MOFA run times
get_targets_running_time(
  target_patterns = c(
    "MOFA supervised" = "^mofa_(?!=unsupervised)",
    "MOFA unsupervised" = "^mofa_unsupervised"
  ),
  patterns_to_methods = c("MOFA", "MOFA")
)
plot_running_time(
  target_patterns = c(
    "MOFA supervised" = "^mofa_(?!=unsupervised)",
    "MOFA unsupervised" = "^mofa_unsupervised"
  ),
  patterns_to_methods = c("MOFA", "MOFA")
)
} # }
```
