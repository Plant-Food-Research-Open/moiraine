# Plots the correlation between factors and covariates for MOFA

Plots the Pearson correlation between MOFA latent factors and covariates
(obtained from the samples metadata). This function provides the
`ggplot2` version of the plot created with
[`correlate_factors_with_covariates`](https://rdrr.io/pkg/MOFA2/man/correlate_factors_with_covariates.html)
(with the `plot` parameter set to `'r'`).

## Usage

``` r
mofa_plot_cor_covariates(
  mofa_output,
  covariates = NULL,
  show_cor = TRUE,
  min_show_cor = 0.2,
  round_cor = 2,
  factor_as_col = TRUE
)
```

## Arguments

- mofa_output:

  The output from
  [`run_mofa`](https://rdrr.io/pkg/MOFA2/man/run_mofa.html).

- covariates:

  Character vector, the covariates to use for the plot. If `NULL`, all
  covariates retrieved via
  `colnames(MOFA2::samples_metadata(mofa_output))` (except `group`, `id`
  and `sample`) will be used. Default value is `NULL`.

- show_cor:

  Logical, should the correlation values be added to the plot? Default
  value is `TRUE`.

- min_show_cor:

  Numeric, the minimum value below which correlation coefficients values
  are not added to the plot (i.e. only a circle will appear for these
  values but no text). Ignored if `show_cor` is `FALSE`. Default value
  is 0.2.

- round_cor:

  Integer, how many decimal places to show for the correlation
  coefficients. Ignored if `show_cor` is `FALSE`. Default value is 2.

- factor_as_col:

  Logical, should the factors be represented in columns? If `FALSE`,
  will be represented in rows instead. Default value is `TRUE`.

## Value

A ggplot.
