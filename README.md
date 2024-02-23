
<!-- README.md is generated from README.Rmd. Please edit that file -->

# moiraine

<!-- badges: start -->

[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
<!-- badges: end -->

`moiraine` is a package for facilitating the construction of a
reproducible analysis pipeline for multi-omics data integration. It
provides functions to automate data import, pre-processing,
transformation, and integration through several tools. It relies on the
[targets](https://books.ropensci.org/targets/) package to generate
reproducible workflows. `moiraine` currently supports multi-omics data
integration through:

- sPLS and DIABLO from the mixOmics package;

- sO2PLS from the omicsPLS package;

- MOFA and MEFISTO from the MOFA2 package.

## Installation

You can install the development version of moiraine from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("PlantandFoodResearch/moiraine")
```

## Example

To get started, create a new analysis pipeline in your working directory
with:

``` r
library(moiraine)

create_moiraine_pipeline()
```

Before using `moiraine`, we encourage you to get familiar with the
`targets` package; the [manual](https://books.ropensci.org/targets/) is
a great place to start.