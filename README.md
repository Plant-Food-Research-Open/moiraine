
<!-- README.md is generated from README.Rmd. Please edit that file -->

# moiraine <img src="man/figures/logo.png" align="right" alt="" width="120" />

<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
<!-- badges: end -->

`moiraine` is a package for facilitating the construction of a
reproducible analysis pipeline for multi-omics data integration. It
provides functions to automate data import, pre-processing,
transformation, integration through several tools, as well as
interpretation and comparison of the integration results. It relies on
the [targets](https://books.ropensci.org/targets/) package to generate
reproducible workflows.

## Overview

The workflow for a typical multi-omics integration analysis handled with
`moiraine` includes the following steps:

- Data import: this covers the import of omics measurements as well as
  associated metadata (i.e. information about the omics features and
  samples) – moiraine relies on the [`MultiDataSet`
  package](https://bioconductor.org/packages/release/bioc/html/MultiDataSet.html)
  to store this information in a consistent format;

- Inspection of the omics datasets: including checking values density
  distribution, samples overlap between omics datasets, or presence of
  missing values;

- Preprocessing of the omics datasets: missing values imputation,
  transformation, and pre-filtering of samples and omics features;

- Integration of the omics datasets by one or more of the supported
  tools; currently, the following integration methods are covered in
  `moiraine`:

  - sPLS and DIABLO from the `mixOmics` package

  - sO2PLS from the `OmicsPLS` package

  - MOFA and MEFISTO from the `MOFA2` package

- Interpretation of the integration results using standardised
  visualisations enriched with features and samples metadata;

- Comparison of the integration results obtained by different methods or
  pre-processing approaches.

An overview of the capabilities of the package is available
[here](https://solid-lamp-kq546rq.pages.github.io/overview.html).

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

The [user manual](https://solid-lamp-kq546rq.pages.github.io/) provides
an in-depth walk-through of a multi-omics integration analysis with the
package.

Before using `moiraine`, we encourage you to get familiar with the
`targets` package; the [`targets`
manual](https://books.ropensci.org/targets/) is a great place to start.
