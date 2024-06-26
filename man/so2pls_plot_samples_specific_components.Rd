% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/so2pls.R
\name{so2pls_plot_samples_specific_components}
\alias{so2pls_plot_samples_specific_components}
\title{Plots sO2PLS specific components samples scores}
\usage{
so2pls_plot_samples_specific_components(so2pls_res, dataset = NULL, ...)
}
\arguments{
\item{so2pls_res}{The output from the \code{\link[OmicsPLS]{o2m}} function.}

\item{dataset}{Character, the name of the dataset for which the specific
components should be plotted. Default is \code{NULL}, i.e. the specific components
of both datasets are plotted.}

\item{...}{Further arguments passed to \code{\link[=plot_samples_score]{plot_samples_score()}}.}
}
\value{
A list of ggmatrix plots (one per dataset), or one plot if \code{dataset}
was used to specify a dataset.
}
\description{
Plots the samples scores for the datasets specific components from an
sO2PLS run (from the \code{OmicsPLS} package).
}
