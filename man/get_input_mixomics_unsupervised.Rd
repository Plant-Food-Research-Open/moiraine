% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mixomics.R
\name{get_input_mixomics_unsupervised}
\alias{get_input_mixomics_unsupervised}
\title{Generate mixomics input data for unsupervised methods}
\usage{
get_input_mixomics_unsupervised(mo_data, datasets = names(mo_data))
}
\arguments{
\item{mo_data}{A \code{\link[MultiDataSet]{MultiDataSet-class}} object.}

\item{datasets}{Character vector, the names of the datasets from \code{mo_data} to include in the analysis.}
}
\value{
A list, in which each element corresponds to one omics dataset, with samples as rows and features as columns.
}
\description{
Creates an object that can be used as input for the MixOmics package.
It contains the omics datasets restricted to common samples.
}
