% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multidataset_utils.R
\name{add_samples_metadata}
\alias{add_samples_metadata}
\title{Adding data-frame to samples metadata}
\usage{
add_samples_metadata(mo_data, df, datasets = NULL)
}
\arguments{
\item{mo_data}{A \link[MultiDataSet:MultiDataSet-class]{MultiDataSet::MultiDataSet} object.}

\item{df}{A tibble or data-frame of samples information, with at least column
\code{id} (giving the sample IDs), and one row per sample ID.}

\item{datasets}{Character vector, name of the datasets to which the samples
information should be added. If \code{NULL} (default value), the information
will be added to the samples metadata of all datasets.}
}
\value{
A MultiDataSet object, with info from \code{df} added to its
corresponding samples metadata.
}
\description{
Adds information from a data-frame to the samples metadata of a MultiDataSet
object.
}
