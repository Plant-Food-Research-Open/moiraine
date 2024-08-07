% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multidataset_pca.R
\name{run_pca}
\alias{run_pca}
\title{Run PCA on MultiDataSet}
\usage{
run_pca(
  mo_data,
  dataset_name,
  n_pcs = 10,
  scale = "none",
  center = TRUE,
  method = NULL
)
}
\arguments{
\item{mo_data}{A \link[MultiDataSet:MultiDataSet-class]{MultiDataSet::MultiDataSet} object.}

\item{dataset_name}{Character, name of the omics dataset on which a PCA
should be run.}

\item{n_pcs}{numeric, number of Principal Components to compute. Default
value is 10.}

\item{scale}{character, type of scaling that should be applied to the dataset
before running the PCA. Should be one of \code{'none'}, \code{'pareto'}, \code{'vector'},
\code{'uv'} (see \code{\link[pcaMethods:pca]{pcaMethods::pca()}}). Default value is \code{'none'}.}

\item{center}{boolean, should the dataset be centred prior to running the
PCA? Default value is \code{TRUE}.}

\item{method}{character, type of PCA that should be applied to the dataset.
See \code{\link[pcaMethods:listPcaMethods]{pcaMethods::listPcaMethods()}}. for a list of available methods.
Default value is \code{'svd'} for datasets with no missing value, and \code{'nipals'}
for datasets with missing values.}
}
\value{
A \link[pcaMethods:pcaRes]{pcaMethods::pcaRes} object containing the result from the
PCA analysis. The attribute \code{dataset_name} specifies the name of the
dataset analysed.
}
\description{
Runs a Principal Component Analysis on an omics dataset from a MultiDataSet
object. This is a wrapper function around the \code{\link[=get_dataset_matrix]{get_dataset_matrix()}} and
\code{\link[=run_pca_matrix]{run_pca_matrix()}} functions.
}
\details{
To facilitate the use of dynamic branching with the \code{targets} package, the
\code{dataset_name} attribute of the resulting object is set as the value of the
\code{dataset_name} parameter, and can be accessed via \code{attr(res_pca, "dataset_name")}.
}
