% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multidataset_utils.R
\name{get_dataset_matrix}
\alias{get_dataset_matrix}
\title{Get multi-omics dataset as matrix}
\usage{
get_dataset_matrix(mo_data, dataset_name, keep_dataset_name = FALSE)
}
\arguments{
\item{mo_data}{A \link[MultiDataSet:MultiDataSet-class]{MultiDataSet::MultiDataSet} object.}

\item{dataset_name}{Character, name of the omics dataset to extract.}

\item{keep_dataset_name}{Logical, should the dataset name be stored in the
\code{'dataset_name'} attribute of the resulting matrix? Default value is
\code{FALSE}.}
}
\value{
A matrix of measurements with features as rows and samples as
columns. The name of the dataset is stored in the \code{'dataset_name'}
attribute if \code{keep_dataset_name} is \code{TRUE}.
}
\description{
Extracts an omics dataset as a matrix of measurements from a MultiDataSet
object.
}
\examples{
\dontrun{
## mo_data is a MultiDataSet object with a dataset called "rnaseq"
mat <- get_dataset_matrix(mo_data, "rnaseq", keep_dataset_name = TRUE)
## with keep_dataset_name = TRUE, can recover dataset name as follows:
attr(mat, "dataset_name")
}
}
