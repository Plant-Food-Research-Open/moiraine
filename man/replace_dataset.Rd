% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multidataset_utils.R
\name{replace_dataset}
\alias{replace_dataset}
\title{Replace matrix dataset within a MultiDataSet object}
\usage{
replace_dataset(mo_data, dataset_name, new_data)
}
\arguments{
\item{mo_data}{A \link[MultiDataSet:MultiDataSet-class]{MultiDataSet::MultiDataSet} object.}

\item{dataset_name}{Character, the name of the dataset for which the matrix
data should be changed.}

\item{new_data}{Matrix, the new data. Should have features as rows and
samples as columns. Rownames should match the corresponding feature IDs,
colnames should match the corresponding sample IDs.}
}
\value{
A MultiDataSet object.
}
\description{
Replaces the matrix for an omics dataset by a new matrix in a MultiDataSet
object.
}
