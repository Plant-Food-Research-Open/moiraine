% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prefiltering.R
\name{run_splsda}
\alias{run_splsda}
\title{Performs sPLS-DA on omics dataset from MultiDataSet object}
\usage{
run_splsda(
  splsda_input,
  perf_res,
  to_keep_n = NULL,
  to_keep_prop = NULL,
  ncomp = NULL,
  seed = NULL
)
}
\arguments{
\item{splsda_input}{Input for the sPLS-DA functions from mixOmics, created
with \code{\link[=get_input_splsda]{get_input_splsda()}}.}

\item{perf_res}{Result of the \code{\link[=perf_splsda]{perf_splsda()}} function. If not supplied,
sPLS-DA will be run on dataset specified by argument \code{dataset_name} with
number of latent components specified by argument \code{comp}.}

\item{to_keep_n}{Integer, the number of features to retain in the dataset.
Should be less than the number of features in the dataset. If \code{NULL} or
\code{NA}, \code{to_keep_prop} will be used instead.}

\item{to_keep_prop}{Numeric, the proportion of features to retain in the
dataset. Will be ignored if \code{to_keep_n} is supplied. Value should be > 0
and < 1.}

\item{ncomp}{Integer, number of latent components to construct. Ignored if
\code{perf_res} is supplied. Default value is \code{NULL}.}

\item{seed}{Integer, seed to use. Default is \code{NULL}, i.e. no seed is set
inside the function.}
}
\value{
A list as per the output of the \code{\link[mixOmics:splsda]{mixOmics::splsda()}} function.
}
\description{
Performs a sPLS-DA (implemented in the \code{mixOmics}) package on a omics dataset
from a MultiDataSet object. This is intended for feature preselection in the
omics dataset (see \code{\link[=get_filtered_dataset_splsda]{get_filtered_dataset_splsda()}}).
}
\details{
This function uses the \code{\link[mixOmics:plsda]{mixOmics::plsda()}} function from the \code{mixOmics}
package. Note that the sPLS-DA method can select the same feature for several
latent components, so the number of features retained for a dataset might be
less than the number specified in the \code{to_keep_n} argument.
}
