% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diablo.R
\name{diablo_tune}
\alias{diablo_tune}
\title{Tunes keepX arg for DIABLO}
\usage{
diablo_tune(
  mixomics_data,
  design_matrix,
  keepX_list = NULL,
  cpus = NULL,
  seed = NULL,
  ...
)
}
\arguments{
\item{mixomics_data}{A \code{mixOmics} input object created with
\code{get_input_mixomics_supervised()}.}

\item{design_matrix}{Either numeric matrix created through
\code{diablo_generate_design_matrix}, or character (accepted values are
\code{'null'}, \code{'weighted_full'}, \code{'full'}). See Details.}

\item{keepX_list}{Named list, gives for each omics dataset in the mixOmics
input (i.e. excluding the response Y) a vector of values to test (i.e.
number of features to return from this dataset). If \code{NULL} (default), a
standard grid will be applied for each dataset and latent component,
testing values: \code{seq(5, 30, 5)}.}

\item{cpus}{Integer, the number of CPUs to use when running the code in
parallel. For advanced users, see the \code{BPPARAM} argument of
\code{\link[mixOmics:tune.block.splsda]{mixOmics::tune.block.splsda()}}.}

\item{seed}{Integer, seed to use. Default is \code{NULL}, i.e. no seed is set
inside the function.}

\item{...}{Arguments to be passed to the \code{\link[mixOmics:tune.block.splsda]{mixOmics::tune.block.splsda()}}
function.}
}
\value{
A list, see \code{\link[mixOmics:tune.block.splsda]{mixOmics::tune.block.splsda()}}.
}
\description{
Performs cross-validation to estimate the optimal number of features to
retain from each dataset for a DIABLO run.
}
\details{
The \verb{design_matrix`` argument can either be a custom design matrix (for example as constructed via the }diablo_generate_design_matrix` function); or
a character indicating the type of design matrix to generate. Possible values
include:
\itemize{
\item \code{'null'}: Off-diagonal elements of the design matrix are set to 0;
\item \code{'weighted_full'}: Off-diagonal elements of the design matrix are set to
0.1;
\item \code{'full'}: Off-diagonal elements of the design matrix are set to 1.
}
}
