% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prefiltering.R
\name{feature_preselection_mad_factory}
\alias{feature_preselection_mad_factory}
\title{Target factory for feature preselection based on Median Absolute Deviation}
\usage{
feature_preselection_mad_factory(
  mo_data_target,
  to_keep_ns,
  to_keep_props = NULL,
  with_ties = TRUE,
  target_name_prefix = "",
  filtered_set_target_name = NULL
)
}
\arguments{
\item{mo_data_target}{Symbol, the name of the target containing the \code{MultiDataSet} object.}

\item{to_keep_ns}{Named integer vector, the number of feature to retain in each dataset to be prefiltered
(names should correspond to a dataset name). Value should be less than the number of features in the
corresponding dataset. Set to \code{NULL} in order to use \code{to_keep_props} instead.}

\item{to_keep_props}{Named numeric vector, the proportion of features to retain in each dataset
to be prefiltered (names should correspond to a dataset name). Value should be > 0 and < 1.
Will be ignored if \code{to_keep_ns} is not \code{NULL}.}

\item{with_ties}{Should ties be kept together? If \code{TRUE}, may return more features than requested. Default value is \code{TRUE}.}

\item{target_name_prefix}{Character, a prefix to add to the name of the targets created by this target factory.
Default value is \code{""}.}

\item{filtered_set_target_name}{Character, the name of the final target containing the filtered \code{MultiDataSet} object.
If NULL, a name will automatically be supplied. Default value is \code{NULL}.}
}
\value{
A list of target objects. With \code{target_name_prefix = ""} and \code{filtered_set_target_name = NULL},
the following targets are created:
\itemize{
\item \code{mad_spec}: a target that generates a grouped tibble where each row corresponds to one dataset to be filtered,
with the columns specifying each dataset name, and associated values from \code{to_keep_ns}, \code{to_keep_props} and \code{with_ties}.
\item \code{mad_mat}: a dynamic branching target that run the \code{\link[=get_dataset_matrix]{get_dataset_matrix()}} function for each dataset.
\item \code{individual_mad_values}: a dynamic branching target that runs the \code{\link[=select_features_mad_matrix]{select_features_mad_matrix()}} function for each dataset.
\item \code{filtered_set_mad}: a target to retain from the original \code{MultiDataSet} object only features selected based on their MAD values.
}
}
\description{
Creates a list of targets to perform feature preselection on datasets from a \code{MultiDataSet}
object by retaining features with the highest Median Absolute Deviation (MAD).
}
\examples{
\dontrun{
## in the _targets.R
library(moiraine)

list(
  ## add code here to load the different datasets

  ## the following target creates a MultiDataSet object from previously
  ## created omics sets (geno_set, trans_set, etc)
  tar_target(
    mo_set,
    create_multiomics_set(geno_set, trans_set, metabo_set, pheno_set)
  ),
  feature_preselection_mad_factory(
    mo_set,
    to_keep_ns = c("rnaseq" = 1000, "metabolome" = 500),
    filtered_set_target_name = "mo_set_filtered"
  ),

  ## Another example using to_keep_props
  feature_preselection_mad_factory(
    mo_set,
    to_keep_ns = NULL,
    to_keep_props = c("rnaseq" = 0.3, "metabolome" = 0.5),
    filtered_set_target_name = "mo_set_filtered"
  )
)
}
}
