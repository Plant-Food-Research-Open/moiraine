% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transformation.R
\name{transform_dataset}
\alias{transform_dataset}
\title{Applies a transformation to a dataset from a MultiDataSet object}
\usage{
transform_dataset(
  mo_data,
  dataset,
  transformation,
  return_multidataset = FALSE,
  return_matrix_only = FALSE,
  verbose = TRUE,
  log_base = 2,
  pre_log_function = zero_to_half_min,
  method,
  ...
)
}
\arguments{
\item{mo_data}{A \code{\link[MultiDataSet]{MultiDataSet-class}} object.}

\item{dataset}{Character, name of the dataset to transform.}

\item{transformation}{Character, transformation to be applied. Possible
values are: \code{vsn}, \code{vst-deseq2}, \code{logx} \code{best-normalize-auto}
or \code{best-normalize-manual}. See \code{Details}.}

\item{return_multidataset}{Logical, should a \code{MultiDataSet} object with the
original data replaced by the transformed data returned? If \code{FALSE}, the
output of the function depends on \code{return_matrix_only}. Default value is
\code{FALSE}.}

\item{return_matrix_only}{Logical, should only the transformed matrix be
returned? If \code{TRUE}, the function will return a matrix. If \code{FALSE}, the
function instead returns a list with the transformed data as well as
other information relevant to the transformation. Ignored if
\code{return_multidataset} is \code{TRUE}. Default value is \code{FALSE}.}

\item{verbose}{Logical, should information about the transformation be
printed? Default value is \code{TRUE}.}

\item{log_base}{Numeric, the base with respect to which logarithms are
computed. Default value is \code{2}. Only used if \code{transformation = 'logx'}.}

\item{pre_log_function}{Function that will be applied to the matrix before
the log transformation (e.g. to apply an offset to the values to avoid
issues with zeros). Default value is the \code{\link[=zero_to_half_min]{zero_to_half_min()}} function.
Only used if \code{transformation = 'logx'}.}

\item{method}{Character, if \code{transformation = 'best-normalize-manual'}, which
normalisation method should be applied. See possible values in
\code{\link[=transform_bestNormalise_manual]{transform_bestNormalise_manual()}}. Ignored for other transformations.}

\item{...}{Further arguments passed to the \code{\link[bestNormalize:bestNormalize]{bestNormalize::bestNormalize()}}
function or the \code{method} function from the \code{bestNormalize} package.}
}
\value{
\itemize{
\item if \code{return_multidataset = TRUE}: a \link[MultiDataSet:MultiDataSet-class]{MultiDataSet::MultiDataSet}
object, in which the original data for the transformed dataset has been
replaced.
\item if \code{return_multidataset = FALSE} and \code{return_matrix_only = TRUE}: a matrix
with the transformed data.
\item if \code{return_multidataset = FALSE} and \code{return_matrix_only = FALSE}: a list
with two elements, \code{transformed_data} containing a matrix of transformed
data, and \code{info_transformation} containing information about the
transformation (depends on the transformation applied).
}
}
\description{
Applies a transformation to a dataset from a \code{MultiDataSet} object.
Implemented transformations are: Variance Stabilising Normalisation (from the
\code{vsn} package), Variance Stabilising Transformation (from the \code{DESeq2}
package - only for count data), and appropriate feature-wise normalisation
through the \code{BestNormalise} package.
}
\details{
Currently implemented transformations and recommendations based on dataset
type:
\itemize{
\item \code{vsn}: Variance Stabilising normalisation, implemented in the
\code{\link[vsn:justvsn]{vsn::justvsn()}} function from the \code{vsn} package. This method was originally
developed for microarray intensities. This transformation is recommended for
microarray, metabolome, chemical or other intensity-based datasets. In
practice, applies the \code{\link[=transform_vsn]{transform_vsn()}} function.
\item \code{vst-deseq2}: Variance Stabilising Transformation, implemented in the
\code{\link[DESeq2:varianceStabilizingTransformation]{DESeq2::varianceStabilizingTransformation()}} function from the \code{DESeq2}
package. This method is applicable to count data only. This transformation is
recommended for RNAseq or similar count-based datasets. In practice, applies
the \code{\link[=transform_vst]{transform_vst()}} function.
\item \code{logx}: log-transformation (default to log2, but base can be specified).
In practice, applies the \code{\link[=transform_logx]{transform_logx()}} function.
\item \code{best-normalize-auto}: most appropriate normalisation method automatically
selected from a number of options, implemented in the
\code{\link[bestNormalize:bestNormalize]{bestNormalize::bestNormalize()}} function from the \code{bestNormalize} package.
This transformation is recommended for phenotypes that are each measured on
different scales (since the transformation method selected will potentially
be different across the features), preferably with a reasonable number of
features (less than 100) to avoid large computation times. In practice,
applies the \code{\link[=transform_bestNormalise_auto]{transform_bestNormalise_auto()}} function.
\item \code{best-normalize-manual}: performs the same transformation (specified
through the \code{method} argument) to each feature of a dataset. This
transformation is recommended for phenotypes data in which the different
phenotypes are measured on the same scale. The different normalisation
methods are:
\itemize{
\item \code{"arcsinh_x"}: data is transformed as \code{log(x + sqrt(x^2 + 1))};
\item \code{"boxcox"}: Box Cox transformation;
\item \code{"center_scale"}: data is centered and scaled;
\item \code{"exp_x"}: data is transformed as \code{exp(x)};
\item \code{"log_x"}: data is transformed as \code{log_b(x+a)} (\code{a} and \code{b} either
selected automatically per variable or passed as arguments);
\item \code{"orderNorm"}: Ordered Quantile technique;
\item \code{"sqrt_x"}: data transformed as \code{sqrt(x + a)} (\code{a} selected automatically
per variable or passed as argument),
\item \code{"yeojohnson"}: Yeo-Johnson transformation.
}
}
}
