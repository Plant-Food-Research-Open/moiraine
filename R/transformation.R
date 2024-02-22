#' Applies Variance Stabilising Normalisation (vsn) to matrix
#'
#' Applies the Variance Stabilising Normalisation performed by the `vsn` package
#' via the \code{\link[vsn]{justvsn}} function.
#'
#' @param mat Numeric matrix.
#' @param return_matrix_only Logical, should only the transformed matrix be returned? If `TRUE`,
#' the function will return a matrix. If `FALSE`, the function instead returns a list with the
#' transformed data and potentially other information relevant to the transformation. Default
#' value is `FALSE`.
#' @param ... Further arguments passed to \code{\link[vsn]{vsn2}}.
#' @return Depending on the `return_matrix_only`, either a matrix of transformed data, or
#' a list with the following elements:
#' \itemize{
#' \item `transformed_data`: matrix of the transformed data;
#' \item `info_transformation`: NULL.
#' }
#' @export
transform_vsn <- function(mat, return_matrix_only = FALSE, ...) {
  if (!requireNamespace("vsn", quietly = TRUE)) {
    stop(
      "Package \"vsn\" must be installed to use this function.",
      call. = FALSE
    )
  }

  res <- vsn::justvsn(mat, ...)

  if (return_matrix_only) {
    return(res)
  }

  return(list(
    transformed_data = res,
    info_transformation = NULL,
    transformation = "vsn"
  ))
}

#' Applies Variance Stabilising Transformation (DESeq2) to matrix
#'
#' Applies the Variance Stabilising Transformation (VST) performed by the `DESeq2` package
#' via the \code{\link[DESeq2]{varianceStabilizingTransformation}} function. Includes a
#' size factor normalisation prior to the VST. Only applies to a matrix of count.
#'
#' @param mat Numeric matrix, must contain integers only.
#' @param return_matrix_only Logical, should only the transformed matrix be returned? If `TRUE`,
#' the function will return a matrix. If `FALSE`, the function instead returns a list with the
#' transformed data and potentially other information relevant to the transformation. Default
#' value is `FALSE`.
#' @return Depending on the `return_matrix_only`, either a matrix of transformed data, or
#' a list with the following elements:
#' \itemize{
#' \item `transformed_data`: matrix of the transformed data;
#' \item `info_transformation`: A \code{\link[DESeq2]{DESeqTransform}} object, with details about the transformation.
#' }
#' @export
transform_vst <- function(mat, return_matrix_only = FALSE) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop(
      "Package \"DESeq2\" must be installed to use this function.",
      call. = FALSE
    )
  }

  ## Rather than called the function directly on the matrix, this way ensures that
  ## we keep information about the transformation (e.g. the size factors used)
  object <- DESeq2::DESeqDataSetFromMatrix(mat, data.frame(row.names = colnames(mat)), ~1)
  res <- DESeq2::varianceStabilizingTransformation(object, blind = TRUE)

  res_mat <- res@assays@data[[1]]

  if (return_matrix_only) {
    return(res_mat)
  }

  return(list(
    transformed_data = res_mat,
    info_transformation = res,
    transformation = "vst-deseq2"
  ))
}

#' Applies the bestNormalize function to rows of a matrix
#'
#' Applies an appropriate normalisation method to each feature (row) in a matrix,
#' via the \code{\link[bestNormalize]{bestNormalize}} function from the `bestNormalize`
#' package.
#'
#' @param mat Numeric matrix.
#' @param return_matrix_only Logical, should only the transformed matrix be returned? If `TRUE`,
#' the function will return a matrix. If `FALSE`, the function instead returns a list with the
#' transformed data and potentially other information relevant to the transformation. Default
#' value is `FALSE`.
#' @param ... Further arguments passed to the \code{\link[bestNormalize]{bestNormalize}} function.
#' @return Depending on the `return_matrix_only`, either a matrix of transformed data, or
#' a list with the following elements:
#' \itemize{
#' \item `transformed_data`: matrix of the transformed data;
#' \item `info_transformation`: A named list with one element per feature (row), giving
#' details of the transformation applied to the feature (see output of \code{\link[bestNormalize]{bestNormalize}}).
#' }
#' @export
transform_bestNormalise_auto <- function(mat, return_matrix_only = FALSE, ...) {
  if (!requireNamespace("bestNormalize", quietly = TRUE)) {
    stop(
      "Package \"bestNormalize\" must be installed to use this function.",
      call. = FALSE
    )
  }

  res <- apply(mat, 1, bestNormalize::bestNormalize, ...)

  res_mat <- t(sapply(res, function(x) {
    x$x.t
  }))

  if (return_matrix_only) {
    return(res_mat)
  }

  return(list(
    transformed_data = res_mat,
    info_transformation = res,
    transformation = "best-normalize-auto"
  ))
}

#' Applies a normalisation method from bestNormalize to rows of a matrix
#'
#' Applies an chosen normalisation method to each feature (row) in a matrix,
#' via the `bestNormalize` package.
#'
#' Applies a normalisation method implemented in the `bestNormalize` package.
#' The `method` argument corresponds to the function from the `bestNormalize`
#' package that will be applied to the rows of the matrix. See the
#' \href{https://cran.r-project.org/web/packages/bestNormalize/vignettes/bestNormalize.html}{vignette from the `bestNormalize` package}
#' for more information about the transformations.
#'
#' @param mat Numeric matrix.
#' @param method Character, name of the normalisation method to apply. Possible values are
#' `"arcsinh_x"`, `"boxcox"`, `"center_scale"`, `"exp_x"`, `"log_x"`, `"orderNorm"`,
#' `"sqrt_x"`, `"yeojohnson"`. See Details.
#' @param return_matrix_only Logical, should only the transformed matrix be returned? If `TRUE`,
#' the function will return a matrix. If `FALSE`, the function instead returns a list with the
#' transformed data and potentially other information relevant to the transformation. Default
#' value is `FALSE`.
#' @param ... Further arguments passed to the `method` function from the `bestNormalize` package.
#' @return Depending on the `return_matrix_only`, either a matrix of transformed data, or
#' a list with the following elements:
#' \itemize{
#' \item `transformed_data`: matrix of the transformed data;
#' \item `info_transformation`: A named list with one element per feature (row), giving
#' details of the transformation applied to the feature (see output for the `bestNormalize`
#' function corresponding to `method`).
#' }
#' @export
transform_bestNormalise_manual <- function(mat, method, return_matrix_only = FALSE, ...) {
  if (!requireNamespace("bestNormalize", quietly = TRUE)) {
    stop(
      "Package \"bestNormalize\" must be installed to use this function.",
      call. = FALSE
    )
  }

  poss_methods <- c("arcsinh_x", "boxcox", "log_x", "sqrt_x",
                    "yeojohnson", "center_scale", "exp_x", "orderNorm")
  .check_names(
    method,
    poss_methods,
    "'method' argument: '_W_' not a valid method. Possible methods are: '_C_'."
  )

  res <- apply(mat, 1, function(x) {
    eval(str2expression(paste0("bestNormalize::", method, "(x, ...)")))
  })

  res_mat <- t(sapply(res, function(x) {x$x.t}))

  if (return_matrix_only) {
    return(res_mat)
  }

  return(list(
    transformed_data = res_mat,
    info_transformation = res,
    transformation = paste0("best-normalize-manual-", method)
  ))
}

#' Applies a log-x transformation to matrix
#'
#' Applies a log-x transformation (by default log2) through the [log()]
#' function.
#'
#' @param mat Numeric matrix.
#' @param return_matrix_only Logical, should only the transformed matrix be
#'   returned? If `TRUE`, the function will return a matrix. If `FALSE`, the
#'   function instead returns a list with the transformed data and potentially
#'   other information relevant to the transformation. Default value is `FALSE`.
#' @param log_base Numeric, the base with respect to which logarithms are
#'   computed.
#' @param pre_log_function Function that will be applied to the matrix before
#'   the log transformation (e.g. to apply an offset to the values to avoid
#'   issues with zeros). Default value is the [zero_to_half_min()] function.
#' @returns Depending on the `return_matrix_only`, either a matrix of
#'   transformed data, or a list with the following elements:
#' * `transformed_data`: matrix of the transformed data;
#' * `info_transformation`: a list with the log base used and the function
#'    applied prior to log-transformation.
#' @export
transform_logx <- function(mat,
                           return_matrix_only = FALSE,
                           log_base = 2,
                           pre_log_function = zero_to_half_min) {

  if (is.null(log_base)) {
    stop("`log_base` argument cannot be `NULL`.")
  }
  if (is.null(pre_log_function)) {
    stop("`pre_log_function` argument cannot be `NULL`.")
  }

  mat <- pre_log_function(mat)

  if (any(mat == 0, na.rm = TRUE)) {
    warning("The matrix contains zero values; log-transformation will yield `-Inf`.")
  }

  res_mat <- log(mat, base = log_base)

  if (return_matrix_only) {
    return(res_mat)
  }

  res <- list(
    log_base = log_base,
    pre_log_function = pre_log_function
  )

  return(list(
    transformed_data = res_mat,
    info_transformation = res,
    transformation = paste0("log", log_base)
  ))
}

#' Replace zeros with half-min in matrix
#'
#' Replace zero values in a matrix by half of the minimum non-null value in the
#' matrix.
#'
#' @param mat Numeric matrix.
#' @returns The matrix with zero values replaced.
#' @export
zero_to_half_min <- function(mat) {
  if (!any(mat == 0)) {
    return(mat)
  }

  min_val <- min(mat[mat != 0])
  mat[mat == 0] <- min_val / 2

  return(mat)
}

#' Applies a transformation to a dataset from a MultiDataSet object
#'
#' Applies a transformation to a dataset from a `MultiDataSet` object.
#' Implemented transformations are: Variance Stabilising Normalisation (from the
#' `vsn` package), Variance Stabilising Transformation (from the `DESeq2`
#' package - only for count data), and appropriate feature-wise normalisation
#' through the `BestNormalise` package.
#'
#' Currently implemented transformations and recommendations based on dataset
#' type:
#' * `vsn`: Variance Stabilising normalisation, implemented in the
#' [vsn::justvsn()] function from the `vsn` package. This method was originally
#' developed for microarray intensities. This transformation is recommended for
#' microarray, metabolome, chemical or other intensity-based datasets. In
#' practice, applies the [transform_vsn()] function.
#' * `vst-deseq2`: Variance Stabilising Transformation, implemented in the
#' [DESeq2::varianceStabilizingTransformation()] function from the `DESeq2`
#' package. This method is applicable to count data only. This transformation is
#' recommended for RNAseq or similar count-based datasets. In practice, applies
#' the [transform_vst()] function.
#' * `logx`: log-transformation (default to log2, but base can be specified).
#' In practice, applies the [transform_logx()] function.
#' * `best-normalize-auto`: most appropriate normalisation method automatically
#' selected from a number of options, implemented in the
#' [bestNormalize::bestNormalize()] function from the `bestNormalize` package.
#' This transformation is recommended for phenotypes that are each measured on
#' different scales (since the transformation method selected will potentially
#' be different across the features), preferably with a reasonable number of
#' features (less than 100) to avoid large computation times. In practice,
#' applies the [transform_bestNormalise_auto()] function.
#' * `best-normalize-manual`: performs the same transformation (specified
#' through the `method` argument) to each feature of a dataset. This
#' transformation is recommended for phenotypes data in which the different
#' phenotypes are measured on the same scale. The different normalisation
#' methods are:
#'   * `"arcsinh_x"`: data is transformed as `log(x + sqrt(x^2 + 1))`;
#'   * `"boxcox"`: Box Cox transformation;
#'   * `"center_scale"`: data is centered and scaled;
#'   * `"exp_x"`: data is transformed as `exp(x)`;
#'   * `"log_x"`: data is transformed as `log_b(x+a)` (`a` and `b` either
#'                selected automatically per variable or passed as arguments);
#'   * `"orderNorm"`: Ordered Quantile technique;
#'   * `"sqrt_x"`: data transformed as `sqrt(x + a)` (`a` selected automatically
#'                 per variable or passed as argument),
#'   * `"yeojohnson"`: Yeo-Johnson transformation.
#'
#' @param mo_data A \code{\link[MultiDataSet]{MultiDataSet-class}} object.
#' @param dataset Character, name of the dataset to transform.
#' @param transformation Character, transformation to be applied. Possible
#'   values are: `vsn`, `vst-deseq2`, `logx` `best-normalize-auto`
#'   or `best-normalize-manual`. See `Details`.
#' @param return_multidataset Logical, should a `MultiDataSet` object with the
#'   original data replaced by the transformed data returned? If `FALSE`, the
#'   output of the function depends on `return_matrix_only`. Default value is
#'   `FALSE`.
#' @param return_matrix_only Logical, should only the transformed matrix be
#'   returned? If `TRUE`, the function will return a matrix. If `FALSE`, the
#'   function instead returns a list with the transformed data as well as
#'   other information relevant to the transformation. Ignored if
#'   `return_multidataset` is `TRUE`. Default value is `FALSE`.
#' @param verbose Logical, should information about the transformation be
#'   printed? Default value is `TRUE`.
#' @param log_base Numeric, the base with respect to which logarithms are
#'   computed. Default value is `2`. Only used if `transformation = 'logx'`.
#' @param pre_log_function Function that will be applied to the matrix before
#'   the log transformation (e.g. to apply an offset to the values to avoid
#'   issues with zeros). Default value is the [zero_to_half_min()] function.
#'   Only used if `transformation = 'logx'`.
#' @param method Character, if `transformation = 'best-normalize-manual'`, which
#'   normalisation method should be applied. See possible values in
#'   [transform_bestNormalise_manual()]. Ignored for other transformations.
#' @param ... Further arguments passed to the [bestNormalize::bestNormalize()]
#'   function or the `method` function from the `bestNormalize` package.
#' @returns
#' * if `return_multidataset = TRUE`: a [MultiDataSet::MultiDataSet-class]
#' object, in which the original data for the transformed dataset has been
#' replaced.
#' * if `return_multidataset = FALSE` and `return_matrix_only = TRUE`: a matrix
#' with the transformed data.
#' * if `return_multidataset = FALSE` and `return_matrix_only = FALSE`: a list
#' with two elements, `transformed_data` containing a matrix of transformed
#' data, and `info_transformation` containing information about the
#' transformation (depends on the transformation applied).
#' @export
transform_dataset <- function(mo_data,
                              dataset,
                              transformation,
                              return_multidataset = FALSE,
                              return_matrix_only = FALSE,
                              verbose = TRUE,
                              log_base = 2,
                              pre_log_function = zero_to_half_min,
                              method,
                              ...) {
  ## We don't want to subset the dataset
  check_is_multidataset(mo_data)
  .check_names(
    dataset,
    names(mo_data),
    "'_W_' datasets are not present in mo_data. Possible dataset names are: '_C_'."
  )

  .check_names(
    transformation,
    c("vsn", "vst-deseq2", "logx", "best-normalize-auto", "best-normalize-manual"),
    "'transformation' argument: '_W_' is not a recognised transformation. Possible values are: '_C_'."
  )

  if (transformation == "best-normalize-manual" & missing(method)) {
    stop("'method' argument should be provided for 'best-normalize-manual' transformation.")
  }

  ## If returning a MultiDataSet object, cannot keep additional information about the transformation
  if (return_multidataset) return_matrix_only <- TRUE

  ## get the dataset
  mat <- get_datasets(mo_data)[[dataset]]

  if (verbose) {
    transf_name <- c(
      "vsn" = "Variance Stabilising Normalisation (vsn)",
      "vst-deseq2" = "Variance Stabilising Transformation (DESeq2)",
      "logx" = "Log Transformation",
      "best-normalize-auto" = "automatic normalisation selection (bestNormalize)",
      "best-normalize-manual" = " transformation (bestNormalize)"
    )

    transf_name_i <- transf_name[transformation]
    if (transformation == "best-normalize-manual") transf_name_i <- paste0(method, transf_name_i)

    message("Applying ", transf_name_i, " to ", dataset, " dataset.")
  }

  res <- switch(
    transformation,
    "vsn" = transform_vsn(
      mat,
      return_matrix_only = return_matrix_only
    ),
    "vst-deseq2" = transform_vst(
      mat,
      return_matrix_only = return_matrix_only
    ),
    "logx" = transform_logx(
      mat,
      return_matrix_only = return_matrix_only,
      log_base = log_base,
      pre_log_function = pre_log_function
    ),
    "best-normalize-auto" = transform_bestNormalise_auto(
      mat,
      return_matrix_only = return_matrix_only,
      ...
    ),
    "best-normalize-manual" = transform_bestNormalise_manual(
      mat,
      method,
      return_matrix_only = return_matrix_only,
      ...
    )
  )

  ## Option to return a MultiDataSet object in which the original dataset is replaced by the transformed dataset
  if (return_multidataset) {
    return(replace_dataset(mo_data, dataset, res))
  }

  ## to keep track of which dataset was analysed
  attr(res, "dataset_name") <- dataset
  ## to keep track of which transformation was applied
  attr(res, "transformation") <- transformation

  return(res)
}

#' Get MultiDataSet with transformed data
#'
#' Replace the original datasets with transformed datasets in a MultiDataSet
#' object from the results of transformations applied to the datasets.
#'
#' @param mo_data A \code{\link[MultiDataSet]{MultiDataSet-class}} object.
#' @param transformation_result A list in which each element is the result of
#' a transformation applied to a different dataset, computed with the
#' \code{\link{transform_dataset}} function.
#' @return A \code{\link[MultiDataSet]{MultiDataSet-class}} object, for which the
#' assay of each dataset is the imputed dataset.
#' @export
get_transformed_data <- function(mo_data, transformation_result) {
  check_is_multidataset(mo_data)

  ## Make sure that the names of the transformation_result list are the name of the datasets
  names(transformation_result) <- sapply(transformation_result, attr, "dataset_name")
  .check_names(
    names(transformation_result),
    names(mo_data),
    "The following datasets are not in 'mo_data': '_W_'."
  )

  res <- mo_data

  for (i in names(transformation_result)) {
    ## Getting the matrix with the transformed data
    if (is.matrix(transformation_result[[i]])) {
      new_mat <- transformation_result[[i]]
    } else {
      if (!("transformed_data" %in% names(transformation_result[[i]]))) {
        stop("Error with 'transformation_result': should be a list ",
             "of matrices or a list of named lists, each with a",
             " 'transformed_data' element.")
      }
      new_mat <- transformation_result[[i]][["transformed_data"]]
    }

    res <- replace_dataset(res, i, new_mat)
  }

  return(res)
}

#' Target factory for datasets transformation
#'
#' Create a list of targets to apply some transformation methods to one or more
#' datasets in a `MultiDataSet` object.
#'
#' Currently implemented transformations and recommendations based on dataset
#' type:
#' * `vsn`: Variance Stabilising normalisation, implemented in the
#' [vsn::justvsn()] function from the `vsn` package. This method was originally
#' developed for microarray intensities. This transformation is recommended for
#' microarray, metabolome, chemical or other intensity-based datasets. In
#' practice, applies the [transform_vsn()] function.
#' * `vst-deseq2`: Variance Stabilising Transformation, implemented in the
#' [DESeq2::varianceStabilizingTransformation()] function from the `DESeq2`
#' package. This method is applicable to count data only. This transformation is
#' recommended for RNAseq or similar count-based datasets. In practice, applies
#' the [transform_vst()] function.
#' * `logx`: log-transformation (default to log2, but base can be specified).
#' In practice, applies the [transform_logx()] function.
#' * `best-normalize-auto`: most appropriate normalisation method automatically
#' selected from a number of options, implemented in the
#' [bestNormalize::bestNormalize()] function from the `bestNormalize` package.
#' This transformation is recommended for phenotypes that are each measured on
#' different scales (since the transformation method selected will potentially
#' be different across the features), preferably with a reasonable number of
#' features (less than 100) to avoid large computation times. In practice,
#' applies the [transform_bestNormalise_auto()] function.
#' * `best-normalize-manual`: performs the same transformation (specified
#' through the `method` argument) to each feature of a dataset. This
#' transformation is recommended for phenotypes data in which the different
#' phenotypes are measured on the same scale. The different normalisation
#' methods are:
#'   * `"arcsinh_x"`: data is transformed as `log(x + sqrt(x^2 + 1))`;
#'   * `"boxcox"`: Box Cox transformation;
#'   * `"center_scale"`: data is centered and scaled;
#'   * `"exp_x"`: data is transformed as `exp(x)`;
#'   * `"log_x"`: data is transformed as `log_b(x+a)` (`a` and `b` either
#'                selected automatically per variable or passed as arguments);
#'   * `"orderNorm"`: Ordered Quantile technique;
#'   * `"sqrt_x"`: data transformed as `sqrt(x + a)` (`a` selected automatically
#'                 per variable or passed as argument),
#'   * `"yeojohnson"`: Yeo-Johnson transformation.
#'
#' @param mo_data_target Symbol, the name of the target containing the
#'   `MultiDataSet` object.
#' @param transformations Named character vector, name of each element is the
#'   name of a dataset to transform, corresponding element gives the type of
#'   transformation to apply to the dataset (e.g. `c(rnaseq = 'vst-deseq2',
#'   phenotypes = 'best-normalize-auto')`). See Details for a list of available
#'   transformations. If `'best-normalize-auto'` is selected, need to provide
#'   the `methods` argument as well.
#' @param return_matrix_only Logical, should only the transformed matrix be
#'   returned for each transformation? If `TRUE`, only transformed matrices will
#'   be stored. If `FALSE`, instead for each transformation, a list with the
#'   transformed data and potentially other information relevant to the
#'   transformation will be saved. Default value is `FALSE`.
#' @param target_name_prefix Character, a prefix to add to the name of the
#'   targets created by this target factory. Default value is `""`.
#' @param transformed_data_name Character, the name of the target containing the
#'   `MultiDataSet` with transformed data to be created. If `NULL`, will be
#'   selected automatically. Default value is `NULL`.
#' @param methods Character or named character list, gives for each dataset for
#'   which the `'best-normalize-manual'` transformation is selected the
#'   normalisation method that should be applied. See possible values in
#'   Details. If one value, will be used for all concerned datasets. Otherwise,
#'   can specify a different method for each concerned dataset by passing a
#'   named list.
#' @param log_bases Numeric or named numeric list, gives for each dataset for
#'   which the `'logx'` transformation is selected the log base to use. If one
#'   value, will be used for all concerned datasets. Otherwise, can specify a
#'   different log-base for each concerned dataset by passing a named list.
#' @param pre_log_functions Function or named list of functions, gives for each
#'   dataset for which the `'logx`` transformation is selected the function that
#'   will be applied to the matrix before the log transformation (e.g. to apply
#'   an offset to the values to avoid issues with zeros). Default value is the
#'   [zero_to_half_min()] function. If one value, will be used for all concerned
#'   datasets. Otherwise, can specify a different log-base for each concerned
#'   dataset by passing a named list.
#' @param ... Further arguments passed to the \code{\link{transform_dataset}}
#'   function or the `method` function from the `bestNormalize` package. Only
#'   relevant for `'best-normalize-XX'` transformations.
#' @returns A list of target objects. With `target_name_prefix = ""` and
#'   `transformed_data_name = NULL`, the following targets are created:
#' * `transformations_spec`: generates a grouped tibble where each row
#' corresponds to one dataset to be tranformed, with the columns specifying each
#' dataset name and the transformation to apply.
#' * `transformations_runs_list`: a dynamic branching target that runs the
#' [transform_dataset()] function on each dataset. Returns a list.
#' * `transformed_set`: a target that returns the `MultiDataSet` object with the
#' original data replaced by the transformed data.
#' @examples
#' \dontrun{
#' ## in the _targets.R
#' library(moiraine)
#'
#' list(
#'   ## add code here to load the different datasets
#'
#'   ## the following target creates a MultiDataSet object from previously
#'   ## created omics sets (geno_set, trans_set, etc)
#'   tar_target(
#'     mo_set,
#'     create_multiomics_set(geno_set, trans_set, metabo_set, pheno_set)
#'   ),
#'
#'   ## Example 1
#'   transformation_datasets_factory(mo_set,
#'     c(
#'       rnaseq = "vst-deseq2",
#'       metabolome = "vsn",
#'       phenotypes = "best-normalize-auto"
#'     ),
#'     return_matrix_only = FALSE,
#'     transformed_data_name = "mo_set_transformed"
#'   ),
#'
#'   ## Example 2 - with a log2 transformation for both datasets
#'   transformation_datasets_factory(
#'     mo_set_complete,
#'     c(
#'       "rnaseq" = "logx",
#'       "metabolome" = "logx"
#'     ),
#'     log_bases = 2,
#'     pre_log_functions = zero_to_half_min
#'   ),
#'
#'   ## Example 3 - with different log bases for each dataset and a different
#'   ## preprocessing function to be run before applying the log
#'   transformation_datasets_factory(
#'     mo_set_complete,
#'     c(
#'       "rnaseq" = "logx",
#'       "metabolome" = "logx"
#'     ),
#'     log_bases = list(rnaseq = 10, metabolome = 2),
#'     pre_log_functions = list(
#'       rnaseq = \(x) x + 0.5,
#'       metabolome = zero_to_half_min
#'      )
#'   )
#' )
#' }
#' @export
transformation_datasets_factory <- function(mo_data_target,
                                            transformations,
                                            return_matrix_only = FALSE,
                                            target_name_prefix = "",
                                            transformed_data_name = NULL,
                                            log_bases = 2,
                                            pre_log_functions = zero_to_half_min,
                                            methods,
                                            ...) {
  if (is.null(names(transformations))) {
    stop("'transformations' vector should be named.")
  }

  if (missing(methods)) methods <- NULL

  if (any(transformations == "best-normalize-manual")) {
    methods <- .make_var_list(
      methods,
      names(transformations)[transformations == "best-normalize-manual"]
    )
  }

  if (any(transformations == "logx")) {
    log_bases <- .make_var_list(
      log_bases,
      names(transformations)[transformations == "logx"]
    )

    pre_log_functions <- .make_var_list(
      pre_log_functions,
      names(transformations)[transformations == "logx"]
    )
  } else {
    log_bases <- pre_log_functions <- NULL
  }

  ## Target names
  transf_spec_name <- paste0(target_name_prefix, "transformations_spec")
  transf_run_name <- paste0(target_name_prefix, "transformations_runs_list")
  if (is.null(transformed_data_name)) {
    transformed_data_name <- paste0(target_name_prefix, "transformed_set")
  }

  ## Target symbols
  trans_spec_target <- as.symbol(transf_spec_name)
  transf_run_target <- as.symbol(transf_run_name)

  dsn_vals <- names(transformations)
  transf_vals <- unname(transformations)
  meth_vals <- purrr::map(dsn_vals, \(x) methods[[x]])
  log_b_vals <- purrr::map(dsn_vals, \(x) log_bases[[x]])
  prelog_f_vals <- purrr::map(dsn_vals, \(x) pre_log_functions[[x]])

  list(
    targets::tar_target_raw(
      transf_spec_name,
      substitute(
        tibble::tibble(
          dsn = dsn_vals,
          transf = transf_vals,
          meth = meth_vals,
          log_b = log_b_vals,
          prelog_f = prelog_f_vals
        ) |>
          dplyr::group_by(dsn) |>
          targets::tar_group()),
      iteration = "group"
    ),

    ## Apply the transformation to each dataset
    targets::tar_target_raw(
      transf_run_name,
      substitute(
        transform_dataset(
          mo_data_target,
          dataset = trans_spec_target$dsn,
          transformation = trans_spec_target$transf,
          return_matrix_only = return_matrix_only,
          method = trans_spec_target$meth[[1]],
          log_base = trans_spec_target$log_b[[1]],
          pre_log_function = trans_spec_target$prelog_f[[1]],
          ...
        )
      ),
      pattern = substitute(map(trans_spec_target)),
      iteration = "list"
    ),

    ## Get the MultiDataSet object with the transformed data
    targets::tar_target_raw(
      transformed_data_name,
      substitute(get_transformed_data(mo_data_target, transf_run_target))
    )
  )
}

.bestNormalize_get_transfo_name <- function(bn_list) {
  if ("bestNormalize" %in% class(bn_list)) {
    txt <- utils::capture.output(bn_list$chosen_transform)[1]
    res <- stringr::str_extract(txt, "^.+(?= Transformation)")
  } else {
    txt <- utils::capture.output(bn_list)[1]
    res <- stringr::str_extract(txt, "^.+(?= Transformation)")
  }

  if (is.na(res)) {
    res <- class(bn_list)[1]

    if (res == "log_x") {
      res <- paste0("log", bn_list$b)
      if (bn_list$a != 0) res <- paste0(res, "(x + ", round(bn_list$a, 2), ")")
    }

    if (!is.null(bn_list$standardize)) {
      if (bn_list$standardize) res <- paste0("Standardised ", res)
    }
  } else {
    if (stringr::str_detect(res, "Log_b")) {
      res <- stringr::str_replace(res, "(?<=Log_)b", paste0(bn_list$b))
      res <- stringr::str_replace(res, " \\+ a", dplyr::if_else(bn_list$a == 0, "", paste0(" + ", bn_list$a)))
    }
  }

  return(res)
}

#' Get table with transformation applied to each dataset
#'
#' From the results of transformations on datasets, generates a table giving for
#' each dataset the transformation that was applied to it.
#'
#' @param transformation_result A list in which each element is the result of a
#'   transformation applied to a different dataset, computed with the
#'   [transform_dataset] function.
#' @param best_normalize_details Logical, should information about the
#'   transformations selected by bestNormalize for each feature be displayed?
#'   Default value is `FALSE`.
#' @returns A tibble with columns `'Dataset'` and `'Transformation'`. If
#'   `best_normalize_details = TRUE`, an additional column `'Details'` lists the
#'   chsoen transformation applied to each feature of the corresponding dataset
#'   for a bestNormalize transformation.
#' @export
get_table_transformations <- function(transformation_result,
                                      best_normalize_details = FALSE) {
    ## for devtools::check
    Features <- Feature <- NULL

    transf_name <- c(
      "vsn" = "Variance Stabilising Normalisation (vsn)",
      "vst-deseq2" = "Variance Stabilising Transformation (DESeq2)",
      "logx" = "Log-X transformation",
      "best-normalize-auto" = "automatic normalisation selection (bestNormalize)",
      "best-normalize-manual" = " transformation (bestNormalize)"
    )

    ## for devtools::check()
    transf <- Transformation <- Chosen_transformation <- NULL

    names(transformation_result) <- sapply(
      transformation_result,
      attr,
      "dataset_name"
    )

    res <- tibble::tibble(
      Dataset = sapply(transformation_result, attr, "dataset_name"),
      transf = sapply(transformation_result, attr, "transformation")
    ) |>
      dplyr::mutate(Transformation = transf_name[transf])

    for (i in which(res$transf == "logx")) {
      res$Transformation[i] <- stringr::str_replace(
        res$Transformation[i],
        "X",
        transformation_result$info_transformation$log_base
      )
    }

    for (i in which(res$transf == "best-normalize-manual")) {
      res$Transformation[i] <- paste0(
          .bestNormalize_get_transfo_name(
            transformation_result[[res$Dataset[i]]][["info_transformation"]][[1]]
          ),
          res$Transformation[i]
        )
    }

    if (best_normalize_details & any(res$transf == "best-normalize-auto")) {
      res$Details <- sapply(transformation_result, function(x) {
        if (attr(x, "transformation") == "best-normalize-auto") {
          info_list <- x$info_transformation

          df <- tibble::tibble(
            Feature = names(info_list),
            Chosen_transformation = sapply(info_list, .bestNormalize_get_transfo_name)
          ) |>
            dplyr::group_by(Chosen_transformation) |>
            dplyr::summarise(Features = paste0(Feature, collapse = ", "))

          return(paste0(
            "- ",
            df$Chosen_transformation,
            ": ",
            df$Features,
            collapse = "\n"
          ))
        }
        return("")
      })
    }

    res <- res |>
      dplyr::select(-transf)

    return(res)
  }

.get_transformed_matrix <- function(res, return_matrix_only) {
  if (return_matrix_only) {
    mat_transformed <- res
  } else {
    mat_transformed <- res[["transformed_data"]]
  }

  return(mat_transformed)
}
