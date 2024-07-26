#' Wrapper for OmicsPLS::crossval_o2m_adjR2 function
#'
#' Wrapper function around the \code{\link[OmicsPLS]{crossval_o2m_adjR2}} function.
#' The main purpose of this wrapper is to add to the result the names of the datasets
#' to facilitate plotting.
#'
#' @param omicspls_input A named list of length 2, produced by \code{\link{get_input_omicspls}}.
#' @param a Vector of positive integers, number of joint components to test.
#' @param ax Vector of non-negative integers, number of specific components to test
#' for the first dataset.
#' @param ay Vector of non-negative integers, number of specific components to test
#' for the second dataset.
#' @param nr_folds Positive integer, number of folds to use for the cross-validation.
#' Default value is `10`.
#' @param ... Further arguments passed to the \code{\link[OmicsPLS]{crossval_o2m_adjR2}} function.
#' @returns A data-frame with four columns: `MSE`, `n`, `nx` and `ny.` Each row
#' corresponds to an element in `a`.
#' @export
so2pls_crossval_o2m_adjR2 <- function(omicspls_input,
                                      a = 1:5,
                                      ax = seq(0, 10, by = 2),
                                      ay = seq(0, 10, by = 2),
                                      nr_folds = 10,
                                      ...) {
  if (!is.list(omicspls_input) | length(omicspls_input) != 2) stop("'omicspls_input' should be a list of length 2 produced by get_input_omicspls().")

  res <- OmicsPLS::crossval_o2m_adjR2(omicspls_input[[1]],
                                      omicspls_input[[2]],
                                      a = a,
                                      ax = ax,
                                      ay = ay,
                                      nr_folds = nr_folds,
                                      ...
  )

  attr(res, "datasets_name") <- names(omicspls_input)

  return(res)
}

#' Plot adjusted cross-validation results for sO2PLS
#'
#' Plots the results of an adjusted cross-validation for an sO2PLS
#' run (from the `OmicsPLS` package)
#'
#' @param cv_res Data-frame, output from the \code{\link[OmicsPLS]{crossval_o2m_adjR2}} function.
#' @param with_labels Boolean, whether the optimal values for `nx` and `ny` for each value
#' of `n` should be displayed. Default value is `TRUE`.
#' @return A ggplot
#' @export
so2pls_plot_cv_adj <- function(cv_res, with_labels = TRUE) {
  if (!is.data.frame(cv_res)) stop("Input should be a data-frame with columns 'MSE', 'n', 'nx' and 'ny'. Make sure input is the output from the OmicsPLS::crossval_o2m_adjR2() function.")
  if (length(setdiff(colnames(cv_res), c("MSE", "n", "nx", "ny")))) stop("Input should be a data-frame with columns 'MSE', 'n', 'nx' and 'ny'. Make sure input is the output from the OmicsPLS::crossval_o2m_adjR2() function.")

  ## For devtools::check()
  n <- MSE <- label <- nx <- ny <- NULL

  p <- tibble::as_tibble(cv_res) |>
    dplyr::mutate(label = paste0(
      "nx: ", nx,
      ", ny: ", ny
    )) |>
    ggplot2::ggplot(aes(x = n, y = MSE)) +
    ggplot2::geom_line() +
    ggplot2::geom_point()

  if (with_labels) {
    p <- p + ggrepel::geom_label_repel(aes(label = label))
  }


  p <- p +
    ggplot2::scale_x_continuous(breaks = min(cv_res$n):max(cv_res$n)) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Number of joint components")

  return(p)
}

#' Print adjusted cross-validation results for sO2PLS
#'
#' Prints the results of an adjusted cross-validation for an sO2PLS
#' run (from the `OmicsPLS` package)
#'
#' @param cv_res Data-frame, output from the \code{\link[OmicsPLS]{crossval_o2m_adjR2}} function.
#' @return A tibble.
#' @export
so2pls_print_cv_adj <- function(cv_res) {
  if (!is.data.frame(cv_res)) stop("Input should be a data-frame with columns 'MSE', 'n', 'nx' and 'ny'. Make sure input is the output from the OmicsPLS::crossval_o2m_adjR2() function.")
  if (length(setdiff(colnames(cv_res), c("MSE", "n", "nx", "ny")))) stop("Input should be a data-frame with columns 'MSE', 'n', 'nx' and 'ny'. Make sure input is the output from the OmicsPLS::crossval_o2m_adjR2() function.")

  ## For devtools::check()
  n <- nx <- ny <- MSE <- NULL

  tibble::as_tibble(cv_res) |>
    dplyr::arrange(MSE) |>
    dplyr::rename(
      nb_joint_components = n,
      nb_Xspecific_components_optimal = nx,
      nb_Yspecific_components_optimal = ny
    )
}

#' Extract optimal number of components from adjusted cross-validation results for sO2PLS
#'
#' Extracts the optimal number of components (joint and dataset-specific) estimated
#' via adjusted cross-validation results for sO2PLS.
#'
#' @param cv_res Data-frame, output from the \code{\link[OmicsPLS]{crossval_o2m_adjR2}} function.
#' @return A vector with three integer values:
#' \itemize{
#' \item `n`: optimal number of joint components
#' \item `nx`: optimal number of specific components for dataset X (first dataset)
#' \item `ny`: optimal number of specific components for dataset Y (second dataset)
#' }
#' @export
so2pls_get_optim_ncomp_adj <- function(cv_res) {
  if (!is.data.frame(cv_res)) stop("Input should be a data-frame with columns 'MSE', 'n', 'nx' and 'ny'. Make sure input is the output from the OmicsPLS::crossval_o2m_adjR2() function.")
  if (length(setdiff(colnames(cv_res), c("MSE", "n", "nx", "ny")))) stop("Input should be a data-frame with columns 'MSE', 'n', 'nx' and 'ny'. Make sure input is the output from the OmicsPLS::crossval_o2m_adjR2() function.")

  unlist(cv_res[which(cv_res$MSE == min(cv_res$MSE)), -1])
}


#' Wrapper for OmicsPLS::crossval_o2m function
#'
#' Wrapper function around the \code{\link[OmicsPLS]{crossval_o2m}} function.
#' The main purpose of this wrapper is to add to the result the names of the datasets
#' to facilitate plotting. If the result of a previous call to \code{\link[OmicsPLS]{crossval_o2m_adjR2}}
#' or \code{\link{so2pls_crossval_o2m_adjR2}} is provided, will be used
#' to set values to test for `a`, `ax` and `ay`.
#'
#' If the result of a previous call to \code{\link[OmicsPLS]{crossval_o2m_adjR2}}
#' or \code{\link{so2pls_crossval_o2m_adjR2}} is provided through the
#' `cv_adj_res` parameter, the optimal values for `n`, `nx` and `ny` are
#' extracted from it, and the values of `a`, `ax` and `ay` are set as follows:
#'
#' * `a` = max(`n` - 1, 1):(`n` + 1)
#'
#' * `ax` = max(`nx` - 1, 0):(`nx` + 1)
#'
#' * `ay` = max(`ny` - 1, 0):(`ny` + 1)
#'
#'
#' @param omicspls_input A named list of length 2, produced by \code{\link{get_input_omicspls}}.
#' @param cv_adj_res Data-frame returned by \code{\link[OmicsPLS]{crossval_o2m_adjR2}}
#' or \code{\link{so2pls_crossval_o2m_adjR2}}. Default value is `NULL`.
#' @param a Vector of positive integers, number of joint components to test.
#' Ignored if `cv_adj_res` is not `NULL`.
#' @param ax Vector of non-negative integers, number of specific components to test
#' for the first dataset. Ignored if `cv_adj_res` is not `NULL`.
#' @param ay Vector of non-negative integers, number of specific components to test
#' for the second dataset. Ignored if `cv_adj_res` is not `NULL`.
#' @param nr_folds Positive integer, number of folds to use for the cross-validation.
#' Default value is `10`.
#' @param ... Further arguments passed to the \code{\link[OmicsPLS]{crossval_o2m_adjR2}} function.
#' @returns A list of class "cvo2m" with the original and sorted Prediction
#' errors and the number of folds used.
#' @export
so2pls_crossval_o2m <- function(omicspls_input,
                                cv_adj_res = NULL,
                                a = 1:5,
                                ax = seq(0, 10, by = 2),
                                ay = seq(0, 10, by = 2),
                                nr_folds = 10,
                                ...) {
  if (!is.list(omicspls_input) | length(omicspls_input) != 2) stop("'omicspls_input' should be a list of length 2 produced by get_input_omicspls().")

  if (!is.null(cv_adj_res)) {
    ## in case the results of a previous adjusted cross-validation are passed on,
    ## use the results to specify the values to test for a, ax, and ay
    so2pls_cv_adj_res <- so2pls_get_optim_ncomp_adj(cv_adj_res)

    a <- max(so2pls_cv_adj_res["n"] - 1, 1):(so2pls_cv_adj_res["n"] + 1)
    ax <- max(so2pls_cv_adj_res["nx"] - 1, 0):(so2pls_cv_adj_res["nx"] + 1)
    ay <- max(so2pls_cv_adj_res["ny"] - 1, 0):(so2pls_cv_adj_res["ny"] + 1)
  }

  res <- OmicsPLS::crossval_o2m(omicspls_input[[1]],
                                omicspls_input[[2]],
                                a = a,
                                ax = ax,
                                ay = ay,
                                nr_folds = nr_folds,
                                ...
  )

  attr(res, "datasets_name") <- names(omicspls_input)

  return(res)
}

#' Plots cross-validation results for sO2PLS
#'
#' Plots the results of a cross-validation for an sO2PLS
#' run (from the `OmicsPLS` package)
#'
#' @param cv_res A `cvo2m` object, output from the \code{\link[OmicsPLS]{crossval_o2m}} function.
#' @param nb_col Integer, the number of columns to use for the faceted plot. Default value is
#' `NULL` (the number of columns will be chosen automatically).
#' @return A ggplot.
#' @export
so2pls_plot_cv <- function(cv_res, nb_col = NULL) {
  if (!is(cv_res, "cvo2m")) {
    stop(
      "Expecting a cvo2m object. ",
      "Make sure input is the output from the so2pls_crossval_o2m() ",
      "or OmicsPLS::crossval_o2m() function."
    )
  }

  names_datasets <- attr(cv_res, "datasets_name")
  if (is.null(names_datasets)) names_datasets <- c("X", "Y")

  ## For devtools::check()
  n <- nx <- ny <- MSE <- is_optimal <- NULL

  toplot <- tibble::as_tibble(cv_res$Original) |>
    dplyr::mutate(nx = rownames(cv_res$Original)) |>
    tidyr::pivot_longer(
      cols = -nx,
      names_to = "names",
      values_to = "MSE"
    ) |>
    tidyr::separate(names, c("ny", "n"), sep = "\\.") |>
    dplyr::mutate(
      dplyr::across(
        dplyr::starts_with("n"),
        ~ as.numeric(stringr::str_extract(.x, "\\d+"))
      ),
      n = paste0("Number of joint components: ", n),
      is_optimal = MSE == min(MSE)
    ) |>
    dplyr::arrange(is_optimal)

  toplot |>
    ggplot2::ggplot(aes(x = nx, y = ny, fill = MSE, colour = is_optimal)) +
    ggplot2::geom_tile(linewidth = 1) +
    ggplot2::facet_wrap(~n, ncol = nb_col) +
    ggplot2::scale_fill_viridis_c(option = "plasma") +
    ggplot2::scale_colour_manual(
      values = c("TRUE" = "red", "FALSE" = "white"),
      guide = "none"
    ) +
    ggplot2::scale_x_continuous(
      breaks = min(toplot$nx):max(toplot$nx),
      expand = ggplot2::expansion()
    ) +
    ggplot2::scale_y_continuous(
      breaks = min(toplot$ny):max(toplot$ny),
      expand = ggplot2::expansion()
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(
      x = paste0("Number of ", names_datasets[1], "-specific components"),
      y = paste0("Number of\n", names_datasets[2], "-specific components"),
      fill = "MSE"
    )
}


#' Extract optimal number of components from cross-validation results for sO2PLS
#'
#' Extracts the optimal number of components (joint and dataset-specific) estimated
#' via cross-validation results for sO2PLS.
#'
#' @param cv_res A `cvo2m` object, output from the \code{\link[OmicsPLS]{crossval_o2m}} function.
#' @return A vector with three integer values:
#' \itemize{
#' \item `n`: optimal number of joint components
#' \item `nx`: optimal number of specific components for dataset X (first dataset)
#' \item `ny`: optimal number of specific components for dataset Y (second dataset)
#' }
#' @export
so2pls_get_optim_ncomp <- function(cv_res) {
  if (!is(cv_res, "cvo2m")) stop("Expecting a cvo2m object. Make sure input is the output from the so2pls_crossval_o2m() or OmicsPLS::crossval_o2m() function.")

  ## For devtools::check()
  n <- nx <- ny <- MSE <- is_optimal <- NULL

  df <- tibble::as_tibble(cv_res$Original) |>
    dplyr::mutate(nx = rownames(cv_res$Original)) |>
    tidyr::pivot_longer(
      cols = -nx,
      names_to = "names",
      values_to = "MSE"
    ) |>
    tidyr::separate(names, c("ny", "n"), sep = "\\.") |>
    dplyr::mutate(
      dplyr::across(
        dplyr::starts_with("n"),
        ~ as.integer(stringr::str_extract(.x, "\\d+"))
      )
    ) |>
    dplyr::filter(MSE == min(MSE)) |>
    dplyr::select(n, nx, ny)

  unlist(df)
}

#' Perform cross-validation to find the optimal number of features/groups to
#' keep for each joint component for sO2PLS
#'
#' Computes the optimal number of features/groups to keep for each joint
#' component for an sO2PLS run. Directly copied from the
#' [OmicsPLS::crossval_sparsity()] function, but improved the output for
#' plotting purposes.
#'
#' @param omicspls_input A named list of length 2, produced by
#'   [get_input_omicspls()].
#' @param n Integer, number of joint PLS components. Must be positive.
#' @param nx Integer, number of orthogonal components in `X`. Negative values
#'   are interpreted as 0.
#' @param ny Integer, number of orthogonal components in `Y`. Negative values
#'   are interpreted as 0.
#' @param nr_folds integer, number of folds for the cross-validation. Default
#'   value is 10.
#' @param keepx_seq Numeric vector, how many features/groups to keep for
#'   cross-validation in each of the joint components of `X`. Sparsity of each
#'   joint component will be selected sequentially.
#' @param keepy_seq Numeric vector, how many features/groups to keep for
#'   cross-validation in each of the joint components of `Y`. Sparsity of each
#'   joint component will be selected sequentially.
#' @param groupx Character vector, group name of each `X`-feature. Its length
#'   must be equal to the number of features in `X`. The order of the group
#'   names must corresponds to the order of the features. If `NULL`, no groups
#'   are considered. Default value is `NULL`.
#' @param groupy Character vector, group name of each `Y`-feature. Its length
#'   must be equal to the number of features in `Y`. The order of the group
#'   names must corresponds to the order of the features. If `NULL`, no groups
#'   are considered. Default value is `NULL`.
#' @param tol Numeric, threshold for which the NIPALS method is deemed
#'   converged. Must be positive. Default value is `1e-10`.
#' @param max_iterations Integer, maximum number of iterations for the NIPALS
#'   method.
#' @param seed Integer, seed to use. Default is `NULL`, i.e. no seed is set
#'   inside the function.
#' @returns A list with the following elements:
#' * `Best`: a vector giving for each join component the number of features
#'   to keep from `X` and `Y` that yield the highest covariance between the
#'   joint components of `X` and `Y` (elements `x1`, `y1`, `x2`, `y2`, etc),
#'   and the number of features to keep from `X` and `Y` yielding the highest
#'   covariance under the 1 standard error rule (elements `x_1sd1`, `y_1sd1`,
#'   `x_1sd2`, `y_1sd2`, etc).
#' * `Covs`: a list, with as many elements as number of joint components (`n`).
#'   Each element is a matrix giving the average covariance between the joint
#'   components of `X` and `Y` obtained across the folds, for each tested values
#'   of `keepx` (columns) and of `keepy` (rows).
#' * `SEcov`: a list, with as many elements as number of joint components (`n`).
#'   Each element is a matrix giving the standard error of the covariance
#'   between the joint components of `X` and `Y` obtained across the folds, for
#'   each tested values of `keepx` (columns) and of `keepy` (rows).
#' @export
so2pls_crossval_sparsity <- function(omicspls_input,
                                     n,
                                     nx,
                                     ny,
                                     nr_folds = 10,
                                     keepx_seq = NULL,
                                     keepy_seq = NULL,
                                     groupx = NULL,
                                     groupy = NULL,
                                     tol = 1e-10,
                                     max_iterations = 100,
                                     seed = NULL) {
  if (!is.list(omicspls_input) | length(omicspls_input) != 2) stop("'omicspls_input' should be a list of length 2 produced by get_input_omicspls().")

  X <- omicspls_input[[1]]
  Y <- omicspls_input[[2]]

  if (!is.null(seed)) set.seed(seed)

  if (is.null(groupx) & is.null(groupy)) {
    method <- "SO2PLS"
    message("Group information not provided, CV for number of variables to keep\n")
    OmicsPLS::cv_lambda_checker(X, Y, keepx_seq, keepy_seq)
    if (is.null(keepx_seq)) {
      keepx_seq <- ncol(X)
    }
    if (is.null(keepy_seq)) {
      keepy_seq <- ncol(Y)
    }
  } else {
    method <- "GO2PLS"
    message("Group information provided, CV for number of groups to keep\n")
    if (is.null(groupx)) {
      groupx <- colnames(X)
    }
    if (is.null(groupy)) {
      groupy <- colnames(Y)
    }
    OmicsPLS::cv_lambda_checker_group(groupx, groupy, keepx_seq, keepy_seq)
    if (is.null(keepx_seq)) {
      keepx_seq <- ncol(X)
    }
    if (is.null(keepy_seq)) {
      keepy_seq <- ncol(Y)
    }
  }
  stopifnot(all(n == round(n)), nr_folds == round(nr_folds))
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  OmicsPLS::input_checker(X, Y)
  if (nx + ny > 0) {
    n2 <- n + max(nx, ny)
    W_C <- suppressMessages(OmicsPLS::pow_o2m(X, Y, n2, tol, max_iterations))
    W <- W_C$W
    C <- W_C$C
    Tt <- W_C$Tt
    U <- W_C$U
    if (nx > 0) {
      E_XY <- X - Tt %*% t(W)
      udv <- svd(t(E_XY) %*% Tt, nu = nx, nv = 0)
      W_Yosc <- udv$u
      T_Yosc <- X %*% W_Yosc
      P_Yosc <- t(solve(t(T_Yosc) %*% T_Yosc) %*% t(T_Yosc) %*%
                    X)
      X <- X - T_Yosc %*% t(P_Yosc)
    }
    if (ny > 0) {
      F_XY <- Y - U %*% t(C)
      udv <- svd(t(F_XY) %*% U, nu = ny, nv = 0)
      C_Xosc <- udv$u
      U_Xosc <- Y %*% C_Xosc
      P_Xosc <- t(solve(t(U_Xosc) %*% U_Xosc) %*% t(U_Xosc) %*%
                    Y)
      Y <- Y - U_Xosc %*% t(P_Xosc)
    }
  }
  N <- length(X[, 1])
  if (N != length(Y[, 1])) {
    stop("N not the same")
  }
  mean_covTU <- srr_covTU <- matrix(NA,
                                    nrow = length(keepy_seq),
                                    ncol = length(keepx_seq)
  )
  rownames(mean_covTU) <- rownames(srr_covTU) <- keepy_seq
  colnames(mean_covTU) <- colnames(srr_covTU) <- keepx_seq
  covTU <- NA * 1:nr_folds
  keepxy_x <- keepxy_y <- x_max <- y_max <- vector()
  mean_covTU_list <- list()
  srr_covTU_list <- list()
  if (method == "SO2PLS") {
    for (comp in 1:n) {
      kx <- 0
      blocks <- cut(seq(1:N), breaks = nr_folds, labels = F)
      folds <- sample(N)
      for (lx in keepx_seq) {
        kx <- kx + 1
        ky <- 0
        for (ly in keepy_seq) {
          ky <- ky + 1
          for (i in 1:nr_folds) {
            ii <- which(blocks == i)
            X_tr <- X[-folds[ii], ]
            X_tst <- X[folds[ii], ]
            Y_tr <- Y[-folds[ii], ]
            Y_tst <- Y[folds[ii], ]
            v <- X_tr[1, ] / OmicsPLS::norm_vec(X_tr[1, ])
            for (k in 1:max_iterations) {
              v_old <- v
              u <- t(Y_tr) %*% (X_tr %*% v)
              u <- OmicsPLS::thresh_n(u, ly)
              u <- u / OmicsPLS::norm_vec(u)
              v <- t(X_tr) %*% (Y_tr %*% u)
              v <- OmicsPLS::thresh_n(v, lx)
              v <- v / OmicsPLS::norm_vec(v)
              if (OmicsPLS::mse(v, v_old) < tol) {
                break
              }
            }
            t_tst <- X_tst %*% v
            u_tst <- Y_tst %*% u
            covTU[i] <- abs(drop(stats::cov(t_tst, u_tst)))
          }
          mean_covTU[ky, kx] <- mean(covTU)
          srr_covTU[ky, kx] <- stats::sd(covTU) / sqrt(nr_folds)
        }
      }
      mean_covTU_list[[comp]] <- mean_covTU
      srr_covTU_list[[comp]] <- srr_covTU
      cov_max <- max(mean_covTU)
      cov_1srr <- cov_max - srr_covTU[which.max(mean_covTU)]
      keepxy <- OmicsPLS::err_back(mean_covTU, which(mean_covTU >
                                                       cov_1srr, arr.ind = T), dim(X)[2], dim(Y)[2])
      v <- X[1, ] / OmicsPLS::norm_vec(X[1, ])
      for (k in 1:max_iterations) {
        v_old <- v
        u <- t(Y) %*% (X %*% v)
        u <- OmicsPLS::thresh_n(u, keepxy$y)
        u <- u / OmicsPLS::norm_vec(u)
        v <- t(X) %*% (Y %*% u)
        v <- OmicsPLS::thresh_n(v, keepxy$x)
        v <- v / OmicsPLS::norm_vec(v)
        if (OmicsPLS::mse(v, v_old) < tol) {
          break
        }
      }
      t_tmp <- X %*% v
      u_tmp <- Y %*% u
      p <- t(X) %*% t_tmp / drop(crossprod(t_tmp))
      q <- t(Y) %*% u_tmp / drop(crossprod(u_tmp))
      X <- X - t_tmp %*% t(p)
      Y <- Y - u_tmp %*% t(q)
      keepxy_x[comp] <- keepxy$x
      keepxy_y[comp] <- keepxy$y
      y_max[comp] <- as.numeric(rownames(mean_covTU)[which(mean_covTU ==
                                                             max(mean_covTU), arr.ind = T)[1]])
      x_max[comp] <- as.numeric(colnames(mean_covTU)[which(mean_covTU ==
                                                             max(mean_covTU), arr.ind = T)[2]])
    }
  } else {
    names_grx <- groupx |> unique()
    names_gry <- groupy |> unique()
    nr_grx <- names_grx |> length()
    nr_gry <- names_gry |> length()
    index_grx <- lapply(1:nr_grx, function(j) {
      index <- which(groupx == names_grx[j])
      size <- length(index)
      return(list(index = index, size = size))
    })
    index_gry <- lapply(1:nr_gry, function(j) {
      index <- which(groupy == names_gry[j])
      size <- length(index)
      return(list(index = index, size = size))
    })
    names(index_grx) <- names_grx
    names(index_gry) <- names_gry
    for (comp in 1:n) {
      kx <- 0
      blocks <- cut(seq(1:N), breaks = nr_folds, labels = F)
      folds <- sample(N)
      for (lx in keepx_seq) {
        kx <- kx + 1
        ky <- 0
        for (ly in keepy_seq) {
          ky <- ky + 1
          for (i in 1:nr_folds) {
            ii <- which(blocks == i)
            X_tr <- X[-folds[ii], ]
            X_tst <- X[folds[ii], ]
            Y_tr <- Y[-folds[ii], ]
            Y_tst <- Y[folds[ii], ]
            v <- X_tr[1, ] / OmicsPLS::norm_vec(X_tr[1, ])
            for (k in 1:max_iterations) {
              v_old <- v
              u <- t(Y_tr) %*% (X_tr %*% v)
              ul <- OmicsPLS::thresh_n_gr(u, ly, index_gry)
              u <- ul$w
              u <- u / OmicsPLS::norm_vec(u)
              v <- t(X_tr) %*% (Y_tr %*% u)
              vl <- OmicsPLS::thresh_n_gr(v, lx, index_grx)
              v <- vl$w
              v <- v / OmicsPLS::norm_vec(v)
              if (OmicsPLS::mse(v, v_old) < tol) {
                break
              }
            }
            t_tst <- X_tst %*% v
            u_tst <- Y_tst %*% u
            covTU[i] <- drop(stats::cov(t_tst, u_tst))
            covTU[i] <- abs(drop(stats::cov(t_tst, u_tst)))
          }
          mean_covTU[ky, kx] <- mean(covTU)
          srr_covTU[ky, kx] <- stats::sd(covTU) / sqrt(nr_folds)
        }
      }
      cov_max <- max(mean_covTU)
      cov_1srr <- cov_max - srr_covTU[which.max(mean_covTU)]
      keepxy <- OmicsPLS::err_back(mean_covTU, which(mean_covTU >
                                                       cov_1srr, arr.ind = T), nr_grx, nr_gry)
      v <- X[1, ] / OmicsPLS::norm_vec(X[1, ])
      for (k in 1:max_iterations) {
        v_old <- v
        u <- t(Y) %*% (X %*% v)
        ul <- OmicsPLS::thresh_n_gr(u, keepxy$y, index_gry)
        u <- ul$w
        u <- u / OmicsPLS::norm_vec(u)
        v <- t(X) %*% (Y %*% u)
        vl <- OmicsPLS::thresh_n_gr(v, keepxy$x, index_grx)
        v <- vl$w
        v <- v / OmicsPLS::norm_vec(v)
        if (OmicsPLS::mse(v, v_old) < tol) {
          break
        }
      }
      t_tmp <- X %*% v
      u_tmp <- Y %*% u
      p <- t(X) %*% t_tmp / drop(crossprod(t_tmp))
      q <- t(Y) %*% u_tmp / drop(crossprod(u_tmp))
      X <- X - t_tmp %*% t(p)
      Y <- Y - u_tmp %*% t(q)
      keepxy_x[comp] <- keepxy$x
      keepxy_y[comp] <- keepxy$y
      y_max[comp] <- as.numeric(rownames(mean_covTU)[which(mean_covTU ==
                                                             max(mean_covTU), arr.ind = T)[1]])
      x_max[comp] <- as.numeric(colnames(mean_covTU)[which(mean_covTU ==
                                                             max(mean_covTU), arr.ind = T)[2]])
    }
  }
  bestsp <- list()
  bestsp$x_1sd <- keepxy_x
  bestsp$y_1sd <- keepxy_y
  bestsp$x <- x_max
  bestsp$y <- y_max

  res <- list(Best = unlist(bestsp), Covs = mean_covTU_list, SEcov = srr_covTU_list)
  attr(res, "datasets_name") <- names(omicspls_input)
  return(res)
}

#' Extract optimal number of features to keep from cross-validation results for
#' sO2PLS
#'
#' Extracts the optimal number of features to retain from datasets `X` and `Y`
#' for the joint components.
#'
#' The 1-SD rule means that we are retaining the smallest number of features
#' yielding an average covariance that is within 1SD of the maximum covariance
#' obtained.
#'
#' @param cv_res List, result from a call to the [so2pls_crossval_sparsity()] or
#'   [OmicsPLS::crossval_sparsity()].
#' @param use_1sd_rule Boolean, should the 1 standard deviation rule be used
#'   when selecting the optimal number of features to retain? See Details.
#' @returns A list with elements `keepx` and `keepy`, each a vector of length
#'   equal to the number of joint components, where the ith element giving the
#'   number of features to retain from dataset `X` (`keepx`) or `Y` (`keepy`)
#'   for the i-th joint component.
#' @export
so2pls_get_optim_keep <- function(cv_res, use_1sd_rule = TRUE) {
  x <- cv_res$Best
  if (use_1sd_rule) {
    res <- list(
      keepx = x[stringr::str_detect(names(x), "x_1sd")],
      keepy = x[stringr::str_detect(names(x), "y_1sd")]
    )
  } else {
    res <- list(
      keepx = x[stringr::str_detect(names(x), "x\\d*$")],
      keepy = x[stringr::str_detect(names(x), "y\\d*$")]
    )
  }

  attr(res, "datasets_name") <- attr(cv_res, "datasets_name")

  return(res)
}

#' Print sparsity cross-validation results for sO2PLS
#'
#' Prints the results of a sparsity cross-validation for an sO2PLS
#' run (from the `OmicsPLS` package)
#'
#' @param cv_res_optim Named list, output from the \code{\link{so2pls_get_optim_keep}} function.
#' @returns A tibble, giving for each dataset (`dataset` column) and joint component
#' (other columns) the optimal number of features to retain, as well as the total number
#' of features per dataset to retain.
#' @export
so2pls_print_cv_sparsity <- function(cv_res_optim) {
  if (!(is.list(cv_res_optim) & identical(names(cv_res_optim), c("keepx", "keepy")))) {
    stop("'cv_res_optim' should be a named list of length 2 as returned by so2pls_get_optim_keep() function.")
  }

  names_datasets <- attr(cv_res_optim, "datasets_name")
  if (is.null(names_datasets)) names_datasets <- c("X", "Y")

  ## for devtools::check
  dataset <- joint_component <- n_features <- NULL

  names(names_datasets) <- c("keepx", "keepy")

  purrr::map_dfr(cv_res_optim,
                 tibble::enframe,
                 name = "joint_component",
                 value = "n_features",
                 .id = "dataset"
  ) |>
    dplyr::mutate(
      dataset = names_datasets[dataset],
      joint_component = stringr::str_extract(joint_component, "\\d+$"),
      joint_component = as.numeric(joint_component)
    ) |>
    ## different formatting in CV results when only 1 joint component
    tidyr::replace_na(list(joint_component = 1)) |>
    dplyr::arrange(joint_component) |>
    dplyr::mutate(
      joint_component = paste0("Joint component ", joint_component)
    ) |>
    tidyr::pivot_wider(
      names_from = joint_component,
      values_from = n_features
    ) |>
    dplyr::mutate(Total = rowSums(dplyr::across(tidyselect::where(is.numeric))))
}

#' Plot sparsity cross-validation results for sO2PLS
#'
#' Plots the results of a sparsity cross-validation for an sO2PLS run (from
#' the \code{OmicsPLS} package).
#'
#' The produced plot has one facet for each joint component. The x-axis
#' corresponds to the number of features retained from the X dataset to
#' construct the joint component, and the y-axis to the number of features
#' retained from the Y dataset to construct the joint component. The colour
#' of each point in the ith facet represents the average covariance obtained
#' between the joint ith components of the two datasets over the cross-validation
#' folds. The size of the points' shadow correspond to the covariance standard
#' error across the cross-validation folds. For each joint component, the
#' setting yielding the maximum average covariance is highlighted in orange,
#' and the one yielding the highest average covariance under the 1-SD rule in
#' red.
#'
#' @param cv_res List, result from a call to the \code{\link{so2pls_crossval_sparsity}}.
#' @return A ggplot.
#' @export
so2pls_plot_cv_sparsity <- function(cv_res) {
  names_datasets <- attr(cv_res, "datasets_name")
  if (is.null(names_datasets)) names_datasets <- c("X", "Y")

  ncomps <- length(cv_res$Covs)

  ## For devtools::check()
  value <- name <- dataset <- type <- Comp <- NULL
  keepX <- keepY <- Covariance <- SE <- NULL

  choice_keep <- tibble::tibble(
    value = cv_res$Best,
    name = names(cv_res$Best)
  ) |>
    dplyr::mutate(
      dataset = dplyr::case_when(
        stringr::str_detect(name, "^x") ~ "keepX",
        stringr::str_detect(name, "^y") ~ "keepY"
      ),
      type = dplyr::case_when(
        stringr::str_detect(name, "_1sd") ~ "1SD",
        TRUE ~ "NoSD"
      ),
      Comp = as.numeric(stringr::str_extract(name, "\\d+$"))
    ) |>
    dplyr::select(-name) |>
    tidyr::pivot_wider(
      names_from = dataset,
      values_from = value
    ) |>
    ## If a pair (keepX, keepY) is kept both with the 1SD rule and the noSD rule, only
    ## use the colour of the 1SD rule
    dplyr::group_by(Comp, keepX, keepY) |>
    dplyr::summarise(
      col = dplyr::case_when(
        "1SD" %in% type ~ "1SD",
        "NoSD" %in% type ~ "NoSD"
      ),
      .groups = "drop"
    )

  ## Different formatting of the names in cv_res$Best when only 1 common component
  if (ncomps == 1 & all(is.na(choice_keep$Comp))) {
    choice_keep <- choice_keep |>
      dplyr::mutate(Comp = 1)
  }

  df <- lapply(1:ncomps, function(comp) {
    dplyr::full_join(
      cv_res$Covs[[comp]] |>
        tibble::as_tibble(rownames = "keepY") |>
        tidyr::pivot_longer(
          cols = -keepY,
          names_to = "keepX",
          values_to = "Covariance"
        ),
      cv_res$SEcov[[comp]] |>
        tibble::as_tibble(rownames = "keepY") |>
        tidyr::pivot_longer(
          cols = -keepY,
          names_to = "keepX",
          values_to = "SE"
        ),
      by = c("keepX", "keepY")
    ) |>
      dplyr::mutate(Comp = comp)
  }) |>
    purrr::reduce(dplyr::bind_rows) |>
    dplyr::mutate(dplyr::across(tidyselect::starts_with("keep"), ~ as.numeric(.x))) |>
    dplyr::left_join(choice_keep, by = c("Comp", "keepX", "keepY"))

  df |>
    dplyr::mutate(Comp = paste0("Joint comp. ", Comp)) |>
    ggplot2::ggplot(aes(x = keepX, y = keepY)) +
    ggplot2::geom_point(aes(size = SE, colour = col), shape = 21, fill = "gray70", colour = "gray60") +
    ggplot2::geom_point(aes(fill = Covariance, colour = col), size = 3, shape = 21) +
    ggplot2::scale_colour_manual(
      values = c("1SD" = "red", "NoSD" = "orange"),
      labels = c("1SD" = "1SD rule", "NoSD" = "Max"),
      na.value = "black",
      guide = ggplot2::guide_legend(title.position = "top")
    ) +
    ggplot2::scale_fill_viridis_c(option = "viridis", direction = 1, guide = ggplot2::guide_colourbar(title.position = "top")) +
    ggplot2::scale_size(range = c(3, 6), guide = ggplot2::guide_legend(title.position = "top")) +
    ggplot2::facet_wrap(~Comp) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.title = ggplot2::element_text(hjust = 0.5),
      legend.position = "bottom"
    ) +
    ggplot2::labs(
      title = paste0("Mean covariance between joint scores"),
      x = paste0("Nb of features retained in ", names_datasets[1], " dataset"),
      y = paste0("Nb of features retained in\n", names_datasets[2], " dataset"),
      fill = "Covariance",
      colour = "Optimal combination",
      size = "SE"
    )
}


#' Wrapper for OmicsPLS::o2m function
#'
#' Wrapper function around the \code{\link[OmicsPLS]{o2m}} function.
#' The main purpose of this wrapper is to add to the result the names of the datasets
#' to facilitate plotting.
#'
#' @param omicspls_input A named list of length 2, produced by \code{\link{get_input_omicspls}}.
#' @param cv_res Named integer vector of length 3, with names `n`, `nx`, `ny`.
#' Should be obtained with \code{\link{so2pls_get_optim_ncomp_adj}} or
#' \code{\link{so2pls_get_optim_ncomp}}.
#' @param sparsity_res Named list of length 2, with names `keepx` and `keepy`.
#' Should be obtained with \code{\link{so2pls_get_optim_keep}}.
#' @param n Positive integer, number of joint components to compute.
#' Ignored if `cv_res` is not `NULL`.
#' @param nx Positive integer, number of specific components to compute for the
#' first dataset. Ignored if `cv_res` is not `NULL`.
#' @param ny Positive integer, number of specific components to compute for the
#' second dataset. Ignored if `cv_res` is not `NULL`.
#' @param sparse Logical, should feature selection be performed? Default value
#' is `FALSE`. If `sparsity_res` is not `NULL`, will be set to `TRUE`.
#' @param keepx Integer or integer vector of length `n`, number of features from the first
#' dataset to retain for each joint component. Ignored if `sparsity_res` is not `NULL`.
#' @param keepy Integer or integer vector of length `n`, number of features from the second
#' dataset to retain for each joint component. Ignored if `sparsity_res` is not `NULL`.
#' @param ... Other arguments passed to \code{\link[OmicsPLS]{o2m}}.
#' @returns A list (see \code{\link[OmicsPLS]{o2m}}).
#' @export
so2pls_o2m <- function(omicspls_input,
                       cv_res = NULL,
                       sparsity_res = NULL,
                       n = NULL,
                       nx = NULL,
                       ny = NULL,
                       sparse = FALSE,
                       keepx = NULL,
                       keepy = NULL,
                       ...) {
  if (!is.list(omicspls_input) || length(omicspls_input) != 2) stop("'omicspls_input' should be a list of length 2 produced by get_input_omicspls().")

  if (!is.null(cv_res)) {
    if (is.null(names(cv_res)) || !all(sort(names(cv_res)) == c("n", "nx", "ny"))) stop("'cv_res' argument should be a named vector of length 3 with names: 'n', 'nx', 'ny'.")

    n <- cv_res[["n"]]
    nx <- cv_res[["nx"]]
    ny <- cv_res[["ny"]]
  } else {
    if (is.null(n)) {
      stop("Need to provide either a cross-validation result through 'cv_res' argument or a value for 'n' argument.")
    }
    if (is.null(nx)) {
      stop("Need to provide either a cross-validation result through 'cv_res' argument or a value for 'nx' argument.")
    }
    if (is.null(ny)) {
      stop("Need to provide either a cross-validation result through 'cv_res' argument or a value for 'ny' argument.")
    }
  }

  if (!is.null(sparsity_res)) {
    if (is.null(names(sparsity_res)) || !all(sort(names(sparsity_res)) == c("keepx", "keepy"))) stop("'sparsity_res' argument should be a named list of length 2 with names: 'keepx', 'keepy'.")

    sparse <- TRUE
    keepx <- sparsity_res[["keepx"]]
    keepy <- sparsity_res[["keepy"]]
  } else {
    if (sparse && is.null(keepx)) {
      stop("'sparse' = TRUE: need to provide either a sparsity cross-validation result through 'sparsity_res' argument or a value for 'keepx' argument.")
    }

    if (sparse && is.null(keepy)) {
      stop("'sparse' = TRUE: need to provide either a sparsity cross-validation result through 'sparsity_res' argument or a value for 'keepy' argument.")
    }
  }

  ds_names <- names(omicspls_input)


  ## Do it this way rather than running the command directly because this
  ## way the call is more informative when using the summary function
  cmd <- paste0(
    "OmicsPLS::o2m(omicspls_input[['", ds_names[1], "']], ",
    "omicspls_input[['", ds_names[2], "']], ",
    "n = ", n, ", ",
    "nx = ", nx, ", ",
    "ny = ", ny, ", ",
    "sparse = ", sparse, ", ",
    "keepx = c(", paste0(keepx, collapse = ","), "), ",
    "keepy = c(", paste0(keepy, collapse = ","), "), ",
    "...)"
  )

  res <- eval(str2expression(cmd))

  attr(res, "datasets_name") <- names(omicspls_input)
  return(res)
}

#' Get list of latent components from sO2PLS results
#'
#' Extracts the list of joint and specific latent components sO2PLS results.
#'
#' @param so2pls_res The sO2PLS results generated with [get_output_so2pls()]
#' function.
#' @returns A list with the following elements:
#' * `joint`: character vector with the name of the joint latent components.
#' * `specific`: named list of length 2. Each element corresponds to a dataset
#'  (names of the list are the datasets name), and is a character vector with
#'  the name of the specific latent components for the corresponding dataset.
#' @export
so2pls_get_components <- function(so2pls_res) {
  ## for devtools::check
  ld <- ds <- data <- NULL

  if (!inherits(so2pls_res, "output_dimension_reduction")) {
    stop("Expecting a o2m object (from so2pls_o2m() function).")
  }

  ld_list <- get_latent_dimensions(so2pls_res)
  ld_joint <- ld_list[stringr::str_detect(ld_list, "joint component")]

  ld_specific <- tibble::tibble(
    ld = ld_list[stringr::str_detect(ld_list, "specific component")],
    ds = stringr::str_extract(ld, "\\w+(?= specific)")
  ) |>
    dplyr::group_by(ds) |>
    tidyr::nest() |>
    dplyr::mutate(
      data = purrr::map(
        data,
        ~ dplyr::pull(.x, "ld")
      )
    ) |>
    tibble::deframe()

  res <- list(
    "joint" = ld_joint,
    "specific" = ld_specific
  )

  return(res)
}

#' Plot summary of sO2PLS run
#'
#' Plots a summary of variation from an sO2PLS run (from
#' the \code{OmicsPLS} package).
#'
#' @param so2pls_res The output from the \code{\link[OmicsPLS]{o2m}} function.
#' @param datasets Optional, a character vector with the names of the datasets
#' for which selected features should be extracted. Default is `NULL`, i.e. both
#' datasets are considered.
#' @returns A \href{https://patchwork.data-imaginist.com/index.html}{patchwork} of ggplots.
#' @export
so2pls_plot_summary <- function(so2pls_res, datasets = NULL) {
  ## for devtools::check
  value <- NULL

  if (!is(so2pls_res, "o2m")) stop("Expecting a o2m object. Make sure input is the output from the OmicsPLS::o2m() function.")

  ## Extract datasets name from the o2m results
  names_datasets <- attr(so2pls_res, "datasets_name")
  if (is.null(names_datasets)) names_datasets <- c("X", "Y")

  ## If datasets not given, will consider both datasets
  ## otherwise need to make sure that the names match the dataset names
  if (is.null(datasets)) {
    datasets <- names_datasets
  } else {
    .check_names(datasets, names_datasets, "'datasets' argument: '_W_' are not existing datasets. Possible values are: '_C_'.")
  }

  lname <- c("X", "Y") |>
    rlang::set_names(names_datasets)

  lname <- lname[datasets]

  summary_res <- summary(so2pls_res)

  ## For devtools::check
  dataset <- label <- type <- perc_explained <- NULL

  p1 <- lname |>
    purrr::map_dfr(
      ~ tibble::tibble(
        joint = summary_res[[paste0("R2_", .x, "joint")]],
        orthogonal = summary_res[[paste0("R2_", .x)]] - joint,
        noise = 1 - summary_res[[paste0("R2_", .x)]]
      ),
      .id = "dataset"
    ) |>
    tidyr::pivot_longer(
      cols = -dataset,
      names_to = "type",
      values_to = "value"
    ) |>
    dplyr::mutate(
      label = round(100 * value, 1),
      label = paste0(label, "%"),
      type = factor(type,
                    levels = c("joint", "orthogonal", "noise"),
                    labels = c("Joint variation", "Orthogonal variation", "Noise")
      ),
      dataset = factor(dataset, levels = datasets)
    ) |>
    ggplot2::ggplot(aes(x = type, y = value, fill = dataset)) +
    ggplot2::geom_col(
      position = "dodge",
      colour = "white"
    ) +
    ggplot2::geom_text(
      aes(label = label, group = dataset),
      hjust = 0.5,
      vjust = -0.3,
      position = ggplot2::position_dodge(width = 0.9)
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent,
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::labs(
      title = "Percentage of variation per dataset",
      x = NULL,
      y = NULL,
      fill = "Dataset"
    )

  p2 <- lname |>
    purrr::map_dfr(
      ~ tibble::tibble(
        perc_explained = summary_res[[paste0("R2_", .x, "pred")]]
      ),
      .id = "dataset"
    ) |>
    dplyr::mutate(
      label = round(100 * perc_explained, 1),
      label = paste0(label, "%"),
      dataset = factor(dataset, levels = datasets)
    ) |>
    ggplot2::ggplot(aes(x = dataset, y = perc_explained, fill = dataset)) +
    ggplot2::geom_col() +
    ggplot2::geom_text(
      aes(label = label),
      hjust = 0.5,
      vjust = -0.3
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent,
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::labs(
      title = "Percentage of joint variation\nexplained by joint component\nof other dataset",
      x = "Dataset",
      y = "Percentage of joint variation explained",
      fill = "Dataset"
    )

  patchwork::wrap_plots(
    p1, p2,
    widths = c(0.6, 0.4)
  )
}

#' Percentage of variance explained for sO2PLS
#'
#' Generates a table giving the percentage of variance explained by each
#' component from an sO2PLS in the corresponding dataset.
#'
#' @param so2pls_res The output from the \code{\link[OmicsPLS]{o2m}} function.
#' @returns A tibble with columns `latent_dimension`, `dataset` and `prop_var_expl`.
#' @export
so2pls_get_variance_explained <- function(so2pls_res) {
  if (!is(so2pls_res, "o2m")) stop("Expecting a o2m object. Make sure input is the output from the OmicsPLS::o2m() function.")

  ## Extract datasets name from the o2m results
  names_datasets <- attr(so2pls_res, "datasets_name")
  if (is.null(names_datasets)) names_datasets <- c("X", "Y")

  lname <- c("X", "Y") |>
    rlang::set_names(names_datasets)

  lname <- lname[names_datasets]

  ## Compute the variance explained by each component
  res_summary <- so2pls_res$flags

  n_specific_comps <- c(
    "X" = res_summary[["nx"]],
    "Y" = res_summary[["ny"]]
  )

  ## To deal with the case where there are no specific components for a dataset
  for (.x in c("varXorth", "varYorth")) {
    if (all(res_summary[[.x]] == 0)) res_summary[[.x]] <- NULL
  }

  lname |>
    purrr::imap_dfr(
      ~ tibble::tibble(
        var = c(
          res_summary[[paste0("var", .x, "joint")]],
          res_summary[[paste0("var", .x, "orth")]]
        ),
        ld_int = c(seq_len(res_summary$n), seq_len(n_specific_comps[[.x]])),
        ld_type = c(
          rep(" joint component ", res_summary$n),
          rep(" specific component ", n_specific_comps[[.x]])
        ),
        dataset = .y
      ) |>
        dplyr::mutate(
          var = var / res_summary[[paste0("ssq", .x)]],
          latent_dimension = paste0(dataset, ld_type, ld_int),
        ) |>
        dplyr::select(latent_dimension, dataset, prop_var_expl = var)
    )
}

#' Screeplot sO2PLS run
#'
#' Plots the percentage of variation explained by each latent component in
#' an sO2PLS (from the \code{OmicsPLS} package).
#'
#' Note that the plots are set up so that it is possible to add a custom colour
#' palette to get different colours for each dataset (see example).
#'
#' @param so2pls_res The output from the \code{\link[OmicsPLS]{o2m}} function.
#' @param datasets Optional, a character vector with the names of the datasets
#' for which selected features should be extracted. Default is `NULL`, i.e. both
#' datasets are considered.
#' @returns A \href{https://patchwork.data-imaginist.com/index.html}{patchwork} of ggplots.
#' @examples
#' \dontrun{
#' ## by default, same colour used for both datasets (cannot find a way to fix that cleanly)
#' so2pls_screeplot(so2pls_final_res)
#'
#' ## Add a colour palette to get different colour for each dataset
#' so2pls_screeplot(so2pls_final_res) &
#'  scale_fill_brewer(palette = "Set1", drop = F)
#' }
#' @export
so2pls_screeplot <- function(so2pls_res, datasets = NULL) {

  ## for devtools::check
  dataset <- prop_var_expl <- type <- component <- var <- var_label <- dataset2 <- latent_dimension <- NULL

  var_expl <- so2pls_get_variance_explained(so2pls_res)
  names_datasets <- unique(var_expl$dataset)

  ## If datasets not given, will consider both datasets
  ## otherwise need to make sure that the names match the dataset names
  if (is.null(datasets)) {
    datasets <- names_datasets
  } else {
    .check_names(datasets, names_datasets, "'datasets' argument: '_W_' are not existing datasets. Possible values are: '_C_'.")
  }

  var_expl |>
    dplyr::filter(dataset %in% datasets) |>
    dplyr::rename(var = prop_var_expl) |>
    dplyr::mutate(
      type = stringr::str_extract(latent_dimension, "(joint)|(specific)"),
      type = stringr::str_to_title(type),
      component = stringr::str_extract(latent_dimension, "\\d+$"),
      component = as.integer(component),
      var_label = round(100 * var, 1),
      var_label = paste0(var_label, "%"),
      component = factor(component),
      dataset = factor(dataset, levels = datasets),
      dataset2 = dataset
    ) |>
    dplyr::group_by(dataset2) |>
    tidyr::nest() |>
    tibble::deframe() |>
    purrr::imap(
      ~ ggplot2::ggplot(.x, aes(x = component, y = var)) +
        ggplot2::geom_col(fill = "steelblue") +
        ggplot2::geom_text(aes(label = var_label),
                           hjust = 0.5,
                           vjust = -0.3
        ) +
        ggplot2::facet_grid(. ~ type, space = "free", scale = "free_x") +
        ggplot2::scale_y_continuous(
          labels = scales::percent,
          expand = ggplot2::expansion(mult = c(0, 0.05))
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5),
          legend.position = "none"
        ) +
        ggplot2::labs(
          title = .y,
          x = "Latent components",
          y = "Percentage of variance explained"
        )
    ) |>
    patchwork::wrap_plots()
}


#' Computes average sample coordinates for sO2PLS joint components
#'
#' Computes the average sample coordinates for sO2PLS joint components
#' across the two datasets.
#'
#' @param so2pls_res The output from the \code{\link[OmicsPLS]{o2m}} function.
#' @returns A matrix of samples coordinates, with samples in rows and joint
#' components in columns.
#' @export
so2pls_get_wa_coord <- function(so2pls_res) {
  res <- (so2pls_res$Tt + so2pls_res$U) / 2

  return(res)
}


#' Compares sO2PLS samples joint component scores between the two datasets
#'
#' Plots a comparison of the samples joint component scores obtained for the two
#' datasets in an sO2PLS run (from the \code{OmicsPLS} package).
#'
#' @param so2pls_res The output from the \code{\link[OmicsPLS]{o2m}} function.
#' @param components Optional, an integer vector with the joint components
#' that should be plotted. Default is `NULL`, i.e. all joint components are
#' represented.
#' @param ... Further arguments passed to \code{\link{plot_samples_score_pair}}.
#' @returns A \href{https://patchwork.data-imaginist.com/index.html}{patchwork}
#' of ggplots.
#' @export
so2pls_compare_samples_joint_components <- function(so2pls_res,
                                                    components = NULL,
                                                    ...) {

  if (!is(so2pls_res, "o2m")) {
    stop("Expecting a o2m object. Make sure input is the output from the OmicsPLS::o2m() function.")
  }

  method_output <- get_output_so2pls(so2pls_res, use_average_dimensions = FALSE)
  ds <- attr(so2pls_res, "datasets_name")

  all_comps <- get_latent_dimensions(method_output)
  all_comps <- all_comps[stringr::str_detect(all_comps, "joint component")] |>
    stringr::str_extract("joint component \\d+") |>
    unique()

  if (!is.null(components)) {
    components <- paste0("joint component ", components)
    .check_names(components, all_comps, "'_W_' do not exist. Possible joint components are: '_C_'.")
  } else {
    components <- all_comps
  }

  fct <- function(.x, ...) {
    plot_samples_score_pair(
      method_output,
      latent_dimensions = paste(ds, .x),
      ...
    ) +
      ggplot2::labs(
        title = .x,
        x = ds[[1]],
        y = ds[[2]]
      )
  }

  p <- purrr::map(
    components,
    fct,
    ...
  ) |>
    patchwork::wrap_plots() +
    patchwork::plot_annotation(title = "Joint components samples score - sO2PLS") +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(hjust = 0.5)
    )

  return(p)
}

#' Plots sO2PLS contributions between datasets joint components
#'
#' Plots the regression coefficients that link the joint components of the two
#' datasets, from an SO2PLS run (from the \code{OmicsPLS} package).
#'
#' @param so2pls_res The output from the \code{\link[OmicsPLS]{o2m}} function.
#' @param datasets Optional, a character vector with the names of the datasets
#' that should be plotted. Default is `NULL`, i.e. both datasets are considered.
#' @returns A \href{https://patchwork.data-imaginist.com/index.html}{patchwork} of ggplots.
#' @export
so2pls_plot_joint_components_coefficients <- function(so2pls_res,
                                                      datasets = NULL) {
  if (!is(so2pls_res, "o2m")) stop("Expecting a o2m object. Make sure input is the output from the OmicsPLS::o2m() function.")

  ## Extract datasets name from the o2m results
  names_datasets <- attr(so2pls_res, "datasets_name")
  if (is.null(names_datasets)) names_datasets <- c("X", "Y")

  ## If datasets not given, will consider both datasets
  ## otherwise need to make sure that the names match the dataset names
  if (is.null(datasets)) {
    datasets <- names_datasets
  } else {
    .check_names(datasets, names_datasets, "'datasets' argument: '_W_' are not existing datasets. Possible values are: '_C_'.")
  }

  lname <- c("B_U", "B_T.") |>
    rlang::set_names(names_datasets)

  lname <- lname[datasets]

  ## for devtools::check
  other_comp <- comp <- contribution <- dataset <- NULL

  lname |>
    purrr::imap(function(.x, .y) {
      mat <- so2pls_res[[.x]]
      colnames(mat) <- paste0(seq_len(ncol(mat)))

      mat |>
        tibble::as_tibble() |>
        dplyr::mutate(other_comp = paste0(seq_len(nrow(mat)))) |>
        tidyr::pivot_longer(
          cols = -other_comp,
          names_to = "comp",
          values_to = "contribution"
        ) |>
        ggplot2::ggplot(aes(x = comp, y = contribution, fill = other_comp)) +
        ggplot2::geom_col(
          position = "dodge",
          colour = "black"
        ) +
        ggplot2::scale_fill_viridis_d(
          option = "plasma",
          guide = ggplot2::guide_legend(
            title.position = "left",
            title.hjust = 0.5,
            title.vjust = 0.5
          )
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "bottom") +
        ggplot2::labs(
          x = paste0(.y, " joint component"),
          y = "Contribution (regression coefficients)",
          fill = paste0(
            "Contribution of ", setdiff(names_datasets, .y),
            "\njoint component:"
          )
        )
    }) |>
    patchwork::wrap_plots()
}

#' Plots sO2PLS joint components samples scores
#'
#' Plots the samples scores for the average joint components from an
#' sO2PLS run (from the \code{OmicsPLS} package).
#'
#' @param so2pls_res The output from the \code{\link[OmicsPLS]{o2m}} function.
#' @param ... Further arguments passed to [plot_samples_score()].
#' @returns A ggmatrix plot.
#' @export
so2pls_plot_samples_joint_components <- function(so2pls_res,
                                                 ...) {
  if (!is(so2pls_res, "o2m")) {
    stop("Expecting a o2m object. Make sure input is the output from the OmicsPLS::o2m() function.")
  }

  method_output <- get_output_so2pls(so2pls_res)
  ld_list <- so2pls_get_components(method_output)[["joint"]]

  res <- plot_samples_score(
    method_output,
    ld_list,
    ...
  )

  return(res)
}


#' Plots sO2PLS specific components samples scores
#'
#' Plots the samples scores for the datasets specific components from an
#' sO2PLS run (from the \code{OmicsPLS} package).
#'
#' @param so2pls_res The output from the \code{\link[OmicsPLS]{o2m}} function.
#' @param dataset Character, the name of the dataset for which the specific
#' components should be plotted. Default is `NULL`, i.e. the specific components
#' of both datasets are plotted.
#' @param ... Further arguments passed to [plot_samples_score()].
#' @returns A list of ggmatrix plots (one per dataset), or one plot if `dataset`
#' was used to specify a dataset.
#' @export
so2pls_plot_samples_specific_components <- function(so2pls_res,
                                                    dataset = NULL,
                                                    ...) {
  if (!is(so2pls_res, "o2m")) stop("Expecting a o2m object. Make sure input is the output from the OmicsPLS::o2m() function.")

  ds <- attr(so2pls_res, "datasets_name")

  if (!is.null(dataset)) {
    .check_names(
      dataset,
      ds,
      "'dataset' argument: '_W_' is not a valid dataset name. Valid names are: '_C_'."
    )
  } else {
    dataset <- ds
  }

  method_output <- get_output_so2pls(so2pls_res)
  ld_list <- so2pls_get_components(method_output)[["specific"]]

  res <- dataset |>
    rlang::set_names() |>
    purrr::map(
      function(.x) {
        if (!(.x %in% names(ld_list))) {
          message(.x, " dataset has no specific components.")
          return(NULL)
        }
        plot_samples_score(
          method_output,
          ld_list[[.x]],
          ...
        )
      }
    ) |>
    purrr::discard(is.null)

  if (length(res) == 1) res <- res[[1]]

  return(res)
}
