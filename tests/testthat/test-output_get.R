test_that("get_output_pca works", {
  pca_res <- test_get_pca_run()

  expect_error(
    get_output_pca("TEST"),
    "Expecting a pcaRes object (from run_pca() or pcaMethods::pca() function).",
    fixed = TRUE
  )

  res <- get_output_pca(pca_res)

  ## Testing general format
  test_helper_output(res, "PCA")
  expect_equal(attr(res, "dataset"), "rnaseq")

  ## Testing dimensions
  expect_equal(nrow(res$features_weight), 10 * 35)
  expect_equal(nrow(res$samples_score), 10 * 15)
  expect_equal(nrow(res$variance_explained), 10)

  ## Testing latent dimensions levels
  expect_equal(levels(res$features_weight$latent_dimension), paste0("Principal component ", 1:10))
  expect_equal(levels(res$samples_score$latent_dimension), paste0("Principal component ", 1:10))
  expect_equal(levels(res$variance_explained$latent_dimension), paste0("Principal component ", 1:10))


  ## Testing datasets levels
  expect_equal(levels(res$features_weight$dataset), c("rnaseq"))
  expect_equal(levels(res$variance_explained$dataset), c("rnaseq"))

  ## Testing values
  expect_equal(test_extract_features_weight(res, "Principal component 2", "rnaseq"), pca_res@loadings[, "PC2", drop = TRUE])

  expect_equal(test_extract_samples_score(res, "Principal component 2"), pca_res@scores[, "PC2", drop = TRUE])

  expect_equal(
    res$variance_explained$prop_var_expl,
    unname(pca_res@R2)
  )
})

test_that("get_output_splsda works", {
  splsda_res <- test_get_splsda_run()

  expect_error(
    get_output_splsda("TEST"),
    "Expecting a mixo_splsda object (from run_splsda() function).",
    fixed = TRUE
  )

  res <- get_output_splsda(splsda_res)

  ## Testing general format
  test_helper_output(res, "sPLS-DA")
  expect_equal(attr(res, "dataset"), "rnaseq")

  ## Testing dimensions
  expect_equal(nrow(res$features_weight), 2 * 35)
  expect_equal(nrow(res$samples_score), 2 * 15)
  expect_equal(nrow(res$variance_explained), 2)

  ## Testing latent dimensions levels
  expect_equal(levels(res$features_weight$latent_dimension), paste0("Component ", 1:2))
  expect_equal(levels(res$samples_score$latent_dimension), paste0("Component ", 1:2))
  expect_equal(levels(res$variance_explained$latent_dimension), paste0("Component ", 1:2))


  ## Testing datasets levels
  expect_equal(levels(res$features_weight$dataset), c("rnaseq"))
  expect_equal(levels(res$variance_explained$dataset), c("rnaseq"))

  ## Testing values
  expect_equal(test_extract_features_weight(res, "Component 2", "rnaseq"), splsda_res$loadings$X[, "comp2", drop = TRUE])

  expect_equal(test_extract_samples_score(res, "Component 2"), splsda_res$variates$X[, "comp2", drop = TRUE])

  expect_equal(
    res$variance_explained$prop_var_expl,
    unname(splsda_res$prop_expl_var$X)
  )
})

test_that("get_output_spls works", {
  spls_res <- test_get_spls_run()

  expect_error(
    get_output_spls("TEST"),
    "Expecting a mixo_pls or mixo_spls object (from spls_run() function).",
    fixed = TRUE
  )

  expect_error(
    get_output_spls(test_get_splsda_run()),
    "Expecting output from sPLS and not sPLS-DA. Please use get_output_splsda() instead.",
    fixed = TRUE
  )

  res1 <- get_output_spls(spls_res, use_average_dimensions = TRUE)
  res2 <- get_output_spls(spls_res, use_average_dimensions = FALSE)

  ## Testing general format
  test_helper_output(res1, "sPLS")
  test_helper_output(res2, "sPLS")

  ## Testing dimensions
  expect_equal(nrow(res1$features_weight), 2 * (35 + 40))
  expect_equal(nrow(res2$features_weight), 2 * (35 + 40))

  expect_equal(nrow(res1$samples_score), 2 * 10)
  expect_equal(nrow(res2$samples_score), 2 * 2 * 10)

  expect_equal(nrow(res1$variance_explained), 4)
  expect_equal(nrow(res2$variance_explained), 4)

  ## Testing latent dimensions levels
  for (i in c("features_weight", "samples_score", "variance_explained")) {
    expect_equal(levels(res1[[i]]$latent_dimension), c("Component 1", "Component 2"))
    expect_equal(
      levels(res2[[i]]$latent_dimension),
      c(
        "rnaseq Component 1", "metabolome Component 1",
        "rnaseq Component 2", "metabolome Component 2"
      )
    )
  }

  ## Testing datasets levels
  for (i in c("features_weight", "variance_explained")) {
    expect_equal(levels(res1[[i]]$dataset), c("rnaseq", "metabolome"))
    expect_equal(levels(res2[[i]]$dataset), c("rnaseq", "metabolome"))
  }

  ## Testing values
  expect_equal(test_extract_features_weight(res1, "Component 2", "metabolome"), spls_res$loadings$Y[, "comp2", drop = TRUE])
  expect_equal(test_extract_features_weight(res2, "metabolome Component 2", "metabolome"), spls_res$loadings$Y[, "comp2", drop = TRUE])

  expect_equal(test_extract_samples_score(res1, "Component 2"), spls_get_wa_coord(spls_res)[, "comp2", drop = TRUE])
  expect_equal(test_extract_samples_score(res2, "metabolome Component 2"), spls_res$variates$Y[, "comp2", drop = TRUE])

  expect_equal(
    res1$variance_explained |>
      dplyr::arrange(dataset, latent_dimension) |>
      dplyr::pull(prop_var_expl),
    unname(unlist(spls_res$prop_expl_var))
  )
  expect_equal(
    res2$variance_explained |>
      dplyr::arrange(dataset, latent_dimension) |>
      dplyr::pull(prop_var_expl),
    unname(unlist(spls_res$prop_expl_var))
  )
})

test_that("get_output_diablo works", {
  diablo_res <- test_get_diablo_run()

  expect_error(
    get_output_diablo("TEST"),
    "Expecting a block.plsda or block.splsda object (from diablo_run() function).",
    fixed = TRUE
  )

  res1 <- get_output_diablo(diablo_res, use_average_dimensions = TRUE)
  res2 <- get_output_diablo(diablo_res, use_average_dimensions = FALSE)

  ## Testing general format
  test_helper_output(res1, "DIABLO")
  test_helper_output(res2, "DIABLO")

  ## Testing dimensions
  expect_equal(nrow(res1$features_weight), 2 * (30 + 35 + 40))
  expect_equal(nrow(res2$features_weight), 2 * (30 + 35 + 40))

  expect_equal(nrow(res1$samples_score), 2 * 5)
  expect_equal(nrow(res2$samples_score), 2 * 3 * 5)

  expect_equal(nrow(res1$variance_explained), 6)
  expect_equal(nrow(res2$variance_explained), 6)

  ## Testing latent dimensions levels
  for (i in c("features_weight", "samples_score", "variance_explained")) {
    expect_equal(levels(res1[[i]]$latent_dimension), c("Component 1", "Component 2"))
    expect_equal(
      levels(res2[[i]]$latent_dimension),
      c(
        "snps+A Component 1", "rnaseq Component 1", "metabolome Component 1",
        "snps+A Component 2", "rnaseq Component 2", "metabolome Component 2"
      )
    )
  }

  ## Testing datasets levels
  for(i in c("features_weight", "variance_explained")) {
    expect_equal(levels(res1[[i]]$dataset), c("snps+A", "rnaseq", "metabolome"))
    expect_equal(levels(res2[[i]]$dataset), c("snps+A", "rnaseq", "metabolome"))
  }

  ## Testing values
  expect_equal(test_extract_features_weight(res1, "Component 2", "metabolome"), diablo_res$loadings$metabolome[, "comp2", drop = TRUE])
  expect_equal(test_extract_features_weight(res2, "metabolome Component 2", "metabolome"), diablo_res$loadings$metabolome[, "comp2", drop = TRUE])

  expect_equal(test_extract_samples_score(res1, "Component 2"), diablo_get_wa_coord(diablo_res)[, "comp2", drop = TRUE])
  expect_equal(test_extract_samples_score(res2, "metabolome Component 2"), diablo_res$variates$metabolome[, "comp2", drop = TRUE])

  expect_equal(
    res1$variance_explained |>
      dplyr::arrange(dataset, latent_dimension) |>
      dplyr::pull(prop_var_expl),
    diablo_res$prop_expl_var[setdiff(names(diablo_res$prop_expl_var), "Y")] |>
      unlist() |>
      unname()
  )

  expect_equal(
    res2$variance_explained |>
      dplyr::arrange(dataset, latent_dimension) |>
      dplyr::pull(prop_var_expl),
    diablo_res$prop_expl_var[setdiff(names(diablo_res$prop_expl_var), "Y")] |>
      unlist() |>
      unname()
  )
})

test_that("get_output_mofa2 works", {
  mofa_res <- test_get_mofa_run()

  expect_error(
    get_output_mofa2("TEST"),
    "Expecting a MOFA object (from MOFA2::run_mofa() function).",
    fixed = TRUE
  )

  res <- get_output_mofa2(mofa_res)

  ## Testing general format
  test_helper_output(res, "MOFA")

  ## Testing dimensions
  expect_equal(nrow(res$features_weight), 1 * (30 + 35 + 40))
  expect_equal(nrow(res$samples_score), 1 * 25)
  expect_equal(nrow(res$variance_explained), 3)

  ## Testing latent dimensions levels
  expect_equal(levels(res$features_weight$latent_dimension), "Factor 1")
  expect_equal(levels(res$samples_score$latent_dimension), "Factor 1")
  expect_equal(levels(res$variance_explained$latent_dimension), "Factor 1")


  ## Testing datasets levels
  expect_equal(levels(res$features_weight$dataset), c("snps+A", "rnaseq", "metabolome"))
  expect_equal(levels(res$variance_explained$dataset), c("snps+A", "rnaseq", "metabolome"))

  ## Testing values
  expect_equal(test_extract_features_weight(res, "Factor 1", "rnaseq"), mofa_res@expectations$W$rnaseq[, 1, drop = TRUE])
  expect_equal(test_extract_samples_score(res, "Factor 1"), mofa_res@expectations$Z$group1[, 1, drop = TRUE])

  expect_equal(
    res$variance_explained |>
      dplyr::arrange(dataset, latent_dimension) |>
      dplyr::pull(prop_var_expl),
    MOFA2::get_variance_explained(mofa_res, as.data.frame = TRUE)$r2_per_factor$value / 100
  )
})

test_that("get_output_so2pls works", {
  o2m_res <- test_get_so2pls_run()

  expect_error(
    get_output_so2pls("TEST"),
    "Expecting a o2m object (from so2pls_o2m() function).",
    fixed = TRUE
  )

  res1 <- get_output_so2pls(o2m_res, use_average_dimensions = TRUE)
  res2 <- get_output_so2pls(o2m_res, use_average_dimensions = FALSE)

  ## Testing general format
  test_helper_output(res1, "sO2PLS")
  test_helper_output(res2, "sO2PLS")

  ## Testing dimensions
  expect_equal(nrow(res1$features_weight), 3 * (35 + 40))
  expect_equal(nrow(res2$features_weight), 3 * (35 + 40))

  expect_equal(nrow(res1$samples_score), 4 * 10)
  expect_equal(nrow(res2$samples_score), (2 * 2 + 2) * 10)

  expect_equal(nrow(res1$variance_explained), 6)
  expect_equal(nrow(res2$variance_explained), 6)

  ## Testing latent dimensions levels
  for (i in c("features_weight", "samples_score", "variance_explained")) {
    expect_equal(
      levels(res1[[i]]$latent_dimension),
      c("joint component 1", "joint component 2", "rnaseq specific component 1", "metabolome specific component 1")
    )
    expect_equal(
      levels(res2[[i]]$latent_dimension),
      c(
        "rnaseq joint component 1", "metabolome joint component 1",
        "rnaseq joint component 2", "metabolome joint component 2",
        "rnaseq specific component 1", "metabolome specific component 1"
      )
    )
  }

  ## Testing datasets levels
  for (i in c("features_weight", "variance_explained")) {
    expect_equal(levels(res1[[i]]$dataset), c("rnaseq", "metabolome"))
    expect_equal(levels(res2[[i]]$dataset), c("rnaseq", "metabolome"))
  }


  ## Testing values
  expect_equal(test_extract_features_weight(res1, "joint component 2", "metabolome"), o2m_res$C.[, 2, drop = TRUE])
  expect_equal(test_extract_features_weight(res2, "metabolome joint component 2", "metabolome"), o2m_res$C.[, 2, drop = TRUE])

  expect_equal(test_extract_features_weight(res1, "rnaseq specific component 1", "rnaseq"), o2m_res$P_Yosc.[, 1, drop = TRUE])
  expect_equal(test_extract_features_weight(res2, "rnaseq specific component 1", "rnaseq"), o2m_res$P_Yosc.[, 1, drop = TRUE])

  expect_equal(test_extract_samples_score(res1, "joint component 2"), (o2m_res$Tt[, 2, drop = TRUE] + o2m_res$U[, 2, drop = TRUE]) / 2)
  expect_equal(test_extract_samples_score(res2, "metabolome joint component 2"), o2m_res$U[, 2, drop = TRUE])

  expect_equal(test_extract_samples_score(res1, "rnaseq specific component 1"), o2m_res$T_Yosc[, 1, drop = TRUE])
  expect_equal(test_extract_samples_score(res2, "rnaseq specific component 1"), o2m_res$T_Yosc[, 1, drop = TRUE])

  expect_equal(
    res1$variance_explained |>
      dplyr::arrange(dataset) |>
      dplyr::pull(prop_var_expl),
    so2pls_get_variance_explained(o2m_res)$prop_var_expl
  )
  expect_equal(
    res2$variance_explained |>
      dplyr::arrange(dataset) |>
      dplyr::pull(prop_var_expl),
    so2pls_get_variance_explained(o2m_res)$prop_var_expl
  )
})
