library(testthat)
library(MetricGraph)


make_full_cv_graph <- function(seed = 1L, missing_response = FALSE) {
  edges <- list(
    rbind(c(0, 0), c(1, 0)),
    rbind(c(1, 1), c(1, 0)),
    rbind(c(-1, 0), c(0, 0))
  )
  graph <- metric_graph$new(edges = edges)
  set.seed(seed)
  n_per_edge <- 6L
  edge_number <- rep(1:2, each = n_per_edge)
  distance_on_edge <- unlist(lapply(
    1:2,
    function(edge) sort(stats::runif(n_per_edge, 0.1, 0.9))
  ))
  response <- stats::rnorm(length(edge_number))
  if (missing_response) response[3L] <- NA_real_
  graph$add_observations(
    data = data.frame(
      y = response,
      edge_number = edge_number,
      distance_on_edge = distance_on_edge
    ),
    normalized = TRUE
  )
  X <- cbind(
    `(Intercept)` = 1,
    ELEV = stats::runif(length(response)),
    SLOPE = stats::runif(length(response)),
    PRECIP = stats::runif(length(response))
  )
  list(graph = graph, y = response, X = X)
}


dense_wm1_loo <- function(theta, graph, precomputed_data, BC = 1L) {
  sigma_e <- exp(theta[1L])
  reciprocal_tau <- exp(theta[2L])
  kappa <- exp(theta[3L])
  Q_list <- MetricGraph:::spde_precision(
    kappa = kappa,
    tau = 1 / reciprocal_tau,
    alpha = 1L,
    graph = graph,
    build = FALSE,
    BC = BC
  )
  Q <- as.matrix(Matrix::sparseMatrix(
    i = Q_list$i,
    j = Q_list$j,
    x = Q_list$x,
    dims = Q_list$dims
  ))

  edge_cache <- precomputed_data$edge_cache[[1L]]
  n_cache <- sum(lengths(edge_cache$y))
  Sinv <- matrix(0, n_cache, n_cache)
  C <- matrix(0, nrow(graph$V), n_cache)
  BtSinvB <- matrix(0, nrow(graph$V), nrow(graph$V))
  y <- numeric(0L)
  X <- matrix(0, 0L, precomputed_data$n_cov)
  cache_position <- 0L

  for (j in seq_along(edge_cache$e)) {
    y_i <- edge_cache$y[[j]]
    X_i <- edge_cache$X[[j]]
    D_i <- edge_cache$D[[j]]
    S <- MetricGraph:::r_1(
      D_i,
      kappa = kappa,
      tau = 1 / reciprocal_tau
    )
    endpoints <- 1L:2L
    observations <- -endpoints
    Bt <- solve(
      S[endpoints, endpoints, drop = FALSE],
      S[endpoints, observations, drop = FALSE]
    )
    Sigma_i <- S[observations, observations, drop = FALSE] -
      S[observations, endpoints, drop = FALSE] %*% Bt
    diag(Sigma_i) <- diag(Sigma_i) + sigma_e^2
    Sinv_i <- solve(Sigma_i)
    Sigma_iB <- Sinv_i %*% t(Bt)
    vertices <- edge_cache$E[j, ]
    positions <- cache_position + seq_along(y_i)

    Sinv[positions, positions] <- Sinv_i
    C[vertices, positions] <- C[vertices, positions, drop = FALSE] +
      t(Sigma_iB)
    BtSinvB[vertices, vertices] <-
      BtSinvB[vertices, vertices, drop = FALSE] + Bt %*% Sigma_iB
    y <- c(y, y_i)
    X <- rbind(X, X_i)
    cache_position <- cache_position + length(y_i)
  }

  P <- Sinv - t(C) %*% solve(Q + BtSinvB, C)
  H <- crossprod(X, P %*% X)
  beta_hat <- as.vector(solve(H, crossprod(X, P %*% y)))
  residual <- y - as.vector(X %*% beta_hat)
  Pr <- as.vector(P %*% residual)
  d <- diag(P)
  list(
    mu = y - Pr / d,
    var = 1 / d,
    beta_hat = beta_hat,
    H = H
  )
}


test_that("directional selected inverse matches the dense reference", {
  setup <- make_full_cv_graph(seed = 101L, missing_response = TRUE)
  theta <- c(log(0.4), log(0.8), log(0.9))
  precomputed <- MetricGraph:::precompute_alpha1_directional(
    setup$graph,
    data_name = "y",
    X_cov = setup$X
  )

  selected <- MetricGraph:::cv_core_alpha1_directional(
    theta,
    precomputed,
    method = "selinv"
  )
  reference <- MetricGraph:::cv_core_alpha1_directional(
    theta,
    precomputed,
    method = "reference"
  )

  expect_equal(selected$mu, reference$mu, tolerance = 1e-8)
  expect_equal(selected$var, reference$var, tolerance = 1e-8)
  expect_equal(selected$beta_hat, reference$beta_hat, tolerance = 1e-8)
  expect_equal(unname(selected$H), unname(reference$H), tolerance = 1e-8)
  expect_equal(selected$idx, setdiff(seq_along(setup$y), 3L))
  expect_length(selected$beta_hat, 4L)
  expect_true(all(is.finite(selected$mu)))
  expect_true(all(is.finite(selected$var) & selected$var > 0))
})


test_that("WM1 plug-in LOO matches an explicit dense Gaussian calculation", {
  setup <- make_full_cv_graph(seed = 202L)
  theta <- c(log(0.5), log(0.7), log(1.1))
  precomputed <- MetricGraph:::precompute_alpha1(
    setup$graph,
    data_name = "y",
    X_cov = setup$X,
    repl = NULL
  )

  actual <- MetricGraph:::cv_core_alpha1(
    theta,
    setup$graph,
    precomputed,
    y_resp_full = setup$y,
    BC = 1L
  )
  expected <- dense_wm1_loo(theta, setup$graph, precomputed, BC = 1L)

  expect_equal(actual$mu, expected$mu, tolerance = 1e-8)
  expect_equal(actual$var, expected$var, tolerance = 1e-8)
  expect_equal(actual$beta_hat, expected$beta_hat, tolerance = 1e-8)
  expect_equal(unname(actual$H), unname(expected$H), tolerance = 1e-8)
  expect_identical(actual$idx, seq_along(setup$y))
})


test_that("the packaged Columbia component has audited content", {
  environment <- new.env(parent = globalenv())
  utils::data(
    "columbia_main_component",
    package = "MetricGraph",
    envir = environment
  )
  component <- environment$columbia_main_component

  expect_s3_class(component, "columbia_main_component")
  expect_equal(nrow(component$edges), 18668L)
  expect_equal(nrow(component$observations), 2080L)
  expect_identical(
    component$fingerprint,
    "sha256:eaf68cffcdd7e33dbbbc44744df3054a0fcfef1dac4990021b02fa3f24a6d2f7"
  )
  expect_true(all(c(
    "STREAM_AUG", "ELEV", "SLOPE", "PRECIP", "columbia_obs_id",
    "edge_number", "distance_on_edge"
  ) %in% names(component$observations)))
  expect_false(anyDuplicated(component$observations$columbia_obs_id) > 0L)

  copyright_path <- system.file("COPYRIGHTS", package = "MetricGraph")
  expect_true(nzchar(copyright_path))
  copyright <- paste(readLines(copyright_path, warn = FALSE), collapse = "\n")
  expect_match(copyright, "Creative Commons Attribution 4.0")
  expect_match(copyright, "10.6084/m9.figshare.24132840.v1", fixed = TRUE)
})


test_that("the packaged Columbia component reconstructs in both directions", {
  skip_on_cran()
  environment <- new.env(parent = globalenv())
  utils::data(
    "columbia_main_component",
    package = "MetricGraph",
    envir = environment
  )
  helper <- test_path(
    "..", "..", "examples", "directional", "columbia_full_graph_helpers.R"
  )
  # 'examples' is in .Rbuildignore, so the helper is only available when the
  # tests run from the source repository, not from a built package.
  skip_if_not(
    file.exists(helper),
    "examples/ is not shipped with the built package."
  )
  source(helper, local = environment)

  tampered <- environment$columbia_main_component
  tampered$observations$ELEV[1L] <- tampered$observations$ELEV[1L] + 1
  expect_error(
    environment$.columbia_require_component(tampered),
    "fingerprint"
  )

  original <- environment$columbia_make_graph(
    environment$columbia_main_component,
    reversed = FALSE
  )
  reversed <- environment$columbia_make_graph(
    environment$columbia_main_component,
    reversed = TRUE
  )
  original_data <- original$get_data(format = "tibble", drop_na = FALSE)
  reversed_data <- reversed$get_data(format = "tibble", drop_na = FALSE)

  expect_equal(original$nE, 18668L)
  expect_equal(reversed$nE, 18668L)
  expect_equal(nrow(original_data), 2080L)
  expect_setequal(
    original_data$columbia_obs_id,
    reversed_data$columbia_obs_id
  )
})
