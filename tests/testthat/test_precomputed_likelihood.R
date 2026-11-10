library(testthat)
library(MetricGraph)

test_that("Standard likelihood equals precomputed likelihood", {
  ## Define graph edges
  edge1 <- rbind(c(0, 0), c(1, 0))
  edge2 <- rbind(c(0, 0), c(0, 1))
  edge3 <- rbind(c(0, 1), c(-1, 1))
  theta_seq <- seq(from = pi, to = 3 * pi/2, length.out = 20)
  edge4 <- cbind(sin(theta_seq), 1 + cos(theta_seq))
  edges <- list(edge1, edge2, edge3, edge4)

  ## Create the metric graph object
  graph <- metric_graph$new(edges = edges)

  ## Model parameters
  kappa <- 0.2
  sigma <- 1.3
  sigma_e <- 0.1
  alpha <- 1
  n_obs_per_edge <- 75

  ## Generate observation locations on each edge
  PtE <- NULL
  for (i in 1:graph$nE) {
    PtE <- rbind(PtE, cbind(rep(i, n_obs_per_edge), runif(n_obs_per_edge)))
  }

  ## Simulate SPDE field and noisy observations
  u <- sample_spde(kappa = kappa, tau = 1/sigma^2, alpha = alpha,
                   graph = graph, PtE = PtE)
  y <- 1 + u + sigma_e * rnorm(n_obs_per_edge * graph$nE)

  ## Create observations data frame
  df_data <- data.frame(y = y,
                        edge_number = PtE[, 1],
                        distance_on_edge = PtE[, 2])

  ## Clear previous observations and add new ones to the graph
  graph$clear_observations()
  graph$add_observations(data = df_data, normalized = TRUE)

  ## Set up the parameter vector (named theta)
  theta.1 <- c(`std. dev` = sigma_e,
               tau = 1/sigma^2,
               kappa = kappa,
               `(Intercept)` = 1)

  ## Covariate matrix (intercept-only)
  X_cov <- rep(1, length(y))

  ## Compute likelihood using the standard approach
  lik <- MetricGraph:::likelihood_alpha1(theta = theta.1,
                                         graph = graph,
                                         data_name = NULL,
                                         manual_y = y,
                                         X_cov = as.matrix(X_cov),
                                         repl = 1,
                                         BC = 0,
                                         parameterization = "spde")

  ## Precompute required data for the alternative likelihood computation
  param <- MetricGraph:::precompute_alpha1(graph = graph,
                                           data_name = NULL,
                                           manual_y = y,
                                           X_cov = as.matrix(X_cov),
                                           repl = 1)

  ## Compute likelihood using precomputed data
  lik_pre <- MetricGraph:::likelihood_alpha1_precompute(theta = theta.1,
                                                        graph = graph,
                                                        precomputeddata = param,
                                                        data_name = NULL,
                                                        manual_y = y,
                                                        X_cov = as.matrix(X_cov),
                                                        repl = 1,
                                                        BC = 0,
                                                        parameterization = "spde")

  ## Test that both likelihood calculations match (within a tolerance)
  expect_equal(lik, lik_pre, tolerance = 1e-6)
})
