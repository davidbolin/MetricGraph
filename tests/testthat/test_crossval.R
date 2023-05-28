test_that("Check likelihoods for alternative models", {
  set.seed(1)
  V <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1),
             c(-1, 1), c(-1, 0), c(0, -1))
  E <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 5),
             c(5, 6), c(6, 1), c(4, 1), c(1, 7))

  graph <- metric_graph$new(V = V, E = E)

  kappa <- 10
  tau <- 1/20
  sigma_e <- 0.1
  theta <-  c(sigma_e, 1/tau, kappa)

  n.obs.per.edge <- 5
  PtE <- NULL
  for(i in 1:graph$nE){
    PtE <- rbind(PtE, cbind(rep(i, n.obs.per.edge), (runif(n.obs.per.edge))))
  }

  u <- sample_spde(kappa = kappa, tau = tau, alpha = 2,
                   graph = graph, PtE = PtE)

  y <- u + sigma_e*rnorm(n.obs.per.edge * graph$nE)

  df_temp <- data.frame(y = y, edge_number = PtE[,1], distance_on_edge = PtE[,2])

  graph$add_observations(data=df_temp, normalized = TRUE)
  graph$compute_resdist()
  lik.exp.v1 <- likelihood_graph_covariance(graph, model = "isoCov", cov_function = exp_covariance, log_scale = FALSE, repl = NULL, X_cov = NULL, y_graph = graph$data$y)
  lik.exp.v1 <- lik.exp.v1(theta)


  graph$observation_to_vertex()
  graph$compute_geodist()
  graph$compute_resdist()
  lik.exp.v2 <- likelihood_graph_covariance(graph, model = "isoCov", cov_function = exp_covariance, log_scale = FALSE, y_graph = graph$data$y, X_cov = NULL, repl = NULL)
  lik.exp.v2 <- lik.exp.v2(theta)

  graph$compute_laplacian(full=TRUE)
  lik.gl1.v1 <- likelihood_graph_covariance(graph, model = "GL1", log_scale = FALSE, y_graph = graph$data$y, repl = NULL, X_cov = NULL)
  lik.gl1.v1 <- lik.gl1.v1(theta)
  lik.gl2.v1 <- likelihood_graph_covariance(graph, model = "GL2", log_scale = FALSE, y_graph = graph$data$y, repl = NULL, X_cov = NULL)
  lik.gl2.v1 <- lik.gl2.v1(theta)

  lik.gl1.v2 <- likelihood_graph_laplacian(graph, alpha = 1, y_graph = graph$data$y, repl = NULL, X_cov = NULL, parameterization = "spde")
  lik.gl1.v2 <- lik.gl1.v2(log(theta))
  lik.gl2.v2 <- likelihood_graph_laplacian(graph, alpha = 2, y_graph = graph$data$y, repl = NULL, X_cov = NULL, parameterization = "spde")
  lik.gl2.v2 <- lik.gl2.v2(log(theta))

  expect_equal(lik.exp.v1, lik.exp.v2, tolerance = 1e-10)
  expect_equal(lik.gl1.v1, lik.gl1.v2, tolerance = 1e-10)
  expect_equal(lik.gl2.v1, lik.gl2.v2, tolerance = 1e-10)
})


test_that("Check agrement beteen covariance and precision cross validation", {
  set.seed(1)
  V <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1),
             c(-1, 1), c(-1, 0), c(0, -1))
  E <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 5),
             c(5, 6), c(6, 1), c(4, 1), c(1, 7))

  graph <- metric_graph$new(V = V, E = E)

  kappa <- 10
  tau <- 1/20
  sigma_e <- 0.1
  theta <-  c(sigma_e, 1/tau, kappa)

  n.obs.per.edge <- 5
  PtE <- NULL
  for(i in 1:graph$nE){
    PtE <- rbind(PtE, cbind(rep(i, n.obs.per.edge), runif(n.obs.per.edge)))
  }

  u <- sample_spde(kappa = kappa, tau = tau, alpha = 2,
                   graph = graph, PtE = PtE)

  y <- u + sigma_e*rnorm(n.obs.per.edge * graph$nE)

  df_temp <- data.frame(y = y, edge_number = PtE[,1], distance_on_edge = PtE[,2])

  graph$add_observations(data=df_temp, normalized = TRUE)

  graph$observation_to_vertex()
  cv.alpha1.v1 <- posterior_crossvalidation_manual(theta, graph, model = "alpha1", data_name = "y")
  cv.alpha1.v2 <- posterior_crossvalidation_covariance_manual(theta, graph,
                                                       model = "alpha1", data_name = "y")
  expect_equal(cv.alpha1.v1$mu, cv.alpha1.v2$mu, tolerance = 1e-10)
  expect_equal(cv.alpha1.v1$var, cv.alpha1.v2$var, tolerance = 1e-10)

  cv.alpha2.v1 <- posterior_crossvalidation_manual(theta, graph, model = "alpha2", data_name = "y")
  cv.alpha2.v2 <- posterior_crossvalidation_covariance_manual(theta, graph,
                                                       model = "alpha2", data_name = "y")
  expect_equal(cv.alpha2.v1$mu, cv.alpha2.v2$mu, tolerance = 1e-10)
  expect_equal(cv.alpha2.v1$var, cv.alpha2.v2$var, tolerance = 1e-10)

  graph$compute_laplacian()
  cv.GL1.v1 <- posterior_crossvalidation_manual(theta, graph, model = "GL1", data_name = "y")
  cv.GL1.v2 <- posterior_crossvalidation_covariance_manual(theta, graph, model = "GL1", data_name = "y")
  expect_equal(cv.GL1.v1$mu, cv.GL1.v2$mu, tolerance = 1e-10)
  expect_equal(cv.GL1.v1$var, cv.GL1.v2$var, tolerance = 1e-10)

  cv.GL2.v1 <- posterior_crossvalidation_manual(theta, graph, model = "GL2", data_name = "y")
  cv.GL2.v2 <- posterior_crossvalidation_covariance_manual(theta, graph, model = "GL2", data_name = "y")
  expect_equal(cv.GL2.v1$mu, cv.GL2.v2$mu, tolerance = 1e-10)
  expect_equal(cv.GL2.v1$var, cv.GL2.v2$var, tolerance = 1e-10)
})
