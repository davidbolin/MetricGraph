test_that("Check likelihoods for alternative models", {
  V <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1),
             c(-1, 1), c(-1, 0), c(0, -1))
  E <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 5),
             c(5, 6), c(6, 1), c(4, 1), c(1, 7))

  graph <- metric_graph$new(V = V, E = E)

  kappa <- 10
  sigma <- 20
  sigma_e <- 0.1
  theta <-  c(sigma_e, sigma, kappa)

  n.obs.per.edge <- 5
  PtE <- NULL
  for(i in 1:graph$nE){
    PtE <- rbind(PtE, cbind(rep(i, n.obs.per.edge), (runif(n.obs.per.edge))))
  }

  u <- sample_spde(kappa = kappa, sigma = sigma, alpha = 2,
                   graph = graph, PtE = PtE)

  y <- u + sigma_e*rnorm(n.obs.per.edge * graph$nE)

  graph$add_PtE_observations(y,PtE, normalized = TRUE)
  graph$compute_resdist()
  lik.exp.v1 <- likelihood_graph_covariance(theta, graph, model = "isoExp")


  graph$observation_to_vertex()
  graph$compute_geodist()
  graph$compute_resdist()
  lik.exp.v2 <- likelihood_graph_covariance(theta, graph, model = "isoExp")

  graph$compute_laplacian()
  lik.gl1.v1 <- likelihood_graph_covariance(theta, graph, model = "GL1")
  lik.gl2.v1 <- likelihood_graph_covariance(theta, graph, model = "GL2")

  lik.gl1.v2 <- likelihood_graph_laplacian(theta, graph, alpha = 1)
  lik.gl2.v2 <- likelihood_graph_laplacian(theta, graph, alpha = 2)

  expect_equal(lik.exp.v1, lik.exp.v2, tolerance = 1e-10)
  expect_equal(lik.gl1.v1, lik.gl1.v2, tolerance = 1e-10)
  expect_equal(lik.gl2.v1, lik.gl2.v2, tolerance = 1e-10)
})


test_that("Check agrement beteen covariance and precision cross validation", {
  V <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1),
             c(-1, 1), c(-1, 0), c(0, -1))
  E <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 5),
             c(5, 6), c(6, 1), c(4, 1), c(1, 7))

  graph <- metric_graph$new(V = V, E = E)

  kappa <- 10
  sigma <- 20
  sigma_e <- 0.1
  theta <-  c(sigma_e, sigma, kappa)

  n.obs.per.edge <- 5
  PtE <- NULL
  for(i in 1:graph$nE){
    PtE <- rbind(PtE, cbind(rep(i, n.obs.per.edge), runif(n.obs.per.edge)))
  }

  u <- sample_spde(kappa = kappa, sigma = sigma, alpha = 2,
                   graph = graph, PtE = PtE)

  y <- u + sigma_e*rnorm(n.obs.per.edge * graph$nE)

  graph$add_PtE_observations(y,PtE, normalized = TRUE)

  graph$observation_to_vertex()
  cv.alpha1.v1 <- posterior_crossvalidation(theta, graph, model = "alpha1")
  cv.alpha1.v2 <- posterior_crossvalidation_covariance(theta, graph,
                                                       model = "alpha1")
  expect_equal(cv.alpha1.v1$mu, cv.alpha1.v2$mu, tolerance = 1e-10)
  expect_equal(cv.alpha1.v1$var, cv.alpha1.v2$var, tolerance = 1e-10)

  cv.alpha2.v1 <- posterior_crossvalidation(theta, graph, model = "alpha2")
  cv.alpha2.v2 <- posterior_crossvalidation_covariance(theta, graph,
                                                       model = "alpha2")
  expect_equal(cv.alpha2.v1$mu, cv.alpha2.v2$mu, tolerance = 1e-10)
  expect_equal(cv.alpha2.v1$var, cv.alpha2.v2$var, tolerance = 1e-10)

  graph$compute_laplacian()
  cv.GL1.v1 <- posterior_crossvalidation(theta, graph, model = "GL1")
  cv.GL1.v2 <- posterior_crossvalidation_covariance(theta, graph, model = "GL1")
  expect_equal(cv.GL1.v1$mu, cv.GL1.v2$mu, tolerance = 1e-10)
  expect_equal(cv.GL1.v1$var, cv.GL1.v2$var, tolerance = 1e-10)

  cv.GL2.v1 <- posterior_crossvalidation(theta, graph, model = "GL2")
  cv.GL2.v2 <- posterior_crossvalidation_covariance(theta, graph, model = "GL2")
  expect_equal(cv.GL2.v1$mu, cv.GL2.v2$mu, tolerance = 1e-10)
  expect_equal(cv.GL2.v1$var, cv.GL2.v2$var, tolerance = 1e-10)
})
