#verify claim that one can solve circular edges using same expression as for non circular
# see likelihood inference in article

#' test agreement between Q and R
test_that("Check agrement beteen covariance and precision matrix formulation", {
  kappa <- 1
  sigma <- 1
  t <- 0:3
  Q0  <- precision_exp_line(kappa = kappa, sigma = sigma, t = t)
  R0  <- r_1(as.matrix(dist(t)),kappa = kappa, sigma = sigma)
  R0_ <- solve(Q0)
  expect_equal(c(as.matrix(solve(R0))), c(as.matrix(Q0)), tolerance=1e-10)
  kappa = 1.5
  sigma = 0.5
  t <- 0:3
  Q0  <- precision_exp_line(kappa = kappa, sigma = sigma, t = t)
  R0  <- r_1(as.matrix(dist(t)),kappa = kappa, sigma = sigma)
  R0_ <- solve(Q0)
  expect_equal(c(as.matrix(solve(R0))),c(as.matrix(Q0)), tolerance=1e-10)
})



test_that("Check agrement beteen covariance and precision matrix formulation", {
  nt <- 10
  kappa <- 0.1
  sigma_e <- 0.1
  sigma   <- 2
  line1 <- Line(rbind(c(30, 80), c(120, 80)))
  line2 <- Line(rbind(c(30, 00), c(30, 80)))

  graph <-  metric_graph$new(sp::SpatialLines(list(Lines(list(line1),ID="1"),
                                                   Lines(list(line2),ID="2"))))
  Q <- spde_precision(kappa = kappa, sigma = sigma, alpha = 1, graph = graph)
  R <- Cholesky(Q,LDL = FALSE, perm = TRUE)
  V0 <- as.vector(solve(R, solve(R,rnorm(3), system = 'Lt')
                        , system = 'Pt'))
  X <- c()
  for(i in 1:length(graph$edge_lengths)){
    X <- rbind(X,cbind(sample_alpha1_line(kappa = kappa,
                                       sigma = sigma,
                                       u_e = V0[graph$E[i, ]],
                                       l_e = graph$edge_lengths[i],
                                       nt = nt), i))
  }
  X[,2] <- X[,2] + sigma_e*rnorm(nt)
  graph$add_observations2(y = X[,2], PtE = X[, c(3, 1)])
  graph$observation_to_vertex()

  theta <-  c(sigma_e, kappa, sigma)
  lik <- likelihood_graph_spde(theta,graph, alpha = 1, version = 1)
  lik.v2 <- likelihood_graph_spde(theta, graph, alpha = 1, version = 2)
  lik.cov <- likelihood_graph_covariance(theta, graph, model = "alpha1")

  #version 1
  expect_equal(as.matrix(lik.v2),as.matrix(lik.cov), tolerance=1e-10)
  #version 2
  expect_equal(as.matrix(lik),as.matrix(lik.cov), tolerance=1e-10)
})
