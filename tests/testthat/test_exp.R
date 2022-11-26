# verify claim that one can solve circular edges using same expression as for
# non circular, see likelihood inference in article.

#' test agreement between Q and R
test_that("Check agrement beteen covariance and precision matrix formulation", {
  kappa <- 1
  sigma <- 1
  t <- 0:3
  Q0  <- precision_exp_line(kappa = kappa, sigma = sigma, t = t)
  R0  <- r_1(as.matrix(dist(t)),kappa = kappa, sigma = sigma)
  R0_ <- solve(Q0)
  expect_equal(c(as.matrix(solve(R0))), c(as.matrix(Q0)), tolerance = 1e-10)
  kappa = 1.5
  sigma = 0.5
  t <- 0:3
  Q0  <- precision_exp_line(kappa = kappa, sigma = sigma, t = t)
  R0  <- r_1(as.matrix(dist(t)), kappa = kappa, sigma = sigma)
  R0_ <- solve(Q0)
  expect_equal(c(as.matrix(solve(R0))), c(as.matrix(Q0)), tolerance = 1e-10)
})



test_that("Check agrement beteen covariance and precision likelihoods", {
  nt <- 10
  kappa <- 1
  sigma_e <- 0.1
  sigma   <- 2
  #line1 <- Line(rbind(c(30, 80), c(120, 80)))
  #line2 <- Line(rbind(c(30, 00), c(30, 80)))
  line1 <- Line(rbind(c(0, 0), c(1e-0, 0)))
  line2 <- Line(rbind(c(0, 1e-0), c(0, 0)))

  graph <-  metric_graph$new(sp::SpatialLines(list(Lines(list(line1), ID = "1"),
                                                   Lines(list(line2), ID = "2")
                                                   )))

  n.obs.per.edge <- 10
  PtE <- NULL
  for(i in 1:graph$nE){
    PtE <- rbind(PtE, cbind(rep(i, n.obs.per.edge), (runif(n.obs.per.edge))))
  }


  nt <- graph$nE* n.obs.per.edge
  PtE <- PtE[sample(1:nt),]
  u <- sample_spde(kappa = kappa, sigma = sigma,
                   alpha = 1, graph = graph, PtE = PtE)

  y <- u + sigma_e*rnorm(nt)
  graph$add_PtE_observations(y = y, PtE = PtE, normalized = TRUE)
  theta <-  c(sigma_e, kappa, sigma)
  lik <- likelihood_graph_spde(graph, alpha = 1, version = 1, log_scale = FALSE)
  lik <- lik(theta)

  graph$observation_to_vertex()
  lik.v2 <- likelihood_graph_spde(graph, alpha = 1, version = 2, log_scale = FALSE)
  lik.v2 <- lik.v2(theta)

  lik.cov <- likelihood_graph_covariance(graph, model = "alpha1", log_scale = FALSE)
  lik.cov <- lik.cov(theta)

  #version 1
  expect_equal(as.matrix(lik.v2),as.matrix(lik.cov), tolerance=1e-10)

  #version 2
  expect_equal(as.matrix(lik),as.matrix(lik.cov), tolerance=1e-10)
})

test_that("Test posterior mean", {
  nt <- 10
  kappa <- 0.1
  sigma_e <- 0.2
  sigma   <- 2
  line1 <- Line(rbind(c(30, 80), c(120, 80)))
  line2 <- Line(rbind(c(30, 00), c(30, 80)))

  graph <-  metric_graph$new(sp::SpatialLines(list(Lines(list(line1), ID = "1"),
                                                   Lines(list(line2), ID = "2")
  )))
  PtE <- rbind(cbind(rep(1,nt/2),
                     seq(from = 0,to =1, length.out = nt/2 + 1)[1:(nt/2)]),
               cbind(rep(2,nt/2),
                     seq(from = 0,to =1, length.out = nt/2 + 1)[1:(nt/2)]))

  u <- sample_spde(kappa = kappa, sigma = sigma,
                   alpha = 1, graph = graph, PtE = PtE)


  y <- u + sigma_e*rnorm(nt)
  graph$add_PtE_observations(y = y, PtE = PtE)
  theta <-  c(sigma_e, sigma, kappa)
  pm <- spde_posterior_mean(theta, alpha = 1, type = "obs", graph = graph)

  graph$observation_to_vertex()

  Q <- spde_precision(kappa = kappa, sigma = sigma, alpha = 1, graph = graph)
  Sigma <- solve(Q)[graph$PtV, graph$PtV]
  Sigma.obs <- Sigma
  diag(Sigma.obs) <- diag(Sigma.obs) + sigma_e^2
  pm2 <- Sigma %*% solve(Sigma.obs, graph$y)

  expect_equal(as.vector(pm),as.vector(pm2), tolerance=1e-10)
})
