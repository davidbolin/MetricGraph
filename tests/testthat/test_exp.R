# verify claim that one can solve circular edges using same expression as for
# non circular, see likelihood inference in article.

#' test agreement between Q and R
test_that("Check agrement beteen covariance and precision matrix formulation", {
  kappa <- 1
  tau <- 1
  t <- 0:3
  Q0  <- precision_exp_line(kappa = kappa, tau = tau, t = t)
  R0  <- r_1(as.matrix(dist(t)),kappa = kappa, tau = tau)
  R0_ <- solve(Q0)
  expect_equal(c(as.matrix(solve(R0))), c(as.matrix(Q0)), tolerance = 1e-10)
  kappa = 1.5
  sigma = 0.5
  t <- 0:3
  Q0  <- precision_exp_line(kappa = kappa, tau = tau, t = t)
  R0  <- r_1(as.matrix(dist(t)), kappa = kappa, tau = tau)
  R0_ <- solve(Q0)
  expect_equal(c(as.matrix(solve(R0))), c(as.matrix(Q0)), tolerance = 1e-10)
})



test_that("Check agrement beteen covariance and precision likelihoods", {
  set.seed(1)
  nt <- 10
  kappa <- 1
  sigma_e <- 0.1
  tau   <- 1/2

  #line1 <- Line(rbind(c(30, 80), c(120, 80)))
  #line2 <- Line(rbind(c(30, 00), c(30, 80)))
  edge1 <- rbind(c(0, 0), c(1e-0, 0))
  edge2 <- rbind(c(0, 1e-0), c(0, 0))

  graph <-  metric_graph$new(list(edge1, edge2))

  n.obs.per.edge <- 10
  PtE <- NULL
  for(i in 1:graph$nE){
    PtE <- rbind(PtE, cbind(rep(i, n.obs.per.edge), (runif(n.obs.per.edge))))
  }


  nt <- graph$nE* n.obs.per.edge
  PtE <- PtE[sample(1:nt),]
  u <- sample_spde(kappa = kappa, tau = tau,
                   alpha = 1, graph = graph, PtE = PtE)

  y <- u + sigma_e*rnorm(nt)

  df_temp <- data.frame(y = y, edge_number = PtE[,1], distance_on_edge = PtE[,2])
  graph$add_observations(data=df_temp, normalized = TRUE)
  theta <-  c(sigma_e, kappa, 1/tau)
  lik <- likelihood_alpha1(theta = theta, graph = graph, data_name = "y", 
                             X_cov = NULL , repl = NULL, BC = 1, parameterization = "spde")

  graph$observation_to_vertex()
  lik.v2 <- likelihood_alpha1_v2(theta = theta, graph = graph, 
              X_cov = matrix(ncol=0,nrow=0), y = graph$get_data()$y, repl = NULL, BC = 1, 
              parameterization = "spde")

  lik.cov <- likelihood_graph_covariance(graph, model = "WM1", log_scale = TRUE, y_graph = graph$get_data()$y, repl = NULL, X_cov = NULL, maximize = TRUE)
  lik.cov <- lik.cov(theta)

  #version 1
  expect_equal(as.matrix(lik.v2),as.matrix(lik.cov), tolerance=1e-10)

  #version 2
  expect_equal(as.matrix(lik),as.matrix(lik.cov), tolerance=1e-10)
})

test_that("Test posterior mean", {
  set.seed(1)
  nt <- 100
  range <- 20
  kappa <- sqrt(4)/range
  sigma_e <- 0.2
  tau   <- 0.5
  line1 <- sp::Line(rbind(c(30, 80), c(120, 80)))
  line2 <- sp::Line(rbind(c(30, 00), c(30, 80)))

  graph <-  metric_graph$new(sp::SpatialLines(list(sp::Lines(list(line1), ID = "1"),
                                                   sp::Lines(list(line2), ID = "2")
  )))
  PtE <- rbind(cbind(rep(1,nt/2),
                     seq(from = 0,to =1, length.out = nt/2 + 1)[1:(nt/2)]),
               cbind(rep(2,nt/2),
                     seq(from = 0,to =1, length.out = nt/2 + 1)[1:(nt/2)]))

  u <- sample_spde(kappa = kappa, tau = tau,
                   alpha = 1, graph = graph, PtE = PtE)


  y <- u + sigma_e*rnorm(nt)
  df_temp <- data.frame(y = y, edge_number = PtE[,1], distance_on_edge = PtE[,2])
  graph$add_observations(data = df_temp, normalized = TRUE)

  #test posterior at observation locations
  res <- suppressWarnings(graph_lme(y ~ -1, graph=graph, model="WM1", parallel = FALSE))
  pm <- predict(res, newdata = df_temp, normalized=TRUE)$mean

  kappa_est <- res$coeff$random_effects[2]
  tau_est <- res$coeff$random_effects[1]
  theta_est <- c(res$coeff$measurement_error, tau_est, kappa_est)

  graph$observation_to_vertex()

  Q <- spde_precision(kappa = kappa_est, tau = tau_est, alpha = 1, graph = graph)
  Sigma <- solve(Q)[graph$PtV, graph$PtV]
  Sigma.obs <- Sigma
  diag(Sigma.obs) <- diag(Sigma.obs) + theta_est[1]^2
  pm2 <- Sigma %*% solve(Sigma.obs, graph$get_data()$y)

  expect_equal(sum((sort(pm)-sort(pm2))^2), 0, tolerance=1e-8)

  expect_equal(sum((pm - df_temp$y)^2), sum((pm2 - graph$get_data()$y)^2), tolerance = 1e-8)

})
