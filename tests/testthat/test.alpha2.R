
test_that("Check agreement covariance function agrees", {
  set.seed(1)
  kappa <- 0.1*runif(1)+0.1
  sigma <- runif(1)+1
  c <- 1/(4*kappa^3)
  x <- seq(0,1,length.out=10)
  expect_equal(MetricGraph:::r_2(x, tau = 1/sigma, kappa = kappa),
               rSPDE::matern.covariance(x, kappa, 3/2, sigma)[,1]*c, tol=1e-9)
})


test_that("Check agreement derivative covariance function agrees", {
  set.seed(1)
  kappa <- 0.1 * runif(1) + 0.1
  sigma <- runif(1) + 1
  c <- 1/(4*kappa^3)
  x <- seq(-1, 1, length.out = 20)
  expect_equal(MetricGraph:::r_2(x, tau = 1/sigma, kappa = kappa, deriv = 1),
               matern_derivative(x, kappa, 3/2, sigma)[,1]*c, tol=1e-9)
})
test_that("Check agreement derivative covariance function agrees", {
  set.seed(1)
  kappa <- 0.1*runif(1)+0.1
  sigma <- runif(1)+1
  c <- 1/(4*kappa^3)
  x <- seq(-1,1,length.out=20)
  expect_equal(MetricGraph:::r_2(x, tau = 1/sigma, kappa = kappa,
                             deriv = 2),
               matern_derivative(x, kappa, 3/2, sigma,2)[,1]*c, tol=1e-9)
})

test_that("Check agreement covariance matrix", {
  set.seed(1)
  kappa <- 0.1 * runif(1) + 0.1
  sigma <- runif(1) + 1
  c <- 1/(4 * kappa^3)
  l_e <- runif(1) + 0.5
  x_ <- c(0,l_e)
  D <- outer(x_, x_, "-")
  r <- rSPDE::matern.covariance(D, kappa = kappa, nu = 3/2, sigma = sigma)
  r1 <- -matern_derivative(D, kappa = kappa, sigma = sigma, nu = 3/2, deriv = 1)
  r2 <- -matern_derivative(D, kappa = kappa, sigma = sigma, nu = 3/2, deriv = 2)
  Sigma.0 <- rbind(cbind(r, r1), cbind(t(r1), r2))*c
  r_00 <- MetricGraph:::r_2(D, tau = 1/sigma, kappa = kappa)
  r_01 <- - MetricGraph:::r_2(D, tau = 1/sigma, kappa = kappa, deriv = 1)
  r_11 <- - MetricGraph:::r_2(D, tau = 1/sigma, kappa = kappa, deriv = 2)

  Sigma_ <- rbind(cbind(r_00, r_01), cbind(t(r_01), r_11))
  testthat::expect_equal( c(Sigma.0), c(Sigma_), tol=1e-9)
})



test_that("test agreement precision matrix and article method", {
  set.seed(1)
  kappa <- 0.1* runif(1) + 0.1
  sigma <- runif(1) + 1
  c <- 1 / (4 * kappa^3)
  l_e <- runif(1) + 0.5
  x_ <- c(0, l_e)
  D <- outer(x_, x_, "-")
  r_00 <- MetricGraph:::r_2(D, tau = 1/sigma, kappa = kappa)
  r_01 <- - MetricGraph:::r_2(D, tau = 1/sigma, kappa = kappa, deriv = 1)
  r_11 <- - MetricGraph:::r_2(D, tau = 1/sigma, kappa = kappa, deriv = 2)
  # order by node not derivative
  R_00 <- matrix(c(r_00[1], r_01[1,1], r_01[1,1], r_11[1,1]),2,2)
  R_01 <- matrix(c(r_00[2], r_01[2,1], r_01[1,2], r_11[2,1]),2,2)
  R_node <- rbind(cbind(R_00, R_01), cbind(t(R_01), R_00))
  Q_adj = solve(R_node) - 0.5 * solve(rbind(cbind(R_00, matrix(0, 2, 2)),
                                            cbind(matrix(0, 2, 2), R_00)))


  build_C_beta1 <- function(L, kappa, sigma){
    C_0 <- matern_neumann_free(c(0,L), c(0,L), kappa, sigma = 1,
                               nu = 3/2, L = L,deriv = c(1,1))
    return(sigma^2 * solve(solve(C_0) - 0.5 * diag(2) / kappa^2))
  }
  C <- build_C_beta1(l_e, kappa, sigma)
  r_free.2 <- matern_neumann_free2(x_, x_,C, kappa, sigma, nu=3/2, L = l_e)
  rd1_free.2 <- matern_neumann_free2(x_, x_,C, kappa, sigma,
                                     nu=3/2, L = l_e, deriv = c(0,1))
  rd2_free.2 <- matern_neumann_free2(x_, x_,C, kappa, sigma,
                                     nu=3/2, L = l_e, deriv = c(1,1))
  Sigma.2  <- rbind(cbind(r_free.2, rd1_free.2),
                    cbind(t(rd1_free.2),rd2_free.2)) * c
  Sigma.2 <- Sigma.2[c(1, 3, 2, 4), c(1, 3, 2, 4)] #order by nodes

  testthat::expect_equal(c(Q_adj), c(solve(Sigma.2)), tol = 1e-7)

  #test adjusted covariance matrix against article formula
  R00R0l <-rbind(cbind(R_00, R_01), cbind(t(R_01), R_00))
  Adj <- solve(rbind(cbind(R_00, -R_01), cbind(-t(R_01), R_00)))
  R_adj <- R_node + R00R0l %*% Adj %*% R00R0l
  testthat::expect_equal(c(Sigma.2), c(R_adj), tol = 1e-9)
})
# Commented because there is no test
# test_that("test agreement neumann of the adjusted", {
#   set.seed(1)
#   kappa <- 0.1* runif(1) + 1
#   sigma <- runif(1) + 1
#   c <- 1 / (4 * kappa^3)
#   l_e <- runif(1) + 0.5
#   x_ <- c(0, l_e)
#   eps <- 1e-4

#   D <- outer(x_, x_, "-")
#   r_00 <- MetricGraph:::r_2(D, tau = 1/sigma, kappa = kappa)
#   r_01 <- - MetricGraph:::r_2(D, tau = 1/sigma, kappa = kappa, deriv = 1)
#   r_11 <- - MetricGraph:::r_2(D, tau = 1/sigma, kappa = kappa, deriv = 2)
#   # order by node not derivative
#   R_00 <- matrix(c(r_00[1], r_01[1,1], r_01[1,1], r_11[1,1]),2,2)
#   R_01 <- matrix(c(r_00[2], r_01[2,1], r_01[1,2], r_11[2,1]),2,2)
#   R_node <- rbind(cbind(R_00, R_01), cbind(t(R_01), R_00))
#   Q_adj = solve(R_node) - 0.5 * solve(rbind(cbind(R_00, matrix(0, 2, 2)),
#                                             cbind(matrix(0, 2, 2), R_00)))

#   Adj <- solve(rbind(cbind(R_00, -R_01), cbind(-t(R_01), R_00)))
#   R00R0l <-rbind(cbind(R_00, R_01), cbind(t(R_01), R_00))
#   R_adj <- R_node + R00R0l %*% Adj %*% R00R0l
#   t <- sum(x_)/2
#   D2 <- (outer(c(x_,t), c(x_,t), "-"))
#   r_00_ <- MetricGraph:::r_2(D2, tau = 1/sigma, kappa = kappa)
#   r_01_ <- - MetricGraph:::r_2(D2, tau = 1/sigma, kappa = kappa, deriv = 1)
#   r_11_ <- - MetricGraph:::r_2(D2, tau = 1/sigma, kappa = kappa, deriv = 2)
#   R_t0 <- matrix(c(r_00_[1,3], r_01_[1,3], r_01_[3,1], r_11_[3,1]),2,2)
#   R_t1 <- matrix(c(r_00_[2,3], r_01_[2,3], r_01_[3,2], r_11_[3,2]),2,2)
#   R_corr <- cbind(R_t0,R_t1)
#   R_t <- R_t0 + cbind(R_00, R_01)%*%Adj %*% t(R_corr)
#   x_[0] <- 0
#   t_eps <- eps
#   D2 <- (outer(c(x_,t_eps), c(x_,t_eps), "-"))
#   r_00_ <- MetricGraph:::r_2(D2, tau = 1/sigma, kappa = kappa)
#   r_01_ <- - MetricGraph:::r_2(D2, tau = 1/sigma, kappa = kappa, deriv = 1)
#   r_11_ <- - MetricGraph:::r_2(D2, tau = 1/sigma, kappa = kappa, deriv = 2)
#   R_eps0 <- matrix(c(r_00_[1,3], r_01_[1,3], r_01_[3,1], r_11_[3,1]),2,2)
#   R_eps1 <- matrix(c(r_00_[2,3], r_01_[2,3], r_01_[3,2], r_11_[3,2]),2,2)
#   R_corr_eps <- cbind(R_eps0,R_eps1)
#   R_00_e <- matrix(c(r_00_[1,1], r_01[1,1], r_01[1,1], r_11[1,1]),2,2)
#   R_01_e <- matrix(c(r_00_[2,1], r_01[2,1], r_01[1,2], r_11[2,1]),2,2)

#   R_t2 <- R_t0 + cbind(R_00_e, R_01_e)%*%Adj %*% t(R_corr)
#   (R_t-R_adj[1:2,1:2])/eps
# })
test_that("test likelihood",{
  set.seed(13)
  nt <- 40
  kappa <- 0.3
  sigma_e <- 0.1
  sigma   <- 1
  theta <-  c(sigma_e,sigma,kappa)
  edge2 <- rbind(c(30, 80), c(140, 80))
  edge1 <- rbind(c(30, 00), c(30, 80))
  edges <- list(edge1, edge2)
  graph <- metric_graph$new(edges = edges)
  Q <- spde_precision(kappa = kappa, tau = 1/sigma,
                      alpha = 2, graph = graph, BC = 1)
  graph$buildC(2, FALSE)
  Qmod <- (graph$CoB$T) %*% Q %*% t(graph$CoB$T)
  Qtilde <- Qmod
  Qtilde <- Qtilde[-c(1:2),-c(1:2)]
  R <- Cholesky(Qtilde,LDL = FALSE, perm = TRUE)
  V0 <- as.vector(solve(R, solve(R,rnorm(6), system = 'Lt')
                        , system = 'Pt'))
  u_e <- t(graph$CoB$T) %*% c(0, 0, V0)
  X <- c()
  for(i in 1:length(graph$edge_lengths)){
    X <- rbind(X,cbind(sample_alpha2_line(kappa = kappa,
                                          tau = 1/sigma,
                                          sigma_e = sigma_e,
                                          u_e = u_e[4*(i-1) +1:4],
                                          l_e = graph$edge_lengths[i],
                                          nt = nt),i))
  }
  X[,2] <- X[,2] + sigma_e*rnorm(2*nt)
  
  df_test <- data.frame(y = X[,2], edge_number = X[,3], distance_on_edge = X[,1])

  graph$add_observations(data = df_test, normalized = FALSE)
  graph$buildC(2, FALSE)

  #standard likelihood
  lik <- -likelihood_alpha2(theta = theta, graph = graph, data_name = "y", 
                             X_cov = NULL, repl = NULL, BC = 1, parameterization = "spde")
  
  graph2 <- graph$clone()
  graph2$observation_to_vertex()
  graph2$buildC(2, FALSE)

  #covariance likelihood

  lik2 <-likelihood_graph_covariance(graph = graph2,
                                     model = "WM2", repl = NULL, y_graph = graph2$get_data()[["y"]],
                                     log_scale = FALSE, X_cov = NULL)
  lik2 <- lik2(exp(theta))
  #likelihood with extended graph
  lik3 <- -likelihood_alpha2(theta = theta, graph = graph2, data_name = "y", 
                             X_cov = NULL, repl = NULL, BC = 1, parameterization = "spde")

  expect_equal(as.matrix(lik), as.matrix(lik2), tolerance = 1e-10)
  expect_equal(as.matrix(lik), as.matrix(lik3), tolerance = 1e-10)
})

test_that("test posterior mean",{
  library(Matrix)
  set.seed(13)
  nt <- 90
  kappa <- 0.3
  sigma_e <- 0.1
  sigma   <- 2
  theta <-  c(sigma_e,sigma,kappa)
  edge2 <- rbind(c(30, 80), c(140, 80))
  edge1 <- rbind(c(30, 00), c(30, 80))
  edges <- list(edge1, edge2)
  graph <- metric_graph$new(edges = edges)
  Q <- spde_precision(kappa = kappa, tau = 1/sigma,
                      alpha = 2, graph = graph, BC = 1)
  graph$buildC(2, FALSE)
  Qmod <- (graph$CoB$T) %*% Q %*% t(graph$CoB$T)
  Qtilde <- Qmod
  Qtilde <- Qtilde[-c(1:2),-c(1:2)]
  R <- Cholesky(Qtilde,LDL = FALSE, perm = TRUE)
  V0 <- as.vector(solve(R, solve(R,rnorm(6), system = 'Lt')
                        , system = 'Pt'))
  u_e <- t(graph$CoB$T) %*% c(0, 0, V0)
  X <- c()
  for(i in 1:length(graph$edge_lengths)){
    X <- rbind(X,cbind(sample_alpha2_line(kappa = kappa,
                                          tau = 1/sigma,
                                          sigma_e = sigma_e,
                                          u_e = u_e[4*(i-1) +1:4],
                                          l_e = graph$edge_lengths[i],
                                          nt = nt),i))
  }
  X <- X[-nt,] # There is a repeated location, let us remove it
  X[,2] <- X[,2] + rnorm(nrow(X), sd = sigma_e)

  df_temp <- data.frame(y = X[,2], edge_number = X[,3], distance_on_edge = X[,1])

  graph$add_observations(data = df_temp, normalized = FALSE)

  #test posterior at observation locations
  res <- graph_lme(y ~ -1, graph=graph, model="WM2", parallel = FALSE)
  pm <- predict(res, newdata = df_temp)$mean

  kappa_est <- res$coeff$random_effects[2]
  tau_est <- res$coeff$random_effects[1]
  theta_est <- c(res$coeff$measurement_error, tau_est, kappa_est)

  graph2 <- graph$clone()
  graph2$observation_to_vertex()
  graph2$buildC(2, FALSE)
  n.v <- dim(graph2$V)[1]
  n.c <- 1:length(graph2$CoB$S)
  Q <- spde_precision(kappa = kappa_est, tau = tau_est,
                      alpha = 2, graph = graph2, BC = 1)
  Qtilde <- (graph2$CoB$T) %*% Q %*% t(graph2$CoB$T)
  Qtilde <- Qtilde[-n.c,-n.c]
  Sigma.overdetermined  = t(graph2$CoB$T[-n.c, ]) %*%
    solve(Qtilde) %*% (graph2$CoB$T[-n.c, ])
  PtE <- graph2$get_PtE()
  index.obs <- 4*(PtE[, 1] - 1) +
    (1 * (abs(PtE[, 2]) < 1e-14)) +
    (3 * (abs(PtE[, 2]) > 1e-14))
  Sigma <-  as.matrix(Sigma.overdetermined[index.obs, index.obs])
  Sigma.Y <- Sigma
  diag(Sigma.Y) <- diag(Sigma.Y) + theta_est[1]^2
  pm2 <- Sigma %*% solve(Sigma.Y, graph2$get_data()$y)

  pm2 <- as.vector(pm2)

  ord1 <- order(graph$get_data()[[".coord_x"]], graph$get_data()[[".coord_y"]])
  ord2 <- order(graph2$get_data()[[".coord_x"]], graph2$get_data()[[".coord_y"]])

  expect_equal(sum((pm2[ord2]-pm[ord1])^2),0, tolerance = 1e-5)
})




