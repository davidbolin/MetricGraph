
test_that("test if the over determined covariances are correct",{
  kappa <- 0.07
  sigma   <- 0.1
  edge1 <- rbind(c(120,80),c(30,80))
  edge2 <- rbind(c(30,80),c(30,00))

  graph <-  metric_graph$new(list(edge1,edge2))
  Q <- spde_precision(kappa = kappa, tau = 1/sigma, alpha = 2, graph = graph, BC = 1)
  graph$buildC(2, FALSE)
  Qtilde <- (graph$CoB$T)%*%Q%*%t(graph$CoB$T)
  Qtilde <- Qtilde[-c(1:2),-c(1:2)]
  Sigma.overdetermined  = t(graph$CoB$T[-c(1:2),])%*%solve(Qtilde)%*%(graph$CoB$T[-c(1:2),])

  t <- c(0,graph$edge_lengths[1],graph$edge_lengths[1],graph$edge_lengths[2]+ graph$edge_lengths[1])

  D <-  rep(1,4)%*%t(t)-t(rep(1,4)%*%t(t))
  R00 <- r_2(D, tau = 1/sigma, kappa = kappa)
  R01 <- -r_2(D,tau = 1/sigma, kappa = kappa, deriv = 1)
  R11 <- -r_2(D,tau = 1/sigma, kappa = kappa, deriv = 2)
  R <- cbind(rbind(R00, R01), rbind(t(R01), R11))
  ind <- rep(0:3, each = 2) + rep(c(1, 5), times = 4) #reoder
  R <- R[ind, ind]
  testthat::expect_equal( c(as.matrix(R)), c(as.matrix(Sigma.overdetermined)), tol=1e-9)
})
