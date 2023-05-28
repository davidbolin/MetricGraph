#build graph Figure covariance between points


test_that("Check if Matern(2) has equal variance on the edges on a circle", {
  theta <- c(1,2)
  P <- rbind(c(0,0),
             c(1,0),
             c(1,1),
             c(0,1))

  #specify edges
  E <- rbind(c(1,2),
             c(2,3),
             c(1,4),
             c(3,4))

  graph <-  metric_graph$new(V = P, E = E)


  #if Line is null assume straight line
  #If we have a line keep
  #compute covarians for a given point
  P1 <- c(1, 0.1) #edge 1, 0.5 length in


  kappa <- theta[1]
  sigma <- theta[2]
  #compute covarains of the two edges of EP[1]
  Q <- spde_precision(kappa = kappa, tau = 1/sigma, alpha = 2, graph = graph)
  if(is.null(graph$CoB))
    graph$buildC(2, FALSE)

  n_const <- length(graph$CoB$S)
  ind.const <- c(1:n_const)
  Tc <- graph$CoB$T[-ind.const,]
  Q_mod <- Tc%*%Q%*%t(Tc)
  R <- Cholesky(Q_mod, LDL = FALSE, perm = TRUE)
  Sigma_const <- t(Tc)%*%solve(Q_mod)%*%Tc
  Sigma_unit <- diag(Sigma_const[seq(1,dim(Q)[1],by=2),seq(1,dim(Q)[1],by=2)])

  expect_equal(c(Sigma_unit - Sigma_unit[1]),rep(0,length(Sigma_unit)), tolerance=1e-10)
})
