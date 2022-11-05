#verify claim that one can solve circular edges using same expression as for non circular
# see likelihood inference in article

#' test agreement between Q and R
test_that("Check agrement beteen covariance and precision matrix formulation", {
  theta <- c(1, 1)
  t <- 0:3
  Q0  <- Q.exp.line(theta, t)
  R0  <- r_1(as.matrix(dist(t)),theta)
  R0_ <- solve(Q0)
  expect_equal(c(as.matrix(solve(R0))),c(as.matrix(Q0)), tolerance=1e-10)
  theta <- c(1.5, 0.5)
  t <- 0:3
  Q0  <- Q.exp.line(theta, t)
  R0  <- r_1(as.matrix(dist(t)),theta)
  R0_ <- solve(Q0)
  expect_equal(c(as.matrix(solve(R0))),c(as.matrix(Q0)), tolerance=1e-10)
})



test_that("Check agrement beteen covariance and precision matrix formulation", {
  nt <- 10
  kappa <- 0.1
  sigma_e <- 0.1
  sigma   <- 2
  theta <-  c(sigma_e,kappa,sigma)
  line.line <- Line(rbind(c(30,80),c(120,80)))
  line.line2 <- Line(rbind(c(30,00),c(30,80)))

  graph <-  gpgraph_graph$new(sp::SpatialLines(list(Lines(list(line.line),ID="1"),
                                                Lines(list(line.line2),ID="2"))))
  Q <- Q.exp(theta[2:3], graph$V, graph$EtV, graph$edge_lengths)
  R <- Cholesky(Q,LDL = FALSE, perm = TRUE)
  V0 <- as.vector(solve(R, solve(R,rnorm(3), system = 'Lt')
                        , system = 'Pt'))
  X <- c()
  for(i in 1:length(graph$edge_lengths)){
    X <- rbind(X,cbind(sample.line.expontial(theta,
                                             V0[graph$EtV[i,2:3]],
                                             Line = graph$Lines[i,],
                                             graph$edge_lengths[i],
                                             nt = nt),i))
  }
  X[,2] <- X[,2] + sigma_e*rnorm(nt)
  graph$add_observations2(y = X[,2], PtE = X[,c(3,1)])
  graph$observation_to_vertex()
  lik <- likelihood.exp.graph(theta,graph)
  #lik_diff <- lik - likelihood.exp.graph(theta+c(1,0,0),graph)


  lik.v2 <- likelihood.exp.graph.v2(theta,graph)
  #lik_diff.v2 <- lik.v2 - likelihood.exp.graph.v2(theta+c(1,0,0),graph)

  lik.cov <-  likelihood.graph.covariance(theta,graph, model="alpha1")

  #version 1
  expect_equal(as.matrix(lik.v2),as.matrix(lik.cov), tolerance=1e-10)
  #version 2
  expect_equal(as.matrix(lik),as.matrix(lik.cov), tolerance=1e-10)
})
