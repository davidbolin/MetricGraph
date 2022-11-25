#check resistance metric and reordering
test_that("Test resistance metric", {
  #check resistance metric and reordering
  library(sp)
  line1 <- Line(rbind(c(0,0),c(1,0)))
  line2 <- Line(rbind(c(1,0),c(2,0)))
  Lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                                Lines(list(line2),ID="2")))
  graph <- metric_graph$new(lines = Lines)


  PtE <- cbind(c(1,1,
                 2,2,2),c(0.2, 0.1,
                          0.9, 0.2, 0.4))

  reo <- order(PtE[, 1],PtE[, 2])
  v <- PtE[,2]
  v[PtE[,1]==2] <- v[PtE[,1]==2] + 1
  y <- 1:dim(PtE)[1]

  D <- graph$compute_resdist(PtE)
  distances <- as.matrix(dist(v))
  dimnames(distances) <- NULL
  expect_equal(as.vector(D), as.vector(distances), tolerance = 1e-10)

  graph$add_PtE_observations(y,PtE, normalized = TRUE)
  graph$compute_resdist()
  D1 <- graph$res_dist
  expect_equal(as.vector(D1), as.vector(distances), tolerance = 1e-10)

  theta <- c(1, 2, 3)
  lik1 <- likelihood_graph_covariance(theta, graph, model = "isoExp")

  Sigma <- theta[2]^2 * exp(-theta[3]*distances)
  diag(Sigma) <- diag(Sigma) + theta[1]
  R <- chol(Sigma)
  lik.true <- as.vector(-sum(log(diag(R))) - (dim(PtE)[1]/2)* log(2*pi) -
    0.5*t(y)%*%solve(Sigma, y))


  graph$observation_to_vertex()
  lik2 <- likelihood_graph_covariance(theta, graph, model = "isoExp")

  PtE.order <- PtE[reo, ]
  y.order <- y[reo]
  graph <- metric_graph$new(lines = Lines)
  graph$add_PtE_observations(y.order,PtE.order, normalized = TRUE)
  graph$compute_resdist()
  lik3 <- likelihood_graph_covariance(theta, graph, model = "isoExp")

  graph$observation_to_vertex()
  lik4 <- likelihood_graph_covariance(theta, graph, model = "isoExp")

  expect_equal(lik.true, lik1, tolerance = 1e-10)
  expect_equal(lik.true, lik2, tolerance = 1e-10)
  expect_equal(lik.true, lik3, tolerance = 1e-10)
  expect_equal(lik.true, lik4, tolerance = 1e-10)


})
