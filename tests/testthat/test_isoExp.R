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
                 2,2,2,1),c(0.2, 0.1,
                          0.9, 0.2, 0.4,0))

  reo <- order(PtE[, 1],PtE[, 2])
  v <- PtE[,2]
  v[PtE[,1]==2] <- v[PtE[,1]==2] + 1
  y <- 1:dim(PtE)[1]

  D <- graph$compute_resdist_PtE(PtE)
  distances <- as.matrix(dist(v))
  dimnames(distances) <- NULL
  expect_equal(as.vector(D), as.vector(distances), tolerance = 1e-10)
  
    df_test <- data.frame(y = y, edge_number = PtE[,1], distance_on_edge = PtE[,2])

  graph$add_observations(data = df_test, normalized = TRUE)
  graph$compute_resdist()
  D1 <- graph$res_dist[[1]]
  expect_equal(as.vector(D1), as.vector(distances[reo,reo]), tolerance = 1e-10)

  theta <- c(1, 2, 3)
  lik1 <- likelihood_graph_covariance(graph, model = "isoCov", y_graph = graph$get_data()$y, repl=NULL, cov_function = exp_covariance, log_scale=FALSE, maximize = TRUE)
  lik1 <- lik1(theta)

  Sigma <- theta[2]^2 * exp(-theta[3]*distances)
  diag(Sigma) <- diag(Sigma) + theta[1]
  R <- chol(Sigma)
  lik.true <- as.vector(-sum(log(diag(R))) - (dim(PtE)[1]/2)* log(2*pi) -
    0.5*t(y)%*%solve(Sigma, y))


  graph$observation_to_vertex()
  lik2 <- likelihood_graph_covariance(graph, model = "isoCov", cov_function = exp_covariance, y_graph = graph$get_data()$y, repl=NULL, log_scale = FALSE, maximize = TRUE)
  lik2 <- lik2(theta)

  PtE.order <- PtE[reo, ]
  y.order <- y[reo]
  graph <- metric_graph$new(lines = Lines)
  df_temp <- data.frame(y = y.order, edge_number = PtE.order[,1], distance_on_edge = PtE.order[,2])
  graph$add_observations(data=df_temp, normalized = TRUE)
  graph$compute_resdist()
  lik3 <- likelihood_graph_covariance(graph,  model = "isoCov", repl = NULL, y = graph$get_data()$y, cov_function = exp_covariance, log_scale = FALSE, maximize = TRUE)
  lik3 <- lik3(theta)

  graph$observation_to_vertex()
  lik4 <- likelihood_graph_covariance(graph, model = "isoCov", repl = NULL, y_graph = graph$get_data()$y, cov_function = exp_covariance, log_scale = FALSE, maximize = TRUE)
  lik4 <- lik4(theta)

  expect_equal(lik.true, lik1, tolerance = 1e-10)
  expect_equal(lik.true, lik2, tolerance = 1e-10)
  expect_equal(lik.true, lik3, tolerance = 1e-10)
  expect_equal(lik.true, lik4, tolerance = 1e-10)


})
