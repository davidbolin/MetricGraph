#check resistance metric and reordering
test_that("Test resistance metric", {
  #check resistance metric and reordering
  library(sp)
  line1 <- Line(rbind(c(0,0),c(1,0)))
  line2 <- Line(rbind(c(1,0),c(2,0)))
  Lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                                Lines(list(line2),ID="2")))
  graph <- metric_graph$new(Lines = Lines)


  PtE <- cbind(c(1,1,
                 2,2,2),c(0.1, 0.2,
                          0.2, 0.9, 0.4))
  v <- PtE[,2]
  v[PtE[,1]==2] <- v[PtE[,1]==2] + 1
  y <- 1:dim(PtE)[1]

  D <- graph$compute_resdist(PtE)
  distances <- as.matrix(dist(v))
  dimnames(distances) <- NULL
  expect_equal(as.vector(D), as.vector(distances), tolerance = 1e-10)

  graph$add_observations2(y,PtE, normalized = TRUE)
  graph$compute_resdist()

  expect_equal(as.vector(graph$res.dist), as.vector(distances), tolerance = 1e-10)
  theta <- c(1, 2, 3)
  lik1 <- likelihood_graph_covariance(theta, graph, model = "isoExp")

  reo <- order(PtE[, 1],PtE[, 2])
  PtE.order <- PtE[reo, ]
  y.order <- y[reo]
  graph <- metric_graph$new(Lines = Lines)
  graph$add_observations2(y.order,PtE.order, normalized = TRUE)
  graph$compute_resdist()
  lik2 <- likelihood_graph_covariance(theta, graph, model = "isoExp")

  lik1 - lik2
  expect_equal(lik1, lik2, tolerance = 1e-10)


})
