test_that("Test graph construction from lines", {

  edge1 <-  rbind(c(0,0),c(0,1))
  edge2 <- rbind(c(0,1),c(1,1))

  graph <-  metric_graph$new(edges = list(edge1, edge2))

  expect_equal(graph$V, matrix(c(0, 0, 1, 0, 1, 1),3, 2), tol = 1e-9)
  expect_equal(graph$E, matrix(c(1, 2, 2, 3), 2, 2), tol = 1e-9)
})

test_that("Test adding point", {


  edge1 <-  rbind(c(0,0),c(0,1))
  edge2 <- rbind(c(0,1),c(1,1))

  graph <-  metric_graph$new(edges = list(edge1, edge2))

  P1 <- matrix(c(1,0.5), nrow=1,ncol=2)
  p <- sp::SpatialPoints(matrix(c(0,0.5),1,2))
  
  graph$add_observations(Spoints = p, data = data.frame(y=1))
  graph$observation_to_vertex()

  expect_equal(graph$V, matrix(c(0, 0, 1, 0, 0, 1, 1, 0.5), 4, 2), tol=1e-9)
  expect_equal(graph$E, matrix(c(1, 2, 4, 4, 3, 2), 3, 2), tol=1e-9)
})


