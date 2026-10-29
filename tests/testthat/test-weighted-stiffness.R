test_that("weighted stiffness matches the exact one-edge integral", {
  G <- weighted_stiffness(matrix(c(1L, 2L), nrow = 1L),
                          lengths = 1, weights = c(1, 4))
  expected <- 2.5 * matrix(c(1, -1, -1, 1), 2, 2)
  expect_equal(as.matrix(G), expected)
  expect_equal(as.numeric(G %*% c(1, 1)), c(0, 0))
  expect_equal(as.numeric(c(0, 1) %*% G %*% c(0, 1)), 2.5)
})

test_that("weighted stiffness remains positive semidefinite on a mesh", {
  edges <- rbind(c(1L, 2L), c(2L, 3L), c(3L, 4L))
  G <- weighted_stiffness(edges, c(0.2, 0.4, 0.3),
                          c(1, 4, 0.5, 2))
  expect_equal(as.matrix(G), t(as.matrix(G)))
  expect_equal(as.numeric(G %*% rep(1, 4)), rep(0, 4),
               tolerance = 1e-12)
  expect_gte(min(eigen(as.matrix(G), symmetric = TRUE)$values), -1e-12)
})
