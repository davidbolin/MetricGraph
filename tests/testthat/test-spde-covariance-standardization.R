test_that("standardized covariance inversion validates its inputs", {
  for (bad in list(1, matrix(numeric(0), 0, 0), matrix(1, 2, 3),
                   matrix(NA_real_, 2, 2), matrix(Inf, 2, 2),
                   matrix(1 + 1i, 2, 2))) {
    expect_error(invert_spde_covariance(bad), "finite, nonempty, real square matrix")
  }
  for (variance in c(0, -1)) {
    expect_error(invert_spde_covariance(diag(c(1, variance))),
                 "strictly positive diagonal")
  }
  # Scaling must not suppress solve()'s singularity check.
  expect_error(invert_spde_covariance(matrix(1, 2, 2)), "singular")
})
