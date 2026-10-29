test_that("binomial predictive draws respect observation trial counts", {
  fit <- list(.args = list(family = "binomial"))
  lp <- matrix(rep(0.5, 3000), nrow = 3)
  set.seed(19)
  samples <- .cv_get_response_samples(lp, fit, 1L, 1000L, FALSE,
                                      trial_counts = c(1L, 2L, 10L))
  expect_length(samples, 3L)
  expect_true(all(vapply(seq_along(samples), function(i) {
    all(samples[[i]] >= 0 & samples[[i]] <= c(1, 2, 10)[i])
  }, logical(1))))
  expect_equal(vapply(samples, mean, numeric(1)), c(0.5, 1, 5),
               tolerance = 0.4)
  expect_error(.cv_get_response_samples(lp, fit, 1L, 1000L, FALSE,
                                         trial_counts = c(1, 2, 0)),
               "positive integers")
})

test_that("loss scores and best-model direction agree in both orientations", {
  loss <- c(A = 2, B = 1)
  negative <- .cv_orient_loss(loss, "negative")
  positive <- .cv_orient_loss(loss, "positive")
  expect_identical(names(which.min(negative)), "B")
  expect_identical(names(which.max(positive)), "B")
  expect_equal(positive, -negative)
})
