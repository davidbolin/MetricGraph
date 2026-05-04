# Leave-one-out pseudo-crossvalidation for `graph_lme` models assuming observations at the vertices of metric graphs

This function performs pseudo-crossvalidation by computing leave-one-out
predictions using the posterior distribution from a fitted model. In
pseudo-crossvalidation, the model parameters are kept fixed at the
values estimated from the full dataset (those provided in the object),
rather than re-estimating them for each fold.

## Usage

``` r
posterior_crossvalidation_loo(
  object,
  factor = 1,
  tibble = TRUE,
  which_repl = NULL
)
```

## Arguments

- object:

  A fitted model using the
  [`graph_lme()`](https://davidbolin.github.io/MetricGraph/reference/graph_lme.md)
  function or a named list of fitted objects using the
  [`graph_lme()`](https://davidbolin.github.io/MetricGraph/reference/graph_lme.md)
  function.

- factor:

  Which factor to multiply the scores. The default is 1.

- tibble:

  Return the scores as a
  [`tidyr::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)

- which_repl:

  Which replicates to consider?

## Value

Vector with the posterior expectations and variances as well as mean
absolute error (MAE), root mean squared errors (RMSE), and three
negatively oriented proper scoring rules: log-score, CRPS, and scaled
CRPS.
