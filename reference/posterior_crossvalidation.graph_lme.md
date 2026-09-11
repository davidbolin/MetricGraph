# Cross-validation for `graph_lme` models assuming observations at the vertices of metric graphs

This function performs cross-validation by computing predictions for
test data using either the posterior distribution from a fitted model
(pseudo-CV) or by refitting the model for each fold (true CV).

## Usage

``` r
# S3 method for class 'graph_lme'
posterior_crossvalidation(
  object,
  scores = c("logscore", "crps", "scrps", "mae", "rmse"),
  mode = "k-fold",
  k = 10,
  percentage = 20,
  number_folds = 10,
  train_test_indices = NULL,
  true_CV = FALSE,
  factor = 1,
  tibble = TRUE,
  parallel_folds = FALSE,
  parallel_fitting = FALSE,
  n_cores = parallel::detectCores() - 1,
  print = FALSE,
  seed = NULL,
  return_indices = FALSE,
  use_precomputed = TRUE,
  ...
)
```

## Arguments

- object:

  A fitted model using the
  [`graph_lme()`](https://davidbolin.github.io/MetricGraph/reference/graph_lme.md)
  function. Lists of fitted models are handled by the list method of the
  generic
  [`rSPDE::posterior_crossvalidation()`](https://davidbolin.github.io/rSPDE/reference/posterior_crossvalidation.html),
  which can also mix `graph_lme` and `rspde_lme` fits.

- scores:

  A vector of scores to compute. The options are "logscore", "crps",
  "scrps", "mae", and "rmse". By default, all scores are computed.

- mode:

  Cross-validation mode. Options are "k-fold", "loo" (leave-one-out), or
  "lpo" (leave-percentage-out). Default is "k-fold".

- k:

  Number of folds for k-fold cross-validation. Default is 10.

- percentage:

  The percentage (from 1 to 99) of the data to be used to train the
  model. Will only be used if `mode` is "lpo". Default is 80.

- number_folds:

  Number of folds to be done if `mode` is "lpo". Default is 10.

- train_test_indices:

  Optional list containing train and test indices for each fold. If
  provided, k, mode, and percentage are ignored.

- true_CV:

  Logical indicating whether to refit the model for each fold (TRUE) or
  use the posterior distribution from the fitted model (FALSE). Default
  is FALSE.

- factor:

  Which factor to multiply the scores. The default is 1.

- tibble:

  Return the scores as a
  [`tidyr::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)

- parallel_folds:

  Logical indicating whether to run computations in parallel across
  folds. Default is FALSE.

- parallel_fitting:

  Logical indicating whether to run model fitting in parallel. Default
  is FALSE.

- n_cores:

  Number of cores to use for parallel computation. Default is
  parallel::detectCores() - 1.

- print:

  Logical indicating whether to print progress of which fold is being
  processed. Default is FALSE.

- seed:

  Random seed for reproducibility in fold creation. Default is NULL.

- return_indices:

  Logical indicating whether to return the train/test indices used.
  Default is FALSE.

- use_precomputed:

  Logical indicating whether to use precomputation for faster CV.
  Default is TRUE.

- ...:

  Currently not used.

## Value

Vector with the posterior expectations and variances as well as mean
absolute error (MAE), root mean squared errors (RMSE), and three
negatively oriented proper scoring rules: log-score, CRPS, and scaled
CRPS.
