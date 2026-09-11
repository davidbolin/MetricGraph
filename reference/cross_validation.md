# Perform cross-validation on a list of fitted inlabru models on metric graphs.

Mirrors
[`rSPDE::cross_validation()`](https://davidbolin.github.io/rSPDE/reference/cross_validation.html)
for `bru` fits (output from
[`inlabru::bru()`](https://inlabru-org.github.io/inlabru/reference/bru.html)).
For each fold, the held-out responses are set to NA with
[`inlabru::bru_set_missing()`](https://inlabru-org.github.io/inlabru/reference/bru_set_missing.html)
and the model is refit with
[`inlabru::bru_rerun()`](https://inlabru-org.github.io/inlabru/reference/bru.html).
When `true_CV = FALSE`, the original fit's `mode$theta` is held fixed so
only the latent posterior is updated; when `TRUE`, theta is
re-estimated. Posterior response samples are then drawn at the held-out
locations via
[`inlabru::generate()`](https://generics.r-lib.org/reference/generate.html).
Works with the exact, non-FEM SPDE models in MetricGraph
(`inla_metric_graph_spde`) and FEM-based rSPDE models alike.

## Usage

``` r
cross_validation(
  models,
  model_names = NULL,
  scores = c("mae", "mse", "crps", "scrps", "dss"),
  cv_type = c("k-fold", "loo", "lpo"),
  weight_thr = NULL,
  k = 5,
  percentage = 20,
  number_folds = 10,
  n_samples = 1000,
  return_scores_folds = FALSE,
  orientation_results = c("negative", "positive"),
  include_best = TRUE,
  train_test_indexes = NULL,
  return_train_test = FALSE,
  return_post_samples = FALSE,
  return_true_test_values = FALSE,
  parallelize_RP = FALSE,
  n_cores_RP = parallel::detectCores() - 1,
  true_CV = TRUE,
  save_settings = FALSE,
  model_options_bru = list(),
  print = TRUE,
  fit_verbose = FALSE
)
```

## Arguments

- models:

  A fitted model from
  [`inlabru::bru()`](https://inlabru-org.github.io/inlabru/reference/bru.html)
  or a list of such models. All models must have the same number of
  likelihoods and be fitted to identical datasets.

- model_names:

  Character vector of model names. Defaults to `names(models)` or
  `Model 1`, `Model 2`, ... if absent.

- scores:

  Subset of
  `c("mae", "mse", "crps", "scrps", "dss", "wcrps", "swcrps")`.

- cv_type:

  One of `"k-fold"`, `"loo"`, `"lpo"`.

- weight_thr:

  Threshold for `wcrps` / `swcrps`. Required if either is requested.

- k:

  Number of folds for `k-fold`.

- percentage:

  Train percentage (1-99) for `lpo`.

- number_folds:

  Number of folds for `lpo`.

- n_samples:

  Number of posterior samples for sample-based scoring.

- return_scores_folds:

  If `TRUE`, return per-fold score matrices.

- orientation_results:

  One of `"negative"` (lower is better) or `"positive"` (higher is
  better).

- include_best:

  Add a `Best` row indicating the best model per score.

- train_test_indexes:

  Optional pre-built fold list. If supplied, `cv_type`, `k`,
  `percentage`, `number_folds` are ignored.

- return_train_test:

  Return the train/test indices used.

- return_post_samples:

  Return posterior response samples (forces
  `return_scores_folds = TRUE`).

- return_true_test_values:

  Return the true response values at test points.

- parallelize_RP:

  Parallelize CRPS/SCRPS pairwise integrals.

- n_cores_RP:

  Number of cores for `parallelize_RP`.

- true_CV:

  If `TRUE`, refit the model on each training fold; if `FALSE`, sample
  directly from the supplied fit without refitting.

- save_settings:

  If `TRUE`, return the CV settings used.

- model_options_bru:

  List of options passed to inlabru.

- print:

  Print partial progress.

- fit_verbose:

  Pass `verbose` through to INLA when refitting.

## Value

Either a `data.frame` (default), or a list with `scores_df` plus the
optional components requested by `return_*` / `save_settings` arguments.
