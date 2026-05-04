# Data extraction for 'spde' models

Extracts data from metric graphs to be used by 'INLA' and 'inlabru'.

## Usage

``` r
graph_data_spde(
  graph_spde,
  name = "field",
  repl = NULL,
  repl_col = NULL,
  group = NULL,
  group_col = NULL,
  likelihood_col = NULL,
  resp_col = NULL,
  covariates = NULL,
  only_pred = FALSE,
  loc_name = NULL,
  tibble = FALSE,
  drop_na = FALSE,
  drop_all_na = TRUE,
  loc = deprecated()
)
```

## Arguments

- graph_spde:

  An `inla_metric_graph_spde` object built with the
  [`graph_spde()`](https://davidbolin.github.io/MetricGraph/reference/graph_spde.md)
  function.

- name:

  A character string with the base name of the effect.

- repl:

  Which replicates? If there is no replicates, one can set `repl` to
  `NULL`. If one wants all replicates, then one sets to `repl` to
  `.all`.

- repl_col:

  Column containing the replicates. If the replicate is the internal
  group variable, set the replicates to ".group". If not replicates, set
  to `NULL`.

- group:

  Which groups? If there is no groups, one can set `group` to `NULL`. If
  one wants all groups, then one sets to `group` to `.all`.

- group_col:

  Which "column" of the data contains the group variable?

- likelihood_col:

  If only a single likelihood, this variable should be `NULL`. In case
  of multiple likelihoods, which column contains the variable indicating
  the number of the likelihood to be considered?

- resp_col:

  If only a single likelihood, this variable should be `NULL`. In case
  of multiple likelihoods, column containing the response variable.

- covariates:

  Vector containing the column names of the covariates. If no
  covariates, then it should be `NULL`.

- only_pred:

  Should only return the `data.frame` to the prediction data?

- loc_name:

  Character with the name of the location variable to be used in
  'inlabru' prediction.

- tibble:

  Should the data be returned as a
  [`tidyr::tibble`](https://tibble.tidyverse.org/reference/tibble.html)?

- drop_na:

  Should the rows with at least one NA for one of the columns be
  removed? DEFAULT is `FALSE`. This option is turned to `FALSE` if
  `only_pred` is `TRUE`.

- drop_all_na:

  Should the rows with all variables being NA be removed? DEFAULT is
  `TRUE`. This option is turned to `FALSE` if `only_pred` is `TRUE`.

- loc:

  **\[deprecated\]** Use `loc_name` instead.

## Value

An 'INLA' and 'inlabru' friendly list with the data.
