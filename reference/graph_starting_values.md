# Starting values for random field models on metric graphs

Computes appropriate starting values for optimization of Gaussian random
field models on metric graphs.

## Usage

``` r
graph_starting_values(
  graph,
  model = c("alpha1", "alpha2", "isoExp", "GL1", "GL2"),
  data = TRUE,
  data_name = NULL,
  range_par = FALSE,
  nu = FALSE,
  manual_data = NULL,
  like_format = FALSE,
  log_scale = FALSE,
  model_options = list(),
  rec_tau = TRUE,
  factor_start_range = 0.3,
  type_start_range_bbox = "diag"
)
```

## Arguments

- graph:

  A `metric_graph` object.

- model:

  Type of model, "alpha1", "alpha2", "isoExp", "GL1", and "GL2" are
  supported.

- data:

  Should the data be used to obtain improved starting values?

- data_name:

  The name of the response variable in `graph$data`.

- range_par:

  Should an initial value for range parameter be returned instead of for
  kappa?

- nu:

  Should an initial value for nu be returned?

- manual_data:

  A vector (or matrix) of response variables.

- like_format:

  Should the starting values be returned with sigma.e as the last
  element? This is the format for the likelihood constructor from the
  'rSPDE' package.

- log_scale:

  Should the initial values be returned in log scale?

- model_options:

  List object containing the model options.

- rec_tau:

  Should a starting value for the reciprocal of tau be given?

- factor_start_range:

  Factor to multiply the max/min/diagonal dimension of the bounding box
  to obtain a starting value for range. Default is 0.5.

- type_start_range_bbox:

  Which dimension from the bounding box should be used? The options are
  'diag', the default, 'max' and 'min'.

## Value

A vector, `c(start_sigma_e, start_sigma, start_kappa)`
