# Precompute expensive quantities for efficient LGCP model fitting

This function precomputes the computationally expensive quantities
needed for fitting log-Gaussian Cox process (LGCP) models on metric
graphs. It enables efficient refitting of multiple models with different
formulas while reusing the same spatial structure, integration points,
and covariates. This is particularly beneficial for model selection,
cross-validation, or sensitivity analysis scenarios.

## Usage

``` r
precompute_lgcp_graph(
  graph,
  resp_variable_name,
  model_name,
  spde_model,
  covariates = NULL,
  interpolate = TRUE,
  manual_integration_points = NULL,
  manual_covariates = NULL,
  use_current_mesh = TRUE,
  new_h = NULL,
  new_n = NULL,
  repl = ".all",
  repl_col = ".group",
  clone_graph = TRUE
)
```

## Arguments

- graph:

  A `metric_graph` object containing the network structure and point
  pattern data. Must have observations added via `add_observations()`.

- resp_variable_name:

  Character. Name of the response variable (typically binary: 0 or 1) in
  the graph data that represents point occurrences.

- model_name:

  Character. Name to be used for the spatial field in INLA formulas
  (e.g., "field" for `f(field, model = spde_model)`).

- spde_model:

  An SPDE model object of class `inla_metric_graph_spde` (from
  [`graph_spde()`](https://davidbolin.github.io/MetricGraph/reference/graph_spde.md))
  or `rspde_metric_graph` (from `rspde.metric_graph()`).

- covariates:

  Character vector. Names of covariates to include in the
  precomputation. Only covariates listed here can be used in subsequent
  model fits with the precomputed data.

- interpolate:

  Logical. If `TRUE` (default), interpolate covariate values from graph
  data to integration points. If `FALSE`, use `manual_covariates`.

- manual_integration_points:

  Data frame with columns `edge_number`, `distance_on_edge`, and `E`
  (integration weights). If `NULL`, automatic integration points are
  created from the mesh or specified parameters.

- manual_covariates:

  Data frame or named list containing covariate values at integration
  points. Required when `interpolate = FALSE`. Must include a `.group`
  column for replicates if using replicated data.

- use_current_mesh:

  Logical. If `TRUE` (default), use the existing mesh in
  `graph$mesh$VtE` as integration points. If `FALSE` or no mesh exists,
  create a new mesh using `new_h` or `new_n`.

- new_h:

  Numeric. Mesh resolution for creating a new mesh when
  `use_current_mesh = FALSE`. Smaller values create finer meshes.

- new_n:

  Integer. Alternative to `new_h`, specifies the approximate number of
  mesh nodes for the new mesh.

- repl:

  Character vector or `".all"`. Specifies which replicates to include in
  the model. Use `".all"` to include all available replicates.

- repl_col:

  Character. Name of the column in the graph data that contains
  replicate identifiers. Default is `".group"`.

- clone_graph:

  Logical. If `TRUE` (default), clone the graph to avoid modifying the
  original object. If `FALSE`, work directly on the original graph
  (faster but modifies the input).

## Value

A list of class `"precomputed_lgcp"` containing:

- graph:

  The modified metric_graph object with integration points

- covariates:

  Vector of available covariate names

- stk:

  Precomputed INLA stack object

- resp_variable_name:

  Name of the response variable

- aux_spde_model:

  The SPDE model object prepared for the graph

- type_model:

  Type of SPDE model ("exact" or "rational")

- model_name:

  Name of the spatial field for formulas

## Details

This function is most beneficial when:

- Fitting multiple models with different covariate combinations

- Performing model selection or cross-validation

- Using exact SPDE models (which have higher setup costs)

- Working with large graphs or complex spatial structures

The speedup typically becomes apparent when fitting 3 or more models,
with greater benefits for more complex spatial structures and exact SPDE
models.

## Note

- All covariates you plan to use must be specified in the initial
  precomputation

- The spatial structure (mesh, integration points) is fixed during
  precomputation

- Manual covariate values should correspond to mesh nodes in
  `graph$mesh$VtE` when `interpolate = FALSE`

## See also

[`lgcp_graph`](https://davidbolin.github.io/MetricGraph/reference/lgcp_graph.md)
for fitting LGCP models with precomputed data,
[`graph_lgcp_sim`](https://davidbolin.github.io/MetricGraph/reference/graph_lgcp_sim.md)
for simulating LGCP data,
[`graph_spde`](https://davidbolin.github.io/MetricGraph/reference/graph_spde.md)
for exact SPDE models,
[`spde_metric_graph_result`](https://davidbolin.github.io/MetricGraph/reference/spde_metric_graph_result.md)
for extracting spatial parameter estimates

## Examples

``` r
if (FALSE) { # \dontrun{
# Setup: Create graph with mesh and add point pattern data
graph <- metric_graph$new()
graph$build_mesh(h = 0.1)

# Add point pattern data with covariates
point_data <- data.frame(
  y = 1,  # All observed locations have y = 1
  edge_number = c(1, 2, 3, 1, 2),
  distance_on_edge = c(0.1, 0.3, 0.8, 0.9, 0.2),
  elevation = c(100, 150, 200, 120, 180),
  temperature = c(15, 12, 8, 14, 10)
)
graph$add_observations(point_data, normalized = TRUE)

# Create SPDE model
spde_model <- graph_spde(graph, alpha = 1)

# Precompute for multiple model fitting scenarios
precomputed_data <- precompute_lgcp_graph(
  graph = graph,
  resp_variable_name = "y",
  model_name = "field",
  spde_model = spde_model,
  covariates = c("elevation", "temperature"),
  use_current_mesh = TRUE
)

# Now fit multiple models efficiently
# Model selection: compare different covariate combinations
fit1 <- lgcp_graph(y ~ elevation + f(field, model = spde_model),
                   graph = graph, precomputed_data = precomputed_data)

fit2 <- lgcp_graph(y ~ temperature + f(field, model = spde_model),
                   graph = graph, precomputed_data = precomputed_data)

fit3 <- lgcp_graph(y ~ elevation + temperature + f(field, model = spde_model),
                   graph = graph, precomputed_data = precomputed_data)

# Compare models using log marginal likelihood
c(fit1$mlik[1], fit2$mlik[1], fit3$mlik[1])

# Using manual covariates (exact values at mesh nodes)
manual_covs <- data.frame(
  elevation = graph$mesh$VtE[,1] * 50 + 100,  # Synthetic elevation
  temperature = 20 - graph$mesh$VtE[,1] * 10, # Synthetic temperature
  .group = 1
)

precomputed_manual <- precompute_lgcp_graph(
  graph = graph,
  resp_variable_name = "y", 
  model_name = "field",
  spde_model = spde_model,
  covariates = c("elevation", "temperature"),
  manual_covariates = manual_covs,
  interpolate = FALSE
)

# For maximum performance (modifies original graph)
precomputed_fast <- precompute_lgcp_graph(
  graph = graph,
  resp_variable_name = "y",
  model_name = "field", 
  spde_model = spde_model,
  covariates = c("elevation", "temperature"),
  clone_graph = FALSE
)

# Example with replicates
replicate_data <- data.frame(
  y = 1,
  edge_number = c(1, 2, 1, 3, 2, 3),
  distance_on_edge = c(0.2, 0.4, 0.7, 0.1, 0.8, 0.6),
  elevation = c(110, 160, 130, 190, 170, 210),
  replicate_id = c(1, 1, 1, 2, 2, 2)
)

graph$clear_observations()
graph$add_observations(replicate_data, normalized = TRUE, group = "replicate_id")

precomputed_reps <- precompute_lgcp_graph(
  graph = graph,
  resp_variable_name = "y",
  model_name = "field",
  spde_model = spde_model, 
  covariates = "elevation",
  repl = ".all",
  repl_col = "replicate_id"
)

fit_reps <- lgcp_graph(y ~ elevation + f(field, model = spde_model, 
                                        replicate = field.repl),
                       graph = graph, precomputed_data = precomputed_reps)
} # }
```
