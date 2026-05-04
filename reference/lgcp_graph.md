# Fit log-Gaussian Cox process models on metric graphs

This function fits log-Gaussian Cox process (LGCP) models for point
pattern data on metric graphs using R-INLA. It handles the complex setup
required for LGCP modeling, including creation of integration points,
data preparation, and interface with INLA's Poisson likelihood
framework. The function supports both exact and rational SPDE
approximations, multiple replicates, and efficient refitting using
precomputed quantities.

## Usage

``` r
lgcp_graph(
  formula,
  graph,
  interpolate = TRUE,
  manual_integration_points = NULL,
  manual_covariates = NULL,
  use_current_mesh = TRUE,
  new_h = NULL,
  new_n = NULL,
  repl = ".all",
  repl_col = ".group",
  clone_graph = TRUE,
  precomputed_data = NULL,
  ...
)
```

## Arguments

- formula:

  A formula object specifying the model structure. Should follow INLA
  syntax, e.g., `y ~ covariate + f(field, model = spde_model)` where
  `field` is the spatial random effect and `spde_model` is an SPDE model
  object.

- graph:

  A `metric_graph` object containing the network structure and point
  pattern data. Must have observations added via `add_observations()`.

- interpolate:

  Logical. If `TRUE` (default), interpolate covariate values from graph
  data to integration points. If `FALSE`, use `manual_covariates`.

- manual_integration_points:

  Data frame with columns `edge_number`, `distance_on_edge`, and `E`
  (integration weights). If `NULL`, automatic integration points are
  created.

- manual_covariates:

  Data frame containing covariate values at integration points when
  `interpolate = FALSE`. Must include a `.group` column for replicates
  if using replicated data.

- use_current_mesh:

  Logical. If `TRUE` (default), use the existing mesh in the graph as
  integration points. If `FALSE`, create a new mesh.

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
  (faster but modifies the input). Only used when `precomputed_data` is
  `NULL`.

- precomputed_data:

  Optional object of class `"precomputed_lgcp"` from
  [`precompute_lgcp_graph()`](https://davidbolin.github.io/MetricGraph/reference/precompute_lgcp_graph.md).
  Enables efficient refitting with different formulas using the same
  spatial structure and covariates.

- ...:

  Additional arguments passed to
  [`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html), such as
  `control.inla`, `control.predictor`, `control.compute`, etc.

## Value

An object of class `"inla"` containing the fitted LGCP model. This
includes posterior marginal distributions for model parameters, fitted
values, and other standard INLA output components. Use
[`spde_metric_graph_result()`](https://davidbolin.github.io/MetricGraph/reference/spde_metric_graph_result.md)
to extract spatial parameter estimates in their original scale.

## Details

The function implements LGCP modeling using the approach of Simpson et
al. (2016), where the log-Gaussian Cox process with intensity
\\\lambda(s) = exp(\eta(s))\\ is approximated using a Poisson likelihood
with carefully constructed integration points and weights.

The key steps are:

1.  Create integration points (typically mesh nodes) across the graph

2.  Set up data with observed points (response = 1, weights = 0) and
    integration points (response = 0, weights = integration weights)

3.  Fit using Poisson regression with the constructed weights as
    exposure

The spatial component can be modeled using:

- Exact SPDE models via
  [`graph_spde()`](https://davidbolin.github.io/MetricGraph/reference/graph_spde.md)
  (slower setup, exact likelihood)

- Rational SPDE approximations via `rspde.metric_graph()` (faster,
  approximate)

## Performance

For fitting multiple models with the same spatial structure:

- Use
  [`precompute_lgcp_graph()`](https://davidbolin.github.io/MetricGraph/reference/precompute_lgcp_graph.md)
  first, then `lgcp_graph()` with `precomputed_data` for substantial
  speedups

- Set `clone_graph = FALSE` for additional performance gains when you
  don't need to preserve the original graph

## References

Simpson, D., Illian, J., Lindgren, F., Sørbye, S., & Rue, H. (2016).
Going off grid: Computationally efficient inference for log-Gaussian Cox
processes. Biometrika, 103(1), 49-70.

## See also

[`graph_lgcp_sim`](https://davidbolin.github.io/MetricGraph/reference/graph_lgcp_sim.md)
for simulating LGCP data,
[`precompute_lgcp_graph`](https://davidbolin.github.io/MetricGraph/reference/precompute_lgcp_graph.md)
for efficient model refitting,
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
graph$add_observations(data = your_point_data, 
                       edge_number = "edge_id", 
                       distance_on_edge = "location")

# Create SPDE model
spde_model <- graph_spde(graph, alpha = 1)
# or: rspde_model <- rspde.metric_graph(graph, nu = 1.5)

# Fit basic LGCP model
fit1 <- lgcp_graph(y ~ 1 + f(field, model = spde_model), 
                   graph = graph)

# Fit model with covariates
fit2 <- lgcp_graph(y ~ elevation + temperature + 
                       f(field, model = spde_model), 
                   graph = graph)

# Extract spatial parameter estimates
spde_result <- spde_metric_graph_result(fit2, "field", spde_model)
summary(spde_result)

# Efficient fitting of multiple models
precomputed <- precompute_lgcp_graph(
  graph = graph,
  resp_variable_name = "y",
  model_name = "field", 
  spde_model = spde_model,
  covariates = c("elevation", "temperature", "slope")
)

# Now fit multiple models efficiently
fit_a <- lgcp_graph(y ~ elevation + f(field, model = spde_model), 
                    graph = graph, precomputed_data = precomputed)
fit_b <- lgcp_graph(y ~ elevation + temperature + f(field, model = spde_model), 
                    graph = graph, precomputed_data = precomputed)
fit_c <- lgcp_graph(y ~ slope + f(field, model = spde_model), 
                    graph = graph, precomputed_data = precomputed)

# Model with replicates
fit_rep <- lgcp_graph(y ~ covariate + f(field, model = spde_model, 
                                       replicate = field.repl), 
                      graph = graph)
} # }
```
