# Simulate log-Gaussian Cox processes on metric graphs

This function simulates point patterns from a log-Gaussian Cox process
(LGCP) driven by Whittle-Matérn Gaussian random fields on metric graphs.
The intensity function is modeled as \\\lambda(s) = exp(\beta + u(s))\\,
where \\\beta\\ is an intercept parameter and \\u(s)\\ is a Gaussian
field with Whittle-Matérn covariance.

## Usage

``` r
graph_lgcp_sim(n = 1, intercept = 0, sigma, range, alpha, graph)
```

## Arguments

- n:

  Integer. Number of replicate point patterns to simulate. Default is 1.

- intercept:

  Numeric scalar or vector. Mean value(s) of the log-intensity field.
  Can be a constant or vary spatially if provided as a vector of length
  equal to the number of mesh nodes.

- sigma:

  Numeric. Marginal standard deviation parameter of the Gaussian field.

- range:

  Numeric. Practical correlation range parameter of the Gaussian field.

- alpha:

  Numeric. Smoothness parameter of the Whittle-Matérn field. Currently
  supports values 1 and 2, corresponding to exponential and Matérn-3/2
  covariance functions respectively.

- graph:

  A `metric_graph` object with a built mesh. The graph must have an
  existing mesh (built with `graph$build_mesh()`) and computed FEM
  matrices (computed with `graph$compute_fem()`).

## Value

If `n = 1`, returns a list with components:

- u:

  Numeric vector of the simulated Gaussian field values at mesh nodes

- edge_number:

  Integer vector of edge numbers where points were simulated

- edge_loc:

  Numeric vector of locations along edges (normalized coordinates)

If `n > 1`, returns a list of length `n`, where each element is a list
with the above components.

## Details

The function implements a two-step simulation procedure:

1.  Simulate the Gaussian random field u(s) using finite element methods

2.  Generate point locations using acceptance-rejection sampling from
    the intensity \\\lambda(s) = exp(\beta + u(s))\\

The Gaussian field is characterized by the SPDE: \\(\kappa^2 -
\Delta)^{\alpha/2} \tau u = \mathcal{W}\\ where \\\kappa, \tau\\ are
derived from the range and sigma parameters, and \\\mathcal{W}\\ is
white noise.

## See also

[`lgcp_graph`](https://davidbolin.github.io/MetricGraph/reference/lgcp_graph.md)
for fitting LGCP models,
[`precompute_lgcp_graph`](https://davidbolin.github.io/MetricGraph/reference/precompute_lgcp_graph.md)
for efficient model fitting

## Examples

``` r
if (FALSE) { # \dontrun{
# Create a metric graph and build mesh
graph <- metric_graph$new()
graph$build_mesh(h = 0.1)
graph$compute_fem()

# Simulate a single point pattern
lgcp_data <- graph_lgcp_sim(
  intercept = -1, 
  sigma = 0.5, 
  range = 2, 
  alpha = 2, 
  graph = graph
)

# Simulate multiple replicates
lgcp_replicates <- graph_lgcp_sim(
  n = 10,
  intercept = -1, 
  sigma = 0.5, 
  range = 2, 
  alpha = 2, 
  graph = graph
)

# Plot the simulated intensity
graph$plot_function(X = exp(lgcp_data$u), vertex_size = 0)
} # }
```
