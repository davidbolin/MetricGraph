# Closed-form marginal variance of the directional OU model

`r(x,x)` at each evaluation point – the diagonal of
[`directional_ou_covariance()`](https://davidbolin.github.io/MetricGraph/reference/directional_ou_covariance.md),
computed directly via the marginal variance recursion along each edge,
without any last-common-ancestor or transfer-factor work.

## Usage

``` r
directional_ou_variance(
  graph,
  kappa,
  tau,
  PtE = NULL,
  sigma_source = NULL,
  normalized = TRUE,
  cpp = TRUE
)
```

## Arguments

- graph:

  A `metric_graph` object with directional edge weights set
  (`set_edge_weights(directional_weights = ...)`) and the vertex weight
  functions set (`setDirectionalWeightFunction()`).

- kappa:

  Rate parameter of the OU process (`kappa > 0`).

- tau:

  White-noise scale. The stationary marginal variance is
  `1/(2*kappa*tau^2)` – the same convention as `MetricGraph:::r_1()`.
  Note this is NOT `reciprocal_tau` (used elsewhere in the package as
  `tau = 1/reciprocal_tau`); pass `tau` directly here.

- PtE:

  Evaluation points, a matrix with columns
  `(edge_number, distance_on_edge)`. Defaults to `graph$get_PtE()` (the
  observation locations).

- sigma_source:

  Anchoring variance at the source vertices (indegree 0). `NULL`
  (default) anchors every source at the stationary variance
  `1/(2*kappa*tau^2)`, matching
  `Qalpha1_edges_custom(..., stationary_points = "all")`. Otherwise a
  numeric vector giving `Var(u)` at each source vertex, either named by
  vertex index (as character) or positional in the order
  `which(graph$get_degrees("indegree") == 0)`.

- normalized:

  If `TRUE` (default), `distance_on_edge` in `PtE`/`PtE2` is in
  `[0, 1]`; if `FALSE`, distances are absolute.

- cpp:

  If `TRUE` (default), use the C++ numeric setup and the shared C++
  covariance fill for oriented in-trees and out-trees. Irregular trees
  use the generic R fallback. Set to `FALSE` to use the pure-R
  implementation throughout.

## Value

A numeric vector of length `nrow(PtE)`.

## See also

[`directional_ou_covariance()`](https://davidbolin.github.io/MetricGraph/reference/directional_ou_covariance.md)
for the full closed-form covariance (including the
`\deqn{r(s,t) = r(a,a) A(a,s) A(a,t)}` formula and the tree/acyclicity
requirements).

## Examples

``` r
edge1 <- rbind(c(0, 0), c(1, 0))
edge2 <- rbind(c(1, 0), c(2, 0))
edge3 <- rbind(c(1, 0), c(1, 1))
graph <- metric_graph$new(edges = list(edge1, edge2, edge3))
#> Starting graph creation...
#> LongLat is set to FALSE
#> Creating edges...
#> Setting edge weights...
#> Computing bounding box...
#> Setting up edges
#> Merging close vertices
#> Total construction time: 0.22 secs
#> Creating and updating vertices...
#> Storing the initial graph...
#> Computing the relative positions of the edges...
graph$set_edge_weights(weights = data.frame(w = c(1, 1, 1)),
                       directional_weights = "w")
graph$setDirectionalWeightFunction(f_in = function(x) sqrt(x / sum(x)))
directional_ou_variance(graph, kappa = 1, tau = 1,
                        PtE = rbind(c(1, 0.5), c(2, 0.5)))
#> [1] 0.5 0.5
```
