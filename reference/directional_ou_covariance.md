# Closed-form covariance of the directional OU model

Evaluates the covariance function of the "proper global OU" directional
model (Whittle-Matern alpha = 1 with directional vertex conditions) in
closed form, via the last-common-ancestor (LCA) representation, without
assembling or inverting the precision matrix.

## Usage

``` r
directional_ou_covariance(
  graph,
  kappa,
  tau,
  PtE = NULL,
  PtE2 = NULL,
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

- PtE2:

  Second set of evaluation points for the cross-covariance
  `Cov(u(PtE), u(PtE2))`. Defaults to `PtE`, giving the symmetric
  `n x n` covariance matrix at the observation locations.

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

A dense `matrix` with `nrow(PtE)` rows and `nrow(PtE2)` columns.

## Details

The covariance between two points `s` and `t` is \$\$r(s,t) = r(a,a) \\
A(a,s) \\ A(a,t)\$\$ where `a` is the last common ancestor (LCA) of `s`
and `t` (or covariance `0` if no common ancestor exists), `r(a,a)` is
the marginal variance at `a`, and \$\$A(x,y) = \exp(-\kappa \\ d(x,y))
\prod_v \beta_v\$\$ is the transfer factor for the (unique, forward)
path from `x` to `y`, the product running over the `beta_v` vertex
weights of every vertex crossed strictly between `x` and `y`.

In plain terms: the last common ancestor of two points is the
furthest-downstream point that is upstream of both of them (it does not
exist if the two points lie on disjoint headwater branches, in which
case their covariance is exactly zero). The transfer factor `A(x,y)` is
how much of the value at `x` propagates forward to `y`: exponential
decay over the distance `d(x,y)`, attenuated by the directional weight
`beta_v` at every vertex the path passes through.

This closed form requires `graph` to represent a directed *tree*: the
undirected skeleton must have no cycles, and the directed graph must
itself be acyclic. If `graph` has a directed cycle, is not a tree, or
has no directional weight function set, this errors (via the internal
`directional_ou_setup()` validation) – see
`graph$setDirectionalWeightFunction()` to set one.

## See also

[`directional_ou_variance()`](https://davidbolin.github.io/MetricGraph/reference/directional_ou_variance.md)
for the diagonal only.

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
#> Total construction time: 0.23 secs
#> Creating and updating vertices...
#> Storing the initial graph...
#> Computing the relative positions of the edges...
graph$set_edge_weights(weights = data.frame(w = c(1, 1, 1)),
                       directional_weights = "w")
graph$setDirectionalWeightFunction(f_in = function(x) sqrt(x / sum(x)))
Sigma <- directional_ou_covariance(graph, kappa = 1, tau = 1,
                                   PtE = rbind(c(1, 0.5), c(2, 0.5)))
Sigma
#>           [,1]      [,2]
#> [1,] 0.5000000 0.1839397
#> [2,] 0.1839397 0.5000000
```
