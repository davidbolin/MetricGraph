# Parallel simulation of a Whittle-Matérn field on a metric graph

Calls
[`simulate.metric_graph()`](https://davidbolin.github.io/MetricGraph/reference/simulate.metric_graph.md)
with `parallel = TRUE`. The edge loop is distributed over a
`doParallel`/`foreach` cluster.

## Usage

``` r
simulate_parallel(graph, ..., n_cores = NULL, cluster = NULL)
```

## Arguments

- graph:

  A `metric_graph` object.

- ...:

  Further arguments passed to
  [`simulate.metric_graph()`](https://davidbolin.github.io/MetricGraph/reference/simulate.metric_graph.md),
  including `alpha`, `method`, `kappa`/`tau` or `range`/`sigma`, `PtE`,
  `type`, `BC`, `nsim`, and `seed`.

- n_cores:

  Number of parallel workers.

- cluster:

  An already-registered `doParallel` cluster.

## Value

Numeric vector (nsim = 1) or matrix with nsim columns.

## See also

[`sample_spde()`](https://davidbolin.github.io/MetricGraph/reference/sample_spde.md)

## Examples

``` r
# \donttest{
V <- rbind(c(0,0), c(1,0), c(1,1), c(0,1))
E <- rbind(c(1,2), c(2,3), c(3,4), c(4,1))
g <- metric_graph$new(V = V, E = E)
#> Starting graph creation...
#> LongLat is set to FALSE
#> Creating edges...
#> Setting edge weights...
#> Computing bounding box...
#> Setting up edges
#> Merging close vertices
#> Total construction time: 0.71 secs
#> Creating and updating vertices...
#> Storing the initial graph...
#> Computing the relative positions of the edges...
t_norm <- seq(0.1, 0.9, by = 0.2)
PtE    <- do.call(rbind, lapply(1:4, function(e) cbind(e, t_norm)))
u <- simulate_parallel(g, alpha = 1, method = "kriging",
                       kappa = 1, tau = 1, PtE = PtE, n_cores = 2)
# }
```
