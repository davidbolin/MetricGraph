# Simulate a Whittle-Matérn field on a metric graph

Draws unconditional (prior) samples of a Whittle-Matérn field on a
metric graph using one of two exact algorithms from Section 6.6 of the
paper.

## Usage

``` r
# S3 method for class 'metric_graph'
simulate(
  object,
  nsim = 1,
  seed = NULL,
  alpha = 1,
  method = c("direct", "kriging", "extended"),
  impl = c("cpp", "R"),
  kappa,
  tau,
  range,
  sigma,
  PtE = NULL,
  type = "manual",
  BC = 1,
  parallel = FALSE,
  n_cores = NULL,
  cluster = NULL,
  ...
)
```

## Arguments

- object:

  A `metric_graph` object.

- nsim:

  Number of samples. Returns a vector (nsim = 1) or matrix (nsim \> 1)
  of sampled values.

- seed:

  Integer seed for [`set.seed()`](https://rdrr.io/r/base/Random.html)
  before the vertex draw. Pass `NULL` (default) to use the current RNG
  state.

- alpha:

  Smoothness parameter (1 or 2).

- method:

  Simulation algorithm: `"direct"` (Method A, \\O(m^3)\\), `"kriging"`
  (Method B, \\O(m)\\), or `"extended"` (add simulation locations as
  graph vertices, single sparse Cholesky draw).

- impl:

  Implementation: `"cpp"` (default, Rcpp/Eigen) or `"R"` (pure R
  reference, for testing/comparison). Ignored for `method = "extended"`.

- kappa:

  Range parameter.

- tau:

  Precision parameter.

- range:

  Practical correlation range (alternative to `kappa`/`tau`).

- sigma:

  Marginal standard deviation (alternative to `kappa`/`tau`).

- PtE:

  Matrix with columns `(edge_number, normalised_position)`. Required
  when `type = "manual"` (the default).

- type:

  Location specification: `"manual"` (use `PtE`), `"mesh"` (mesh nodes),
  or `"obs"` (observation locations).

- BC:

  Boundary condition for degree-1 vertices: 1 = stationary (default), 0
  = Neumann.

- parallel:

  Logical. Use a parallel cluster for the edge loop?

- n_cores:

  Number of parallel workers. Ignored if `cluster` is supplied. Defaults
  to `detectCores() - 1`.

- cluster:

  An already-registered `doParallel` cluster. If `NULL` and
  `parallel = TRUE`, a cluster is created and stopped automatically.

- ...:

  Ignored.

## Value

Numeric vector (nsim = 1) or matrix with nsim columns of field values at
the requested locations, ordered to match `PtE`.

## Details

**Algorithm A ("direct")**: interior values on each edge are drawn from
the bridge conditional \\\mathcal{N}(S_e b_e, \Sigma^\*\_e)\\ via a
dense Cholesky of the \\m \times m\\ bridge covariance. Cost: \\O(m^3)\\
per edge. Preferred when \\m\\ is small (roughly \\m \lesssim 100\\).

**Algorithm B ("kriging")**: an unconditioned Matérn process is
simulated on each edge by a linear Markov recursion, then the bridge
boundary data are imposed by a kriging correction (Theorem 2 /
Proposition 4). Cost: \\O(m)\\ per edge. Preferred for large \\m\\.

Both algorithms produce draws from exactly the same finite-dimensional
distribution; they differ only in computational cost.

The two-step procedure is:

1.  Draw the vertex (boundary) state from the vertex precision via a
    sparse Cholesky.

2.  For each edge (independently given the vertex state), draw the
    interior values using the chosen algorithm.

When `parallel = TRUE` the edge loop is distributed over a
`doParallel`/`foreach` cluster. Reproducibility requires passing a
`seed` so that per-edge seeds are derived deterministically.

## See also

[`sample_spde()`](https://davidbolin.github.io/MetricGraph/reference/sample_spde.md)

## Examples

``` r
# Small square graph, alpha = 1
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
#> Total construction time: 0.72 secs
#> Creating and updating vertices...
#> Storing the initial graph...
#> Computing the relative positions of the edges...
t_norm <- seq(0.1, 0.9, by = 0.2)
PtE    <- do.call(rbind, lapply(1:4, function(e) cbind(e, t_norm)))
u <- simulate(g, alpha = 1, method = "direct", kappa = 1, tau = 1, PtE = PtE)
# \donttest{
# Method B on a larger graph, parallel
u2 <- simulate(g, alpha = 2, method = "kriging", range = 1, sigma = 1,
               PtE = PtE, parallel = TRUE, n_cores = 2)
# }
```
