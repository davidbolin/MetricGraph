# Precision matrix for Whittle-Matérn fields

Computes a sparse precision matrix for a Whittle-Matérn field. Order
three returns an unconstrained edge-endpoint precision. Fractional
orders return a list of unconstrained precisions for the
rational-approximation components.

## Usage

``` r
spde_precision(
  kappa,
  tau,
  alpha,
  graph,
  BC = 1,
  build = TRUE,
  m = 2,
  type_rational_approx = "chebfun",
  type_interp = "spline"
)
```

## Arguments

- kappa:

  Positive range parameter; zero is also supported for alpha 1.

- tau:

  Positive precision parameter.

- alpha:

  Smoothness parameter, strictly greater than 0.5 and at most 3.

- graph:

  A connected or disconnected `metric_graph` object. The deprecated
  `graph_components` class is retained for compatibility.

- BC:

  Boundary-condition code, one of 0, 1, 2, or 3, passed unchanged to the
  alpha-one and alpha-two backends. `BC = 1` adds stationary endpoint
  corrections; `BC = 0` adds no additional corrections. Codes 2 and 3
  retain their backend-specific interpretation. Ignored for alpha 3 and
  fractional alpha: these return boundaryless raw precisions, with
  constraints imposed separately.

- build:

  If `TRUE`, sparse precision matrices are returned; otherwise triplet
  lists with entries `i`, `j`, `x`, and `dims` are returned. For
  fractional alpha, either representation is returned separately for
  each rational component, including on a disconnected `metric_graph`.
  `build = FALSE` is not supported for the deprecated `graph_components`
  class.

- m:

  Rational approximation order for fractional alpha. A positive integer;
  ignored for alpha 1, 2, and 3. Defaults to 2.

- type_rational_approx:

  Rational approximation type passed to rSPDE; one of `"chebfun"`,
  `"brasil"`, and `"chebfunLB"`. Used only for fractional alpha.
  Defaults to `"chebfun"`.

- type_interp:

  Coefficient interpolation method passed to rSPDE; one of `"linear"`,
  `"log"`, `"spline"`, and `"logspline"`. Used only for fractional
  alpha. Defaults to `"spline"`.

## Value

For integer alpha, a sparse precision matrix or a triplet list,
retaining the existing representation. For fractional alpha, a named
list of component matrices or triplet lists. For fractional
`graph_components`, a named list of components, each containing a named
list of graph matrices.

## Details

For fractional alpha, the result is a named list `Q0`, `Q1`, ..., `Qm`,
containing one sparse matrix (`build = TRUE`) or one triplet list
(`build = FALSE`) per rational component. This structure is the same for
connected and disconnected `metric_graph` objects: each `Qi` covers the
entire graph in its original vertex or edge-endpoint coordinate order.
Disconnected pieces are not extracted or reordered with
`get_components()`. To assemble the joint raw precision, use
`Matrix::bdiag(spde_precision(..., build = TRUE))`. Its coordinates
match constraint templates built on that same full graph, as in the
alpha-two workflow; disconnected graphs require no additional coordinate
permutation. The zeroth component has `graph$nV` vertex coordinates when
alpha is below one, and `2 * floor(alpha) * graph$nE` endpoint
coordinates otherwise. Each remaining component has
`2 * ceiling(alpha) * graph$nE` coordinates. Endpoint coordinates follow
edge order, with the process and its successive derivatives at the start
followed by those at the end of each edge. Parameter scaling and the
boundaryless edge correction are included, but no constraint basis or
observation map is applied. For alpha 3, the result is a single
`6 * graph$nE` square raw matrix, ordered as process value, first
derivative, and second derivative at each edge start and then end. No
rational approximation is used for alpha 3. The returned raw matrix need
not be positive definite before the graph constraints are imposed.
Neither observations nor `graph$PtV` are required, and `C`, `CoB`, and
the constraint cache are left unchanged. For compatibility with the
deprecated `graph_components` class only, fractional output has an
additional graph level: each `Qi` is a named list `G1`, `G2`, ..., in
`graph$graphs` order, containing one sparse matrix per graph. This
nested structure does not apply to a disconnected `metric_graph`. The
deprecated path performs no block-diagonal assembly across graphs or
rational components.
