# Precision matrix for Whittle-Matérn fields

Computes a sparse precision matrix for a Whittle-Matérn field. Order
three returns an unconstrained edge-endpoint precision. Fractional
orders return the unconstrained joint precision of the
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

  A `metric_graph` object.

- BC:

  Boundary-condition code, one of 0, 1, 2, or 3, passed unchanged to the
  alpha-one and alpha-two backends. `BC = 1` adds stationary endpoint
  corrections; `BC = 0` adds no additional corrections. Codes 2 and 3
  retain their backend-specific interpretation. Ignored for alpha 3 and
  fractional alpha: these return boundaryless raw precisions, with
  constraints imposed separately.

- build:

  If `TRUE`, the precision matrix is returned. Otherwise a list with
  entries `i`, `j`, `x`, and `dims` is returned. As for the integer
  cases, `build = FALSE` is not supported for `graph_components`
  objects.

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

Precision matrix or triplet list. The existing alpha-one and alpha-two
representations are unchanged.

## Details

For fractional alpha, the component order is `Q0`, `Q1`, ..., `Qm`. The
zeroth component has `graph$nV` vertex coordinates when alpha is below
one, and `2 * floor(alpha) * graph$nE` endpoint coordinates otherwise.
Each remaining component has `2 * ceiling(alpha) * graph$nE`
coordinates. Endpoint coordinates follow edge order, with the process
and its successive derivatives at the start followed by those at the end
of each edge. Parameter scaling and the boundaryless edge correction are
included, but no constraint basis or observation map is applied. For
alpha 3, the result is a single `6 * graph$nE` square raw matrix,
ordered as process value, first derivative, and second derivative at
each edge start and then end. No rational approximation is used for
alpha 3. The returned raw matrix need not be positive definite before
the graph constraints are imposed. Neither observations nor `graph$PtV`
are required, and `C`, `CoB`, and the constraint cache are left
unchanged. For `graph_components`, entire joint matrices are assembled
in component order, as in the integer cases.
