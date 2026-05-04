# Variancefor Whittle-Matérn fields

Computes the variance function for a Whittle-Matérn field. Warning is
not feasible for large graph due to matrix inversion

## Usage

``` r
spde_variance(
  kappa,
  tau,
  range,
  sigma,
  alpha,
  graph,
  BC = 1,
  include_vertices = FALSE,
  directional = F
)
```

## Arguments

- kappa:

  Parameter kappa from the SPDE.

- tau:

  Parameter tau from the SPDE.

- range:

  Range parameter.

- sigma:

  Standard deviation parameter.

- alpha:

  Smoothness parameter (1 or 2).

- graph:

  A `metric_graph` object.

- BC:

  boundary conditions

- include_vertices:

  Should the variance at the vertices locations be included in the
  returned vector?

- directional:

  bool is the model a directional or not. directional only works for
  alpha=1

## Value

Vector with the variance function evaluate at the mesh locations.

## Details

Compute the variance \\\rho(s_i,s_i)\\ where \\s_i\\ are all locations
in the mesh of the graph.
