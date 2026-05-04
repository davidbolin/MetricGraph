# Precision matrix for Whittle-Matérn fields

Computes the precision matrix for all vertices for a Whittle-Matérn
field.

## Usage

``` r
spde_precision(kappa, tau, alpha, graph, BC = 1, build = TRUE)
```

## Arguments

- kappa:

  Range parameter.

- tau:

  Precision parameter.

- alpha:

  Smoothness parameter (1 or 2).

- graph:

  A `metric_graph` object.

- BC:

  Set boundary conditions for degree=1 vertices. BC =0 gives Neumann
  boundary conditions and BC=1 gives stationary boundary conditions.

- build:

  If `TRUE`, the precision matrix is returned. Otherwise a list
  list(i,j,x, nv) is returned.

## Value

Precision matrix or list.
