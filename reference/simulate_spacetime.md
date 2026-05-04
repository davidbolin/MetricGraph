# space-time simulation based on implicit Euler discretization in time

Simulation with starting value u0

## Usage

``` r
simulate_spacetime(graph, t, kappa, rho, gamma, alpha, beta, sigma, u0, BC = 0)
```

## Arguments

- graph:

  A `metric_graph` object.

- t:

  Vector of time points.

- kappa:

  Spatial range parameter.

- rho:

  Drift parameter.

- gamma:

  Temporal range parameter.

- alpha:

  Smoothness parameter (integer) for spatial operator.

- beta:

  Smoothness parameter (integer) for Q-Wiener process.

- sigma:

  Variance parameter.

- u0:

  Starting value.

- BC:

  Which boundary condition to use (0,1). Here, 0 is no adjustment on the
  boundary and 1 results in making the boundary condition stationary.

## Value

Precision matrix.
