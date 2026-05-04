# Space-time precision operator Euler discretization

The precision matrix for all vertices for space-time field

## Usage

``` r
make_Q_euler(graph, t, kappa, rho, gamma, alpha, beta, sigma, theta = 1)
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

- theta:

  Parameter theta for the Euler scheme.

## Value

Precision matrix.
