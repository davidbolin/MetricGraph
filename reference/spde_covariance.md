# Covariance function for Whittle-Matérn fields

Computes the covariance function for a Whittle-Matérn field.

## Usage

``` r
spde_covariance(P, kappa, tau, range, sigma, alpha, graph, directional = F)
```

## Arguments

- P:

  Location (edge number and normalized location on the edge) for the
  location to evaluate the covariance function at.

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

- directional:

  bool is the model a directional or not. directional only works for
  alpha=1

## Value

Vector with the covariance function evaluate at the mesh locations.

## Details

Compute the covariance function \\\rho(P,s_i)\\ where P is the provided
location and \\s_i\\ are all locations in the mesh of the graph.
