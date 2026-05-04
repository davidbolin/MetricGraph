# Metric graph SPDE result extraction from 'INLA' estimation results

Extract field and parameter values and distributions for a metric graph
spde effect from an 'INLA' result object.

## Usage

``` r
spde_metric_graph_result(
  inla,
  name,
  metric_graph_spde,
  compute.summary = TRUE,
  n_samples = 5000,
  n_density = 1024
)
```

## Arguments

- inla:

  An 'INLA' object obtained from a call to `inla()`.

- name:

  A character string with the name of the 'rSPDE' effect in the model.

- metric_graph_spde:

  The `inla_metric_graph_spde` object used for the random effect in the
  model.

- compute.summary:

  Should the summary be computed?

- n_samples:

  The number of samples to be used if parameterization is `matern`.

- n_density:

  The number of equally spaced points to estimate the density.

## Value

If the model was fitted with `matern` parameterization (the default), it
returns a list containing:

- marginals.range:

  Marginal densities for the range parameter.

- marginals.log.range:

  Marginal densities for log(range).

- marginals.sigma:

  Marginal densities for std. deviation.

- marginals.log.sigma:

  Marginal densities for log(std. deviation).

- marginals.values:

  Marginal densities for the field values.

- summary.log.range:

  Summary statistics for log(range).

- summary.log.sigma:

  Summary statistics for log(std. deviation).

- summary.values:

  Summary statistics for the field values.

If `compute.summary` is `TRUE`, then the list will also contain

- summary.kappa:

  Summary statistics for kappa.

- summary.tau:

  Summary statistics for tau.

If the model was fitted with the `spde` parameterization, it returns a
list containing:

- marginals.kappa:

  Marginal densities for kappa.

- marginals.log.kappa:

  Marginal densities for log(kappa).

- marginals.log.tau:

  Marginal densities for log(tau).

- marginals.tau:

  Marginal densities for tau.

- marginals.values:

  Marginal densities for the field values.

- summary.log.kappa:

  Summary statistics for log(kappa).

- summary.log.tau:

  Summary statistics for log(tau).

- summary.values:

  Summary statistics for the field values.

If `compute.summary` is `TRUE`, then the list will also contain

- summary.kappa:

  Summary statistics for kappa.

- summary.tau:

  Summary statistics for tau.
