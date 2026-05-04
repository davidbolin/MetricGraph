# Samples a Whittle-Matérn field on a metric graph

Obtains samples of a Whittle-Matérn field on a metric graph.

## Usage

``` r
sample_spde(
  kappa,
  tau,
  range,
  sigma,
  sigma_e = 0,
  alpha = 1,
  directional = FALSE,
  graph,
  PtE = NULL,
  type = "manual",
  posterior = FALSE,
  nsim = 1,
  method = c("conditional", "Q"),
  BC = 1
)
```

## Arguments

- kappa:

  Range parameter.

- tau:

  Precision parameter.

- range:

  Practical correlation range parameter.

- sigma:

  Marginal standard deviation parameter.

- sigma_e:

  Standard deviation of the measurement noise.

- alpha:

  Smoothness parameter.

- directional:

  should we use directional model currently only for alpha=1

- graph:

  A `metric_graph` object.

- PtE:

  Matrix with locations (edge, normalized distance on edge) where the
  samples should be generated.

- type:

  If "manual" is set, then sampling is done at the locations specified
  in `PtE`. Set to "mesh" for simulation at mesh nodes, and to "obs" for
  simulation at observation locations.

- posterior:

  Sample conditionally on the observations?

- nsim:

  Number of samples to be generated.

- method:

  Which method to use for the sampling? The options are "conditional"
  and "Q". Here, "Q" is more stable but takes longer.

- BC:

  Boundary conditions for degree 1 vertices. BC = 0 gives Neumann
  boundary conditions and BC = 1 gives stationary boundary conditions.

## Value

Matrix or vector with the samples.

## Details

Samples a Gaussian Whittle-Matérn field on a metric graph, either from
the prior or conditionally on observations \$\$y_i = u(t_i) + \sigma_e
e_i\$\$ on the graph, where \\e_i\\ are independent standard Gaussian
variables. The parameters for the field can either be specified in terms
of tau and kappa or practical correlation range and marginal standard
deviation.
