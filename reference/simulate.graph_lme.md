# Simulation of models on metric graphs

The function samples a Gaussian random field based on a fitted model
using
[`graph_lme()`](https://davidbolin.github.io/MetricGraph/reference/graph_lme.md).

## Usage

``` r
# S3 method for class 'graph_lme'
simulate(
  object,
  nsim = 1,
  seed = NULL,
  sample_latent = FALSE,
  posterior = FALSE,
  which_repl = NULL,
  ...
)
```

## Arguments

- object:

  A `graph_lme` object

- nsim:

  The number of simulations.

- seed:

  an object specifying if and how the random number generator should be
  initialized ('seeded').

- sample_latent:

  If `FALSE`, samples for the response variable will be generated. If
  `TRUE`, samples for the latent model will be generated. The default is
  `FALSE`.

- posterior:

  Should posterior samples be generated? If `FALSE`, samples will be
  computed based on the estimated prior distribution. The default is
  `FALSE`.

- which_repl:

  Which replicates to generate the samples. If `NULL` samples will be
  generated for all replicates. Default is `NULL`.

- ...:

  Currently not used.

## Value

A list containing elements `samples`, `edge_number` and
`distance_on_edge`. Each of them is a list, whose indexes are the
replicates, and in `samples` a matrix is given with `nsim` columns, each
one being a sample. `edge_number` and `distance_on_edges` contain the
respective edge numbers and distances on edge for each sampled element.
The locations of the samples are the location of the data in which the
model was fitted.
