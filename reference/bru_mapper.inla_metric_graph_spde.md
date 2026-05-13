# Metric graph 'inlabru' mapper

Metric graph 'inlabru' mapper

## Usage

``` r
# S3 method for class 'inla_metric_graph_spde'
bru_get_mapper(model, ...)

# S3 method for class 'bru_mapper_inla_metric_graph_spde'
ibm_n(mapper, ...)

# S3 method for class 'bru_mapper_inla_metric_graph_spde'
ibm_values(mapper, ...)

# S3 method for class 'bru_mapper_inla_metric_graph_spde'
ibm_jacobian(mapper, input, ...)
```

## Arguments

- model:

  An `inla_metric_graph_spde` for which to construct or extract a mapper

- ...:

  Arguments passed on to other methods

- mapper:

  A `bru_mapper.inla_metric_graph_spde` object

- input:

  The values for which to produce a mapping matrix
