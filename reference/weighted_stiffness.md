# Assemble a weighted stiffness matrix on a graph mesh

Uses a piecewise linear positive coefficient on each mesh edge. The
coefficient integral on an edge is its length times the mean of its two
endpoint values.

## Usage

``` r
weighted_stiffness(edges, lengths, weights)
```

## Arguments

- edges:

  Two-column integer matrix of mesh edge endpoints.

- lengths:

  Positive mesh edge lengths, in the same order as rows of edges.

- weights:

  Positive coefficient values at mesh nodes.

## Value

A symmetric sparse positive-semidefinite stiffness matrix.
