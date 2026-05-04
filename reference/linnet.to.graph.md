# Convert a `linnet` object to a metric graph object

This function converts a `linnet` object (from the `spatstat` package)
into a metric graph object.

## Usage

``` r
linnet.to.graph(linnet.object, crs, ...)
```

## Arguments

- linnet.object:

  A `linnet` object to be converted.

- crs:

  The coordinate reference system of the graph.

- ...:

  Additional arguments to be passed to the `metric_graph` constructor.

## Value

A metric graph object with edges defined by the network.
