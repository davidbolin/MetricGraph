# Match Data Frame Rows to Graph Mesh Order

Reorders the rows of a data frame to align with the specific ordering of
points along the edges of a metric graph's mesh. This ensures that data
associated with locations on the graph are correctly aligned with the
graph's internal mesh structure, which is essential for spatial
operations and visualizations.

## Usage

``` r
match_mesh_data(
  graph,
  data,
  edge_col = ".edge_number",
  dist_col = ".distance_on_edge"
)
```

## Arguments

- graph:

  A metric graph object with a mesh component containing VtE.

- data:

  A data.frame to be reordered. It must contain columns for edge numbers
  and distances along those edges.

- edge_col:

  Character. Name of the column in `data` that contains the edge
  numbers. Default: ".edge_number".

- dist_col:

  Character. Name of the column in `data` that contains the distance
  along the edge for each point. Default: ".distance_on_edge".

## Value

A reordered data.frame where rows correspond to the order of points in
`graph$mesh$VtE`. If the input `data` was an `sf` object, the returned
object will also be an `sf` object with its geometry reordered
accordingly.
