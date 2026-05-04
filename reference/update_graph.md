# Update an older version metric graph to the current package version

This function creates a new metric_graph object from an existing one,
preserving all data, edge weights, mesh, distances, and other components
while ensuring compatibility with the current package version. The new
graph is created with minimal changes by disabling merges and other
transformations.

## Usage

``` r
update_graph(old_graph, verbose = TRUE)
```

## Arguments

- old_graph:

  A metric_graph object from an older version of the package

- verbose:

  Logical. Print progress messages. Default is TRUE.

## Value

A new metric_graph object compatible with the current package version

## Details

This function is useful when loading metric graphs created with older
versions of the MetricGraph package. It ensures that all components
(vertices, edges, data, mesh, distances, etc.) are preserved while
updating the internal structure to be compatible with the current
package version.

The function works by:

- Extracting all components from the old graph

- Creating a new graph with settings that minimize structural changes

- Restoring all computed components (mesh, distances, etc.)

- Re-adding observation data if present

## Examples

``` r
if (FALSE) { # \dontrun{
# Load an old graph object
old_graph <- readRDS("old_graph.rds")

# Update to current version
new_graph <- update_graph(old_graph)

# Or use the method directly
new_graph <- old_graph$update_graph()
} # }
```
