# MetricGraph (development version)
* Changed argument `data` to `newdata` in `predict` methods. The argument `data` was deprecated.
* Bugfixes on sample_spde and when adding observations based on Euclidean positions.
* Added options `vertex_unit` and `length_unit` on graph creations. Units are given in edge lengths with the `get_edge_lengths()` method.
* Added a method to check if the graph is a tree.
* The graph construction was thoroughly refactored. The resulting construction is faster and cleaner.
* The graph constructions now accepts list of coordinates (where the coordinates are given as either matrices or data frames), `SpatialLines`, `SpatialLinesDataFrames` or `MULTILINESTRING`.
* Adding two options (`sf` package or `sp` package) for handling `longlat` by using the `which_longlat` option.
* Adding `crs` (if using `sf`) and `proj4string` (if using `sp`) for handling general coordinate reference systems.
* Moving `data` to the `private` environment.
* Several data manipulation helper tools and methods were introduced, together with a vignette with a brief tutorial on these tools.
* The method `mesh_A()` has been deprecated, use `fem_basis()` instead.
* Several quality of life improvements.
* Improved the `plot()` method with the option `plotly=TRUE`.
* Improved the `plot_function()` method to accept `data` and `newdata`.
* Included a `process_data()` method for metric graphs

# MetricGraph 1.1.2
* Adjusts to ensure compatibility with future releases of the Matrix package.

# MetricGraph 1.1.1
* Adjusts on documentation for CRAN.

# MetricGraph 1.1.0
* Improved the documentation.
* Reorganized some functions.

# MetricGraph 1.0.0
* First version of the package.
