library(sp)
library(MetricGraph)

graph <- metric_graph$new(tolerance = list(vertex_vertex = 1e-1, vertex_line = 1e-3, line_line = 1e-3),
                          remove_deg2 = FALSE)
graph$plot(degree = TRUE)
graph$build_mesh(h = 0.01)
