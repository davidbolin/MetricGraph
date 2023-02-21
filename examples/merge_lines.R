library(sp)
library(MetricGraph)

graph <- metric_graph$new(lines = logo_lines(),tolerance = list(vertex_vertex = 1e-1, vertex_line = 1e-3, line_line = 1e-3),
                          remove_deg2 = TRUE)
graph$plot(degree = TRUE)

