library(sp)
library(MetricGraph)

n <- 100
line1 <- Line(rbind(c(0,0),c(1,0)))
line2 <- Line(rbind(c(1,0),c(2,0)))
lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                              Lines(list(line2),ID="2")))
graph <- metric_graph$new(lines = lines,tolerance = list(vertex_vertex = 0, vertex_line = 0, line_line = 0),
                          remove_deg2 = TRUE)
graph$plot(degree = TRUE)

plot(graph$lines)
lines(graph$lines[1],col=2)
