library(sp)
library(MetricGraph)


line1 <- Line(rbind(c(0,0),c(1,0)))
line2 <- Line(rbind(c(0,0),c(0,1)))
line3 <- Line(rbind(c(0,1),c(-1,1)))
theta <- seq(from=pi,to=3*pi/2,length.out = 20)
line4 <- Line(cbind(sin(theta),1+ cos(theta)))
lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                              Lines(list(line2),ID="2"),
                              Lines(list(line3),ID="3"),
                              Lines(list(line4),ID="4")))

graph <- metric_graph$new(lines = lines)
graph$plot(degree = TRUE)

graph <- metric_graph$new(lines = lines, remove_deg2 = TRUE)
graph$plot(degree = TRUE)




lines <- logo_lines()
graph <- metric_graph$new(lines = lines, tolerance = list(vertex_vertex = 1e-2, vertex_line = 1e-2, line_line = 0),
                          remove_deg2 = TRUE, remove_circles = 0.1)
graph$plot(degree = TRUE)

