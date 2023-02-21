library(sp)
library(MetricGraph)

n <- 100
line1 <- Line(rbind(c(8,4),c(8,8)))
theta <- seq(from=pi,to=3*pi/2,length.out = n)
line2 <- Line(cbind(10+2*sin(theta),6+2*cos(theta)))
theta <- seq(from=3*pi/2,to=2*pi,length.out = n)
line3 <- Line(cbind(10+2*sin(theta),6+2*cos(theta)))
theta <- seq(from=0,to=pi,length.out = n)
line6 <- Line(cbind(7+sin(theta),7+cos(theta)))
lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                              Lines(list(line2),ID="2"),
                              Lines(list(line3),ID="3"),
                              Lines(list(line6),ID="4")))
plot(lines)
graph <- metric_graph$new(lines = lines,tolerance = list(vertex_vertex = 1e-2, vertex_line = 1e-3, line_line = 1e-3))
graph$plot(degree = TRUE)

plot(graph$lines)
lines(graph$lines[1],col=2)
