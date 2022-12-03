rm(list=ls())
# test vertex vertex
line1 <- Line(rbind(c(0,0),c(1,0)))
line2 <- Line(rbind(c(0,0.03),c(0,0.06)))
line3 <- Line(rbind(c(0,0.09),c(0,0.97)))
line4 <- Line(cbind(seq(from=0,to=1, length.out = 100), rep(1,100)))
line5 <- Line(rbind(c(1,0.03),c(1,0.97)))
line6 <- Line(rbind(c(0.97,0.03),c(0.03,0.97)))

theta <- seq(from=0,to=2*pi,length.out = 50)
line7 <- Line(cbind(0.3*sin(theta),-0.33+0.3*cos(theta)))

lines = SpatialLines(list(Lines(list(line1),ID="1"),
                          Lines(list(line2),ID="2"),
                          Lines(list(line3),ID="3"),
                          Lines(list(line4),ID="4"),
                          Lines(list(line5),ID="5"),
                          Lines(list(line6),ID="6"),
                          Lines(list(line7),ID="7")))
graph <- metric_graph$new(lines, tolerance = list(vertex_vertex = 0.2))
graph$plot()

#test line vertex

line1 <- Line(rbind(c(0,0),c(1,0)))
line2 <- Line(rbind(c(0.5,0.01),c(0.5,0.97)))
line3 <- Line(cbind(seq(from=0,to=1, length.out = 100), rep(1,100)))
line4 <- Line(rbind(c(1,0.03),c(1,0.97)))
theta <- seq(from=0,to=2*pi,length.out = 50)
line5 <- Line(cbind(0.5+0.3*sin(theta),-0.31+0.3*cos(theta)))

lines = SpatialLines(list(Lines(list(line1),ID="1"),
                          Lines(list(line2),ID="2"),
                          Lines(list(line3),ID="3"),
                          Lines(list(line4),ID="4"),
                          Lines(list(line5),ID="5")))
graph <- metric_graph$new(lines, tolerance = list(vertex_vertex = 0.2,
                                                  vertex_line = 0.2))
graph$plot()

# test line line

#half circle line
line1 <- Line(rbind(c(1,0),c(2,0)))
line2 <- Line(rbind(c(1.5,-1),c(1.5,1)))
theta <- seq(from=pi,to=2*pi,length.out = 50)
line3 <- Line(cbind(1.05+sin(theta),cos(theta)))
line4 <- Line(rbind(c(0,1),c(0,-1)))

lines = SpatialLines(list(Lines(list(line1),ID="1"),
                          Lines(list(line2),ID="2"),
                          Lines(list(line3),ID="3"),
                          Lines(list(line4),ID="4")))
plot(lines)
graph <- metric_graph$new(lines, tolerance = list(vertex_vertex = 0.05,
                                                  vertex_line = 0.05,
                                                  line_line = 0.15))
graph$plot()


# circle line
line1 <- Line(rbind(c(-0.5,0),c(1,0)))
theta <- seq(from=0,to=2*pi,length.out = 50)
line2 <- Line(cbind(1.05+sin(theta),cos(theta)))
line3 <- Line(rbind(c(0,1),c(0,-1)))

lines = SpatialLines(list(Lines(list(line1),ID="1"),
                          Lines(list(line2),ID="2"),
                          Lines(list(line3),ID="3")))
plot(lines)
graph <- metric_graph$new(lines, tolerance = list(vertex_vertex = 0.2,
                                                  vertex_line = 0.2,
                                                  line_line = 0.2))
graph$plot()


#complicated with multiple circles
theta <- seq(from=0,to=2*pi,length.out = 50)
line1 <- Line(cbind(sin(theta),cos(theta)))
line2 <- Line(cbind(2.1+sin(theta),cos(theta)))
line3 <- Line(cbind(4.3+sin(theta),cos(theta)))
line4 <- Line(rbind(c(3.2,1),c(3.2,-1)))
line5 <- Line(rbind(c(-0.75,1),c(-0.75,-1)))
lines = SpatialLines(list(Lines(list(line1),ID="1"),
                          Lines(list(line2),ID="2"),
                          Lines(list(line3),ID="3"),
                          Lines(list(line4),ID="4"),
                          Lines(list(line5),ID="5")))
plot(lines)
graph <- metric_graph$new(lines, tolerance = list(line_line = 0.3))
graph$plot()
