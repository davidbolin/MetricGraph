library(GPGraph)
library(sp)
#create a simple graph
line1 <- Line(rbind(c(0,0),c(1,0)))
line2 <- Line(rbind(c(0,0),c(0,1)))
line3 <- Line(rbind(c(0,1),c(-1,1)))
theta <- seq(from=pi,to=3*pi/2,length.out = 20)
line4 <- Line(cbind(sin(theta),1+ cos(theta)))
Lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                              Lines(list(line2),ID="2"),
                              Lines(list(line3),ID="3"),
                              Lines(list(line4),ID="4")))
graph <- metric_graph$new(Lines = Lines)
graph$plot()
theta <- c(2, 5)
P1 <-  #edge 1, 0.5 length in

C <- covariance.point.to.graph.exp(EP = c(1,0.1), theta = c(2,5), graph = graph, n.p = 50)
gg <- graph$plot_function(C, flat = FALSE)
gg <- graph$plot_function(C)

y <- c(1,2,3)
PtE <- matrix(c(1, 2, 3, 0.5, 0.5, 0.7),3,2)
graph$add_observations2(y, PtE)
p <- graph$plot(plot_data=TRUE)


line1 <- Line(rbind(c(0,0),c(10,0)))
line2 <- Line(rbind(c(0,0),c(0,1)))
line3 <- Line(rbind(c(0,1),c(-1,1)))
Lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                              Lines(list(line2),ID="2"),
                              Lines(list(line3),ID="3")))
graph <- metric_graph$new(Lines = Lines)
p <- graph$plot()

y <- c(1,2,3)
PtE <- matrix(c(1, 2, 3, 0.5, 0.5, 0.7),3,2)
graph$add_observations2(y, PtE)
p <- graph$plot(plot_data=TRUE)
