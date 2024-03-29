library(MetricGraph)
library(sp)

# Example 1: plot of graph, observations, and covariance for alpha = 1
line1 <- Line(rbind(c(0,0),c(1,0)))
line2 <- Line(rbind(c(0,0),c(0,1)))
line3 <- Line(rbind(c(0,1),c(-1,1)))
theta <- seq(from=pi,to=3*pi/2,length.out = 20)
line4 <- Line(cbind(sin(theta),1+ cos(theta)))
Lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                              Lines(list(line2),ID="2"),
                              Lines(list(line3),ID="3"),
                              Lines(list(line4),ID="4")))
graph <- metric_graph$new(lines = Lines)
graph$plot()

y <- c(1,2,3)
PtE <- matrix(c(1, 2, 3, 0.5, 0.5, 0.7),3,2)
graph$add_PtE_observations(y, PtE)
p <- graph$plot(data=TRUE)

kappa <- 10
sigma <- 2

graph$build_mesh(h = 0.1)
C <- covariance_alpha1_mesh(P = c(1,0.1), kappa = kappa, sigma = sigma, graph = graph)
graph$plot_function_mesh(C, plotly = TRUE)

C <- covariance_alpha1(P = c(1,0.1), kappa = kappa, sigma = sigma, graph = graph, n.p = 50)
gg <- graph$plot_function(C, plotly = FALSE)

gg <- graph$plot_function(C)


C <- covariance_alpha2(P = c(1,0.1), kappa = kappa, sigma = sigma, graph = graph, n.p = 50)
gg <- graph$plot_function(C, plotly = FALSE)


# Example 2: test with different graph
line1 <- Line(rbind(c(0,0),c(10,0)))
line2 <- Line(rbind(c(0,0),c(0,1)))
line3 <- Line(rbind(c(0,1),c(-1,1)))
Lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                              Lines(list(line2),ID="2"),
                              Lines(list(line3),ID="3")))
graph <- metric_graph$new(lines = Lines)
p <- graph$plot()

y <- c(1,2,3)
PtE <- matrix(c(1, 2, 3, 0.5, 0.5, 0.7),3,2)
graph$add_PtE_observations(y, PtE)
p <- graph$plot(data=TRUE)


# Example 3: test with alpha = 2
line1 <- Line(rbind(c(0,0),c(1,0)))
line2 <- Line(rbind(c(0,0),c(0,1)))
line3 <- Line(rbind(c(0,1),c(-1,1)))
theta <- seq(from=pi,to=3*pi/2,length.out = 20)
line4 <- Line(cbind(sin(theta),1+ cos(theta)))
Lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                              Lines(list(line2),ID="2"),
                              Lines(list(line3),ID="3"),
                              Lines(list(line4),ID="4")))
graph <- metric_graph$new(lines = Lines)
graph$plot()
kappa <- 10
sigma <- 2
C <- covariance_alpha2(P = c(1,0.1), kappa = kappa, sigma = sigma, graph = graph, n.p = 50)
gg <- graph$plot_function(C, plotly = FALSE)

library(Matrix)
#fix plot function
