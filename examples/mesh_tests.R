library(GPGraph)
library(sp)

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
graph$build_mesh(h = 0.01)
graph$plot(mesh=TRUE)
C <- covariance_alpha1_mesh(P = c(1,0.1), kappa = 10, sigma = 2, graph = graph)
graph$plot_function_mesh(C, flat = FALSE)

C <- covariance_alpha2_mesh(P = c(1,0.1), kappa = 10, sigma = 2, graph = graph)
graph$plot_function_mesh(C, flat = FALSE)

u <- sample_spde_mesh(kappa = 10, sigma = 2, graph = graph)
graph$plot_function_mesh(u, flat = FALSE)
graph$plot_function_mesh(u)

u <- sample_spde_mesh(kappa = 10, sigma = 2, alpha = 2, graph = graph)
graph$plot_function_mesh(u, flat = FALSE)
