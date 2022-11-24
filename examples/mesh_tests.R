library(MetricGraph)
library(sp)
library(INLA)

line1 <- Line(rbind(c(0,0),c(1,0)))
line2 <- Line(rbind(c(0,0),c(0,1)))
line3 <- Line(rbind(c(0,1),c(-1,1)))
theta <- seq(from=pi,to=3*pi/2,length.out = 20)
line4 <- Line(cbind(sin(theta),1+ cos(theta)))
Lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                              Lines(list(line2),ID="2"),
                              Lines(list(line4),ID="3"),
                              Lines(list(line3),ID="4")))
graph <- metric_graph$new(lines = Lines)
graph$build_mesh(h = 0.01)
graph$plot(mesh=TRUE)
C <- covariance_alpha1_mesh(P = c(1,0.1), kappa = 10, sigma = 2, graph = graph)
X <- cbind(graph$mesh$VtE, C)
graph$plot_function(X, plotly = TRUE, vertex_size = 5)

C <- covariance_alpha2_mesh(P = c(1,0.1), kappa = 10, sigma = 2, graph = graph)
graph$plot_function_mesh(C, plotly = FALSE)

u <- sample_spde(kappa = 2/5, sigma = 2, graph = graph, type="mesh", alpha=1)
graph$plot_function_mesh(u, plotly = TRUE)
graph$plot_function_mesh(u)

u <- sample_spde(kappa = 10, sigma = 2, alpha = 2, graph = graph, type = "mesh")
graph$plot_function_mesh(u, plotly = FALSE)

#Test FEM
line1 <- Line(rbind(c(0,0),c(0,1)))
line2 <- Line(rbind(c(0,1),c(0,2.5)))
Lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                              Lines(list(line2),ID="2")))
graph <- metric_graph$new(lines = Lines)

graph$build_mesh(h = 1)
graph$plot(mesh=TRUE)
graph$compute_fem()

ind <- c(1,2,4,3)
mesh <- inla.mesh.1d(c(0,1,1.75, 2.5))
tmp <- inla.mesh.fem(mesh)

graph$mesh$G - tmp$g1[ind,ind]
graph$mesh$C - tmp$c1[ind,ind]


# test FEM A

line1 <- Line(rbind(c(0,0),c(1,0)))
line2 <- Line(rbind(c(0,0),c(0,1)))
line3 <- Line(rbind(c(0,1),c(-1,1)))
theta <- seq(from=pi,to=3*pi/2,length.out = 20)
line4 <- Line(cbind(sin(theta),1+ cos(theta)))
Lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                              Lines(list(line2),ID="2"),
                              Lines(list(line4),ID="3"),
                              Lines(list(line3),ID="4")))
graph <- metric_graph$new(lines = Lines)
graph$build_mesh(h = 0.2)
graph$plot(mesh=TRUE)

PtE <- rbind(c(1,0.99),c(4,0.7))

A <- graph$mesh_A(PtE)

graph$plot_function_mesh(A[1,],plotly=FALSE)
