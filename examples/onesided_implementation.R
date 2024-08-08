library(ggplot2)
library(viridis)

#'
#' @param theta (sigma_e, reciprocal_tau, kappa)
#'
Sigma.vertices.dalpha1 <- function(theta, graph, parameterization="matern"){


  sigma_e <- exp(theta[1])
  if(parameterization == "matern"){
    kappa = sqrt(8 * 0.5) / exp(theta[3])
  } else{
    kappa = exp(theta[3])
  }

  reciprocal_tau <- exp(theta[2])
  tau <- 1/reciprocal_tau
  if(is.null(graph$C)){
    graph$buildDirectionalConstraints(1)
  } else if(graph$CoB$alpha == 2){
    graph$buildDirectionalConstraints(1)
  }

  Q_edges <- MetricGraph:::Qalpha1_edges(c(tau,kappa), graph, w = 0,BC=1, build=TRUE)
  n_const <- length(graph$CoB$S)
  ind.const <- c(1:n_const)
  Tc <- graph$CoB$T[-ind.const, ]
  Q_T <- Matrix::forceSymmetric(Tc%*%Q_edges%*%t(Tc))
  Sigma.overdetermined <- as.matrix(t(Tc)%*%solve(Q_T)%*%Tc)
  PtE = graph$get_PtE()
  index.obs <- 2*(PtE[, 1] - 1) + (PtE[, 2]> 1-0.001) + 1
  Sigma <-  as.matrix(Sigma.overdetermined[index.obs, index.obs])
  Sigma.o <- Sigma
  diag(Sigma.o) <- diag(Sigma.o) + sigma_e^2
  return(list(Sigma = Sigma,Sigma.o =  Sigma.overdetermined, Q= Q_edges))
}

#simple example
edge2 <- rbind(c(2,1), c(1,0))
edge3 <- rbind(c(2,-1), c(1,0))
edge4 <- rbind(c(3,0),c(2,-1))
edge7 <- rbind(c(3,0), c(2,1))
edges = list(edge2, edge3, edge4, edge7)


graph <- metric_graph$new(edges = edges)
graph$build_mesh(h=0.1, continuous = FALSE, continuous.outs = TRUE, continuous.deg2 = TRUE)
graph$plot(mesh = TRUE, direction = TRUE)
v <- rep(0, dim(graph$mesh$V)[1]); v[4] = 1; graph$plot_function(v, plotly = TRUE)

graph$compute_fem(petrov = TRUE)
kappa <- 1
hfull <- c(rep(1/(1+kappa), graph$mesh$n.bc),graph$mesh$h_e)
hfull[graph$mesh$h0] = 0

L <- kappa*graph$mesh$Cpet + graph$mesh$Gpet
Sigma <- sqrt(8 * 0.5)*solve(L,diag(hfull)%*%t(solve(L)))
ind_ <- (graph$get_degrees("indegree") == 1 & graph$get_degrees("outdegree") == 1) | graph$get_degrees()==1
Sigma_V <- Sigma[1:graph$nV,1:graph$nV]
Sigma_V <- Sigma_V[ind_, ind_]
I <- 1:graph$nV
data <- data.frame(edge_number      = graph$mesh$PtE[I[ind_],1],
                   distance_on_edge = graph$mesh$PtE[I[ind_],2],
                   y                = NA)
graph$add_observations(data = data, normalized = TRUE)
#graph$observation_to_vertex()
Sigma.vert <- Sigma.vertices.dalpha1(c(-Inf, log(2)/2, 2*log(kappa)),graph,"")$Sigma
print(Sigma_V - Sigma.vert)

p <- graph$plot(direction = TRUE)
graph$plot_function(diag(Sigma), p = p)

#rooted tree
edge1 <- rbind(c(0,0),c(1,0))
edge2 <- rbind(c(1,0),c(2,1))
edge3 <- rbind(c(1,0),c(2,-1))
edge4 <- rbind(c(2,-1),c(3,0))
edge5 <- rbind(c(2,-1),c(3,-2))
edges = list(edge1, edge2, edge3, edge4, edge5)
graph <- metric_graph$new(edges = edges)
graph$build_mesh(h=0.1, continuous = FALSE, continuous.outs = TRUE, continuous.deg2 = TRUE)
graph$plot(mesh=TRUE, direction = TRUE)

graph$compute_fem(petrov = TRUE)
kappa <- 1
hfull <- c(rep(1/(1+kappa), graph$mesh$n.bc),graph$mesh$h_e)
hfull[graph$mesh$h0] = 0

L <- kappa*graph$mesh$Cpet + graph$mesh$Gpet
Sigma <- sqrt(8 * 0.5)*solve(L,diag(hfull)%*%t(solve(L)))
ind_ <- (graph$get_degrees("indegree") == 1 & graph$get_degrees("outdegree") == 1) | graph$get_degrees()==1
Sigma_V <- Sigma[1:graph$nV,1:graph$nV]
Sigma_V <- Sigma_V[ind_, ind_]
I <- 1:graph$nV
data <- data.frame(edge_number      = graph$mesh$PtE[I[ind_],1],
                   distance_on_edge = graph$mesh$PtE[I[ind_],2],
                   y                = NA)
graph$add_observations(data = data, normalized = TRUE)
#graph$observation_to_vertex()
Sigma.vert <- Sigma.vertices.dalpha1(c(-Inf, log(2)/2, 2*log(kappa)),graph,"")$Sigma
print(Sigma_V - Sigma.vert)

p <- graph$plot(direction = TRUE)
p1 <- graph$plot_function(diag(Sigma), p = p)
p1

W <- rnorm(n=length(hfull),mean=0,sd = sqrt(hfull))
u <- solve(L,W)
graph$plot_function(u,plotly=TRUE)


#non-rooted tree
edge1 <- rbind(c(1,0), c(0,0))
edge2 <- rbind(c(2,1), c(1,0))
edge3 <- rbind(c(2,-1), c(1,0))
edge4 <- rbind(c(3,0),c(2,-1))
edge5 <- rbind(c(3,-2),c(2,-1))
edge6 <- rbind(c(2,-1), c(1,-2))
edges = list(edge1, edge2, edge3, edge4, edge5)

graph <- metric_graph$new(edges = edges)

graph$build_mesh(h=0.1, continuous = FALSE, continuous.outs = TRUE, continuous.deg2 = TRUE)
graph$plot(mesh = TRUE, direction = TRUE)

graph$compute_fem(petrov=TRUE)

kappa <- 1
hfull <- c(rep(1/(1+kappa), graph$mesh$n.bc),graph$mesh$h_e)
hfull[graph$mesh$h0] = 0

L <- kappa*graph$mesh$Cpet + graph$mesh$Gpet
Sigma <- sqrt(8*0.5)*solve(L,diag(hfull)%*%t(solve(L)))
ind_ <- (graph$get_degrees("indegree") == 1 & graph$get_degrees("outdegree") == 1) | graph$get_degrees()==1
Sigma_V <- Sigma[1:graph$nV,1:graph$nV]
Sigma_V <- Sigma_V[ind_, ind_]
I <- 1:graph$nV
data <- data.frame(edge_number      = graph$mesh$PtE[I[ind_],1],
                   distance_on_edge = graph$mesh$PtE[I[ind_],2],
                   y                = NA)
graph$add_observations(data = data, normalized = TRUE)
#graph$observation_to_vertex()
Sigma.vert <- Sigma.vertices.dalpha1(c(-Inf, log(2)/2, 2*log(kappa)),graph,"")$Sigma
print(Sigma_V - Sigma.vert)
p2 <- graph$plot_function(diag(Sigma), p = p)
p2

W <- rnorm(n=length(hfull),mean=0,sd = sqrt(hfull))
u <- solve(L,W)
graph$plot_function(u,plotly=TRUE)


# graph that is not a tree
edge1 <- rbind(c(1,0), c(0,0))
edge2 <- rbind(c(2,1), c(1,0))
edge3 <- rbind(c(2,-1), c(1,0))
edge4 <- rbind(c(3,0),c(2,-1))
edge5 <- rbind(c(3,-2),c(2,-1))
edge6 <- rbind(c(2,-1), c(1,-2))
edge7 <- rbind(c(3,0), c(2,1))
edges = list(edge1, edge2, edge3, edge4, edge5, edge6, edge7)


graph <- metric_graph$new(edges = edges)
graph$build_mesh(h=0.1, continuous = FALSE, continuous.outs = TRUE, continuous.deg2 = TRUE)
graph$plot(mesh = TRUE, direction = TRUE)

graph$compute_fem(petrov=TRUE)

kappa <- 1
hfull <- c(rep(1/(1+kappa), graph$mesh$n.bc),graph$mesh$h_e)
hfull[graph$mesh$h0] = 0

L <- kappa*graph$mesh$Cpet + graph$mesh$Gpet
Sigma <- sqrt(8*0.5)*solve(L,diag(hfull)%*%t(solve(L)))
ind_ <- (graph$get_degrees("indegree") == 1 & graph$get_degrees("outdegree") == 1) | graph$get_degrees()==1
Sigma_V <- Sigma[1:graph$nV,1:graph$nV]
Sigma_V <- Sigma_V[ind_, ind_]
I <- 1:graph$nV
data <- data.frame(edge_number      = graph$mesh$PtE[I[ind_],1],
                   distance_on_edge = graph$mesh$PtE[I[ind_],2],
                   y                = NA)
graph$add_observations(data = data, normalized = TRUE)
#graph$observation_to_vertex()
Sigma.vert <- Sigma.vertices.dalpha1(c(-Inf, log(2)/2, 2*log(kappa)),graph,"")$Sigma
print(Sigma_V - Sigma.vert)
p <- graph$plot(direction = TRUE)
p3 <- graph$plot_function(diag(Sigma), p = p)
p3

library("gridExtra")
p1 <- p1 + theme(legend.position = "none") + scale_colour_gradientn(colours = viridis(100), limits=c(0.25, 0.51))
p2 <- p2 + theme(legend.position = "none") + scale_colour_gradientn(colours = viridis(100), limits=c(0.25, 0.51))
p3 <- p3 + scale_colour_gradientn(colours = viridis(100), limits=c(0.25, 0.51))
grid.arrange(p1, p2, p3, ncol = 3, nrow = 1, widths = c(1,1,1.16))



## plot fem
edge1 <- rbind(c(1,0), c(0,0))
edge2 <- rbind(c(2,1), c(1,0))
edge3 <- rbind(c(2,-1), c(1,0))
edges = list(edge1, edge2, edge3)
graph <- metric_graph$new(edges = edges)
graph$build_mesh(h=0.5, continuous = FALSE, continuous.outs = FALSE, continuous.deg2 = TRUE)
p1 <- graph$plot(mesh = TRUE, direction = TRUE)
p1
u <- rep(0,dim(graph$mesh$V)[1]); u[1] = 1
p2 <- graph$plot_function(u, plotly = TRUE, line_color = "blue", line_width = 2)
u <- rep(0,dim(graph$mesh$V)[1]); u[8] = 1
p2 <- graph$plot_function(u, plotly = TRUE, line_color = "red", p = p2, line_width = 2)
u <- rep(0,dim(graph$mesh$V)[1]); u[11] = 1
p2 <- graph$plot_function(u, plotly = TRUE, line_color = "green", p = p2, line_width = 2)

for(i in setdiff(1:11, c(1,8,11))) {
  u <- rep(0,dim(graph$mesh$V)[1]); u[i] = 1
  p2 <- graph$plot_function(u, plotly = TRUE, line_color = "gray", p = p2, line_width = 2)
}
p2 <- plotly::layout(p2, scene = list(camera = list(eye = list(x = -3, y = 1.25, z = 1.25))))
p2


build_mesh(graph,h=0.5, continuous = TRUE)
u <- rep(0,dim(graph$mesh$V)[1]); u[1] = 1
p3 <- graph$plot_function(u, plotly = TRUE, line_color = "blue", line_width = 2)
for(i in 2:9) {
  u <- rep(0,dim(graph$mesh$V)[1]); u[i] = 1
  p3 <- graph$plot_function(u, plotly = TRUE, line_color = "gray", p = p3, line_width = 2)
}
p3 <- plotly::layout(p3, scene = list(camera = list(eye = list(x = -3, y = 1.25, z = 1.25))))
p3
