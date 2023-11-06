library(MetricGraph)
library(Matrix)
library(rSPDE)

#eigenfunctions of the Laplacian on the tadpole graph with lengths 1 and 2
tadpole.eig <- function(k,graph){
  x1 <- c(0,graph$get_edge_lengths()[1]*graph$mesh$PtE[graph$mesh$PtE[,1]==1,2])
  x2 <- c(0,graph$get_edge_lengths()[2]*graph$mesh$PtE[graph$mesh$PtE[,1]==2,2])

  if(k==0){
    f.e1 <- rep(1,length(x1))
    f.e2 <- rep(1,length(x2))
    f1 = c(f.e1[1],f.e2[1],f.e1[-1], f.e2[-1])
    f = list(phi=f1/sqrt(3))

  } else {
    f.e1 <- -2*sin(pi*k*1/2)*cos(pi*k*x1/2)
    f.e2 <- sin(pi*k*x2/2)

    f1 = c(f.e1[1],f.e2[1],f.e1[-1], f.e2[-1])

    if((k %% 2)==1){
      f = list(phi=f1/sqrt(3))
    } else {
      f.e1 <- (-1)^{k/2}*cos(pi*k*x1/2)
      f.e2 <- cos(pi*k*x2/2)
      f2 = c(f.e1[1],f.e2[1],f.e1[-1],f.e2[-1])
      f <- list(phi=f1,psi=f2/sqrt(3/2))
    }
  }

  return(f)
}

circ <- TRUE

if(circ){
  edge1 <- rbind(c(0,0),c(1,0))
  theta <- seq(from=-pi,to=pi,length.out = 100)
  edge2 <- cbind(1+1/pi+cos(theta)/pi,sin(theta)/pi)
  edges = list(edge1, edge2)
} else {
  edge1 <- rbind(c(1,0),c(0,0))
  edge2 <- rbind(c(0,1/(1+pi/4)),c(0,0))
  edge3 <- rbind(c(-1/(1+pi/4),1/(1+pi/4)),c(0,1/(1+pi/4)))
  theta <- seq(from=pi,to=3*pi/2,length.out = 50)
  edge4 <- cbind(sin(theta)/(1+pi/4),(1+ cos(theta))/(1+pi/4))
  edges = list(edge1, edge2, edge3, edge4)
}

graph <- metric_graph$new(edges = edges)
graph$prune_vertices()
graph$build_mesh(h=0.01)
graph$plot(mesh=TRUE)


kappa = 3
nu = 0.8
tau <- 1
alpha <- nu + 1/2
#check KL expansion
Sigma.kl <- matrix(0,nrow = dim(graph$mesh$V)[1],ncol = dim(graph$mesh$V)[1])
for(i in 0:100){
  phi <- tadpole.eig(i,graph)$phi
  Sigma.kl <- Sigma.kl + (1/(kappa^2 + (i*pi/2)^2)^(alpha))*phi%*%t(phi)
  if(i>0 && (i %% 2)==0){
    psi <- tadpole.eig(i,graph)$psi
    Sigma.kl <- Sigma.kl + (1/(kappa^2 + (i*pi/2)^2)^(alpha))*psi%*%t(psi)
  }

}
Sigma.kl <- Sigma.kl/tau^2
graph$compute_fem()

rspde.order <- 4
op <- matern.operators(alpha = nu+1/2, kappa = kappa, tau = tau,
                       m = rspde.order, graph = graph)

point.PtE <- c(2,0.1)
point <- graph$coordinates(PtE = point.PtE)
c_cov <- op$cov_function_mesh(matrix(point.PtE,1,2))

m1 <- which.min((graph$mesh$V[,1]-point[1])^2 + (graph$mesh$V[,2]-point[2])^2)
p <- graph$plot_function(Sigma.kl[m1,], plotly = TRUE)
graph$plot_function(c_cov,p=p, line_color = "red", plotly = TRUE)
