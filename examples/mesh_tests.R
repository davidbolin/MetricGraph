library(GPGraph)
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
graph <- metric_graph$new(Lines = Lines)
graph$build_mesh(h = 1)
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

#Test FEM
line1 <- Line(rbind(c(0,0),c(0,1)))
line2 <- Line(rbind(c(0,1),c(0,2.5)))
Lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                              Lines(list(line2),ID="2")))
graph <- metric_graph$new(Lines = Lines)

graph$build_mesh(h = 1)
graph$plot(mesh=TRUE)
graph$compute_fem()

ind <- c(1,2,4,3)
mesh <- inla.mesh.1d(c(0,1,1.75, 2.5))
tmp <- inla.mesh.fem(mesh)

graph$mesh$G - tmp$g1[ind,ind]
graph$mesh$C - tmp$c1[ind,ind]

#' convert positions (edge, position) from graph$V to mesh$V
PtE_to_mesh <- function(graph,PtE){
  #First create joint PtE for the mesh (that includes the vertices)
  VtE_mesh <- matrix(0,dim(graph$mesh$V)[1],2)
  for(i in 1:dim(graph$V)[1]) {
    tmp <- which(graph$E == i,arr.ind = TRUE)[1,]
    colnames(tmp) <- NULL
    VtE_mesh[i,1] <- tmp[1] #edge number
    VtE_mesh[i,1] <- ifelse(tmp[2]==1,0,1) #start or end of edge?
  }
  VtE_mesh <- rbind(VtE_mesh, graph$mesh$PtE)

  #now go through and update PtE
  n <- dim(PtE)[1]
  PtE_update <- PtE
  for(i in 1:n){
    ind <- which(VtE_mesh[,1]==PtE[i,1])
    v2 <- which((PtE[i,2]-VtE_mesh[ind,2])>0)[1] #upper edge
    v1 <- which((VtE_mesh[ind,2]-PtE[i,2])>0)[1] #lower edge
    ei <- which(graph$mesh$E[,1]==v1 && graph$mesh$E[,2]==v2)
    d <- (PtE[i,2]-VtE_mesh[v1,2])/(VtE_mesh[v2,2]-VtE_mesh[v1,2])
    PtE_update[i,] <- c(ei,d)
  }
  return(PtE_update)
}

#' x is a matrix with locations (edge, position) on mesh.
#' The function assumes straight edges between mesh nodes.
mesh_A <- function(mesh,x){
  n <- dim(x)[1]
  A <- Matrix(0,nrow=n,ncol=dim(mesh$V)[1])
  p.u <- unique(mesh$PtE[,1])
  for(i in 1:n){
    A[i,mesh$V[mesh$E[x[i],1]]] = 1-x[i,2]
    A[i,mesh$V[mesh$E[x[i],2]]] = x[i,2]
  }
  return(A)
}
