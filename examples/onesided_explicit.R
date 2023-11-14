library(ggplot2)
library(MetricGraph)
library(viridis)
library(Matrix)

graphics.off()
#simple example
edge1 <- rbind(c(0,0),c(1,0))
edge2 <- rbind(c(1,0),c(1,1))
edge3 <- rbind(c(1,0),c(2,0))
#edge4 <- rbind(c(2,0),c(2,1))
#edge5 <- rbind(c(2,0),c(3,0))
edges = list(edge1, edge2, edge3)#, edge4, edge5)
graph <- metric_graph$new(edges = edges)
graph$build_mesh(h=1, continuous = FALSE, continuous.outs = TRUE, continuous.deg2 = TRUE)
print(graph$plot(mesh=TRUE, direction = TRUE))
tau = 0.5
kappa <- 1
Q= MetricGraph:::Qalpha1_v2(c(tau,kappa), graph, w = 0 ,BC = 3, build = TRUE)
Q.sym= MetricGraph:::Qalpha1_v2(c(tau,kappa), graph, w = 0.5 ,BC = 1, build = TRUE)
Q_alpha2 <- MetricGraph:::Qalpha2(c(tau,kappa), graph, w = 0 ,BC = 3, build = TRUE)
# build graph with independence then create constraints matrix through the boundary conditions
Q_edges <- MetricGraph:::Qalpha1_edges(c(tau,kappa), graph, w = 0,BC=1, build=TRUE)


#add constraint that if a graph has multiple outwards they are the same node!
# and no inward
# also only one of them should be stationary!

#' Build K matrix,
#' assumes that there are no single vertex (i.e. no out nor in)
#' constraints are
#' u_output = \sum w_i u_input
#' the number of constraints equal the number output edges where indput vertix >0 x number of deriv
buildC = function(graph, alpha = 1) {

  V_indegree = graph$get_degrees("indegree")
  V_outdegree = graph$get_degrees("outdegree")
  index_outdegree <- V_outdegree > 0 & V_indegree >0
  index_in0      <- V_indegree == 0
  nC = sum(V_outdegree[index_outdegree] *(1 + V_indegree[index_outdegree]) + sum(V_outdegree[index_in0]-1)) * alpha
  i_  =  rep(0, nC)
  j_  =  rep(0, nC)
  x_  =  rep(0, nC)
  Vs <- which(index_outdegree)
  count_constraint <- 0
  count <- 0
  for (v in Vs) {
    out_edges   <- which(graph$E[, 1] %in% v)
    in_edges    <- which(graph$E[, 2] %in% v)
    #for each out edge
    n_in <- length(in_edges)
    for(i in 1:length(out_edges)){
      for(der in 1:alpha){
        i_[count + 1:(n_in+1)] <- count_constraint + 1
        j_[count + 1:(n_in+1)] <- c(2 * alpha * (out_edges[i]-1) + der,
                                    2 * alpha * (in_edges-1)  + alpha + der)

        x_[count + 1:(n_in+1)] <- c(-1,
                                    rep(1/n_in,n_in))
        count <- count + (n_in+1)
        count_constraint <- count_constraint + 1
      }
    }
  }
  Vs0 <- which(index_in0)
  for (v in Vs0) {
    out_edges   <- which(graph$E[, 1] %in% v)
    #for each out edge
    if(length(out_edges)>1){
      for(i in 2:length(out_edges)){
        for(der in 1:alpha){
          i_[count + 1:2] <- count_constraint + 1
          j_[count + 1:2] <- c(2 * alpha * (out_edges[i]-1) + der,
                               2 * alpha * (out_edges[i-1]-1)  + alpha + der)

          x_[count + 1:2] <- c(-1,
                                1)
          count <- count + 2
          count_constraint <- count_constraint + 1
        }
      }
    }
  }
  C <- Matrix::sparseMatrix(i = i_[1:count],
                            j = j_[1:count],
                            x = x_[1:count],
                            dims = c(count_constraint, 2*alpha*graph$nE))
  graph$C = C

  graph$CoB <- MetricGraph:::c_basis2(graph$C)
  graph$CoB$T <- t(graph$CoB$T)
  return(graph)
}
graph <- buildC(graph,alpha = 1)
n_const <- length(graph$CoB$S)
ind.const <- c(1:n_const)
Tc <- graph$CoB$T[-ind.const, ]
#add soft constraint
Q_mod <- forceSymmetric(Tc%*%Q_edges%*%t(Tc))

Sigma_E <- t(Tc)%*%solve(Q_mod)%*%(Tc)


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
graph$build_mesh(h=1, continuous = FALSE, continuous.outs = TRUE, continuous.deg2 = TRUE)
graph$plot(mesh = TRUE, direction = TRUE)

graph$compute_fem(petrov=TRUE)


Q_edges <- MetricGraph:::Qalpha1_edges(c(tau,kappa), graph, w = 0,BC=1, build=TRUE)

graph <- buildC(graph,alpha = 1)
n_const <- length(graph$CoB$S)
ind.const <- c(1:n_const)
Tc <- graph$CoB$T[-ind.const, ]
#add soft constraint
Q_mod <- forceSymmetric(Tc%*%Q_edges%*%t(Tc))

Sigma_E <- t(Tc)%*%solve(Q_mod)%*%(Tc)

hfull <- c(rep(1/(1+kappa), graph$mesh$n.bc),graph$mesh$h_e)
hfull[graph$mesh$h0] = 0

L <- kappa*graph$mesh$Cpet + graph$mesh$Gpet
Sigma <- solve(L,diag(hfull)%*%t(solve(L)))
p <- graph$plot(direction = TRUE)
p3 <- graph$plot_function(diag(Sigma), p = p)
print(p3)

