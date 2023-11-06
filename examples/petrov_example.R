edge1 <- rbind(c(0,0),c(1,0))
edge2 <- rbind(c(1,0),c(2,0))
edge3 <- rbind(c(2,0),c(3,0))
edge4 <- rbind(c(1,0),c(1,1))
edge5 <- rbind(c(2,0),c(2,1))
edges = list(edge1, edge2, edge3, edge4, edge5)
graph <- metric_graph$new(edges = edges)
graph$build_mesh(h=0.05)
graph$plot(mesh=TRUE)
graph$compute_fem(petrov=TRUE)

#set starting value
u.start <- 1
kappa <- 1
L = kappa*graph$mesh$Cpet + graph$mesh$Gpet
L <- rbind(sparseMatrix(i=1,j=1,x=1,dims=c(1,dim(graph$mesh$Cpet)[1])),
           t(L))
W <- c(u.start,rnorm(n=length(graph$mesh$h_e),mean=0,sd = sqrt(graph$mesh$h_e)))
u <- solve(L,W)
graph$plot_function(u,plotly=TRUE)
