library(sp)
library(Matrix)
library(MetricGraph)


line1 <- Line(rbind(c(0,0),c(1,0)))
line2 <- Line(rbind(c(1,0),c(1,1)))
lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                              Lines(list(line2),ID="2")))

graph <- metric_graph$new(lines = lines)
graph$plot(direction = TRUE)
graph$build_mesh(h = 0.01)
graph$compute_fem()

kappa <- 10
rho <- 1
n <- dim(graph$mesh$C)[1]
C <- graph$mesh$C
C <- Diagonal(rowSums(C),n = n)
L <- graph$mesh$G + kappa^2*C + rho*graph$mesh$B
Q <- t(L)%*%solve(C, L)

A <- graph$mesh_A(matrix(c(1,0.5),1,2))

r <- solve(Q,t(A))
vars <- diag(solve(Q))
graph$plot_function(r,plotly = TRUE)

dt <- 0.1
I <- Diagonal(n,1)

u0 <- rep(0,n)
u0[1] <- 1
T <- 10
U <- matrix(0,nrow=n,ncol = T)
U[,1] <- u0

for(i in 1:(T-1)){
  U[,i+1] <- as.vector((I - dt*L)%*%U[,i] )
}

graph$plot_movie(U)