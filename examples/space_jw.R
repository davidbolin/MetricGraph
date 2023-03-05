library(sp)
library(Matrix)
library(MetricGraph)


line1 <- Line(rbind(c(0,-4),c(1,0)))
line2 <- Line(rbind(c(1,0),c(1,1)))
line3 <- Line(rbind(c(1,1),c(2,1)))
lines = sp::SpatialLines(list(Lines(list(line1),ID="1")))

graph <- metric_graph$new(lines = lines)
graph$plot(direction = TRUE)
graph$build_mesh(h = 0.05)
graph$compute_fem()



###########
kappa <- 10
rho <- 0.1
sigma <- 1.0
dt <- 0.01
n <- dim(graph$mesh$C)[1]
R <- chol(graph$mesh$C)
L <- graph$mesh$G + kappa^2*graph$mesh$C + rho*graph$mesh$B

u0 <- rep(0,n)
u0[1] <- 1
T <- 20
U <- matrix(0,nrow=n,ncol = T)
U[,1] <- u0

for(i in 1:(T-1)){
  U[,i+1] <- as.vector(solve((graph$mesh$C + dt*L), graph$mesh$C%*%U[,i] + sqrt(dt)*sigma*R%*%rnorm(n)))
}
fig <- graph$plot_movie(U)
fig
