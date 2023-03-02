library(sp)
library(Matrix)
library(MetricGraph)
library(htmlwidgets)

save.plot <- FALSE

line1 <- Line(rbind(c(1,0),c(0,0)))
line2 <- Line(rbind(c(0,0),c(0,1)))
line3 <- Line(rbind(c(0,0),c(0,-1)))
line4 <- Line(rbind(c(0,1),c(1,0)))
lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                              Lines(list(line2),ID="2"),
                              Lines(list(line3),ID="3"),
                              Lines(list(line4),ID="4")))

graph <- metric_graph$new(lines = lines)
graph$plot(direction = TRUE)
h = 0.01
graph$build_mesh(h = h)
graph$compute_fem()

make.Q <- function(L0,L,dt,C, T) {
  n <- dim(L)[1]
  Lbar = C + dt*L
  CC = kronecker(Diagonal(n=T), dt*C)
  LL <- Matrix(0,nrow = T*n, ncol = T*n)
  LL[1:n,1:n] <- L0
  for(t in 1:(T-1)) {
    LL[(t*n+1):((t+1)*n),(t*n+1):((t+1)*n)] <- Lbar
    LL[(t*n+1):((t+1)*n),((t-1)*n+1):(t*n)] <- -C
  }
  return(list(L = LL, C = CC))
}

kappa <- 15
rho <- 0
sigma <- 75
dt <- 0.2
n <- dim(graph$mesh$C)[1]
L <- graph$mesh$G + kappa^2*graph$mesh$C + rho*graph$mesh$B
L0 <- graph$mesh$C
T = 200
Q <- make.Q(L0, L, dt, graph$mesh$C, T)

k <- 50
v1 <- rep(0,dim(Q$L)[1]); v1[k] <- 1
v2 <- rep(0,dim(Q$L)[1]); v2[n+k] <- 1
v3 <- rep(0,dim(Q$L)[1]); v3[floor(T/2)*n+k] <- 1
v4 <- rep(0,dim(Q$L)[1]); v4[(T-1)*n+k] <- 1

c1 <- sigma^2*solve(Q$L,Q$C%*%solve(t(Q$L),v1))[1:n]
c2 <- sigma^2*solve(Q$L,Q$C%*%solve(t(Q$L),v2))[(n+1):(2*n)]
c3 <- sigma^2*solve(Q$L,Q$C%*%solve(t(Q$L),v3))[(floor(T/2)*n+1):((floor(T/2)+1)*n)]
c4 <- sigma^2*solve(Q$L,Q$C%*%solve(t(Q$L),v4))[((T-1)*n+1):(T*n)]

L <- graph$mesh$G + kappa^2*graph$mesh$C + rho*graph$mesh$B

v1 <- rep(0,dim(L)[1]); v1[k] <- 1
ct <- sigma^2*solve(L,graph$mesh$C%*%solve(t(L),v1))/dt
p <- graph$plot_function(as.vector(c1), plotly = TRUE, support_width = 0, line_color = "red")
p <- graph$plot_function(as.vector(c2), plotly = TRUE, p = p, support_width = 0, line_color = "green")
p <- graph$plot_function(as.vector(c3), plotly = TRUE, support_width = 0, line_color = "blue", line_width = 2)
p <- graph$plot_function(as.vector(c4), plotly = TRUE, p = p,support_width = 0, line_color = "gray", line_width = 2)
p <- graph$plot_function(as.vector(ct*max(c4)/max(ct)), plotly = TRUE, p = p, support_width = 1, line_color = "black")
p
cat(max(ct)/max(c4))





###########
kappa <- 10
rho <- 0
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
  U[,i+1] <- as.vector(solve((C + dt*L), C%*%U[,i] + sqrt(dt)*sigma*R%*%rnorm(n)))
}
fig <- graph$plot_movie(U)
fig
if(save.plot){
  htmlwidgets::saveWidget(fig, file = "spacetime.HTML", selfcontained = TRUE)
}



