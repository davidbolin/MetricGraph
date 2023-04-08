library(sp)
library(Matrix)
library(MetricGraph)
library(htmlwidgets)
library(plotly)


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

nt = 200
T = 0.5
t <- seq(from=0, to = T, length.out = nt)
kappa <- 10
rho <- 40
gamma <- 0.05
sigma <- 2
alpha <- 1
beta <- 0

#show spacetime covariances
Q <- make.Q.spacetime(graph,t,kappa,rho,gamma, alpha, beta, sigma)
ct.dir <- MetricGraph:::plot.spacetime.covariances(graph,Q = Q,
                                                   t.ind = c(nt/2),
                                                   s.ind = 235,
                                                   t = t)


#Simulate field
nt = 30
T = 0.2
kappa <- 0.1
rho <- -100
gamma <- 0.05
sigma <- 0.001
alpha <- 1
beta <- 0
t <- seq(from=0, to = T, length.out = nt)
u0 <- rep(0,dim(graph$mesh$C)[1])
u0[50] <- 1
U <- MetricGraph:::simulate.spacetime(graph, t, kappa, rho, gamma, alpha, beta, sigma, u0)
fig <- graph$plot_movie(U)
fig

if(save.plot){
  htmlwidgets::saveWidget(fig, file = "spacetime.HTML", selfcontained = TRUE)
}



## TESTS
if(alpha==1){
  Qs <- 2*gamma*make.L(beta+alpha,kappa,graph$mesh$C,graph$mesh$G)/sigma^2
}else if (alpha == 2) {
  C <- graph$mesh$C
  G <- graph$mesh$G
  Cd <- Diagonal(rowSums(C),n=n)
  Ci <- Diagonal(1/rowSums(C),n=n)
  Qs <- 2*gamma*(kappa^4*Cd + (2*kappa^2+rho^2)*G + G%*%Ci%*%G)/sigma^2
}

n <- dim(graph$mesh$C)[1]
Cb <- Diagonal(0,n=dim(graph$mesh$C)[1])
diag(Cb)[1:2] <- graph$get_degrees()
Cd <- Diagonal(rowSums(graph$mesh$C),n=n)
Cd[1,1]<- Cd[1,1]
Cd[n,n] <- Cd[n,n]
Ci <- Diagonal(1/diag(Cd),n=n)
B <- Diagonal(0,n=n)
B[1,1]<- 1
B[2,2] <- 1
Qs <- 2*gamma*(graph$mesh$G - 0.00*rho*B + kappa^2*graph$mesh$C)/sigma^2
ct.dir <- plot.covariances(graph,Q = Q/sigma^2,Qs = Qs,
                           t.ind = c(nt/2), s.ind = 5,t.shift = c(0),t = t)







kappa <- 20
rho <- -10
gamma <- 0.01
sigma <- 50
alpha <- 1
beta <- 0
nt = 200
T = 0.25
t <- seq(from=0, to = T, length.out = nt)

Q <- make.Q.euler(graph, t, kappa, rho, gamma, alpha, beta,sigma,theta=1)
ct.eul <- plot.covariances(graph,Q, Qs = 2*gamma*make.L(alpha+beta,kappa,graph$mesh$C,
                                                        graph$mesh$G)/sigma^2,
                 t.ind = c(1, nt/2, nt-2), s.ind = 50,t.shift = 0,t = t)



Q <- make.Q.direct(graph,t,kappa,rho,gamma)
ct.dir <- plot.spacetime.covariances(graph,Q = Q/sigma^2,Qs = 2*gamma*L0/sigma^2,
                                     t.ind = c(12,nt/2,nt-2), s.ind = 50,t.shift = c(0),t = t)

###########




