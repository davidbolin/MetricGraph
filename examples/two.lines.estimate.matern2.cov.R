#'
#' Building simple graph (two lines) simulating and estimating parameters
#' for matern alpha=2
#'
#rm(list=ls())
library(GPGraph)
library(Matrix)
library(sp)
set.seed(15)
graphics.off()
nt <- 300
kappa <- 0.2
sigma_e <- 0.1
sigma   <- 0.1
theta <-  c(sigma_e,kappa,sigma)
line.line2 <- Line(rbind(c(30,80),c(140,80)))
line.line <- Line(rbind(c(30,00),c(30,80)))

graph <-  graph.obj$new(sp::SpatialLines(list(Lines(list(line.line),ID="1"),
                                              Lines(list(line.line2),ID="2"))))
Q <- Q.matern2(theta[2:3], graph$V, graph$EtV, graph$El, BC = 1)
graph$buildA(2, F)
Qtilde <- (graph$CBobj$T)%*%Q%*%t(graph$CBobj$T)
Qtilde <- Qtilde[-c(1:2),-c(1:2)]
R <- Cholesky(Qtilde,LDL = FALSE, perm = TRUE)
V0 <- as.vector(solve(R, solve(R,rnorm(6), system = 'Lt')
                      , system = 'Pt'))
print(round(t(graph$CBobj$T[-c(1:2),])%*%solve(Qtilde)%*%(graph$CBobj$T[-c(1:2),]),2))
u_e <- t(graph$CBobj$T)%*%c(0,0,V0)
X <- c()
for(i in 1:length(graph$El)){
  X <- rbind(X,cbind(sample.line.matern2(theta,
                                           u_e[4*(i-1) +1:4],
                                           Line = graph$Lines[i,],
                                           graph$El[i],
                                           nt = nt),i))
}
X[,2] <- X[,2] #+ sigma_e*rnorm(nt)


graph$add_observations2(y = X[,2], PtE = X[,c(3,1)])

lik <- likelihood.matern2.graph(theta,graph)
res <- optim(log(theta), function(x) -likelihood.matern2.graph(exp(x),graph) )
print(exp(res$par))
lik2 <- likelihood.matern2.graph(exp(res$par),graph)
M <- posterior.mean.matern2(exp(res$par),graph)
r_p <- 4.743859/theta[2]
t <- seq(0, r_p,length.out = 100)
plot(t, r_2(t,theta[2:3]),type='l', ylim=c(0,max(r_2(0,theta[2:3]),
                                                 r_2(0,exp(res$par[2:3])
                                                 ))))
lines(t, r_2(t,exp(res$par[2:3])),col='red')
plot(c((1:nt),(nt):(2*nt-1)),X[,2],type='l')
lines(X[X[,3]==1,2],col='red')
lines((nt):(2*nt-1),X[X[,3]==2,2],col='blue')
#plot(c(X1[,1],X2[,1]+l_e),c(X1[,2],X2[,2]),type='l')
