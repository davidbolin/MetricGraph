#'
#' Building simple graph (two lines) simulating and estimating parameters
#' for matern alpha=2
#'
#rm(list=ls())
library(GPGraph)
library(Matrix)
library(sp)
set.seed(13)
graphics.off()
nt <- 100
kappa <- 0.3
sigma_e <- 0.1
sigma   <- 1
theta <-  c(sigma_e,kappa,sigma)
line.line2 <- Line(rbind(c(30,80),c(140,80)))
line.line <- Line(rbind(c(30,00),c(30,80)))
line.line3 <- Line(rbind(c(30,80),c(30,-20)))

graph <-  graph.obj$new(sp::SpatialLines(list(Lines(list(line.line),ID="1"),
                                              Lines(list(line.line2),ID="2"),
                                              Lines(list(line.line3),ID="3"))))
Q <- Q.matern2(theta[2:3], graph$V, graph$EtV, graph$El, BC = 1)
graph$buildA(2, F)
n.c <- length(graph$CBobj$S)
Qmod <- (graph$CBobj$T)%*%Q%*%t(graph$CBobj$T)
Qtilde <- Qmod
Qtilde <- Qtilde[-c(1:n.c),-c(1:n.c)]
R <- Cholesky(Qtilde,LDL = FALSE, perm = TRUE)
V0 <- as.vector(solve(R, solve(R,rnorm(dim(Q)[1]-n.c), system = 'Lt')
                      , system = 'Pt'))
print(round(t(graph$CBobj$T[-c(1:n.c),])%*%solve(Qtilde)%*%(graph$CBobj$T[-c(1:n.c),]),2))
u_e <- t(graph$CBobj$T)%*%c(rep(0,n.c),V0)
X <- c()
for(i in 1:length(graph$El)){
  X <- rbind(X,cbind(sample.line.matern2(theta,
                                           u_e[4*(i-1) +1:4],
                                           Line = graph$Lines[i,],
                                           graph$El[i],
                                           nt = nt),i))
}
X[,2] <- X[,2] + sigma_e*rnorm(nt)
X1 <- X[X[,3]==1,]
X2 <- X[X[,3]==2,]
X3 <- X[X[,3]==3,]
par(mfrow=c(3,1))
plot(c((1:nt),(nt):(2*nt-1)),c(X1[,2],X2[,2]),type='l')
lines(X1[,2],col='red')
lines((nt):(2*nt-1),X2[,2],col='blue')
plot(c((1:nt),(nt):(2*nt-1)),c(X1[,2],X3[,2]),type='l')
lines(X1[,2],col='red')
lines((nt):(2*nt-1),X3[,2],col='blue')
plot(c((1:nt),(nt):(2*nt-1)),c(X3[nt:1,2],X2[,2]),type='l')
lines(X3[nt:1,2],col='red')
lines((nt):(2*nt-1),X2[,2],col='blue')

print(X1[(nt-1):nt,])
print(X2[1:2,])
graph$add_PtE_observations(y = X[,2], PtE = X[,c(3,1)])
graph$buildA(2, F)
lik <- likelihood.matern2.graph(theta,graph)
res <- optim(log(theta), function(x) -likelihood.matern2.graph(exp(x),graph) )
print(exp(res$par))
graph$observation_to_vertex()
graph$buildA(2, F)
lik2 <-likelihood.graph.covariance(theta, graph, alpha=2)
res2 <- optim(log(theta), function(x) -likelihood.graph.covariance(exp(x),graph, 2) )
print(exp(res2$par))
lik2 <- likelihood.matern2.graph(exp(res2$par),graph)
M <- posterior.mean.matern2(exp(res$par),graph)
r_p <- 4.743859/theta[2]
t <- seq(0, r_p,length.out = 100)
plot(t, r_2(t,theta[2:3]),type='l', ylim=c(0,max(r_2(0,theta[2:3]),
                                                 r_2(0,exp(res$par[2:3])
                                                 ))))
lines(t, r_2(t,exp(res$par[2:3])),col='red')
points(t, r_2(t,exp(res2$par[2:3])),col='blue')

#plot(c(X1[,1],X2[,1]+l_e),c(X1[,2],X2[,2]),type='l')
