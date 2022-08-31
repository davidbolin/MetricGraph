#'
#' Building simple graph (line) simulating and estimating parameters
#'
#'
#rm(list=ls())
library(GPGraph)
library(Matrix)
library(sp)
set.seed(2)
graphics.off()
nt <- 200
kappa <- 5
sigma_e <- 0.01
sigma   <- 1
theta <-  c(sigma_e,kappa,sigma)
line.line <- Line(rbind(c(30,80),c(80,80)))
graph <- graph.obj$new(sp::SpatialLines(list(Lines(list(line.line),ID="1"))))
Q <- Q.exp(theta[2:3], graph$V, graph$EtV, graph$El)
R <- chol(Q)
X0 <- as.vector(solve(R, solve(R,rnorm(2), system = 'Lt')
                     , system = 'Pt'))
X <- sample.line.expontial(theta,
                           X0,
                           Line = graph$Lines[1,],
                           graph$El[1],
                           nt = nt)

X[,2] <- X[,2] + sigma_e*rnorm(nt)
points <- rgeos::gInterpolate(graph$Lines[1,], X[,1], normalized = T)
graph$add_observations2(y = X[,2], PtE = cbind(1,X[,1]), Spoints = points)

lik <- likelihood.exp.graph(theta,graph)
res <- optim(log(theta), function(x) -likelihood.exp.graph(exp(x),graph) )
print(exp(res$par))
plot(X[,1],X[,2])

X_ <-sample.line.expontial(theta,
                          X0,
                          Line = graph$Lines[1,],
                          graph$El[1],
                          nt = 200)
lines(X_[,1],X_[,2],col='green')
X_ <-sample.line.expontial(exp(res$par),
                           X0,
                           Line = graph$Lines[1,],
                           graph$El[1],
                           nt = 200)

lines(X_[,1],X_[,2],col='red')

d_ <- abs(X[,1]-X[1,1])
cov_true <- r_1(d_, theta[2:3])
cov_true[1] <- cov_true[1]  + theta[1]^2
cov_est <- r_1(d_, exp(res$par)[2:3])
cov_est[1] <- cov_est[1]  + exp(res$par)[1]^2
plot(d_, cov_true,type='l', ylim=c(0,max(cov_est,cov_true)))

lines(d_, cov_est,col='red')
