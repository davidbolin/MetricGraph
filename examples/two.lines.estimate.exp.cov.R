#'
#' Building simple graph (two lines) simulating and estimating parameters
#'
#'
#rm(list=ls())
library(GPGraph)
library(Matrix)
library(sp)
set.seed(2)
graphics.off()
nt <- 10
kappa <- 0.1
sigma_e <- 0.1
sigma   <- 2
theta <-  c(sigma_e,kappa,sigma)
line.line <- Line(rbind(c(30,80),c(120,80)))
line.line2 <- Line(rbind(c(30,00),c(30,80)))

graph <-  graph.obj$new(sp::SpatialLines(list(Lines(list(line.line),ID="1"),
                                               Lines(list(line.line2),ID="2"))))
Q <- Q.exp(theta[2:3], graph$V, graph$EtV, graph$El)
R <- Cholesky(Q,LDL = FALSE, perm = TRUE)
V0 <- as.vector(solve(R, solve(R,rnorm(3), system = 'Lt')
                      , system = 'Pt'))
X <- c()
for(i in 1:length(graph$El)){
  X <- rbind(X,cbind(sample.line.expontial(theta,
                             V0[graph$EtV[i,2:3]],
                             Line = graph$Lines[i,],
                             graph$El[i],
                             nt = nt),i))
}

X[,2] <- X[,2] + sigma_e*rnorm(nt)


graph$add_PtE_observations(y = X[,2], PtE = X[,c(3,1)])

lik <- likelihood.exp.graph(theta,graph)
res <- optim(log(theta), function(x) -likelihood.exp.graph(exp(x),graph) )
print(exp(res$par))

ind <- X[,3]==1
d_ <- abs(X[ind,1]-X[1,1])
cov_true <- r_1(d_, theta[2:3])
cov_true[1] <- cov_true[1]  + theta[1]^2
cov_est <- r_1(d_, exp(res$par)[2:3])
cov_est[1] <- cov_est[1]  + exp(res$par)[1]^2
plot(d_, cov_true,type='l', ylim=c(0,max(cov_est,cov_true)))

lines(d_, cov_est,col='red')

i <- 1
plot(X[X[,3]==i,1],X[X[,3]==i,2])

X_ <-sample.line.expontial(theta,
                           V0[graph$EtV[i,2:3]],
                           Line = graph$Lines[i,],
                           graph$El[i],
                           nt = nt,
                           sample=F,
                           py = graph$PtE[graph$PtE[,1]==i,2],
                           y = graph$y[graph$PtE[,1]==i])
lines(X_[,1],X_[,2],col='green')
X_ <-sample.line.expontial(exp(res$par),
                           V0[graph$EtV[i,2:3]],
                           Line = graph$Lines[i,],
                           graph$El[i],
                           nt = nt,
                           py = graph$PtE[graph$PtE[,1]==i,2],
                           y = graph$y[graph$PtE[,1]==i],
                           sample=F)

lines(X_[,1],X_[,2],col='red')
