

nt <- 4
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
graph$add_observations2(y = X[,2], PtE = X[,c(3,1)])
graph$observation_to_vertex()

lik <- likelihood.exp.graph(theta,graph)
#lik_diff <- lik - likelihood.exp.graph(theta+c(1,0,0),graph)


lik.v2 <- likelihood.exp.graph.v2(theta,graph)
#lik_diff.v2 <- lik.v2 - likelihood.exp.graph.v2(theta+c(1,0,0),graph)
#cat('lik_diff- lik_diff.v2 =',lik_diff-lik_diff.v2,'\n')
