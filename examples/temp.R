nt <- 30
kappa <- 0.1
sigma_e <- 0.1
sigma   <- 0.1
theta <-  c(sigma_e,kappa,sigma)
line.line2 <- Line(rbind(c(30,80),c(140,80)))
line.line <- Line(rbind(c(30,00),c(30,80)))
line.line3 <- Line(rbind(c(30,-10),c(30,80)))

graph <-  graph.obj$new(sp::SpatialLines(list(Lines(list(line.line),ID="1"),
                                              Lines(list(line.line2),ID="2"),
                                              Lines(list(line.line2),ID="3"))))
Q <- Q.matern2(theta[2:3], graph$V, graph$EtV, graph$El, BC = 1)
graph$buildA(2, F)
Qtilde <- (graph$CBobj$T)%*%Q%*%t(graph$CBobj$T)
nc <- 1:length(graph$CBobj$S)
Qtilde <- Qtilde[-nc,-nc]
R <- Cholesky(Qtilde,LDL = FALSE, perm = TRUE)
V0 <- as.vector(solve(R, solve(R,rnorm(7), system = 'Lt')
                      , system = 'Pt'))
#print(round(t(graph$CBobj$T[-nc,])%*%solve(Qtilde)%*%(graph$CBobj$T[-nc,]),2))
u_e <- t(graph$CBobj$T)%*%c(rep(0,length(graph$CBobj$S)),V0)
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
graph$observation_to_vertex()
graph$buildA(2, F)
n.c <- 1:length(graph$CBobj$S)
Q <- Q.matern2(c(kappa,sigma), graph$V, graph$EtV, graph$El, BC = 1)
Qtilde <- (graph$CBobj$T)%*%Q%*%t(graph$CBobj$T)
Qtilde <- Qtilde[-n.c,-n.c]
Sigma.overdetermined  = t(graph$CBobj$T[-n.c,])%*%solve(Qtilde)%*%(graph$CBobj$T[-n.c,])
index.obs <-  4*(graph$PtE[,1]-1) + (1 * (abs(graph$PtE[,2])<1e-14)) + (3 * (abs(graph$PtE[,2])>1e-14))
Sigma <-  as.matrix(Sigma.overdetermined[index.obs, index.obs])
