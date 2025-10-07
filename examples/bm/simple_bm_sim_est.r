# very simple brownian moition
sigma2 <- 1
sigma2.eps <- 0.1
edge1 <- rbind(c(0,-1),c(0,0))
edge2 <- rbind(c(0,0),c(0,1))
edge3 <- rbind(c(0,1),c(0,2))
edge4 <- rbind(c(0,2),c(0,3))
edges = list(edge1, edge2, edge3,edge4)
graph <- metric_graph$new(edges = edges)
print(graph$plot())
graph$build_mesh(h = 0.1)
graph$mesh$edge_lengths <- graph$mesh$h
graph$mesh$nE <- dim(graph$mesh$E)[1]
graph$mesh$nV <- dim(graph$mesh$V)[1]
Q <- MetricGraph:::Qrandomwalk(1/sigma2, graph$mesh)
# simulate from the intrinsic model using KL , sums to zero
eQ <- eigen(Q)
d <- sum(eQ$values > 1e-10)
u <- rep(0, graph$mesh$nV)
for(i in 1:d){
  u <- u + rnorm(1,0,1/sqrt(eQ$values[i])) * eQ$vectors[,i]
}

print(graph$plot_function(X = u, type = "plotly"))

#remove observations on the vertices for simplisity
y <- u[-c(1:graph$nV)] + sqrt(sigma2.eps)*rnorm(length(u)-graph$nV)

df_data <- data.frame(y = y, edge_number = graph$mesh$PtE[,1],
                      distance_on_edge = graph$mesh$PtE[,2])
graph$add_observations(data = df_data, normalized = TRUE)

negLogLik <- function(theta){
  return(-MetricGraph:::likelihood_randomwalk(c(theta[1],theta[2]), graph, data_name = 'y'))
}
res <- optim(c(0,0),negLogLik)
cat('sigma.eps = ',exp(res$par[1]),' true = ',sqrt(sigma2.eps),'\n')
cat('sigma = ',exp(0.5*res$par[2]),' true = ',sqrt(sigma2),'\n')
cat('diff(u) = ',var(diff(u[-c(1:graph$nV)] ))/0.01,' true = ',sigma2,'\n')
cat('diff(y) = ',var(diff(y))/0.01,'\n')
PtE.obs <- graph$mesh$PtE
meanU <- posterior_mean_random_walk( c(sqrt(sigma2.eps), 1/sigma2),
                                       graph,
                                       df_data$y,
                                     PtE.obs)
cat('E[U|Y] = ',meanU,'\n')
cat('U = ',u[1:graph$nV],'\n')

graph$build_mesh(h = 0.01)
EU_Y <- MetricGraph:::posterior_mean_obs_random_walk( c(sqrt(sigma2.eps), 1/sigma2),
                                           graph,
                                           df_data$y,
                                           PtE.obs,
                                        graph$mesh$PtE,
                                           type = "PtE")


p <- graph$plot_function(X = EU_Y, type = "plotly")
print(graph$plot(data = "y", p=p, type = "plotly"))
