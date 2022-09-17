library('sf')
library(GPGraph)

# Load data
Lines <- read_sf('data.pems/lines.shp')
V <- read.csv('data.pems/V.csv',header=T, row.names = NULL)
EtV <- read.csv('data.pems/E.csv',header=T, row.names = NULL)
PtE <- read.csv('data.pems/PtE.csv',header=T, row.names = NULL)
PtE[,1] <- PtE[,1] + 1
Y_loc <- read.csv('data.pems/graph_sensor_locations_bay.csv',header=F)
Y <- read.csv('data.pems/Y.csv',header=T, row.names = NULL)
Y <- as.matrix(Y[,-1])
TE <- table(c(EtV$V1,EtV$V2))

# Create graph object
graph <-  graph.obj$new()
graph$Lines <- as_Spatial(Lines)
graph$V <- as.matrix(V)
graph$El <- as.matrix(EtV[,4])
graph$EtV <- as.matrix(EtV[,1:3])
graph$PtE <- as.matrix(PtE)
graph$y <- Y[1,]#-colMeans(Y)#as.matrix(Y[,-1])[1,]
graph$y <- graph$y - mean(graph$y ) #temporary
fig <- plot_obs(graph,graph$y, y_loc = Y_loc) + scale_colour_gradientn(colours = heat.colors(10))
print(fig)
graph$observation_to_vertex()


# Fit isotropic model
graph$compute_resdist()
like.exp <- function(theta,graph){
  Sigma <- theta[1] ^2*exp(-theta[2]*graph$res.dist[graph$PtV,graph$PtV]) +
    theta[3]^2*diag(rep(1,length(graph$y)))
  R <- chol(Sigma)
  return(-sum(log(diag(R))) - 0.5*t(graph$y)%*%solve(Sigma,graph$y))
}
theta <- c(1,1,1)
res <- optim(log(theta), function(x) -like.exp(exp(x),graph))
theta <- exp(res$par)
Sigma <- theta[1]^2*exp(-theta[2]*graph$res.dist)
Sigma.i <- theta[1]^2*exp(-theta[2]*graph$res.dist[graph$PtV,graph$PtV])
Sigma.o <- theta[1]^2*exp(-theta[2]*graph$res.dist[graph$PtV,graph$PtV]) + theta[3]^2*diag(rep(1,length(graph$y)))

mu.p <- Sigma[,graph$PtV]%*%solve(Sigma.o,graph$y)
y_f <- mu.p[graph$PtV]
y_L_iso <- rep(0, length(graph$y))
for(i in 1:length(graph$PtV)){
  mu.p_i <- Sigma[,graph$PtV[-i]]%*%solve(Sigma.o[-i,-i],graph$y[-i])
  y_L_iso[i] <- mu.p_i[graph$PtV[i]]
}
cat('var res= ',var(graph$y-y_L_iso),'\n')
fig <- plot_obs(graph,y_L_iso, y_loc = Y_loc) + scale_colour_gradientn(colours = heat.colors(10))
print(fig)

# Fit alpha=1 model
theta.exp <- c(9.4726902736, 0.0001559032, 0.3745814561)
res <- optim(log(theta.exp), function(x) -likelihood.exp.graph.v2(exp(x),graph) )
theta.exp <- exp(res$par)
Y_L_exp <- posterior.leave.v2(theta.exp,graph)
cat('var exp= ',var(Y_L_exp-graph$y),'\n')
fig <- plot_obs(graph,Y_L_exp, y_loc = Y_loc) + scale_colour_gradientn(colours = heat.colors(10))
print(fig)

# Fit graph Laplace model
graph$compute_laplacian()
theta.graph <- c(9.478923, 7.154005, 1.160011e-05)
res <- optim(log(theta.graph), function(x) -likelihood.graph_laplacian(exp(x),graph) )
theta.graph <- exp(res$par)
Y_L_graph <- posterior.leave.graph_laplacian(theta.graph,graph)
cat('var graph= ',var(Y_L_graph-graph$y),'\n')
fig <- plot_obs(graph,Y_L_graph, y_loc = Y_loc) + scale_colour_gradientn(colours = heat.colors(10))
print(fig)
