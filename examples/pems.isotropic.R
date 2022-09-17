#'
#'
#'

# Load library
library('sf')
library(GPGraph)
# Load shapefile
Lines <- read_sf('data.pems/lines.shp')
V <- read.csv('data.pems/V.csv',header=T, row.names = NULL)
EtV <- read.csv('data.pems/E.csv',header=T, row.names = NULL)
PtE <- read.csv('data.pems/PtE.csv',header=T, row.names = NULL)
PtE[,1] <- PtE[,1] + 1
Y_loc <- read.csv('data.pems/graph_sensor_locations_bay.csv',header=F)
Y <- read.csv('data.pems/Y.csv',header=T, row.names = NULL)
#clean and merge lines
Y <- as.matrix(Y[,-1])
TE <- table(c(EtV$V1,EtV$V2))
#length(names(TE)[TE1=2]) 320 namn

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
y_L <- rep(0, length(graph$y))
for(i in 1:length(graph$PtV)){
  mu.p_i <- Sigma[,graph$PtV[-i]]%*%solve(Sigma.o[-i,-i],graph$y[-i])
  y_L[i] <- mu.p_i[graph$PtV[i]]
}
cat('var res= ',var(graph$y-y_L),'\n')
fig <- plot_obs(graph,y_L, y_loc = Y_loc) + scale_colour_gradientn(colours = heat.colors(10))
print(fig)

theta.exp <- c(9.4726902736, 0.0001559032, 0.3745814561)
res <- optim(log(theta.exp), function(x) -likelihood.exp.graph.stupid(exp(x),graph) )
theta.exp <- exp(res$par)
Y_1_leave2 <- posterior.leave.stupid(theta.exp,graph)
cat('var exp= ',var(Y_1_leave2-graph$y),'\n')
Q <- Q.exp(theta.exp[2:3], graph$V, graph$EtV, graph$El)
Sigma <- as.matrix(solve(Q))[graph$PtV,graph$PtV]
SigmaO <- Sigma
diag(SigmaO) <- diag(SigmaO)  +  theta.exp[1]^2
i <- 1
fig <- plot_obs(graph,Sigma[i,]/Sigma[i,i], y_loc = Y_loc) + scale_colour_gradientn(colours = heat.colors(10))
print(fig)
fig <- plot_obs(graph,Sigma.i[i,]/Sigma.i[i,i], y_loc = Y_loc) + scale_colour_gradientn(colours = heat.colors(10))
print(fig)

Eige <- eigen(Sigma.o)
fig <- plot_obs(graph,Eige$vectors[,4], y_loc = Y_loc) + scale_colour_gradientn(colours = heat.colors(10))
print(fig)
