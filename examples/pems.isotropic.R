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

graph$observation_to_vertex()

graph$compute_resdist()

like.exp <- function(theta,graph){
  Sigma <- theta[1]*exp(-theta[2]*graph$res.dist[graph$PtV,graph$PtV]) +
    theta[3]*diag(rep(1,length(graph$y)))
  R <- chol(Sigma)
  return(-sum(log(diag(R))) - 0.5*t(graph$y)%*%solve(Sigma,graph$y))
}
theta <- c(1,200,1)
res <- optim(log(theta), function(x) -like.exp(exp(x),graph))

theta <- exp(res$par)

Sigma <- theta[1]*exp(-theta[2]*graph$res.dist)
Sigma.o <- theta[1]*exp(-theta[2]*graph$res.dist[graph$PtV,graph$PtV]) + theta[3]*diag(rep(1,length(graph$y)))

mu.p <- Sigma[,graph$PtV]%*%solve(Sigma.o,graph$y)

fig <- plot_obs(graph, mu.p[graph$PtV]) + scale_colour_gradientn(colours = heat.colors(10))
print(fig)
