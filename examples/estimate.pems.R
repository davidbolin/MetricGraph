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
graph$buildA(2, F)
theta <- c(9.4726902736, 0.0001559032, 0.3745814561)
#plot covariance for the parameters
lik <- likelihood.exp.graph(theta,graph)
res <- optim(log(theta), function(x) -likelihood.exp.graph(exp(x),graph) )
res <- optim(log(theta), function(x) -likelihood.exp.graph.stupid(exp(x),graph) )
#res <- optim(log(theta[c(1,3)]), function(x) -likelihood.exp.graph(c(exp(x[1]),theta[2],exp(x[2])),graph) )
#theta <- c(exp(res$par[1]),theta[2],exp(res$par[2]))

theta2 <- c(6.531191638, 0.006921311, 0.020689930)
res <- optim(log(theta2), function(x) -likelihood.matern2.graph( exp(x),graph) )
lik2 <- likelihood.matern2.graph(theta2, graph)
#theta2 <- exp(res$par)

#plot_r_1(theta)

#plot points

fig <- plot_obs(graph) +  scale_colour_gradientn(colours = heat.colors(10))
print(fig)
gg <- plot_posterior_mean(graph, posterior.mean.matern2, sample.line.matern2, theta2, byVertex=F, size=0.3)
gg <- plot_posterior_mean(graph, posterior.mean.exp, sample.line.expontial, theta, size=0.3)
gg <- gg +  scale_colour_gradientn(colours = heat.colors(10))
print(gg)
#leavoe one line out cross val
Y_2 <- posterior.mean.stupid(theta,graph)[graph$PtV]
Y_1 <- posterior.mean.obs.exp(theta,graph)
Y_2 <- posterior.mean.obs.matern2(theta2,graph)
fig <- plot_obs(graph, Y_2-graph$y ) + scale_colour_gradientn(colours = heat.colors(10))
print(fig)


Y_1_leave <- posterior.mean.obs.exp(theta,graph,leave.edge.out =T)
Y_1_leave2 <- posterior.leave.stupid(theta,graph)
Y_2_leave <- posterior.mean.obs.matern2(theta2,graph,leave.edge.out =T)

fig <- plot_obs(graph, Y_1-graph$y) + scale_colour_gradientn(colours = heat.colors(10))
print(fig)
