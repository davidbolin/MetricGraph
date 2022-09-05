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

TE <- table(c(EtV$V1,EtV$V2))
#length(names(TE)[TE1=2]) 320 namn

graph <-  graph.obj$new()
graph$Lines <- as_Spatial(Lines)
graph$V <- as.matrix(V)
graph$El <- as.matrix(EtV[,4])
graph$EtV <- as.matrix(EtV[,1:3])
graph$PtE <- as.matrix(PtE)
graph$y <- as.matrix(Y[,-1])[1,]
graph$y <- graph$y - mean(graph$y ) #temporary

theta <-  c(7.043030174, 0.000139784, 0.537848369)
#plot covariance for the parameters
lik <- likelihood.exp.graph(theta,graph)
#res <- optim(log(theta), function(x) -likelihood.exp.graph(exp(x),graph) )
#res <- optim(log(theta[c(1,3)]), function(x) -likelihood.exp.graph(c(exp(x[1]),theta[2],exp(x[2])),graph) )
#theta <- c(exp(res$par[1]),theta[2],exp(res$par[2]))

plot_r_1(theta)

#plot points

fig <- plot_obs(graph) +  scale_colour_gradientn(colours = heat.colors(10))
print(fig)


gg <- plot_posterior_mean(graph, poster.mean.exp, sample.line.expontial, theta, size=0.3)
gg <- gg +  scale_colour_gradientn(colours = heat.colors(10))
print(gg)
#leavoe one line out cross val
Y_ <- poster.mean.obs.exp(theta,graph)

fig <- plot_obs(graph, Y_) + scale_colour_gradientn(colours = heat.colors(10))
print(fig)


Y_leave <- poster.mean.obs.exp(theta,graph,leave.edge.out =T)

fig <- plot_obs(graph, Y_leave) + scale_colour_gradientn(colours = heat.colors(10))
print(fig)
