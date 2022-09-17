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
  Sigma <- theta[3]^2*exp(-theta[2]*graph$res.dist[graph$PtV,graph$PtV]) +
    theta[1]^2*diag(rep(1,length(graph$y)))
  R <- chol(Sigma)
  return(-sum(log(diag(R))) - 0.5*t(graph$y)%*%solve(Sigma,graph$y))
}
theta <- c(1,1,1)
res <- optim(log(theta), function(x) -like.exp(exp(x),graph))
theta.res <- exp(res$par)

# Fit alpha=1 model
theta.exp <- c(9.4726902736, 0.0001559032, 0.3745814561)
res <- optim(log(theta.exp), function(x) -likelihood.exp.graph.v2(exp(x),graph) )
theta.exp <- exp(res$par)

# Fit graph Laplace model
graph$compute_laplacian()
theta.graph <- c(9.478923, 7.154005, 1.160011e-05)
res <- optim(log(theta.graph), function(x) -likelihood.graph_laplacian(exp(x),graph) )
theta.graph <- exp(res$par)

#cross validation
K <- 5 #number of groups in cross validation
n.y <- length(graph$y)
n.g <- floor(n.y/K)
ind <- rep(1:K,n.g+1)[1:n.y]
ind <- ind[sample(1:n.y,n.y)]

y_L_iso <- posterior.crossvalidation(theta.res, graph, model = "isoExp",ind = ind)
Y_L_exp <- posterior.crossvalidation(theta.exp,graph, model = "alpha1",ind = ind)
Y_L_graph <- posterior.crossvalidation(theta.graph,graph, model ="GL",ind = ind)
result <- data.frame(MSE = c(var(graph$y-y_L_iso),
                            var(Y_L_exp-graph$y),
                            var(Y_L_graph-graph$y)),
                    row.names = c("isoExp","alpha1","GL"))
print(result)


fig <- plot_obs(graph,Y_L_graph, y_loc = Y_loc) + scale_colour_gradientn(colours = heat.colors(10))
print(fig)
