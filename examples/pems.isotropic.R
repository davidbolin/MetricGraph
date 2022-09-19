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

theta <- c(1,1,1)
res <- optim(log(theta), function(x) -likelihood.graph.covariance(exp(x),graph, model = "isoExp"))
theta.res <- exp(res$par)
like.res <- -res$value

# Fit alpha=1 model
theta.alpha1 <- c(9.4726902736, 0.0001559032, 0.3745814561)
res <- optim(log(theta.alpha1), function(x) -likelihood.graph.covariance(exp(x),graph, model = "alpha1") )
theta.alpha1 <- exp(res$par)
like.alpha1 <- -res$value

# Fit graph Laplace model
graph$compute_laplacian()
theta.graph <- c(9.478923, 7.154005, 1.160011e-05)
res <- optim(log(theta.graph), function(x) -likelihood.graph.covariance(exp(x),graph, model = "GL") )
theta.graph <- exp(res$par)
like.graph <- -res$value

# Fit alpha=2 model
graph$buildA(2, F)
theta.alpha2 <- c(1.710348e+01, 1.084454e-04, 0.0001)
res <- optim(log(theta.alpha2), function(x) -likelihood.graph.covariance(exp(x),graph,model = "alpha2") )
theta.alpha2 <- exp(res$par)
like.alpha2 <- -res$value

#cross validation
K <- 5 #number of groups in cross validation
n.y <- length(graph$y)
n.g <- floor(n.y/K)
ind <- rep(1:K,n.g+1)[1:n.y]
ind <- ind[sample(1:n.y,n.y)]

cv.res <- posterior.crossvalidation.covariance(theta.res, graph, model = "isoExp",ind = ind)
cv.alpha1 <- posterior.crossvalidation.covariance(theta.exp,graph, model = "alpha1",ind = ind)
cv.gl <- posterior.crossvalidation.covariance(theta.graph,graph, model ="GL",ind = ind)
cv.alpha2 <- posterior.crossvalidation.covariance(theta.alpha2,graph, model = "alpha2",ind = ind)
result <- data.frame(RMSE = c(cv.res$rmse, cv.alpha1$rmse, cv.gl$rmse, cv.alpha2$rmse),
                     MAE = c(cv.res$mae, cv.alpha1$mae, cv.gl$mae, cv.alpha2$mae),
                     LS = -c(cv.res$logscore, cv.alpha1$logscore, cv.gl$logscore, cv.alpha2$logscore),
                     CRPS = -c(cv.res$crps, cv.alpha1$crps, cv.gl$crps, cv.alpha2$crps),
                     SCRPS = -c(cv.res$scrps, cv.alpha1$scrps, cv.gl$scrps, cv.alpha2$scrps),
                     nlike = -c(like.res, like.alpha1, like.graph, like.alpha2),
                    row.names = c("isoExp","alpha1","GL", "alpha2"))
print(result)


fig <- plot_obs(graph,Y_L_graph, y_loc = Y_loc) + scale_colour_gradientn(colours = heat.colors(10))
print(fig)
