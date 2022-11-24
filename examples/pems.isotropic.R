library('sf')
library(GPGraph)
library(xtable)
library(scales)
save.fig=F
set.seed(1)
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
graph$y <- colMeans(Y)
graph$y <- graph$y - mean(graph$y ) #temporary
graph$nV <- dim(graph$V)[1]
graph$nE <- length(graph$El)
fig <- plot_obs(graph,graph$y, y_loc = Y_loc) + scale_colour_gradientn(colours = heat.colors(10))
print(fig)



graph$buildA(2, F)
theta.alpha2 <- c(1.710348e+01, 1.084454e-04, 0.0001)
res <- optim(log(theta.alpha2), function(x) -likelihood.matern2.graph(exp(x),graph) )
theta.alpha2 <- exp(res$par)
like.alpha2 <- -res$value

graph$observation_to_vertex()

# Fit isotropic model
graph$compute_resdist()

theta.exp <- c(9.4726902736, 0.0001559032, 0.3745814561)
res.exp <- optim(log(theta.exp), function(x) -likelihood.graph.covariance(exp(x),graph, model = "isoExp"))
theta.exp <- exp(res.exp$par)
like.res <- -res.exp$value

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
#theta.alpha2 <- c(1.048223e+01, 4.787329e-04, 4.046737e-04)
#res <- optim(log(theta.alpha2), function(x) -likelihood.graph.covariance(exp(x),graph,model = "alpha2") )
#theta.alpha2 <- exp(res$par)
#like.alpha2 <- -res$value

# Fit graph Laplace model 2
graph$compute_laplacian()
theta.graph2 <- c(9.478923, 7.154005, 1.160011e-05)
res <- optim(log(theta.graph2), function(x) -likelihood.graph.covariance(exp(x),graph, model = "GL2") )
theta.graph2 <- exp(res$par)
like.graph2 <- -res$value

#cross validation
K <- 5 #number of groups in cross validation
n.y <- length(graph$y)
n.g <- floor(n.y/K)
ind <- rep(1:K,n.g+1)[1:n.y]
ind <- ind[sample(1:n.y,n.y)]

cv.res <- posterior_crossvalidation_covariance(theta.exp, graph, model = "isoExp",ind = ind)
cv.alpha1 <- posterior_crossvalidation_covariance(theta.alpha1,graph, model = "alpha1",ind = ind)
cv.gl <- posterior_crossvalidation_covariance(theta.graph,graph, model ="GL",ind = ind)
cv.gl2 <- posterior_crossvalidation_covariance(theta.graph2,graph, model ="GL2",ind = ind)
cv.alpha2 <- posterior_crossvalidation_covariance(theta.alpha2,graph, model = "alpha2",ind = ind)
result <- data.frame(RMSE  = sqrt(c(cv.res$rmse, cv.alpha1$rmse, cv.gl$rmse,cv.gl2$rmse, cv.alpha2$rmse)),
                     MAE   = c(cv.res$mae, cv.alpha1$mae, cv.gl$mae, cv.gl2$mae, cv.alpha2$mae),
                     LS    = -c(cv.res$logscore, cv.alpha1$logscore, cv.gl$logscore, cv.gl2$logscore, cv.alpha2$logscore),
                     CRPS  = -c(cv.res$crps, cv.alpha1$crps, cv.gl$crps, cv.gl2$crps, cv.alpha2$crps),
                     SCRPS = -c(cv.res$scrps, cv.alpha1$scrps, cv.gl$scrps, cv.gl2$scrps, cv.alpha2$scrps),
                     nlike = -c(like.res, like.alpha1, like.graph, like.graph2, like.alpha2),
                    row.names = c("isoExp","alpha1","GL", "GL2","alpha2"))
print(xtable(result))

##
# print observation
##
fig <- plot_obs(graph,colMeans(Y), y_loc = Y_loc, size_path=0.1, size_obs=3) + scale_colour_viridis_c(breaks =c(30,40,60))+
       xlab("long") + ylab("lat") + theme_classic()+
  labs(colour="mph")+theme(axis.text=element_text(size=20),
                           axis.title=element_text(size=20,face="bold"),
                           legend.text = element_text(size=20),
                           legend.title = element_text(size=20,face="bold"))
print(fig)
if(save.fig)
  ggsave("observation.pdf",fig)

###
# print posterior mean
##
X <- graph_posterior_mean_matern2(graph,  theta.alpha2)
X[,2] <- X[,2] + mean(colMeans(Y))
gg <- plot_curve(X[X[,3]==1,1:2]  , graph$Lines[1,], normalized=T,size=0.8)
for(i in 1:dim(graph$EtV)[1]){
  gg <- plot_curve(X[X[,3]==i,1:2] , graph$Lines[i,], normalized=T, gg = gg,size=0.8)
}

gg <- gg +  scale_colour_viridis_c(breaks =c(30,40,60))+
  xlab("long") + ylab("lat") + theme_classic()+
  labs(colour="mph")+theme(axis.text=element_text(size=20),
                           axis.title=element_text(size=20,face="bold"),
                           legend.text = element_text(size=20),
                           legend.title = element_text(size=20,face="bold"))
print(gg)
if(save.fig)
  ggsave("interpolation.pdf",gg)
gg_zoom <- gg +  xlim(-121.9, -121.885) + ylim(37.32,37.3265)

fig_zoom <- fig + xlim(-121.9, -121.885) + ylim(37.32,37.327)
X.sample <- graph_posterior_mean_matern2(graph,  theta.alpha2,T)
