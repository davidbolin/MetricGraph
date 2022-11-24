library('sf')
library(xtable)
library(scales)
library(gridExtra)

save.fig=FALSE
set.seed(1)

# Load data
Lines <- read_sf('examples/data.pems/lines.shp')
V <- read.csv('examples/data.pems/V.csv',header=T, row.names = NULL)
EtV <- read.csv('examples/data.pems/E.csv',header=T, row.names = NULL)
PtE <- read.csv('examples/data.pems/PtE.csv',header=T, row.names = NULL)
PtE[,1] <- PtE[,1] + 1
Y_loc <- read.csv('examples/data.pems/graph_sensor_locations_bay.csv',header=F)
Y <- read.csv('examples/data.pems/Y.csv',header=T, row.names = NULL)
Y <- as.matrix(Y[,-1])
TE <- table(c(EtV$V1,EtV$V2))

# Create graph object
if(0){
  graph <-  metric_graph$new(lines = as_Spatial(Lines))
  y <- colMeans(Y)
  graph$add_PtE_observations(y,PtE)
} else {
  graph <-  metric_graph$new(V = as.matrix(V),
                             E = as.matrix(EtV[,2:3]),
                             edge_lengths = as.matrix(EtV[,4]))
  graph$add_PtE_observations(colMeans(Y),as.matrix(PtE))
  graph$Lines <- as_Spatial(Lines)
}


p <- graph$plot(marker_size = 0, line_width = 0.1, data = TRUE, data_size = 1)
p <- p + xlab(NULL) + ylab(NULL)
p1 <- graph$plot(marker_size = 0, line_width = 0.1, data = TRUE, data_size = 1)
p1 <- p1 + coord_cartesian(xlim =c(-121.91,-121.88), ylim = c(37.313, 37.328)) +
  theme(legend.position = "none")
p1 <- p1 + xlab(NULL) + ylab(NULL)

p2 <- graph$plot(marker_size = 0, line_width = 0.1, data = TRUE, data_size = 1)
p2 <- p2 + coord_cartesian(xlim =c(-121.94,-121.88), ylim = c(37.35, 37.375)) +
  theme(legend.position = "none")
p2 <- p2 + xlab(NULL) + ylab(NULL)
pp <- grid.arrange(p1, p2, p, ncol=2, layout_matrix = cbind(c(1,2), c(3,3)),
             widths=c(3, 6), heights=c(4, 4),
             bottom = "Longitude", left = "Latitude")

pp
if(save.fig){
  ggsave("data.png",pp, device = "png")
}

#center data to assume zero mean
graph$y <- y - mean(y)

# Fit alpha=1 model
theta.alpha1 <- c(0.1624661, 10.5330059,  0.2525539)
res <- optim(log(theta.alpha1), function(x) -likelihood_graph_spde(exp(x),
                                                                   graph,
                                                                   alpha = 1))
theta.alpha1 <- exp(res$par)
like.alpha1 <- -res$value


graph$build_mesh(h=1e-05)
p1 <- graph$plot_function_mesh(mu_alpha1)
p1 <- p1 + coord_cartesian(xlim =c(-121.91,-121.88), ylim = c(37.313, 37.328))
p1 <- p1 + xlab(NULL) + ylab(NULL)
p1 <- p1 + scale_colour_continuous(type = "viridis", limits = c(47,51.5))
p1
p2 <- graph$plot_function_mesh(mu_alpha1)
p2 <- p2 + coord_cartesian(xlim =c(-121.94,-121.88), ylim = c(37.35, 37.375))
p2 <- p2 + xlab(NULL) + ylab(NULL)
p2 <- p2 + scale_colour_continuous(type = "viridis", limits = c(44,53))
p2

pp <- grid.arrange(p1, p2, ncol=2,
                   bottom = "Longitude", left = "Latitude")

pp
if(save.fig){
  ggsave("data_krig.png",pp, device = "png")
}

#fit alpha = 2 model
graph$buildC(2)
theta.alpha2 <- c(0.3747406, 6.7792918, 0.3729290)
res <- optim(log(theta.alpha2), function(x) -likelihood_graph_spde(exp(x),
                                                                   graph,
                                                                   alpha = 2))
theta.alpha2 <- exp(res$par)
like.alpha2 <- -res$value

mu_alpha2 <- mean(y) + spde_posterior_mean(theta.alpha2, graph, alpha = 2, type = "mesh")


graph$observation_to_vertex()

# Fit isotropic model
graph$compute_resdist()

theta.exp <- c(9.4726902736, 0.0001559032, 0.3745814561)
res.exp <- optim(log(theta.exp), function(x) -likelihood.graph.covariance(exp(x),
                                                                          graph,
                                                                          model = "isoExp"))
theta.exp <- exp(res.exp$par)
like.res <- -res.exp$value



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


###
# print posterior mean
##
theta.alpha2 <- c(0.0001, 1.084454e-04, 1.710348e+01)
X <- spde_posterior_mean(theta.alpha2, graph, alpha = 1, type = "mesh")

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
