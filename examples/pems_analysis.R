rm(list=ls())
library(MetricGraph)

library('sf')
library(xtable)
library(scales)
library(gridExtra)
library(ggmap)


show.plot=FALSE
set.seed(1)

# Load data
Lines <- read_sf('examples/data.pems/lines.shp')
EtV <- read.csv('examples/data.pems/E.csv',header=T, row.names = NULL)
PtE <- read.csv('examples/data.pems/PtE.csv',header=T, row.names = NULL)
PtE[,1] <- PtE[,1] + 1
Y <- read.csv('examples/data.pems/Y.csv',header=T, row.names = NULL)
Y <- colMeans(as.matrix(Y[,-1]))

lines <- as_Spatial(Lines)
graph <-  metric_graph$new(lines = as_Spatial(Lines), longlat = TRUE)

#convert PtE to relative distances and add observations
edge_length_m <- EtV[,4]
PtE[,2] = PtE[,2]/edge_length_m[PtE[,1]]

data.PtE <- data.frame(y = Y - mean(Y),
                      edge_number = PtE[,1],
                      distance_on_edge = PtE[,2])
graph$add_observations(data = data.PtE,
                            normalized = TRUE)


#fit alpha = 2 model
graph$buildC(2)
theta.alpha2 <- graph_starting_values(graph, model = "alpha2",data_name="y")
loglik_alpha2 <- likelihood_graph_spde(graph, alpha = 2,data_name="y")
res <- optim(log(theta.alpha2), loglik_alpha2)
theta.alpha2 <- exp(res$par)
like.alpha2<- -res$value

graph$buildC(2,TRUE)
loglik_alpha2_BC <- likelihood_graph_spde(graph, alpha = 2,data_name="y", BC=0)
res <- optim(log(theta.alpha2), loglik_alpha2_BC)
theta.alpha2_BC <- exp(res$par)
like.alpha2_BC<- -res$value


# Fit alpha=1 model
theta.alpha1 <- graph_starting_values(graph, model = "alpha1",data_name="y")
loglik_alpha1 <- likelihood_graph_spde(graph, alpha = 1, data_name="y")
res <- optim(log(theta.alpha1), loglik_alpha1)
theta.alpha1 <- exp(res$par)
like.alpha1<- -res$value

loglik_alpha1_BC <- likelihood_graph_spde(graph, alpha = 1, data_name="y", BC=0)
res <- optim(log(theta.alpha1), loglik_alpha1_BC)
theta.alpha1_BC <- exp(res$par)
like.alpha1_BC<- -res$value


# Fit isotropic model
graph$compute_resdist()
theta.exp <- graph_starting_values(graph, model = "isoExp",data_name="y" )

loglik_isoExp <- likelihood_graph_covariance(graph, model = "isoCov",
                                             cov_function = exp_covariance)
res.exp <- optim(log(theta.exp), loglik_isoExp)
theta.exp <- exp(res.exp$par)
like.exp <- -res.exp$value

# Fit graph Laplace model
graph$observation_to_vertex()
graph$compute_laplacian()
theta.GL1 <- graph_starting_values(graph, model = "GL1")
loglik_GL1 <- likelihood_graph_laplacian(graph, alpha = 1)
res <- optim(log(theta.GL1), loglik_GL1)
theta.GL1 <- exp(res$par)
like.GL1 <- -res$value

# Fit graph Laplace model 2
theta.GL2 <- graph_starting_values(graph, model = "GL2")
loglik_GL2 <- likelihood_graph_laplacian(graph, alpha = 2)
res <- optim(log(theta.GL2), loglik_GL2)
theta.GL2 <- exp(res$par)
like.GL2 <- -res$value

#cross validation
K <- 5 #number of groups in cross validation
n.y <- length(graph$data[['y']])
n.g <- floor(n.y/K)
ind <- rep(1:K,n.g+1)[1:n.y]
ind <- ind[sample(1:n.y,n.y)]

cv.exp <- posterior_crossvalidation(theta.exp, graph,
                                    model = "isoExp", ind = ind)
cv.alpha1 <- posterior_crossvalidation(theta.alpha1, graph,
                                       model = "alpha1", ind = ind)
cv.GL1 <- posterior_crossvalidation(theta.GL1, graph, model ="GL1", ind = ind)
cv.GL2 <- posterior_crossvalidation(theta.GL2, graph, model ="GL2", ind = ind)
cv.alpha2 <- posterior_crossvalidation_covariance(theta.alpha2, graph,
                                                  model = "alpha2", ind = ind)
result <- data.frame(RMSE  = sqrt(c(cv.exp$rmse, cv.alpha1$rmse, cv.GL1$rmse,
                                    cv.GL2$rmse, cv.alpha2$rmse)),
                     MAE   = c(cv.exp$mae, cv.alpha1$mae, cv.GL1$mae,
                               cv.GL2$mae, cv.alpha2$mae),
                     LS    = c(cv.exp$logscore, cv.alpha1$logscore,
                                cv.GL1$logscore, cv.GL2$logscore,
                                cv.alpha2$logscore),
                     CRPS  = c(cv.exp$crps, cv.alpha1$crps, cv.GL1$crps,
                                cv.GL2$crps, cv.alpha2$crps),
                     SCRPS = c(cv.exp$scrps, cv.alpha1$scrps, cv.GL1$scrps,
                                cv.GL2$scrps, cv.alpha2$scrps),
                     nlike = -c(like.exp, like.alpha1, like.GL1, like.GL2,
                                like.alpha2),
                     row.names = c("isoExp","alpha1","GL1", "GL2","alpha2"))
print(result, digits = 3)

#print(xtable(result))

if(show.plot) {
  map2=get_map(c(left = -121.91, bottom = 37.3,
                 right = -121.85, top = 37.35),
               maptype = "toner-lite")
  map3=get_map(c(left = -121.96, bottom = 37.34,
                 right = -121.85, top = 37.38),
               maptype = "toner-lite")
  map1 <- get_map(c(left = -122.2, bottom = 37.05,
                    right = -121.7, top = 37.6),
                  maptype = "toner-lite")

  p <- ggmap(map1)
  p <- p + xlab(NULL) + ylab(NULL)
  p <- graph$plot(vertex_size = 0,
                  line_width = 0.6, data = TRUE, data_size = 2,
                  p = p)
  p <- p + coord_cartesian(xlim =c(-122.1,-121.8), ylim = c(37.2, 37.45))
  p1 <- ggmap(map2)
  p1 <- p1 + xlab(NULL) + ylab(NULL)
  p1 <- graph$plot(vertex_size = 0, line_width = 0.4, data = TRUE, data_size = 2,
                   p = p1)
  p1 <- p1 + coord_cartesian(xlim =c(-121.905,-121.88), ylim = c(37.313, 37.328)) +
    theme(legend.position = "none")
  p2 <- ggmap(map3)
  p2 <- p2 + xlab(NULL) + ylab(NULL)
  p2 <- graph$plot(vertex_size = 0, line_width = 0.4, data = TRUE,
                   data_size = 2, p = p2)
  p2 <- p2 + coord_cartesian(xlim =c(-121.93,-121.88), ylim = c(37.35, 37.375)) +
    theme(legend.position = "none")
  pp <- grid.arrange(p1, p2, p, ncol=2, layout_matrix = cbind(c(1,2), c(3,3)),
                     widths=c(3, 6), heights=c(4, 4),
                     bottom = "Longitude", left = "Latitude")

  pp
  if(save.fig){
    ggsave("data.png",pp, device = "png")
  }

  graph$y <- Y

  p1 <- ggmap(map2)
  p1 <- p1 + xlab(NULL) + ylab(NULL)
  p1 <- graph$plot_function(mu_alpha2, p = p1, line_width = 0.75)
  p1 <- graph$plot(data = TRUE, p = p1, vertex_size = 0, line_width = 0,
                   data_size = 2)
  p1 <- p1 + coord_cartesian(xlim =c(-121.905,-121.88), ylim = c(37.316, 37.328))
  p1 <- p1 + xlab(NULL) + ylab(NULL)
  p1 <- p1 + scale_colour_continuous(type = "viridis") +
    theme(legend.position = "none")

  p2 <- ggmap(map3)
  p2 <- p2 + xlab(NULL) + ylab(NULL)
  p2 <- graph$plot_function_mesh(mu_alpha2, p = p2, line_width = 0.5)
  p2 <- graph$plot(data = TRUE, p = p2, vertex_size = 0, line_width = 0,
                   data_size = 2)
  p2 <- p2 + coord_cartesian(xlim =c(-121.94,-121.88), ylim = c(37.35, 37.375))
  p2 <- p2 + xlab(NULL) + ylab(NULL)
  p2 <- p2 + scale_colour_continuous(type = "viridis")#, limits = c(44,53))
  p2

  pp <- grid.arrange(p1, p2, ncol=2,
                     bottom = "Longitude", left = "Latitude",
                     widths=c(4, 5))

  pp
  if(save.fig){
    ggsave("data_krig.png",pp, device = "png")
  }

  graph$build_mesh(h=1e-04)
  mu_alpha2 <- mean(Y) + spde_posterior_mean(theta.alpha2, graph, alpha = 2,
                                             type = "mesh")

}
