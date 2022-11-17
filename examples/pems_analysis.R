library('sf')
library(xtable)
library(scales)

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
graph <-  metric_graph$new(Lines = as_Spatial(Lines))

plot_2d = function(graph, show = TRUE,
                line_width = 1,
                marker_size = 1,
                vertex_color = 'rgb(0,0,0)',
                edge_color = 'rgb(0,0,0)',
                data = FALSE,
                data_size = 1,
                mesh = FALSE,
                fix_layout = TRUE,
                zoom = 1,
                ...){
  if(is.null(graph$Lines)){
    data.plot <- data.frame(x = c(graph$V[E[,1],1],graph$V[E[,2],1]),
                            y = c(graph$V[E[,1],2],graph$V[E[,2],2]),
                            z = rep(0,2 * graph$nE),
                            i = c(1:graph$nE, 1:graph$nE))
  } else {
    x <- y <- ei <- NULL
    for(i in 1:graph$nE) {
      xi <- graph$Lines@lines[[i]]@Lines[[1]]@coords[,1]
      yi <- graph$Lines@lines[[i]]@Lines[[1]]@coords[,2]
      ii <- rep(i,length(xi))
      x <- c(x,xi)
      y <- c(y,yi)
      ei <- c(ei,ii)
    }
    data.plot <- data.frame(x = x, y = y,
                            z = rep(0,length(x)), i = ei)
  }
  p <- plot_ly(data=data.plot, x = ~y, y=~x)
  p <- p %>% add_trace(data=data.plot, x = ~y, y=~x,
                       mode="lines",type="scatter",
                       line = list(width = line_width,
                                   color = edge_color, ...),
                       split=~i, showlegend=FALSE)

  if(marker_size > 0) {
    data.plot2 <- data.frame(x=graph$V[,1],y=graph$V[,2],z=rep(0,graph$nV))
    p <- p %>% add_trace(data=data.plot2, x = ~y,y = ~x, z = ~z,
                         type="scatter", mode = "markers",
                         marker = list(size = marker_size,
                                       color = vertex_color, ...))
  }

  if(data){
    x <- y <- NULL
    for(i in 1:length(graph$y)){
      Line <- graph$Lines[PtE[i,1],]
      val_line <- gProject(Line, as(Line, "SpatialPoints"), normalized=TRUE)
      Point <- gInterpolate(Line, PtE[i,2], normalized=TRUE)
      x <- c(x, Point@coords[1])
      y <- c(y, Point@coords[2])
    }
    data.plot <- data.frame(x = x, y = y,
                            z = rep(0,length(x)),
                            val = graph$y)
    p <- p %>% add_trace(data=data.plot, x = ~y, y = ~x, z = ~z,
                         type="scatter3d", mode = "markers",
                         marker = list(size = marker_size,
                                       color = ~val,
                                       colorbar=list(title='', len = 0.5),
                                       colorscale='Viridis'),
                         showlegend=FALSE)
  }
  if (mesh) {
    data.plot <- data.frame(x = graph$mesh$V[,1],
                            y = graph$mesh$V[,2],
                            z = rep(0,dim(graph$mesh$V)[1]))
    p <- p %>% add_trace(data=data.plot, x = ~y, y = ~x, z = ~z,
                         type="scatter3d", mode = "markers",
                         marker = list(size = marker_size/2,
                                       color = 'rgb(100,100,100)'),
                         showlegend=FALSE)
  }
  xr <- (diff(range(graph$V[,1])) + diff(range(graph$V[,2])))/zoom
  if(fix_layout){
    ax <- list(title = '',
               zeroline = FALSE,
               showgrid = FALSE,
               showticklabels=FALSE)
    p <- p %>% layout(title = '',
                      scene = list(xaxis = ax, yaxis = ax, zaxis = ax,
                                   camera = list(eye = list(x = 0, y = 0, z = -xr),
                                                 up = list(x=1,y=0,z=0)),
                                   aspectmode='data'))

  }

  if(show){
    print(p)
  }
  return(p)
}
#graph$plot(marker_size = 0)

#graph$V <- as.matrix(V)
#graph$edge_lengths <- as.matrix(EtV[,4])
#graph$E <- as.matrix(EtV[,1:3])
#graph$PtE <- as.matrix(PtE)
y <- colMeans(Y)
y <- y - mean(y) #temporary
#graph$nV <- dim(graph$V)[1]
#graph$nE <- length(graph$edge_lengths)
graph$add_observations2(y,PtE)


p <- graph$plot(data = TRUE, marker_size = 0, data_size = 5, zoom = 0.6)
p %>% config(toImageButtonOptions = list(width=300,height=200))
print(p)

line1 <- Line(rbind(c(0,0),c(1,0)))
line2 <- Line(rbind(c(0,0),c(0,1)))
line3 <- Line(rbind(c(0,1),c(-1,1)))
theta <- seq(from=pi,to=3*pi/2,length.out = 20)
line4 <- Line(cbind(sin(theta),1+ cos(theta)))
Lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                              Lines(list(line2),ID="2"),
                              Lines(list(line4),ID="3"),
                              Lines(list(line3),ID="4")))
graph <- metric_graph$new(Lines = Lines)
p <- graph$plot()

library(plotly)

us_cities = read.csv("https://raw.githubusercontent.com/plotly/datasets/master/us-cities-top-1k.csv")
dat <- data.frame(lat = graph$V[,1], lon = graph$V[,2])
fig <- dat
fig <- fig %>%
  plot_ly(
    lat = ~lat,
    lon = ~lon,
    marker = list(color = "fuchsia"),
    type = 'scattermapbox',
    hovertext = dat[,1])
fig <- fig %>%
  layout(
    mapbox = list(
      style = 'open-street-map',
      zoom =5.5,
      center = list(lon = -121.9376, lat = 37.3479)))

fig


graph$buildC(2)
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

cv.res <- posterior.crossvalidation.covariance(theta.exp, graph, model = "isoExp",ind = ind)
cv.alpha1 <- posterior.crossvalidation.covariance(theta.alpha1,graph, model = "alpha1",ind = ind)
cv.gl <- posterior.crossvalidation.covariance(theta.graph,graph, model ="GL",ind = ind)
cv.gl2 <- posterior.crossvalidation.covariance(theta.graph2,graph, model ="GL2",ind = ind)
cv.alpha2 <- posterior.crossvalidation.covariance(theta.alpha2,graph, model = "alpha2",ind = ind)
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
