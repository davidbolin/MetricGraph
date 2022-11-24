library(MetricGraph)

line1 <- Line(rbind(c(0,0),c(1,0)))
line2 <- Line(rbind(c(0,0),c(0,1)))
line3 <- Line(rbind(c(0,1),c(-1,1)))
theta <- seq(from=pi,to=3*pi/2,length.out = 20)
line4 <- Line(cbind(sin(theta),1+ cos(theta)))
Lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                              Lines(list(line2),ID="2"),
                              Lines(list(line4),ID="3"),
                              Lines(list(line3),ID="4")))
graph <- metric_graph$new(lines = Lines)
graph$build_mesh(h = 0.01)
graph$plot(line_width = 0.3, mesh = TRUE)

obs.per.edge <- 2
obs.loc <- NULL
for(i in 1:(graph$nE-1)) {
  obs.loc <- rbind(obs.loc,
                   cbind(rep(i,obs.per.edge), runif(obs.per.edge)))
}
obs.loc <- rbind(c(1,0.5),c(2,0.5),c(3,0.5))

x <- y <- NULL
for(i in 1:dim(obs.loc)[1]){
  Line <- graph$Lines[obs.loc[i,1],]
  val_line <- gProject(Line, as(Line, "SpatialPoints"), normalized=TRUE)
  Point <- gInterpolate(Line, obs.loc[i,2], normalized=TRUE)
  x <- c(x, Point@coords[1])
  y <- c(y, Point@coords[2])
}
Vobs <- cbind(x,y)

y <- rnorm(dim(obs.loc)[1])
graph$add_PtE_observations(y,obs.loc)
alpha = 2
#check alpha = 1
if(alpha==1){
  theta = c(1,1,1)
  mu_alpha1 <- spde_posterior_mean(theta, graph, alpha = 1, type = "mesh")
  mu_alpha1_obs <- spde_posterior_mean(theta, graph, alpha = 1, type = "obs")

  p <- graph$plot_function_mesh(mu_alpha1, marker_size = 0)

  p <- p + geom_point(data=data.frame(x=Vobs[,1],y=Vobs[,2],
                                      val = mu_alpha1_obs),
                      mapping= aes(x, y, color = val), size = 2) + coord_fixed()
  p

  graph$plot_function_mesh(mu_alpha1,plotly = TRUE, marker_size = 0)

} else {
  #check alpha = 2
  theta = c(1,1,1)
  graph$buildC(2)
  mu_alpha2 <- spde_posterior_mean(theta, graph, alpha = 2, type = "mesh")
  p <- graph$plot_function_mesh(mu_alpha2, marker_size = 0, plotly = TRUE)

  mu_alpha2_obs <- spde_posterior_mean(theta, graph, alpha = 2, type = "obs")
  p <- graph$plot_function_mesh(mu_alpha2, marker_size = 0)
  p <- p + geom_point(data=data.frame(x=Vobs[,1],y=Vobs[,2],
                                      val = mu_alpha2_obs),
                      mapping= aes(x, y, color = val), size = 2) + coord_fixed()
  print(p)



}

