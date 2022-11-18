library(GPGraph)
library(sp)

# Example 1: plot of graph, observations, and covariance for alpha = 1
line1 <- Line(rbind(c(0,0),c(1,0)))
line2 <- Line(rbind(c(0,0),c(0,1)))
line3 <- Line(rbind(c(0,1),c(-1,1)))
theta <- seq(from=pi,to=3*pi/2,length.out = 20)
line4 <- Line(cbind(sin(theta),1+ cos(theta)))
Lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                              Lines(list(line2),ID="2"),
                              Lines(list(line3),ID="3"),
                              Lines(list(line4),ID="4")))
graph <- metric_graph$new(Lines = Lines)
graph$plot()
kappa <- 10
sigma <- 2
C <- covariance_alpha1(P = c(1,0.1), kappa = kappa, sigma = sigma, graph = graph, n.p = 50)
gg <- graph$plot_function(C, flat = FALSE)
gg <- graph$plot_function(C)

C <- covariance_alpha2(P = c(1,0.1), kappa = kappa, sigma = sigma, graph = graph, n.p = 50)
gg <- graph$plot_function(C, flat = FALSE)

y <- c(1,2,3)
PtE <- matrix(c(1, 2, 3, 0.5, 0.5, 0.7),3,2)
graph$add_observations2(y, PtE)
p <- graph$plot(data=TRUE)

# Example 2: test with different graph
line1 <- Line(rbind(c(0,0),c(10,0)))
line2 <- Line(rbind(c(0,0),c(0,1)))
line3 <- Line(rbind(c(0,1),c(-1,1)))
Lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                              Lines(list(line2),ID="2"),
                              Lines(list(line3),ID="3")))
graph <- metric_graph$new(Lines = Lines)
p <- graph$plot()

y <- c(1,2,3)
PtE <- matrix(c(1, 2, 3, 0.5, 0.5, 0.7),3,2)
graph$add_observations2(y, PtE)
p <- graph$plot(data=TRUE)


# Example 3: test with alpha = 2
line1 <- Line(rbind(c(0,0),c(1,0)))
line2 <- Line(rbind(c(0,0),c(0,1)))
line3 <- Line(rbind(c(0,1),c(-1,1)))
theta <- seq(from=pi,to=3*pi/2,length.out = 20)
line4 <- Line(cbind(sin(theta),1+ cos(theta)))
Lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                              Lines(list(line2),ID="2"),
                              Lines(list(line3),ID="3"),
                              Lines(list(line4),ID="4")))
graph <- metric_graph$new(Lines = Lines)
graph$plot()
kappa <- 10
sigma <- 2
C <- covariance_alpha2(P = c(1,0.1), kappa = kappa, sigma = sigma, graph = graph, n.p = 50)
gg <- graph$plot_function(C, flat = FALSE)

library(Matrix)
EtL <-Matrix::sparseMatrix(i    = 1:dim(graph$E)[1],
                           j    = c(1:length(graph$Lines)),
                           x    = rep(1,dim(graph$E)[1]),
                           dims = c(dim(graph$E)[1], length(graph$Lines)) )
#fix plot function

#' @param data.to.plot
plot_curve = function(data.to.plot,
                      Line_edge,
                      flat = TRUE,
                      normalized = TRUE,
                      color = 'rgb(0,0,200)',
                      p = NULL, ...){

  data.to.plot.order <- data.to.plot[order(data.to.plot[, 1]), ,drop=F]
  p2 <- rgeos::gInterpolate(Line_edge, data.to.plot.order[, 1,drop=F],
                            normalized = normalized)
  coords <-p2@coords
  data <- data.frame(x = coords[,1], y = coords[,2],
                     z = data.to.plot.order[,2])
  if(flat){
    if(is.null(p)){
      p <- ggplot2::ggplot(data = data,
                           ggplot2::aes(x = x, y = y,
                                        colour = z)) +
        ggplot2::geom_path(...)
    } else {
      p <- p + ggplot2::geom_path(data = data,...)
    }
  } else {
    if(is.null(p)){
      p <- plot_ly(data=data, x = ~y, y = ~x, z = ~z)
      p <- p %>% add_trace(mode="lines", type="scatter3d",
                           line = list(color = color, ...),
                           showlegend = FALSE)
    } else {
      p <- p %>% add_trace(data=data, x = ~y, y = ~x, z = ~z,
                           mode = "lines", type = "scatter3d",
                           line = list(color = color, ...),
                           showlegend = FALSE)
    }
  }
  return(p)
}
