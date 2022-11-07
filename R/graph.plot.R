#' plot a metric graph
#'
#' @param graph metric_graph object
#' @param show show the plot?
#' @return a plotly plot object
#' @export
metric_graph.plot <- function(graph, show = TRUE, ...){

  data <- data.frame(x = c(graph$V[E[,1],1],graph$V[E[,2],1]),
                     y = c(graph$V[E[,1],2],graph$V[E[,2],2]),
                     z = rep(0,2 * graph$nE),
                     i = c(1:graph$nE, 1:graph$nE))
  p <- plot_ly(data=data, x = ~y, y=~x,z=~z)
  p <- p %>% add_trace(data=data, x = ~y, y=~x, z=~z, mode="lines",type="scatter3d",
                       line = list(width = 1,color='rgb(0,0,0)'),
                       split=~i, showlegend=FALSE)

  data2 <- data.frame(x=graph$V[,1],y=graph$V[,2],z=rep(0,graph$nV))
  p <- p %>% add_trace(data=data2,x=~y,y=~x,z=~z,type="scatter3d", mode = "markers",
                       marker = list(width = 10,color='rgb(0,0,0)'))
  if(show){
    print(p)
  }
  return(p)
}


#' plot a continuous curve
#'
#' @param data.to.plot (m x 2) [,1] position on curve (in length)
#'                             [,2] value
#' @param Line         (SpatialLinesDataFrame,SpatialLines)
#' @param normalized  (bool) is the position in relative length
#' @param gg           (ggplot pbj) add to gg obj if not null
#' @export
#'
plot_curve <- function(data.to.plot, Line_edge, normalized=F ,gg = NULL,...){

  data.to.plot.order <- data.to.plot[order(data.to.plot[,1]),]
  p2 <- rgeos::gInterpolate(Line_edge, data.to.plot.order[,1], normalized = normalized)
  coords <-p2@coords
  if(is.null(gg)){
    gg <- ggplot2::ggplot(data.frame(x = coords[,1],
                                     y = coords[,2],
                                     val = data.to.plot.order[,2]),
                          ggplot2::aes(x=x, y=y, colour=val ) ) +
      ggplot2::geom_path(...)
  }else{
    gg <- gg + ggplot2::geom_path(data = data.frame(x = coords[,1],
                                                    y = coords[,2],
                                                    val = data.to.plot.order[,2]),...)
  }
  return(gg)
}


#' plot a straight curve curve
#'
#' @param data.to.plot (m x 2) [,1] position on curve (in length)
#'                             [,2] value
#' @param V (2 x 2) the position of the two vertices in 2D
#' @param p plot object, add to p if not null
#' @export
#'
plot_straight_curve <- function(data.to.plot, V, flat = FALSE, p = NULL,...){

  data.to.plot.order <- data.to.plot[order(data.to.plot[,1]),]
  l <- sqrt(sum( (V[1,] - V[2,])^2) )
  alpha <- data.to.plot.order[,1]/l
  coords <- cbind((1-alpha) * V[1,1] + alpha * V[2,1],
                  (1-alpha) * V[1,2] + alpha * V[2,2])
  data <- data.frame(x = coords[,1], y = coords[,2],
                     z = data.to.plot.order[,2])
  if(flat){
    if(is.null(p)){
      p <- ggplot2::ggplot(data = data,
                            ggplot2::aes(x=x, y=y, colour=z)) +
        ggplot2::geom_path(...)
    }else{
      p <- p + ggplot2::geom_path(data = data, ...)
    }
  } else {
    if(is.null(p)){
      p <- plot_ly(data=data, x = ~y, y=~x,z=~z)
      p <- p %>% add_trace(mode="lines", type="scatter3d",
                           line = list(...), showlegend = FALSE)
    } else {
      p <- p %>% add_trace(data=data,x=~y,y=~x,z=~z,mode="lines",type="scatter3d",
                           line = list(...), showlegend = FALSE)
    }
  }

  return(p)
}

#' plot generic data object X on the graph
#'
#' @param X (m x 3) [, 1] edge number
#'                  [, 2] position on curve (in length)
#'                  [, 3] value
#' @param graph metric_graph object
#' @param flat plot in 2D or 3D?
#' @export
plot_X_to_graph<- function(X, graph, flat = TRUE,
                           graph_color = 'rgb(0,0,0)',
                           graph_width = 1,
                           marker_size = 10,
                           ...){
  if(flat == FALSE){
    p <- graph$plot(color = graph_color, line_width = graph_width,
                    marker_size = marker_size)
  } else {
    p <- NULL
  }
  for(i in 1:graph$nE){
    if(is.null(graph$Lines) == TRUE)
      p <- plot_straight_curve(X[X[,1]==i,2:3],
                                graph$V[graph$E[i,],],
                                flat = flat,
                                p = p, ...)
  }
  return(p)
}

#' (depricated)
#' @param graph          - graph.obj
#' @param posterior.mean - (function) takes theta, and graph computes posterior mean of the verteces
#' @param sample.line    - (function) takes theta, and graph computes posterior mean
#' @param theta          - params
#' @param byVertex       - if byVertex V.i is by vertex ortherwise by edge
#' @export
plot_posterior_mean <- function(graph, posterior.mean, sample.line, theta, byVertex=T, X.intercept=0, ...){

  V.post.mean <- posterior.mean(theta, graph)
  for(i in 1:dim(graph$EtV)[1]){

    if(byVertex)
      V.i <-   V.post.mean[graph$E[i,]]
    else{
      V.i <-   V.post.mean[4*(i-1) + 1:4]
    }
    ind <- which(graph$PtE[,1] == i)
    if(length(ind)==0){
      X.i   <- sample.line(theta,
                           V.i,
                           graph$edge_lengths[i],
                           nt = 100,
                           sample=FALSE)
    }else{
      X.i   <- sample.line(theta,
                           V.i,
                           graph$edge_lengths[i],
                           nt = 100,
                           y = graph$y[ind],
                           py =graph$PtE[ind,2],
                           sample=FALSE)

    }
    print(X.i)
    X.i[,2] = X.i[,2] +X.intercept
    if(i ==1){
      gg <- plot_curve(X.i  , graph$Lines[i,], normalized=FALSE,  ...)
    }else{
      gg <- plot_curve(X.i, graph$Lines[i,], normalized=FALSE, gg = gg,  ...)
    }


  }
  return(gg)

}

#' plot the observation
#' @param graph (graph.obj)
#' @param yhat - (n x 1) prediction of the obsveration
#' @export
plot_obs<- function(graph, y, y_loc, size_path=0.1, size_obs=1){
  xyl <- c()
  ind <- 1:dim(graph$EtV)[1]
  for(i in ind){

    xyl <- rbind(xyl, cbind(fortify(graph$Lines[i,])[,1:2],i))
  }
  obs <- y
  fig <- ggplot()+
    geom_path(data= data.frame(x=xyl[,1], y=xyl[,2], group=xyl[,3]), mapping = aes(x=x, y=y, group=group), size = size_path)+
    geom_point(data=data.frame(x=y_loc[,2],y=y_loc[,3],obs=obs),mapping= aes(y, x, color = obs,), size=size_obs)

    return(fig)
}
