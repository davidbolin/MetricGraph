
library(ggplot2)

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
#' @param V            (2 x 2) the position of the two vertices in 2D
#' @param gg           (ggplot pbj) add to gg obj if not null
#' @export
#'
plot_straight_curve <- function(data.to.plot, V ,gg = NULL,...){

  data.to.plot.order <- data.to.plot[order(data.to.plot[,1]),]
  l <- sqrt(sum( (V[1,] - V[2,])^2) )
  alpha <- data.to.plot.order[,1]/l
  coords <- cbind((1-alpha) * V[1,1] + alpha * V[2,1],
                  (1-alpha) * V[1,2] + alpha * V[2,2])
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

#' plot generic data object X on the graph
#'
#' @param X            (m x 3)  [, 1] edge number
#'                              [, 2] position on curve (in length)
#'                              [, 3] value
#' @param graph        (graph.object)
#' @param gg           (ggplot pbj) add to gg obj if not null
#' @export
plot_X_to_graph<- function(X, graph,Lines=NULL ,gg = NULL,...){

  for(i in 1:dim(graph$EtV)[1]){
    if(is.null(Lines) == TRUE)
      gg <- plot_straight_curve(X[X[,1]==graph$EtV[i,1],2:3], graph$V[graph$EtV[i,2:3],],gg = gg, ...)
  }

  return(gg)
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
      V.i <-   V.post.mean[graph$EtV[i,2:3]]
    else{
      V.i <-   V.post.mean[4*(i-1) + 1:4]
    }
    ind <- which(graph$PtE[,1] == i)
    if(length(ind)==0){
      X.i   <- sample.line(theta,
                           V.i,
                           graph$El[i],
                           nt = 100,
                           sample=F)
    }else{
      X.i   <- sample.line(theta,
                           V.i,
                           graph$El[i],
                           nt = 100,
                           y = graph$y[ind],
                           py =graph$PtE[ind,2],
                           sample=F)

    }
    print(X.i)
    X.i[,2] = X.i[,2] +X.intercept
    if(i ==1){
      gg <- plot_curve(X.i  , graph$Lines[i,], normalized=F,  ...)
    }else{
      gg <- plot_curve(X.i, graph$Lines[i,], normalized=F, gg = gg,  ...)
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
