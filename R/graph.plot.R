
library(ggplot2)

#' plot a continous curve
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

#'
#' @param graph          - graph.obj
#' @param posterior.mean - (function) takes theta, and graph computees posterior mean
#' @param sample.line    - (function) takes theta, and graph computees posterior mean
#' @param theta          - params
#' @param byVertex       - if byVertex V.i is by vertex ortherwise by edge
#' @export
plot_posterior_mean <- function(graph, posterior.mean, sample.line, theta, byVertex=T, ...){

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
    if(i ==1){
      gg <- plot_curve(X.i, graph$Lines[i,], normalized=F,  ...)
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
plot_obs<- function(graph, y, y_loc){
  xyl <- c()
  ind <- 1:dim(graph$EtV)[1]
  for(i in ind){

    xyl <- rbind(xyl, cbind(fortify(graph$Lines[i,])[,1:2],i))
  }
  obs <- y
  fig <- ggplot()+
    geom_path(data= data.frame(x=xyl[,1], y=xyl[,2], group=xyl[,3]), mapping = aes(x=x, y=y, group=group), size = 0.1)+
    geom_point(data=data.frame(x=y_loc[,2],y=y_loc[,3],obs=obs),mapping= aes(y, x, color = obs))

    return(fig)
}
