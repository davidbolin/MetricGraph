
library(ggplot2)

#' plot a continous curve
#'
#' @param data.to.plot (m x 2) [,1] position on curve (in length)
#'                             [,2] value
#' @param Line         (SpatialLinesDataFrame,SpatialLines)
#'
#' @export
#'
plot_curve <- function(data.to.plot, Line_edge, gg = NULL){

  data.to.plot.order <- data.to.plot[order(data.to.plot[,1]),]
  p2 <- rgeos::gInterpolate(Line_edge, data.to.plot.order[,1])
  coords <-p2@coords
  if(is.null(gg)){
    gg <- ggplot2::ggplot(data.frame(x = coords[,1],
                                     y = coords[,2],
                                     val = data.to.plot.order[,2]),
                          ggplot2::aes(x=x, y=y, colour=val ) ) +
      ggplot2::geom_path(size=2)
  }else{
    gg <- gg + ggplot2::geom_path(data = data.frame(x = coords[,1],
                                                    y = coords[,2],
                                                    val = data.to.plot.order[,2]),
                                  size=2)
  }
  return(gg)
}
