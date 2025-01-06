


#'
#' function to export a stlpp object, from stlnpp, to metric graph object
#'
#'@export
stlpp.to.graph <- function(stlpp.obj){
  graph <- linnet.to.graph(stlpp.obj$domain)
  data <- as.data.frame(stlpp.obj$data)
  graph$add_observations(data = data.frame(y = rep(1, dim(data)[1]),
                                           time = data[,3],
                                           coord_x = data[,1],
                                           coord_y = data[,2]),
                         coord_x= "coord_x",
                         coord_y = "coord_y",
                          data_coords ="spatial")
  return(graph)
}

#' from
#' function to export a linnet object, from spatstat package, to metric graph object
#'
#'@export
linnet.to.graph <- function(linnet.object){
  n <- length(linnet.object$from)
  vertices <- as.data.frame(linnet.object$vertices)
  lines <- vector("list", n)
  for(i in 1:n){
    lines[[i]] <- rbind(vertices[linnet.object$from[i],1:2],
                        vertices[linnet.object$to[i],1:2])
  }
  return(metric_graph$new(edges = lines))
}

#'
#' function to export a psp object, from spatstat package, to metric graph object
#'
#'@export
psp.to.graph <- function(psp.object){

  n <- dim(psp.object$ends)[1]
  lines <- vector("list", n)
  for(i in 1:n){
    lines[[i]] <- rbind(as.matrix(psp.object$ends[i,1:2]),as.matrix(psp.object$ends[i,3:4]))
  }
  graph <- metric_graph$new(edges = lines)
  return(graph)
}
