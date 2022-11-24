#internal function for checking metric_graph inputs
gpgraph_check_graph <- function(graph)
{
  if (!inherits(graph, "GPGraph::graph")) {
    stop("The graph object is not a metric graph")
  }
  out <- list(has.mesh = FALSE,
              has.obs = FALSE)
  if(!is.null(graph$mesh)){
    out$has.mesh = TRUE
  }
  if(!is.null(graph$y) && !is.null(graph$PtE))
    out$has.data = TRUE
  return(out)
}


#' Create metric graphs for connected components of a SpatialLines object
#'
#' @param lines Object of class `SpatialLines`
#' @param by_length Sort the components by total edge length? If FALSE,
#' the components are sorted by the number of vertices.
#' @param only_largest if TRUE, only return the largest component.
#' Otherwise return an ordered list with the components (largest first)
#'
#' @return A `metric_graph` object created from the largest component, or a
#' list of `metric_graph` objects for all connected components
#' @export
#'
#' @examples
#' library(sp)
#' line1 <- Line(rbind(c(0, 0), c(1, 0)))
#' line2 <- Line(rbind(c(1, 0), c(2, 0)))
#' line3 <- Line(rbind(c(1, 1), c(2, 1)))

#' Lines <-  SpatialLines(list(Lines(list(line1), ID = "1"),
#'                            Lines(list(line2), ID = "2"),
#'                            Lines(list(line3), ID = "3")))
#' graphs <- graph_components(Lines, only_largest = FALSE)
#' p <- graphs[[1]]$plot(edge_color = "red")
#' graphs[[2]]$plot(p = p, edge_color = "blue")
graph_components <- function(lines, by_length = TRUE, only_largest = TRUE) {

  graph <- metric_graph$new(lines = lines)
  g <- graph(edges = c(t(graph$E)), directed = FALSE)
  igraph::E(g)$weight <- graph$edge_lengths
  components <- igraph::clusters(g, mode="weak")
  Graphs <- list()
  nc <- components$no
  for(k in 1:nc) {
    vert_ids <- igraph::V(g)[components$membership == k]
    edge_rem <- NULL
    for (i in 1:graph$nE) {
      if(!(graph$E[i, 1] %in% vert_ids) && !(graph$E[i, 2] %in% vert_ids))
        edge_rem <- c(edge_rem, i)
    }
    edge_keep <- setdiff(1:graph$nE, edge_rem)
    Graphs[[k]] = metric_graph$new(lines = Lines[edge_keep])
  }
  sizes <- components$csize
  lengths <- unlist(lapply(1:nc, function(x) sum(Graphs[[x]]$edge_lengths)))
  if(by_length) {
    reo <- sort(lengths, decreasing = TRUE, index.return = TRUE)$ix
  } else {
    reo <- sort(sizes, decreasing = TRUE, index.return = TRUE)$ix
  }
  g <- Graphs[reo]
  if(only_largest) {
    return(g[[1]])
  } else {
    return(g)
  }
}
