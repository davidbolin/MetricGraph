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
