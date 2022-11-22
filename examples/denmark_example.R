
largest_component <- function(Lines) {
  graph <- metric_graph$new(Lines = Lines)
  g <- graph(edges = c(t(graph$E)), directed = FALSE)
  components <- igraph::clusters(g, mode="weak")
  biggest_cluster_id <- which.max(components$csize)

  vert_ids <- V(g)[components$membership == biggest_cluster_id]
  edge_rem <- NULL
  for(i in 1:graph$nE){
    if(!(graph$E[i,1] %in% vert_ids) && !(graph$E[i,2] %in% vert_ids))
      edge_rem <- c(edge_rem, i)
  }
  edge_keep <- setdiff(1:graph$nE, edge_rem)
  return(Lines[edge_keep])
}

library(GPGraph)
library(ggplot2)
library(maptools)
library(Matrix)
library(igraph)
library(ggplot2)
library(osmdata)
set_overpass_url("https://maps.mail.ru/osm/tools/overpass/api/interpreter")

call <- opq(bbox = c(12.5,55.7,12.6,55.8))
call <- add_osm_feature(call, key = "highway",value=c("motorway",
                                                      "primary",
                                                      "secondary",
                                                      "tertiary"))

data <- osmdata_sp(call)

Lines <- largest_component(SpatialLines(data$osm_lines@lines))

graph <- metric_graph$new(Lines = Lines)

#split circular edges (only needed for alpha = 2)
if(0){
  circ.ind <- which(graph$E[,1]==graph$E[,2])
  PtE <- cbind(circ.ind, rep(0.5,length(circ.ind)))
  graph$add_observations2(y = rep(0, length(circ.ind)), PtE, normalized = TRUE)
  graph$observation_to_vertex()
  graph$clear_observations()
}

p <- graph$plot(marker_size = 0)

p <- p + theme_void() + theme(legend.position="none")
p

graph$build_mesh(h = 0.001)

type = "isoExp"
if(type == "isoExp") { #isotropic exponential
  graph$compute_resdist_mesh()
  Sigma <- exp(-50*graph$mesh$res.dist)
  u <- as.vector(t(chol(forceSymmetric(Sigma)))%*%rnorm(dim(Sigma)[1]))
} else if (type == "alpha1") { #alpha =1 model looks strange
  u <- sample_spde(kappa = 50, sigma = 1, alpha = 1,
                   graph = graph, type = "mesh")
} else if (type == "alpha2") { #not working
  u <- sample_spde(kappa = 50, sigma = 1, alpha = 2,
                   graph = graph, type = "mesh")
}

graph$plot_function_mesh(X = as.vector(u))
