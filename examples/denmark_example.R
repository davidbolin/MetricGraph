library(hexSticker)
library(GPGraph)
library(ggplot2)
library(maptools)
library(Matrix)
library(igraph)
library(ggplot2)

library(osmdata)
set_overpass_url("https://maps.mail.ru/osm/tools/overpass/api/interpreter")


call <- opq(bbox = c(12.4,55.6,12.6,55.8))
call <- add_osm_feature(call, key = "highway",value=c("motorway",
                                                      "primary",
                                                      "secondary",
                                                      "tertiary"))

data <- osmdata_sp(call)
graph <- metric_graph$new(SpatialLines(data$osm_lines@lines))

#find largest connected component:
g <- graph(edges = c(t(graph$E)), directed = FALSE)
components <- igraph::clusters(g, mode="weak")
biggest_cluster_id <- which.max(components$csize)

vert_ids <- V(g)[components$membership == biggest_cluster_id]
g_up <- igraph::induced_subgraph(g, vert_ids)
E <- as_adj(g_up)
V <- graph$V[vert_ids,]
E <- get.edgelist(g_up)

#create new graph
graph <- metric_graph$new(P = V, E = E)

p <- graph$plot(marker_size = 0)

p <- p + theme_void() + theme(legend.position="none")
p

#sample function on mesh
graph$build_mesh(h = 0.02)

#plot


# plot_function_mesh_2d = function(X, graph, ...){
#   if(is.null(graph$mesh)) {
#     stop("no mesh provided")
#   }
#   if(length(X) != dim(graph$V)[1] + dim(graph$mesh$PtE)[1]){
#     stop("X does not have the correct size")
#   }

#   n.v <- dim(graph$V)[1]
#   XV <- X[1:n.v]
#   x.loc <- y.loc <- z.loc <- i.loc <- NULL
#   kk = 1
#   for(i in 1:graph$nE){
#     ind <- graph$mesh$PtE[,1] == i
#     if(sum(ind)==0){
#       vals <- rbind(c(0, XV[graph$E[i,1]]),
#                     c(1, XV[graph$E[i,2]]))

#     }else{
#       vals <- rbind(c(0, XV[graph$E[i,1]]),
#                     cbind(graph$mesh$PtE[ind,2],X[n.v + which(ind)]),
#                     c(1, XV[graph$E[i,2]]))

#     }
#     if(is.null(graph$Lines) == TRUE) {
#       data.to.plot.order <- vals[order(vals[, 1]), ]
#       V <- graph$V[graph$E[i,],]
#       l <- 1#sqrt(sum((V[1,] - V[2,])^2))
#       if(l > 0){
#         alpha <- data.to.plot.order[,1]/l
#         coords <- cbind((1 - alpha) * V[1, 1] + alpha * V[2, 1],
#                         (1 - alpha) * V[1, 2] + alpha * V[2, 2])

#         x.loc = c(x.loc, coords[, 1])
#         y.loc = c(y.loc, coords[, 2])
#         z.loc = c(z.loc, data.to.plot.order[, 2])
#         i.loc = c(i.loc, rep(kk, length(coords[, 1])))
#         kk = kk+1
#       }
#     } else {
#       index <- (graph$LtE@p[i]+1):(graph$LtE@p[i+1])
#       LinesPos <- cbind(graph$LtE@i[index] + 1, graph$LtE@x[index])
#       LinesPos <- LinesPos[order(LinesPos[,2]),,drop = FALSE]
#       for(j in 1:length(index)){
#         if(j==1){
#           index_j <- vals[,1] <= LinesPos[j,2]
#         }else{
#           index_j <- (vals[,1] <= LinesPos[j,2]) &  (vals[,1] > LinesPos[j-1,2])
#         }
#         if(sum(index_j) == 0)
#           next
#         rel.pos = vals[index_j,1]
#         if(j == 1){
#           rel.pos <- rel.pos/LinesPos[j,2]
#         }else{
#           rel.pos <- (rel.pos-LinesPos[j-1,2])/(LinesPos[j,2]-LinesPos[j-1,2])
#         }

#         if(j== dim(LinesPos)[1] )
#           rel.pos = graph$ELend[i]*rel.pos
#         if(j==1)
#           rel.pos = rel.pos + graph$ELstart[i]

#         data.to.plot <- cbind(rel.pos,vals[index_j,2])
#         Line_edge <- SpatialLines(list(graph$Lines@lines[[LinesPos[j,1]]]))

#         data.to.plot.order <- data.to.plot[order(data.to.plot[, 1]), , drop = FALSE]
#         p2 <- rgeos::gInterpolate(Line_edge, data.to.plot.order[, 1, drop = FALSE],
#                                   normalized = TRUE)
#         coords <-p2@coords
#         x.loc <- c(x.loc, coords[,1])
#         y.loc <- c(y.loc, coords[,2])
#         z.loc <- c(z.loc, data.to.plot.order[,2])
#         i.loc <- c(i.loc, rep(kk, length(coords[,1])))
#         kk = kk+1
#       }
#     }
#   }
#   data <- data.frame(x = x.loc, y = y.loc, z = z.loc, i = i.loc)

#   p <- ggplot(data = data, aes(x = x, y = y, group = i, colour = z)) + geom_path() + scale_color_viridis()

#   return(p)
# }
u <- sample_spde(kappa = 1, sigma = 1, alpha = 1,
                 graph = graph, type = "mesh")

# plot_function_mesh_2d(u, graph)

graph$plot_function_mesh(X = u)
