rm(list=ls())
# test vertex vertex
line1 <- Line(rbind(c(0,0),c(1,0)))
line2 <- Line(rbind(c(0,0.03),c(0,0.06)))
line3 <- Line(rbind(c(0,0.09),c(0,0.97)))
line4 <- Line(cbind(seq(from=0,to=1, length.out = 100), rep(1,100)))
line5 <- Line(rbind(c(1,0.03),c(1,0.97)))
line6 <- Line(rbind(c(0.97,0.03),c(0.03,0.97)))

theta <- seq(from=0,to=2*pi,length.out = 50)
line7 <- Line(cbind(0.3*sin(theta),-0.33+0.3*cos(theta)))

lines = SpatialLines(list(Lines(list(line1),ID="1"),
                          Lines(list(line2),ID="2"),
                          Lines(list(line3),ID="3"),
                          Lines(list(line4),ID="4"),
                          Lines(list(line5),ID="5"),
                          Lines(list(line6),ID="6"),
                          Lines(list(line7),ID="7")))
graph <- metric_graph$new(lines, tolerance = list(vertex_vertex = 0.2))
graph$plot()

#test line vertex

line1 <- Line(rbind(c(0,0),c(1,0)))
line2 <- Line(rbind(c(0.5,0.01),c(0.5,0.97)))
line3 <- Line(cbind(seq(from=0,to=1, length.out = 100), rep(1,100)))
line4 <- Line(rbind(c(1,0.03),c(1,0.97)))
theta <- seq(from=0,to=2*pi,length.out = 50)
line5 <- Line(cbind(0.5+0.3*sin(theta),-0.31+0.3*cos(theta)))

lines = SpatialLines(list(Lines(list(line1),ID="1"),
                          Lines(list(line2),ID="2"),
                          Lines(list(line3),ID="3"),
                          Lines(list(line4),ID="4"),
                          Lines(list(line5),ID="5")))
graph <- metric_graph$new(lines, tolerance = list(vertex_vertex = 0.2,
                                                  vertex_line = 0.2))
graph$plot()

# test line line

find_line_line_points <- function(graph, tol) {
  dists <- gWithinDistance(lines, dist = tol, byid = TRUE)
  points_add <- NULL
  points_add_PtE <- NULL
  for(i in 1:(length(lines)-1)) {
    #lines within tol of line i
    inds <- i+which(as.vector(dists[i, (i+1):length(graph$lines)]))
    if(length(inds)>0) {
      for(j in inds) {
        #first check if there are intersections
        intersect_tmp <- rgeos::gIntersection(graph$lines[i], graph$lines[j])
        p_cur <- NULL
        if(!is.null(intersect_tmp)) {
          coord_tmp <- coordinates(intersect_tmp)
          for(k in 1:length(intersect_tmp)) {
            p <- matrix(coord_tmp[k,],1,2)
            #add points if they are not close to V or previous points
            if(min(spDists(graph$V, p))>tol) {

              p_cur <- rbind(p_cur,p)
              p2 <- snapPointsToLines(SpatialPoints(p),graph$lines[i])
              points_add <- rbind(points_add, p, coordinates(p2))
              points_add_PtE <- rbind(points_add_PtE,
                                      c(i,gProject(graph$lines[i],
                                                   SpatialPoints(p))),
                                      c(j,gProject(graph$lines[j],SpatialPoints(p))))

            }
          }
        }
        #now check if there are intersections with buffer
        tmp_line <- gBuffer(graph$lines[i], width = tol)
        intersect_tmp <- rgeos::gIntersection(tmp_line, graph$lines[j])
        if(!is.null(intersect_tmp)) {
          for(k in 1:length(intersect_tmp)) {
            if(class(intersect_tmp) == "SpatialLines") {
              coord_tmp <-gInterpolate(intersect_tmp[k], d=0.5, normalized = TRUE)
              p <- matrix(coordinates(coord_tmp),1,2)
            } else {
              p <- matrix(coordinates(intersect_tmp[k]),1,2)
            }
            #add points if they are not close to V or previous points
            if(min(spDists(graph$V, p))>tol) {
              if(is.null(p_cur) || gDistance(SpatialPoints(p_cur), intersect_tmp[k])>tol) {
                p2 <- snapPointsToLines(SpatialPoints(p),graph$lines[i])
                points_add <- rbind(points_add, p, coordinates(p2))
                points_add_PtE <- rbind(points_add_PtE,
                                        c(i,gProject(graph$lines[i],SpatialPoints(p))),
                                        c(j,gProject(graph$lines[j],SpatialPoints(p))))

              }
            }
          }
        }
      }
    }
  }
  return(list(points = points_add, PtE = points_add_PtE))
}

add_vertices = function(self, PtE, tolerance = 1e-10) {

  e.u <- unique(PtE[,1])
  for (i in 1:length(e.u)) {
    dists <- sort(PtE[which(PtE[,1]==e.u[i]),2])
    self$split_edge(e.u[i], dists[1], tolerance)
    if(length(dists)>1) {
      dists_up <- dists
      for(j in 2:length(dists)){
        dists_up[j] <- (dists[j] - dists[j-1])/(1 - dists[j-1])
        self$split_edge(self$nE, dists_up[j], tolerance)
      }
    }
  }
  return(self)
}

line1 <- Line(rbind(c(1,0),c(2,0)))
line2 <- Line(rbind(c(1.5,-1),c(1.5,1)))
theta <- seq(from=pi,to=2*pi,length.out = 50)
line3 <- Line(cbind(1.05+sin(theta),cos(theta)))
line4 <- Line(rbind(c(0,1),c(0,-1)))

lines = SpatialLines(list(Lines(list(line1),ID="1"),
                          Lines(list(line2),ID="2"),
                          Lines(list(line3),ID="3"),
                          Lines(list(line4),ID="4")))
plot(lines)
graph <- metric_graph$new(lines, tolerance = list(vertex_vertex = 0.05,
                                                  vertex_line = 0.05))
graph$plot()

points_add <- find_line_line_points(graph, tol = 0.1)
PtE <- points_add$PtE
PtE[,2] <- PtE[,2]/graph$edge_lengths[PtE[,1]]
graph <- add_vertices(graph, PtE, tolerance = 0.1)
graph$plot(degree=TRUE)

#half circle line
line1 <- Line(rbind(c(1,0),c(2,0)))
line2 <- Line(rbind(c(1.5,-1),c(1.5,1)))
theta <- seq(from=pi,to=2*pi,length.out = 50)
line3 <- Line(cbind(1.05+sin(theta),cos(theta)))
line4 <- Line(rbind(c(0,1),c(0,-1)))

lines = SpatialLines(list(Lines(list(line1),ID="1"),
                          Lines(list(line2),ID="2"),
                          Lines(list(line3),ID="3"),
                          Lines(list(line4),ID="4")))
plot(lines)
graph <- metric_graph$new(lines, tolerance = list(vertex_vertex = 0.05))
points_add <- find_line_line_points(graph, tol = 0.1)
PtE <- points_add$PtE
PtE[,2] <- PtE[,2]/graph$edge_lengths[PtE[,1]]
graph <- add_vertices(graph, PtE, tolerance = 0.1)
graph$plot(degree=TRUE)


# circle line
line1 <- Line(rbind(c(-0.5,0),c(1,0)))
theta <- seq(from=0,to=2*pi,length.out = 50)
line2 <- Line(cbind(1.05+sin(theta),cos(theta)))
line3 <- Line(rbind(c(0,1),c(0,-1)))

lines = SpatialLines(list(Lines(list(line1),ID="1"),
                          Lines(list(line2),ID="2"),
                          Lines(list(line3),ID="3")))
plot(lines)
graph <- metric_graph$new(lines, tolerance = list(vertex_vertex = 0.2))
graph$plot()
points_add <- find_line_line_points(graph, tol = 0.1)
PtE <- points_add$PtE
PtE[,2] <- PtE[,2]/graph$edge_lengths[PtE[,1]]
graph <- add_vertices(graph, PtE, tolerance = 0.1)
graph$plot(degree=TRUE)


#circle line intersect
theta <- seq(from=0,to=2*pi,length.out = 50)
line1 <- Line(cbind(sin(theta),cos(theta)))
line2 <- Line(rbind(c(-0.75,1),c(-0.75,-1)))
lines = SpatialLines(list(Lines(list(line1),ID="1"),
                          Lines(list(line2),ID="2")))
plot(lines)
graph <- metric_graph$new(lines, tolerance = list(vertex_vertex = 0.2))
graph$plot()
points_add <- find_line_line_points(graph, tol = 0.1)
PtE <- points_add$PtE
PtE[,2] <- PtE[,2]/graph$edge_lengths[PtE[,1]]
graph <- add_vertices(graph, PtE, tolerance = 0.1)
graph$plot(degree=TRUE)


#two circles
theta <- seq(from=0,to=2*pi,length.out = 50)
line1 <- Line(cbind(sin(theta),cos(theta)))
line2 <- Line(cbind(2.01+sin(theta),cos(theta)))
line3 <- Line(rbind(c(-0.75,1),c(-0.75,-1)))
lines = SpatialLines(list(Lines(list(line1),ID="1"),
                          Lines(list(line2),ID="2"),
                          Lines(list(line3),ID="3")))
plot(lines)
graph <- metric_graph$new(lines, tolerance = list(vertex_vertex = 0.1))
graph$plot()
points_add <- find_line_line_points(graph, tol = 0.1)
PtE <- points_add$PtE
PtE[,2] <- PtE[,2]/graph$edge_lengths[PtE[,1]]
graph <- add_vertices(graph, PtE, tolerance = 0.1)
graph$plot(degree=TRUE)


#complicated with multiple circles
theta <- seq(from=0,to=2*pi,length.out = 50)
line1 <- Line(cbind(sin(theta),cos(theta)))
line2 <- Line(cbind(2.05+sin(theta),cos(theta)))
line3 <- Line(cbind(4.15+sin(theta),cos(theta)))
line4 <- Line(rbind(c(3.1,1),c(3.1,-1)))
line5 <- Line(rbind(c(-0.75,1),c(-0.75,-1)))
lines = SpatialLines(list(Lines(list(line1),ID="1"),
                          Lines(list(line2),ID="2"),
                          Lines(list(line3),ID="3"),
                          Lines(list(line4),ID="4"),
                          Lines(list(line5),ID="5")))
plot(lines)
graph <- metric_graph$new(lines, tolerance = list(vertex_vertex = 0.1))
graph$plot()
points_add <- find_line_line_points(graph, tol = 0.2)
PtE <- points_add$PtE
PtE[,2] <- PtE[,2]/graph$edge_lengths[PtE[,1]]
graph <- add_vertices(graph, PtE, tolerance = 0.2)
graph$plot(degree=TRUE)
