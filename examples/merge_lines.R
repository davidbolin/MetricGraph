library(sp)
library(MetricGraph)
line1 <- Line(rbind(c(0,0),c(1,0)))
line2 <- Line(rbind(c(0,0),c(0,1)))
line3 <- Line(rbind(c(0,1),c(-1,1)))
theta <- seq(from=pi,to=3*pi/2,length.out = 20)
line4 <- Line(cbind(sin(theta),1+ cos(theta)))
lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                              Lines(list(line2),ID="2"),
                              Lines(list(line3),ID="3"),
                              Lines(list(line4),ID="4")))

graph <- metric_graph$new(lines = lines)
graph$plot()


graph <- metric_graph$new(lines = lines, remove_deg2 = TRUE)
graph$plot()


remove.deg2 <- function(graph, j = 1) {
  ind <- which(graph$get_degrees()==2)[j]
  e1 <- which(graph$E[,2]==ind)
  e2 <- which(graph$E[,1]==ind)
  e_rem <- sort(c(e1,e2))
  v1 <- setdiff(graph$E[e_rem[1],],ind)
  v2 <- setdiff(graph$E[e_rem[2],],ind)
  if(v1 > ind) {
    v1 <- v1-1
  }
  if(v2 > ind) {
    v2 <- v2 - 1
  }

  line_keep1 <- line_keep2 <- NULL
  if(e_rem[1]>1) {
    line_keep1 <- graph$lines[1:(e_rem[1]-1)]
  }
  if(e_rem[2]<length(graph$lines)) {
    line_keep2 <- graph$lines[(e_rem[2]+1):length(graph$lines)]
  }

  line_merge <- list()
  coords <- graph$lines@lines[[e_rem[1]]]@Lines[[1]]@coords
  tmp <- graph$lines@lines[[e_rem[2]]]@Lines[[1]]@coords
  diff_ss <- norm(as.matrix(coords[1,] - tmp[1,]))
  diff_se <- norm(as.matrix(coords[1,] - tmp[dim(tmp)[1],]))
  diff_es <- norm(as.matrix(coords[dim(coords)[1],] - tmp[1,]))
  diff_ee <- norm(as.matrix(coords[dim(coords)[1],] - tmp[dim(tmp)[1],]))
  diffs <- c(diff_ss, diff_se, diff_es, diff_ee)
  if(which.min(diffs) == 1) {
    coords <- rbind(coords[rev(1:dim(coords)[1]),], tmp)
    E_new <- c(v2,v1)
  } else if(which.min(diffs)==2){
    coords <- rbind(tmp,coords)
    E_new <- c(v2,v1)
  } else if(which.min(diffs)==3) {
    coords <- rbind(coords, tmp)
    E_new <- c(v1,v2)
  } else {
    coords <- rbind(coords, tmp[rev(1:dim(tmp)[1]),])
    E_new <- c(v1,v2)
  }
  line_merge <-  Lines(list(Line(coords)), ID = sprintf("new%d",i))

  if(!is.null(line_keep1) && !is.null(line_keep2)) {
    line_new <- SpatialLines(c(line_keep1@lines, line_merge, line_keep2@lines))
  } else if (is.null(line_keep1)) {
    line_new <- SpatialLines(c(line_merge, line_keep2@lines))
  } else if (is.null(line_keep2)) {
    line_new <- SpatialLines(c(line_keep1@lines, line_merge))
  } else {
    line_new <- SpatialLines(c(line_merge))
  }
  for(i in 1:length(line_new)) {
    slot(line_new@lines[[i]],"ID") <- sprintf("%d",i)
  }

  #update lines
  graph$lines <- line_new

  #update vertices
  graph$V <- graph$V[-ind,]
  graph$nV <- graph$nV - 1

  #update edges
  graph$E[graph$E >= ind] <- graph$E[graph$E >= ind] - 1
  graph$E <- graph$E[-e_rem[2],]
  graph$E[e_rem[1],] <- E_new
  graph$EID <- graph$EID[-ind]
  graph$edge_lengths[e_rem[1]] <- graph$edge_lengths[e_rem[1]] + graph$edge_lengths[e_rem[2]]
  graph$edge_lengths <- graph$edge_lengths[-e_rem[2]]
  graph$ELend <- graph$ELend[-e_rem[2]]
  graph$ELstart <- graph$ELstart[-e_rem[2]]
  graph$nE <- graph$nE - 1
  graph$LtE <- graph$LtE[-e_rem[2],-e_rem[2]]
  return(graph)
}

lines <- logo_lines()
graph <- metric_graph$new(lines = lines)
graph$plot(degree = TRUE)

graph <- remove.deg2(graph)
graph$plot(degree = TRUE)

graph <- remove.deg2(graph)
graph$plot(degree = TRUE)

graph <- remove.deg2(graph)
graph$plot(degree = TRUE)

graph <- remove.deg2(graph)
graph$plot(degree = TRUE)

graph <- remove.deg2(graph)
graph$plot(degree = TRUE)
