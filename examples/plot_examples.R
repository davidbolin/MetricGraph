library(GPGraph)
library(sp)
#create a simple graph
line1 <- Line(rbind(c(0,0),c(1,0)))
line2 <- Line(rbind(c(0,0),c(0,1)))
line3 <- Line(rbind(c(0,1),c(-1,1)))
theta <- seq(from=pi,to=3*pi/2,length.out = 20)
x <- sin(theta)
y <- 1+ cos(theta)
line4 <- Line(cbind(x,y))
Lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                              Lines(list(line2),ID="2"),
                              Lines(list(line3),ID="3"),
                              Lines(list(line4),ID="4")))
plot(Lines)
graph <- metric_graph$new(Lines = Lines)

lines <- c()
  for(i in 1:length(Lines)){
    points <- Lines@lines[[i]]@Lines[[1]]@coords
    n <- dim(points)[1]
    lines <- rbind(lines,c(i, points[1,], sp::LineLength( Lines@lines[[i]]@Lines[[1]])),
                         c(i, points[n,], sp::LineLength( Lines@lines[[i]]@Lines[[1]])))
  }

  index.dub <- duplicated(lines[,2:3])
  vertex <- cbind( 1:sum(c(!index.dub)), lines[!index.dub,2:3,drop=F])

  lvl <- matrix(0, nrow= max(lines[,1]), 4)
  for(i in 1:max(lines[,1])){
    which.line <- sort(which(lines[,1]==i))
    line <- lines[which.line,]
    ind1 <- (abs( vertex[,2] - line[1,2] )< 1e-10)* (abs( vertex[,3] - line[1,3] )< 1e-10)==1
    ind2 <-  (abs( vertex[,2] - line[2,2] )< 1e-10)* (abs( vertex[,3] - line[2,3] )< 1e-10)==1

    lvl[i,1] <- i
    lvl[i,2] <- which(ind1)
    lvl[i,3] <- which(ind2)
    lvl[i,4] <- line[1,4]
  }
