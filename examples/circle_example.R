#' Circular `sp` SpatialPolygons
#'
#'   This routine will create an "sp" circular SpatialPolygons object based on
# ''  the arguments...
#'
#' @param radius the radius of the circle
#' @param spUnits the CRS units
#' @param centerPoint the circle's cente
#' @param nptsPerimeter the number of points forming the perimeter of the polygon
#' @param spID the spatial ID for the object
#' @param ... additional arguments to be passed.
#'
#' @return a list with...
#'       spCircle = the SpatialPolygons circular object
#'       location =  the SpatialPoints center point
#' @export
spCircle = function(radius,
                    spUnits = CRS(projargs=as.character(NA)),
                    centerPoint = c(x=0, y=0),   #centerPoint
                    nptsPerimeter = 100,
                    spID = paste('circle',round(1000*stats::runif(1)),sep=':'),
                    ...
                   )
{
#---------------------------------------------------------------------------
#
#   This routine will create an "sp" circular SpatialPolygons object based on
#   the arguments...
#
#   Arguments...
#     radius = the radius of the circle
#     spUnits = the CRS units
#     centerPoint = the circle's center
#     nptsPerimeter = the number of points forming the perimeter of the polygon
#     spID = the spatial ID for the object
#
#   Returns...
#     a list with...
#       spCircle = the SpatialPolygons circular object
#       location =  the SpatialPoints center point
#
#   Please note that you may want to rename the components of, e.g., the
#   spCircle object to make more sense. An example where it is used for
#   dbh is in the standingTree constructor, where spDBH == spCircle...
#
#       names(spDBH@polygons) = 'pgsDBH'
#       names(spDBH@polygons$pgsDBH@Polygons) = 'pgDBH'
#
#
#Author...									Date: 24-Oct-2011
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   make sure the center is a named vector of length 2...
#
    if(any(is.na(match(c('x','y'),names(centerPoint)))))
      stop('Please use names x and y for circle centerPoint vector')
    if(length(centerPoint) != 2)
      stop('Please supply one set of (x,y) coordinates for the plot center location.')


#
#   some other checks...
#
    if(radius <= 0)
      stop('radius must be positive!')
    if(nptsPerimeter < 20) {
      warning('Using less than 20 points for the circle perimeter is not recommended--set to 20')
      nptsPerimeter = 20
    }

    area = pi*radius*radius

    ##location = centerPoint

#
#   left half of the circle, then right...
#
    circ =  seq(0, 2*pi, len=nptsPerimeter)

#
#   make the circle outline...
#
    circle = matrix(c(centerPoint['x'] + radius*cos(circ),
                    centerPoint['y'] + radius*sin(circ), rep(1,nptsPerimeter) ), nrow=nptsPerimeter)

#
#   any little difference between start & end pts with identical() can mess up the
#   the sp package Polygon routine, so set the end point exactly to start, then transform...
#
    circle = rbind(circle, circle[1,])
    dimnames(circle) = list(NULL,c('x','y','hc'))

#
#   and make a SpatialPolygons object...
#
    pgCircle = Polygon(circle[,-3])                               #sans hc
    pgsCircle = Polygons(list(circPlot = pgCircle), ID = spID)
    spCircle = SpatialPolygons(list(pgsCircle = pgsCircle))       #takes a list of Polygons objects

#
#   no id for center point, but it can be added to be the same as spID when
#   we make a container class for the center points elsewhere...
#
    loc = matrix(centerPoint, nrow=1)
    colnames(loc) = names(centerPoint)
    location = SpatialPoints(loc, proj4string = spUnits)


    return( list(spCircle=spCircle, location=location) )
}   #spCircle








graphics.off()
library(GPGraph)
library(sp)
library(maptools)
dbh = 20
kappa <- 0.5
sp.dbh = spCircle(dbh/2, centerPoint=c(x=30,y=80), spID='tree.1')
line.cirlce <- Line(sp.dbh$spCircle@polygons[[1]]@Polygons$circPlot@coords)
line.line <- Line(rbind(c(30,80),c(40,80)))
graph <- graph.obj$new(sp::SpatialLines(list(Lines(list(line.cirlce),ID="1"))))
graph2 <-  graph.obj$new(sp::SpatialLines(list(Lines(list(line.cirlce),ID="1"),
                                               Lines(list(line.line),ID="2"))))


Q <- Q.exp(kappa, graph$V, graph$EtV, graph$El)
Q2 <- Q.exp(kappa, graph2$V, graph2$EtV, graph2$El)

#https://tmieno2.github.io/R-as-GIS-for-Economists/

#simulate Expontial on a Line given end points
#simulate Expontial on a Line given end points and observations

Q_ <- Q.exp.line(kappa, c(0,30,64))

xc = c(1.2,31)
yc = c(1.5,80)
Spoints = SpatialPoints(cbind(xc, yc))
Spoints = SpatialPointsDataFrame(Spoints,  data.frame(a=1:2))
SP <- snapPointsToLines(Spoints, graph2$Lines)

graph2$add_observations(Spoints, c(1,2))

X <- sample.line.expontial(c(1,1,kappa),
                           c(0,0),
                           Line = graph2$Lines[1,],
                           graph2$El[1])
fig <- plot_curve(X, graph2$Lines[1,])
print(fig)

X2 <- sample.line.expontial(c(0.0001,1,kappa),
                           c(0,0),
                           Line = graph2$Lines[1,],
                           graph2$El[1],
                           py = graph2$PtE[1,2],
                           y  = 10
                           )

fig <- plot_curve(X2, graph2$Lines[1,])
print(fig)
