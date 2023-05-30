

## Functions snapPointsToLines, nearestPointOnLine and nearestPointOnSegment are
## from maptools package - under GPL-2 license
##
## We are providing them here since the maptools package will be deprecated in 2023
## However, if they are merged to sp, we will remove these functions and use the
## corresponding ones in sp.
##
## Authors: Roger Bivand, Nicholas Lewin-Koh, Edzer Pebesma, Eric Archer,
## Adrian Baddeley, Nick Bearman, Hans-Jörg Bibiko, Steven Brey, Jonathan Callahan,
## German Carrillo , Stéphane Dray , David Forrest , Michael Friendly , Patrick Giraudoux ,
## Duncan Golicher , Virgilio Gómez Rubio , Patrick Hausmann , Karl Ove Hufthammer , Thomas Jagger ,
## Kent Johnson , Matthew Lewis ORCID iD , Sebastian Luque , Don MacQueen , Andrew Niccolai ,
## Edzer Pebesma , Oscar Perpiñán Lamigueiro , Ethan Plunkett , Ege Rubak ORCID iD , Tom Short ,
## Greg Snow , Ben Stabler , Murray Stokely , Rolf Turner

## https://cran.r-project.org/src/contrib/Archive/maptools/


#' Snap a set of points to a set of lines
#'
#' This function snaps a set of points to a set of lines based on the minimum
#' distance of each point to any of the lines. This function does not work with
#' geographic coordinates. (Function from maptools package - under GPL-2)
#'
#' @param points An object of the class SpatialPoints or SpatialPointsDataFrame.
#' @param lines An object of the class SpatialLines or SpatialLinesDataFrame.
#' @param maxDist Numeric value for establishing a maximum distance to avoid snapping points that
#' are farther apart; its default value is NA.
#' @param withAttrs Boolean value for preserving (TRUE) or getting rid (FALSE) of the original point
#' attributes. Default: TRUE. This parameter is optional.
#' @param idField A string specifying the field which contains each line's id. This id will
#' be transferred to the snapped points data set to distinguish the line which each
#' point was snapped to.
#'
#' @noRd
snapPointsToLines <- function( points, lines, maxDist=NA, withAttrs=TRUE, idField=NA) {

    if (is(points, "SpatialPoints") && missing(withAttrs))
        withAttrs = FALSE

    if (is(points, "SpatialPoints") && withAttrs==TRUE)
        stop("A SpatialPoints object has no attributes! Please set withAttrs as FALSE.")

    # d = rgeos::gDistance(points, lines, byid=TRUE)
    d = distance2(points, lines, byid=TRUE)
    if(!is.na(maxDist)){
      distToLine <- apply(d, 2, min, na.rm = TRUE)
      validPoints <- distToLine <= maxDist  # indicates which points are within maxDist of a line
      distToPoint <- apply(d, 1, min, na.rm = TRUE)
      validLines <- distToPoint <= maxDist

      # Drop elements beyond maxdist
      points <- points[validPoints, ]
      lines = lines[validLines, ]
      d  = d[ validLines,  validPoints, drop = FALSE]
      distToLine <- distToLine[validPoints]

      # If no points are within maxDist return an empty SpatialPointsDataFrame object
      if(!any(validPoints)){
        if(is.na(idField)){
          idCol = character(0)
        } else {
          idCol = lines@data[,idField][0]
        }
        newCols = data.frame(nearest_line_id  = idCol, snap_dist = numeric(0))
        if(withAttrs) df <- cbind(points@data, newCols) else df <- newCols
        res <- SpatialPointsDataFrame(points, data=df,
                               proj4string=CRS(sp::proj4string(points)), match.ID = FALSE)
        return(res)
      }

    } else { # If no maxDist arg still calculate distToLine so it can be returned
      distToLine = apply(d, 2, min, na.rm = TRUE)
    }

    nearest_line_index = apply(d, 2, which.min) # Position of each nearest line in lines object

    coordsLines = coordinates(lines)
    coordsPoints = coordinates(points)

    # Get coordinates of nearest points lying on nearest lines
    mNewCoords = vapply(1:length(points),
        function(x)
            nearestPointOnLine(coordsLines[[nearest_line_index[x]]][[1]],
                coordsPoints[x,]), FUN.VALUE=c(0,0))

    # Recover lines' Ids (If no id field has been specified, take the sp-lines id)
    if (!is.na(idField)) {
      nearest_line_id = lines@data[,idField][nearest_line_index]
    }  else {
      nearest_line_id = sapply(slot(lines, "lines"), function(i) slot(i, "ID"))[nearest_line_index]
    }
    # Create data frame and sp points
    if (withAttrs) df = cbind(points@data, data.frame(nearest_line_id, snap_dist = distToLine))
    else df = data.frame(nearest_line_id, snap_dist = distToLine, row.names=names(nearest_line_index))

    SpatialPointsDataFrame(coords=t(mNewCoords), data=df,
        proj4string=sp::CRS(sp::proj4string(points)))
}

#' @noRd
nearestPointOnLine <- function (coordsLine, coordsPoint)
{
    nearest_points = vapply(2:nrow(coordsLine), function(x) nearestPointOnSegment(coordsLine[(x -
        1):x, ], coordsPoint), FUN.VALUE = c(0, 0, 0))
    nearest_points[1:2, which.min(nearest_points[3, ])]
}

#' @noRd
nearestPointOnSegment <- function(s, p){
    # Adapted from http://pastebin.com/n9rUuGRh
    ap = c(p[1] - s[1,1], p[2] - s[1,2])
    ab = c(s[2,1] - s[1,1], s[2,2] - s[1,2])
    t = sum(ap*ab) / sum(ab*ab)
    t = ifelse(t<0,0,ifelse(t>1,1,t))
    t = ifelse(is.na(t), 0, t) # when start and end of segment are identical t is NA
    x = s[1,1] + ab[1] * t
    y = s[1,2] + ab[2] * t
    result = c(x, y, sqrt((x-p[1])^2 + (y-p[2])^2))  # Return nearest point and distance
    names(result) = c("X","Y","distance")
    result
}