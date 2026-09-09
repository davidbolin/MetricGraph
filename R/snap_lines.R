

## Functions snapPointsToLines, nearestPointOnLine and nearestPointOnSegment are either obtained or modifications
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
snapPointsToLines <- function( points, lines, longlat, crs, idx = NULL) {
    if(!is.list(lines)){
      lines <- list(lines)
    }
    if(!is.null(points)){
      class(points) <- setdiff(class(points), "metric_graph_edge")
    }
    points <- as.matrix(points)
    storage.mode(points) <- "double"
    points <- points[, 1:2, drop = FALSE]

    # Not using longlat, since all the projections (to get the "closest" point
    # on the graph) are considering the Euclidean distances.
    #
    # The dense point-by-line distance matrix that used to drive this function
    # costs O(nPoints * nLines) time and memory, which is prohibitive for
    # city-sized graphs. `nearest_edge_cpp()` indexes the lines in a uniform
    # grid instead and returns the same nearest line, snapped coordinates and
    # distance.
    res <- nearest_edge_cpp(lines, points)

    nearest_line_index <- res[["index"]]
    distToLine <- res[["dist"]]
    mNewCoords <- res[["coords"]]
    # The previous implementation built the coordinates with `vapply()` over
    # `nearestPointOnLine()`, which named the rows; keep that for callers that
    # index the result by name.
    rownames(mNewCoords) <- c("X", "Y")

    if(!is.null(idx)){
      nearest_line_index <- idx
    }
    # Create data frame and sp points
    df = data.frame(nearest_line_index = nearest_line_index, snap_dist = distToLine)

    list(coords = mNewCoords, df = df)
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