\name{weighted_stiffness}
\alias{weighted_stiffness}
\title{Assemble a weighted stiffness matrix on a graph mesh}
\usage{
weighted_stiffness(edges, lengths, weights)
}
\arguments{
\item{edges}{Two-column integer matrix of mesh edge endpoints.}
\item{lengths}{Positive mesh edge lengths, in the same order as rows of edges.}
\item{weights}{Positive coefficient values at mesh nodes.}
}
\value{A symmetric sparse positive-semidefinite stiffness matrix.}
\description{
Uses a piecewise linear positive coefficient on each mesh edge. The
coefficient integral on an edge is its length times the mean of its two
endpoint values.
}
