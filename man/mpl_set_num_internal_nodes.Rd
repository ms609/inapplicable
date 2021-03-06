% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/morphyex.R
\name{mpl_set_num_internal_nodes}
\alias{mpl_set_num_internal_nodes}
\title{Sets the number of internal nodes in the dataset}
\usage{
mpl_set_num_internal_nodes(numnodes, morphyobj)
}
\arguments{
\item{morphyobj}{An instance of the Morphy object.}

\item{nnodes}{The desired number of internal nodes.}
}
\value{
A Morphy error code.
}
\description{
This specifies the number of internal nodes over which
reconstruction sets need to be made. It is up to the caller to ensure the 
correct number of nodes and the relationships between them.
}
\author{
Martin Brazeau
}
\keyword{internal}
