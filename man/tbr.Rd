\name{tbr}
\alias{tbr}

\title{Tree bisection and reconnection}
\description{This function performs a single random \acronym{TBR} iteration.}
\usage{
tbr(tree, edge.to.break = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{a fully resolved tree in \code{\link{phyDat}} format;}
  \item{edge.to.break}{the index of an edge to bisect, generated randomly if not specified.}
}
\details{
Branch lengths are not (yet) supported.
}
\value{This function returns a tree in \code{phyDat} format that has undergone one \acronym{TBR} iteration.}
\references{
The \acronym{TBR} algorithm is summarized in

Felsenstein, J. 2004. \cite{Inferring Phylogenies.} Sinauer Associates, Sunderland, Massachusetts.
}
\author{
Martin R. Smith
}

\seealso{
\code{\link{rooted.tbr}}, useful when the position of the root node should be retained.
}
\examples{
tree <- rtree(20, br=NULL)
tbr(tree)
}