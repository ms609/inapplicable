\name{tbr}
\alias{tbr}

\title{Tree bisection and reconnection}
\description{This function performs a single random \acronym{TBR} iteration.}
\usage{
tbr(tree, edge.to.break = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
A fully resolved tree in phyDat format
}
  \item{edge.to.break}{Edge to bisect, generated randomly if not specified}
}
\details{
Branch lengths are not (yet) supported.
}
\value{a phyDat object}
\references{
The TBR algorithm is summarized in
Felsenstein, J. 2004. \cite{Inferring Phylogenies.} Sinauer Associates, Sunderland, Massachusetts.
}
\author{
Martin R. Smith
}

\seealso{
\code{\link{rooted.tbr}}, when the position of the root node should be retained
}
\examples{
tree <- rtree(20, br=NULL)
tbr(tree)
}