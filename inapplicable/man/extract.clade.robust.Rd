\name{extract.clade.robust}
\alias{extract.clade.robust}

\title{
Extract clade
}
\description{
Extracts a clade from a phylogenetic tree
}
\usage{
extract.clade.robust(phy, node)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{
A phylogenetic tree in phyDat format
}
  \item{node}{
The number of the node at the base of the clade to be extracted
}
}
\details{
Modified from ape:::extract.clade.  Added support for extracting tips.  Removed features unneccesary for TBR.
}
\value{
Returns a phylogenetic tree that represents a clade extracted from \code{phy}
}
\author{
Martin R. Smith
}
\seealso{
\code{\link{extract.clade}}
}
\examples{
tree <- rtree(20, br=NULL)
plot(tree); nodelabels(); nodelabels(33, 33, bg='yellow'); dev.new()
plot(extract.clade.robust(tree, 33))
}
