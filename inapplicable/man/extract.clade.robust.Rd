\name{extract.clade.robust}
\alias{extract.clade.robust}
\title{Extract a clade}
\description{Safely extracts a clade from a phylogenetic tree.}
\usage{
extract.clade.robust(phy, node)
}
\arguments{
  \item{phy}{A phylogenetic tree in phyDat format;}
  \item{node}{The number of the node at the base of the clade to be extracted.}
}
\details{
Modified from \code{ape:::\link{extract.clade}}, which sometimes behaves erratically.  
The function is intended for use with the function \code{\link{tbr}}, and features present in 
\code{extract.clade} but unnecessary for this goal have been removed. Unlike extract.clade, 
this function supports the extraction of 'clades' that constitute a single tip.
}
\value{
Returns a phylogenetic tree that represents a clade extracted from the original tree.
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
