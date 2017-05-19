\name{ExtractClade}
\alias{ExtractClade}
\title{Extract a clade}
\description{Safely extracts a clade from a phylogenetic tree.}
\usage{
ExtractClade(phy, node)
}
\arguments{
  \item{phy}{A phylogenetic tree in \code{phylo} format;}
  \item{node}{The number of the node at the base of the clade to be extracted.}
}
\details{
Modified from the \pkg{ape} function \code{\link{extract.clade}}, which sometimes behaves erratically.  
The function is intended for use with the function \code{\link{TBR}}, and features present in 
\code{extract.clade} but unnecessary for this goal have been removed. Unlike extract.clade, 
this function supports the extraction of 'clades' that constitute a single tip.
}
\value{
This function returns a phylogenetic tree that represents a clade extracted from the original tree.
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
plot(ExtractClade(tree, 33))
}
