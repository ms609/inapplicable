\name{add.tip}
\alias{add.tip}
\title{Add a tip to a phylogenetic tree}
\description{This function adds a tip to a phylogenetic tree at a specified location.}
\usage{
add.tip(tree, where, label)
}
\arguments{
  \item{tree}{A phylogenetic tree of class phyDat;}
  \item{where}{The node or tip that should form the sister taxon to the new node.  To add a new tip at the root, use "where = 0";}
  \item{label}{A character string providing the label the new tip.}
}
\details{
Extends \code{ape:::\link{bind.tree}}, which cannot handle single-taxon trees.
}
\value{
A tree of class phyDat with an additional tip at the desired location.
}
\author{Martin R. Smith}

\seealso{
\code{\link{bind.tree}}

\code{\link{nodelabels}}
}
\examples{
  plot(tree <- rtree(10, br=NULL)); nodelabels(); nodelabels(15, 15, bg='green'); dev.new()
  plot(add.tip(tree, 15, 'NEW_TIP'))
}
\keyword{ tree }
