\name{add.tip}
\alias{add.tip}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add a tip to a phylogenetic tree
}
\description{
This function adds a tip to a phylogenetic tree at a specified location
}
\usage{
add.tip(tree, where, label)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{A phylogenetic tree of class phyDat}
  \item{where}{The node or tip that should form the sister taxon to the new node.  To add a new tip at the root, use "where = 0"}
  \item{label}{The character string to entitle the new tip}
}
\details{
Can be seen as a special case of ape:::bind.tree.fast.
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
