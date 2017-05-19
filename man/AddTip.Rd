\name{AddTip}
\alias{AddTip}
\title{Add a tip to a phylogenetic tree}
\description{This function adds a tip to a phylogenetic tree at a specified location.}
\usage{
AddTip(tree, where, label)
}
\arguments{
  \item{tree}{A phylogenetic tree of class \code{\link{phylo}};}
  \item{where}{The node or tip that should form the sister taxon to the new node.  To add a new tip at the root, use "where = 0";}
  \item{label}{A character string providing the label the new tip.}
}
\details{
Extends \code{\link{bind.tree}}, which cannot handle single-taxon trees.
}
\value{
This function returns a tree of class \code{phylo} with an additional tip at the desired location.
}
\author{Martin R. Smith}

\seealso{
\itemize{
\item \code{\link{bind.tree}}
\item \code{\link{nodelabels}}
}
}
\examples{
  plot(tree <- rtree(10, br=NULL)); nodelabels(); nodelabels(15, 15, bg='green'); dev.new()
  plot(AddTip(tree, 15, 'NEW_TIP'))
}
\keyword{ tree }
