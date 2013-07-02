\name{pratchet.inapp}
\alias{pratchet.inapp}
\title{
Parsimony ratchet
}
\description{
This function uses the Parsimony Ratchet to search for a more parsimonious tree
}
\usage{
pratchet.inapp(start.tree, data, outgroup = NULL, maxit = 5000, maxiter = 500, maxhits = 20, k = 10, trace = 0, rearrangements = "NNI", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{start.tree}{A phyDat object denoting the topology to begin the search from}
  \item{data}{The output of \code{\link{prepare.data}} (or a matrix in phyDat format)}
  \item{outgroup}{a vector specifying all tips in the outgroup, if \code{start.tree} is not rooted appropriately}
  \item{maxit}{maximum pratchet iterations to perform}
  \item{maxiter}{maximum pratchet iterations to perform}
  \item{maxhits}{ maximum times to hit best score before terminating each pratchet iteration}
  \item{k}{stop when k ratchet iterations have found the same best score}
  \item{trace}{Larger number provide more verbose reporting in stdout}
  \item{rearrangements}{method to use when rearranging trees: NNI, SPR or TBR}
}
\details{}
\value{Returns a tree modified by parsimony ratchet iteration, retaining the position of the root}
\references{Nixon, K. C. (1999). The Parsimony Ratchet, a New Method for Rapid Parsimony Analysis. Cladistics, 15(4), 407–414. doi:10.1111/j.1096-0031.1999.tb00277.x}
\author{
Martin R. Smith
Adapted from phangorn:::pratchet
}
\seealso{
\code{\link{pratchet}}
\code{\link{tree.search}}
\code{\link{sectorial.search}}
}
\examples{
load(SigSut)
outgroup <- c('Lingula', 'Mickwitzia', 'Neocrania')
njtree <- root(nj(dist.hamming(SigSut.phy)), outgroup, resolve.root=TRUE)
njtree$edge.length <- NULL; njtree<-set.outgroup(njtree, outgroup)
\dontrun{pratchet.inapp(njtree, SigSut.preparedata, outgroup, maxit=1, maxiter=50)}
}
\keyword{ tree }