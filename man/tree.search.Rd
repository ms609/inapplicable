\name{tree.search}
\alias{tree.search}
\title{Search for Most Parsimonious Trees}
\description{Run standard search algorithms (NNI, SPR or TBR) to search for a more parsimonious tree.}
\usage{
tree.search(start.tree, data, outgroup, method = "NNI", maxiter = 100, maxhits = 20, forest.size = 1, cores = 4, trace = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{start.tree}{A fully-resolved tree in phyDat format, with the desired outgroup. Edge lengths are not supported and will be deleted.}
  \item{data}{A matrix, ideally in the format produced by \code{\link{prepare.data}}, 
  alternatively as a phyDat object.  May contain inapplicable data.}
  \item{outgroup}{A vector listing the taxa in the outgroup}
  \item{method}{Rearrangements to perform; one of 'NNI', 'SPR', or 'TBR'}
  \item{maxiter}{The maximum number of iterations to perform before abandoning the search}
  \item{maxhits}{The maximum times to hit the best pscore before abandoning the search}
  \item{forest.size}{The maximum number of trees to return - useful in concert with \code{\link{consensus}}}
  \item{cores}{Number of trees to examine in parallel on each iteration; may speed up search on multicore machines}
  \item{trace}{Higher values provide more verbose user feedback in stdout}
  \item{\dots}{Other arguments to pass to subsequent functions}
}

\value{
A tree, with an attribute \code{pscore} conveying its parsimony score.
Note that the pscore will be inherited from the tree's attributes, which is only valid if it 
was generated using the same \code{data} that is passed here.
}
\author{
Martin R. Smith
}
\seealso{
\code{\link{parsimony.inapp}}, calculates parsimony score, supports inapplicable tokens.
\code{\link{rooted.nni}}, conducts tree rearrangements.
\code{\link{sectorial.search}}, alternative heuristic, useful for larger trees.
\code{\link{pratchet.inapp}}, alternative heuristic, useful to escape local optima.
}
\examples{
data('SigSut')
outgroup <- c('Lingula', 'Mickwitzia', 'Neocrania')
njtree <- root(nj(dist.hamming(SigSut.phy)), outgroup, resolve.root=TRUE)
njtree$edge.length <- NULL; njtree<-set.outgroup(njtree, outgroup)
\dontrun{
tree.search(njtree, SigSut.preparedata, outgroup, maxiter=20, method='NNI')
tree.search(njtree, SigSut.preparedata, outgroup, maxiter=20, method='SPR')
tree.search(njtree, SigSut.preparedata, outgroup, maxiter=20, method='TBR')}
}

\keyword{ tree }
