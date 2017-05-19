\name{TreeSearch}
\alias{TreeSearch}
\title{Search for most parsimonious trees}
\description{Run standard search algorithms (\acronym{NNI}, \acronym{SPR} or \acronym{TBR}) to search for a more parsimonious tree.}
\usage{
TreeSearch(tree, data, outgroup, concavity = NULL, method = "NNI", maxiter = 100, 
  maxhits = 20, forest.size = 1, cluster = NULL, trace = 1, ...)
}
\arguments{
  \item{tree}{a fully-resolved starting tree in \code{\link{phylo}} format, with the desired outgroup; edge lengths are not supported and will be deleted;}
  \item{data}{a matrix, ideally in the format produced by \code{\link{PrepareData}}, 
  alternatively as a \code{\link{phyDat}} object.  May contain inapplicable data;}
  \item{outgroup}{a vector listing the taxa in the outgroup;}
  \item{concavity}{concavity constant for implied weighting (not currently implemented!); 
    see \code{\link{InapplicableParsimony}};}
  \item{method}{rearrangements to perform; one of \kbd{NNI}, \kbd{SPR}, or \kbd{TBR};}
  \item{maxiter}{the maximum number of iterations to perform before abandoning the search;}
  \item{maxhits}{the maximum times to hit the best pscore before abandoning the search;}
  \item{forest.size}{the maximum number of trees to return - useful in concert with \code{\link{consensus}};}
  \item{cluster}{a cluster prepared using \code{\link{PrepareCluster}}; may speed up search on multicore machines;}
  \item{trace}{higher values provide more verbose user feedback in stdout;}
  \item{\dots}{other arguments to pass to subsequent functions.}
}

\value{
This function returns a tree, with an attribute \code{pscore} conveying its parsimony score.
Note that the parsimony score will be inherited from the tree's attributes, which is only valid if it 
was generated using the same \code{data} that is passed here.
}
\author{
Martin R. Smith
}
\seealso{
\itemize{
\item \code{\link{InapplicableParsimony}}, calculates parsimony score, supports inapplicable tokens;
\item \code{\link{RootedNNI}}, conducts tree rearrangements;
\item \code{\link{SectorialSearch}}, alternative heuristic, useful for larger trees;
\item \code{\link{InapplicablePratchet}}, alternative heuristic, useful to escape local optima.
}}
\examples{
data('SigSut')
outgroup <- c('Lingula', 'Mickwitzia', 'Neocrania')
njtree <- root(nj(dist.hamming(SigSut.phy)), outgroup, resolve.root=TRUE)
njtree$edge.length <- NULL; njtree<-SetOutgroup(njtree, outgroup)
\dontrun{
TreeSearch(njtree, SigSut.preparedata, outgroup, maxiter=20, method='NNI')
TreeSearch(njtree, SigSut.preparedata, outgroup, maxiter=20, method='SPR')
TreeSearch(njtree, SigSut.preparedata, outgroup, maxiter=20, method='TBR')}
}

\keyword{ tree }
