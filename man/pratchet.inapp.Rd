\name{pratchet.inapp}
\alias{pratchet.inapp}
\title{Parsimony ratchet}
\description{This function uses the parsimony ratchet (Nixon 1999) to search for a more parsimonious tree.}
\usage{
pratchet.inapp(tree, data, outgroup = NULL, concavity = NULL, all.trees = FALSE, maxit = 100, 
  maxiter = 5000, maxhits = 40, k = 10, trace = 0, rearrangements = "NNI", ...)
}
\arguments{
  \item{tree}{An object of class \code{phyDat} denoting the topology to begin the search from;}
  \item{data}{the output of \code{\link{prepare.data}}, or a matrix in phyDat format
    (which the function will itself pass through \code{prepare.data});}
  \item{outgroup}{a vector specifying all tips in the outgroup, used if \code{start.tree} is not rooted appropriately;}
  \item{concavity}{concavity constant for implied weighting (not currently implemented!); 
    see \code{\link{parsimony.inapp}};}
  \item{all.trees}{Set to TRUE to report all MPTs encountered during the search, perhaps to analyze consensus}
  \item{maxit}{maximum ratchet iterations to perform;}
  \item{maxiter}{maximum rearrangements to perform on each bootstrap or ratchet iteration;}
  \item{maxhits}{maximum times to hit best score before terminating a tree search within a pratchet iteration;}
  \item{k}{stop when k ratchet iterations have found the same best score;}
  \item{trace}{larger numbers provides more verbose feedback to the user;}
  \item{rearrangements}{method to use when rearranging trees: NNI, SPR or TBR;}
  \item{\dots}{other arguments to pass to subsequent functions.}
}
\value{This function returns a tree modified by parsimony ratchet iteration, retaining the position of the root.}
\references{Nixon, K. C. (1999). \cite{The Parsimony Ratchet, a new method for rapid parsimony analysis.} Cladistics, 15(4), 407-414. doi:\href{http://dx.doi.org/10.1111/j.1096-0031.1999.tb00277.x}{10.1111/j.1096-0031.1999.tb00277.x}}
\author{
Martin R. Smith

Adapted from \code{\link{pratchet}} in the \pkg{phangorn} package, which does not preserve the position of the root.
}
\seealso{
\itemize{
\item \code{\link{pratchet}}
\item \code{\link{tree.search}}
\item \code{\link{sectorial.search}}
}
}
\examples{
data('SigSut')
outgroup <- c('Lingula', 'Mickwitzia', 'Neocrania')
njtree <- root(nj(dist.hamming(SigSut.phy)), outgroup, resolve.root=TRUE)
njtree$edge.length <- NULL; njtree<-set.outgroup(njtree, outgroup)
pratchet.inapp(njtree, SigSut.preparedata, outgroup, maxit=1, maxiter=50)
}
\keyword{ tree }