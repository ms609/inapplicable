\name{sectorial.search}
\alias{sectorial.search}
\alias{sectorial.inapp}
\title{Sectorial search}
\description{Performs a sectorial search on a tree, preserving the position of the root.}
\usage{
sectorial.search(start.tree, data, outgroup, rearrangements = "NNI")
sectorial.inapp(start.tree, data, outgroup = NULL, maxit = 100, maxiter = 500, k = 5,
  trace = 0, smallest.sector = 4, largest.sector = 1e+06, rearrangements = "NNI", ...)
}
\arguments{
  \item{start.tree}{a rooted, resolved tree in \code{\link{phylo}} format from which to start the search;}
  \item{data}{a dataset, in the format produced by \code{\link{prepare.data}} (or in \code{phyDat} format);}
  \item{outgroup}{a vector listing the taxa that form the outgroup;}
  \item{maxit}{maximum number of sectorial iterations to perform;}
  \item{maxiter}{maximum number of rearrangements to perform on each sectorial iteration;}
  \item{k}{stop when \code{k} searches have improved their sectorial score;}
  \item{trace}{integer determining how verbose the reporting to stdout will be;}
  \item{smallest.sector}{sectors with fewer than \code{smallest.sector} taxa will not be selected; \kbd{4} is the smallest sensible value;}
  \item{largest.sector}{sectors with more than \code{largest.sector} taxa will not be selected;}
  \item{rearrangements}{method to use when rearranging subtrees: NNI, SPR or TBR;}
  \item{\dots}{other arguments to pass to subsequent functions.}
}
\details{
\code{sectorial.inapp} performs a sectorial search on the tree specified. A sectorial search 
detaches a random part of the tree, performs rearrangments on this subtree, then reattaches it 
to the main tree (Goloboff, 1999).
The improvement to local \var{pscore} hopefully (but not necessarily) improves the overall \var{pscore}.
As such, the output of \code{sectorial.inapp} should be treated by further \acronym{TBR (/SPR/NNI)}
rearrangements and only retained if the ultimate parsimony score is better than 
that of the original tree.

\code{sectorial.search} is a basic recipe that runs \code{sectorial.inapp} followed by a few rounds
of tree rearrangement, returning a tree whose \var{pscore} is no worse than that of \code{start.tree}.
}
\value{
These functions return a rooted tree of class \code{phylo}.
}
\references{
Goloboff, P. (1999). \cite{Analyzing large data sets in reasonable times: solutions for composite optima.} Cladistics, 15(4), 415-428. doi:\href{http://dx.doi.org/10.1006/clad.1999.0122}{10.1006/clad.1999.0122}
}
\author{
Martin R. Smith
}

\seealso{
\itemize{ 
\item \code{\link{tree.search}}
\item \code{\link{pratchet.inapp}}
}
}

\examples{
data('SigSut')
outgroup <- c('Lingula', 'Mickwitzia', 'Neocrania')
njtree <- root(nj(dist.hamming(SigSut.phy)), outgroup, resolve.root=TRUE)
njtree$edge.length <- NULL; njtree<-set.outgroup(njtree, outgroup)
sectorial.inapp(njtree, SigSut.preparedata, outgroup, maxit=1, maxiter=50, largest.sector=7)
\dontrun{sectorial.search(njtree, SigSut.preparedata, outgroup, 'SPR')} # Will be time-consuming

## sectorial.search is currently defined as
function (start.tree, data, outgroup, rearrangements='NNI') {
  best.score <- attr(start.tree, 'pscore')
  if (length(best.score) == 0) best.score <- parsimony.inapp(start.tree, data)
  sect <- sectorial.inapp(start.tree, data, outgroup=outgroup, trace=0, maxit=30, maxiter=200, maxhits=15, smallest.sector=6, largest.sector=length(start.tree$edge[,2])*0.25, rearrangements=rearrangements)
  sect <- tree.search(sect, data, outgroup, method='NNI', maxiter=2000, maxhits=20, trace=3)
  sect <- tree.search(sect, data, outgroup, method='TBR', maxiter=2000, maxhits=25, trace=3)
  sect <- tree.search(sect, data, outgroup, method='SPR', maxiter=2000, maxhits=50, trace=3)
  sect <- tree.search(sect, data, outgroup, method='NNI', maxiter=2000, maxhits=50, trace=3)
  if (attr(sect, 'pscore') <= best.score) {
    return (sect)
  } else return (set.outgroup(start.tree, outgroup))
}
}
\keyword{ tree }