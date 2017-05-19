\name{SectorialSearch}
\alias{SectorialSearch}
\alias{InapplicableSectorial}
\title{InapplicableSectorial search}
\description{Performs a sectorial search on a tree, preserving the position of the root.}
\usage{
SectorialSearch(tree, data, outgroup, concavity = NULL, rearrangements = "NNI",
  maxiter = 2000, cluster = NULL, trace = 3)
InapplicableSectorial(tree, data, outgroup = NULL, concavity = NULL, maxit = 100, maxiter = 500, k = 5,
  trace = 0, smallest.sector = 4, largest.sector = 1e+06, rearrangements = "NNI", ...)
}
\arguments{
  \item{tree}{a rooted, resolved tree in \code{\link{phylo}} format from which to start the search;}
  \item{data}{a dataset, in the format produced by \code{\link{PrepareData}} (or in \code{phyDat} format);}
  \item{outgroup}{a vector listing the taxa that form the outgroup;}
  \item{concavity}{concavity constant for implied weighting (not currently implemented!); 
    see \code{\link{InapplicableParsimony}};}
  \item{maxit}{maximum number of sectorial iterations to perform;}
  \item{maxiter}{maximum number of rearrangements to perform on each sectorial iteration;}
  \item{cluster}{a cluster prepared using \code{\link{PrepareCluster}}; may speed up search on multicore machines;}
  \item{k}{stop when \code{k} searches have improved their sectorial score;}
  \item{trace}{integer determining how verbose the reporting to stdout will be;}
  \item{smallest.sector}{sectors with fewer than \code{smallest.sector} taxa will not be selected; \kbd{4} is the smallest sensible value;}
  \item{largest.sector}{sectors with more than \code{largest.sector} taxa will not be selected;}
  \item{rearrangements}{method to use when rearranging subtrees: NNI, SPR or TBR;}
  \item{\dots}{other arguments to pass to subsequent functions.}
}
\details{
\code{InapplicableSectorial} performs a sectorial search on the tree specified. A sectorial search 
detaches a random part of the tree, performs rearrangments on this subtree, then reattaches it 
to the main tree (Goloboff, 1999).
The improvement to local \var{pscore} hopefully (but not necessarily) improves the overall \var{pscore}.
As such, the output of \code{InapplicableSectorial} should be treated by further \acronym{TBR (/SPR/NNI)}
rearrangements and only retained if the ultimate parsimony score is better than 
that of the original tree.

\code{SectorialSearch} is a basic recipe that runs \code{InapplicableSectorial} followed by a few rounds
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
\item \code{\link{TreeSearch}}
\item \code{\link{InapplicablePratchet}}
}
}

\examples{
data('SigSut')
outgroup <- c('Lingula', 'Mickwitzia', 'Neocrania')
njtree <- root(nj(dist.hamming(SigSut.phy)), outgroup, resolve.root=TRUE)
njtree$edge.length <- NULL; njtree<-SetOutgroup(njtree, outgroup)
InapplicableSectorial(njtree, SigSut.preparedata, outgroup, maxit=1, maxiter=50, largest.sector=7)
\dontrun{SectorialSearch(njtree, SigSut.preparedata, outgroup, 'SPR')} # Will be time-consuming

## SectorialSearch is currently defined as
function (start.tree, data, outgroup, rearrangements='NNI') {
  best.score <- attr(start.tree, 'pscore')
  if (length(best.score) == 0) best.score <- InapplicableParsimony(start.tree, data)
  sect <- InapplicableSectorial(start.tree, data, outgroup=outgroup, trace=0, maxit=30, maxiter=200, maxhits=15, smallest.sector=6, largest.sector=length(start.tree$edge[,2])*0.25, rearrangements=rearrangements)
  sect <- TreeSearch(sect, data, outgroup, method='NNI', maxiter=2000, maxhits=20, trace=3)
  sect <- TreeSearch(sect, data, outgroup, method='TBR', maxiter=2000, maxhits=25, trace=3)
  sect <- TreeSearch(sect, data, outgroup, method='SPR', maxiter=2000, maxhits=50, trace=3)
  sect <- TreeSearch(sect, data, outgroup, method='NNI', maxiter=2000, maxhits=50, trace=3)
  if (attr(sect, 'pscore') <= best.score) {
    return (sect)
  } else return (SetOutgroup(start.tree, outgroup))
}
}
\keyword{ tree }