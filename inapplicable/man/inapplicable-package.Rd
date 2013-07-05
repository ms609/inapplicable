\name{inapplicable-package}
\alias{inapplicable-package}
\alias{inapplicable}
\docType{package}
\title{Inapplicable data in parsimony analysis}
\description{
This package correctly calculates tree parsimony scores when datasets contain inapplicable data, and
provides functions to rearrange trees whilst preserving the outgroup. The package also contains heuristic
search methods to locate the most parsimonious tree.
}
\details{
This package calculates the parsimony score on phylogenetic trees.  It can also be used to find
the most parsimonious tree.  
}
\section{Usage}{
\itemize{
\item If necessary, install the package using \code{install.packages('inapplicable')}.
\item Load the package using \code{library('inapplicable')}.
\item Load your data: perhaps from a nexus file using \code{my.data <- \link{read.nexus.data}('filename.nex')}.  
Cells that are ambiguous (but definitely not inapplicable) should be denoted with a different token (perhaps \code{+})
from cells that are ambiguous (and may or may not be applicable), which are usually denoted \code{?}.
\item Define the tokens used in your dataset using a contrast matrix: see the 'phangorn-specials' 
vignette  in the \pkg{phangorn} package, accesible from \code{help('phangorn')}\var{>index>package vignettes}, for details.
\item Convert the data to \code{\link{phyDat}} format using \code{my.phyDat <- phyDat(my.data, type='USER', contrast=my.contrast)}
\item Optimize the data for analysis using my.prepdata <- prepare.data(my.phyDat)
\item Define the outgroup or root taxon using my.outgroup <- c('taxon1', 'taxon2')
\item Load a tree, perhaps using \code{tree <- \link{read.nexus}(treename.nex)} --
 or generate a random starting tree, perhaps by \code{tree <- rtree(length(my.prepdata), names(my.prepdata), br=NULL)}
 or \code{tree <- root(\link{nj}(dist.hamming(my.phyDat)), outgroup, resolve.root=TRUE); tree$edge.length <- NULL;}
\item Calculate the tree's parsimony score with \code{\link{parsimony.inapp}(tree, my.prepdata)}
\item Search for a better tree using \code{better.tree <- \link{tree.search}(tree, my.prepdata, my.outgroup, method='SPR')}
\item Try a sectorial search: \code{better.tree <- \link{sectorial.search}(better.tree, my.prepdata, my.outgroup)}
\item Try the parsimony ratchet: \code{better.tree <- \link{pratchet.inapp}(better.tree, my.prepdata, my.outgroup)}
\item View the results: \code{\link{plot}(better.tree)}
\item Once you have reached the most parsimonious tree, retain multiple trees to determine the consensus: 
\code{best.trees <- \link{tree.search}(better.tree, my.prepdata, my.outgroup, method='TBR', forest.size=100, maxiter=10000, maxhits=100)}
\item Calculate and display the consensus tree: \code{plot(my.consensus <- \link{consensus}(best.trees))}
}}

\author{
Martin R. Smith
}
\references{
M. R. Smith (submitted). \cite{Inapplicable data in parsimony analysis.} Systematic Biology XXX, doi: \href{http://dx.doi.org/10.1093/sysbio/syXXXX}{10.1093/sysbio/syXXXX}.

This R package is deposited at Data Dryad, doi:\href{http://dx.doi.org/10.5061/dryad.n2k53/3}{10.5061/dryad.n2k53/3}
}
\keyword{ package }
\keyword{ tree }
\seealso{
\itemize{
\item \code{\link[ape:ape-package]{ape}}
\item \code{\link[phangorn:phangorn-package]{phangorn}}
}}
\examples{
## Load data and generate some example trees:
data('SigSut')
my.contrast <- matrix(data=c(
1,0,0,
0,1,0,
0,0,1,
1,1,0,
1,1,1), ncol=3, byrow=TRUE); 
dimnames(my.contrast) <- list(
  c(0, 1, '-', '+', '?'), # rows in the contrast matrix correspond to tokens in the data matrix
  c(0, 1, '-') # columns correspond to the possible states of transformation series
 )
data <- SigSut.preparedata; outgroup <- c('Lingula', 'Mickwitzia', 'Neocrania')
njtree <- root(nj(dist.hamming(SigSut.phy)), outgroup, resolve.root=TRUE)
njtree$edge.length <- NULL; njtree<-set.outgroup(njtree, outgroup)

parsimony.inapp(njtree, data)
\dontrun{ (time consuming)
tree.search(tree, data, outgroup, method='NNI')
sectorial.inapp(njtree, data, outgroup, maxit=1, maxiter=50, largest.sector=7)
sectorial.search(njtree, data, outgroup)
pratchet.inapp(njtree, data, outgroup)}
}
