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
\tabular{ll}{
Package: \tab inapplicable\cr
Type: \tab Package\cr
Version: \tab 0.01.153\cr
Date: \tab 2013-07-02\cr
License: \tab CC-0\cr
}
}
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
