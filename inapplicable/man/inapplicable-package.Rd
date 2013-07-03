\name{inapplicable-package}
\alias{inapplicable-package}
\alias{inapplicable}
\docType{package}
\title{Inapplicable data in parsimony analysis}
\description{
Performs parsimony calculations and tree-search on datasets that contain inapplicable data.
An implementation of the algorithm proposed in 
M. R. Smith (submitted) Systematic Biology XXX, doi: 10.1093/sysbio/syXXXX
Deposited at Data Dryad, doi: 10.5061/dryad.n2k53/1
}
\details{
\tabular{ll}{
Package: \tab inapplicable\cr
Type: \tab Package\cr
Version: \tab 0.01.153\cr
Date: \tab 2013-07-02\cr
License: \tab CC-0\cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
Martin R. Smith

Maintainer: Martin R. Smith <martins@gmail.com>
}
\references{
M. R. Smith (submitted) Systematic Biology XXX, doi: 10.1093/sysbio/syXXXX
}
\keyword{ package }
\keyword{ tree }
\seealso{
\code{\link[ape:ape-package]{ape}} ~~
\code{\link[phangorn:phangorn-package]{phangorn}} ~~
}
\examples{
# Set up examples:
data('SigSut')
data <- SigSut.preparedata; outgroup <- c('Lingula', 'Mickwitzia', 'Neocrania')
njtree <- root(nj(dist.hamming(SigSut.phy)), outgroup, resolve.root=TRUE)
njtree$edge.length <- NULL; njtree<-set.outgroup(njtree, outgroup)

parsimony.inapp(njtree, data)
\dontrun{tree.search(tree, data, outgroup, method='NNI')
sectorial.inapp(njtree, data, outgroup, maxit=1, maxiter=50, largest.sector=7)
sectorial.search(njtree, data, outgroup)
pratchet.inapp(njtree, data, outgroup)}
}
