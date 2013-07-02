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
load('SigSut.RData')
data <- SigSut.preparedata; tree <- SigSut.parstree; outgroup<-c('Lingula', 'Mickwitzia', 'Neocrania')
parsimony.inapp(tree, data)
\dontrun{tree.search(tree, data, outgroup, method='NNI')
pratchet.inapp(tree, data, outgroup)
sectorial.search(tree, data, outgroup)}
}
