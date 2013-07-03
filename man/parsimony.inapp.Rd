\name{parsimony.inapp}
\alias{parsimony.inapp}
\alias{fitch.inapp}
\title{Calculate Parsimony Score}
\description{Determines the true parsimony score for a rooted tree.}
\usage{
  parsimony.inapp(tree, data)
  fitch.inapp(tree, data)
}

\arguments{
  \item{tree}{ A rooted tree in \code{phylo} format }
  \item{data}{ Data to use to score the tree, preferably formatted using \code{\link{prepare.data}} }
}
\value{
parsimony.inapp returns the parsimony score of the tree, given the data.
fitch.inapp returns a list comprising two elements: 
  \item{ 1 }{Number of evolutionary events (i.e. parsimony score) implied by the tree topology and the data}
  \item{ 2 }{A vector listing the number of evolutionary events implied for each transformation series in \code{'data'}}
}
\references{
M. R. Smith (submitted) Systematic Biology XXX, doi: 10.1093/sysbio/syXXXX
}
\author{Martin R. Smith}
\examples{
data('SigSut')
random.tree <- rtree(34, tip.label=names(SigSut.data), br=NULL)
parsimony.inapp(random.tree, SigSut.preparedata)
}
\keyword{ tree }
