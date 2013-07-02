\name{parsimony.inapp}
\alias{parsimony.inapp}
\alias{fitch.inapp}
\title{Parsimony score}
\description{Determines the true parsimony score for a rooted tree}
\usage{
  parsimony.inapp(tree, data)
  fitch.inapp(tree, data)
}

\arguments{
  \item{tree}{ A rooted tree in phylo format }
  \item{data}{ The output of prepare.data() }
}
\value{
parsimony.inapp returns the parsimony score of the tree, given the data.
fitch.inapp returns
  \item{ comp1 }{Number of evolutionary events (i.e. parsimony score) implied by the tree topology and the data}
  \item{ comp2 }{A vector listing the number of evolutionary events implied for each transformation series in \code{'data'}}
}
\references{
M. R. Smith (submitted) Systematic Biology XXX, doi: 10.1093/sysbio/syXXXX
}
\author{Martin R. Smith}
\examples{
data('SigSut')
parsimony.inapp(SigSut.parstree, SigSut.preparedata)
parsimony.inapp(SigSut.besttree, SigSut.preparedata)
}
\keyword{ tree }
