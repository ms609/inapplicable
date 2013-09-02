\name{parsimony.inapp}
\alias{parsimony.inapp}
\alias{fitch.inapp}
\title{Calculate parsimony score}
\description{Calculates the parsimony score for a rooted tree, correctly handling the inapplicable token.}
\usage{
  parsimony.inapp(tree, data)
  parsimony.inapp(tree, data, k=3)
  fitch.inapp(tree, data)
}
\arguments{
  \item{tree}{ A rooted tree in \code{phylo} format;}
  \item{data}{ Data to use to score the tree, preferably formatted using \code{\link{prepare.data}}.}
  \item{k}{ Concavity constant for implied weighting }
}
\value{
\code{parsimony.inapp} returns the parsimony score of the tree, calculated from \code{data} using Algorithm III in Smith (submitted).
The parsimony score is calculated under equal weights unless a value is provided for \code{k}.  
If \code{k} is specified, the returned score will be calculated using the method of Goloboff (1993) 
with a concavity constant of \code{k}.
 
\code{fitch.inapp} returns a list comprising two elements: 
  \item{ 1 }{Number of evolutionary events (i.e. unweighted parsimony score) 
             implied by the tree topology and the data}
  \item{ 2 }{A vector listing the number of evolutionary events implied
             for each transformation series in \code{data}}
}
\references{
Goloboff, P. A. (1993). \cite{Estimating Character Weights During Tree Search.} Cladistics 9: 83-91. doi: 10.1111/j.1096-0031.1993.tb00209.x
M. R. Smith (submitted). \cite{Incorporating inapplicable data in parsimony analyses.} Systematic Biology XXX, doi: 10.1093/sysbio/syXXXX
}
\author{Martin R. Smith}
\examples{
data('SigSut')
random.tree <- rtree(34, tip.label=names(SigSut.data), br=NULL)
parsimony.inapp(random.tree, SigSut.preparedata)
}
\keyword{ tree }
