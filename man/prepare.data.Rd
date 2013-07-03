\name{prepare.data}
\alias{prepare.data}
\title{Prepare phyDat object}
\description{Prepares a phyDat object for analysis with modified Fitch algorithm.}
\usage{prepare.data(data)}
\arguments{\item{data}{A phyDat object}}
\details{Optimizes the phyDat object provided for faster processing.}
\value{
   A list of class '*phyDat', with each list item representing a tip and housing a matrix.
     Each column of the matrix corresponds to a transformation series
     Each row of the matrix corresponds to a token, as specified by colnames(attr(data, 'contrast'))
     Each cell contains a logical denoting whether the transformation series could take the value indicated by the row
}
\author{Martin R. Smith}
\seealso{\code{\link{parsimony.inapp}}}
\examples{
data('SigSut')
prepare.data(SigSut.phy)
}