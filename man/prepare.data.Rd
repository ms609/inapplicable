\name{prepare.data}
\alias{prepare.data}
\title{Prepare phyDat object}
\description{Prepares a \code{phyDat} object for analysis with \code{\link{fitch.inapp}}.}
\usage{prepare.data(data)}
\arguments{\item{data}{An object of class \code{\link{phyDat}}}}
\details{Optimizes the \code{data} for efficient processing.}
\value{
   This function returns a list of class \code{*phyDat}, with each list item representing a tip and housing a matrix.
   \itemize{
\item     Each column of the matrix corresponds to a transformation series;
\item     Each row of the matrix corresponds to a token, as specified by colnames(attr(data, 'contrast'));
\item     Each cell contains a logical denoting whether the transformation series could take the value indicated by the row.
}
}
\author{Martin R. Smith}
\seealso{\code{\link{parsimony.inapp}}}
\examples{
data('SigSut')
prepare.data(SigSut.phy)
}
