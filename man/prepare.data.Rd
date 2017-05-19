\name{PrepareData}
\alias{PrepareData}
\title{Prepare phyDat object}
\description{Prepares a \code{phyDat} object for analysis with \code{\link{InapplicableFitch}}.}
\usage{PrepareData(data)}
\arguments{\item{data}{An object of class \code{\link{phyDat}}}}
\details{Optimizes the \code{data} for efficient processing.}
\value{
   This function returns a matrix of class \code{*phyDat}.
   \itemize{
\item     Each column of the matrix corresponds to an operational taxonomic unit;
\item     Each row of the matrix represents a transformation series;
\item     Each cell in the matrix contains an integer that denotes the possible tokens.
          This integer is the decimal representation of the binary figure represented by a row
          in the contrast matrix (in reverse order), after
          phangorn:::prepareFitchData.          
  }
  In addition to the standard \link{phyDat} attributes, the returned object bears an attribute
  \code{inapp.level}, specifying the level that corresponds to the inapplicable token, and 
  \code{min.steps}, specifying the minimum number of steps necessary to obtain the observed
  distribution of tokens (necessary when calculating implied weights).
}
\author{Martin R. Smith}
\seealso{\code{\link{InapplicableParsimony}}}
\examples{
data('SigSut')
PrepareData(SigSut.phy)
}
