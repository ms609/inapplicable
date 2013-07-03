\name{single.taxon.tree}
\alias{single.taxon.tree}
\title{Single taxon tree}
\description{Create a phylogenetic 'tree' that comprises a single taxon.}
\usage{single.taxon.tree(label)}
\arguments{  \item{label}{a character vector specifying the label of the tip.}}
\value{This function returns a phyDat object containing a single tip with the specified label.}
\seealso{\code{\link{two.tip.tree}}}
\examples{single.taxon.tree('Homo_sapiens')}
\keyword{ tree }