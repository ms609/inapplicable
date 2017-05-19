\name{SingleTaxonTree}
\alias{SingleTaxonTree}
\title{Single taxon tree}
\description{Create a phylogenetic 'tree' that comprises a single taxon.}
\usage{SingleTaxonTree(label)}
\arguments{  \item{label}{a character vector specifying the label of the tip.}}
\value{This function returns a \code{phylo} object containing a single tip with the specified label.}
\seealso{\code{\link{TwoTipTree}}}
\examples{SingleTaxonTree('Homo_sapiens')}
\keyword{ tree }