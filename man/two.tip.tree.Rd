\name{TwoTipTree}
\alias{TwoTipTree}
\title{Two-tipped tree}
\description{This function generates a tree of class \code{\link{phylo}} containing two tips.}
\usage{TwoTipTree(tip1, tip2)}
\arguments{
  \item{tip1}{A character string representing the label for tip 1}
  \item{tip2}{A character string representing the label for tip 2}
}
\value{This function returns a \code{phylo} object with a single root node and two tips with
the specified labels.}
\author{Martin R. Smith}
\seealso{\code{\link{SingleTaxonTree}}
\code{\link{bind.tree}}}
\examples{TwoTipTree('Homo', 'Pan')}
\keyword{ tree }

