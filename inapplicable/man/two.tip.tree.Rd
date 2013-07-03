\name{two.tip.tree}
\alias{two.tip.tree}
\title{Two-tipped tree}
\description{This function generates a tree of class \code{\link{phylo}} containing two tips.}
\usage{two.tip.tree(tip1, tip2)}
\arguments{
  \item{tip1}{A character string representing the label for tip 1}
  \item{tip2}{A character string representing the label for tip 2}
}
\value{This function returns a \code{phylo} object with a single root node and two tips with
the specified labels.}
\author{Martin R. Smith}
\seealso{\code{\link{single.taxon.tree}}
\code{\link{bind.tree}}}
\examples{two.tip.tree('Homo', 'Pan')}
\keyword{ tree }

