\name{two.tip.tree}
\alias{two.tip.tree}
\title{Two-tipped phylogenetic tree}
\description{This function generates a phylogenetic 'tree' containing two tips.}
\usage{two.tip.tree(tip1, tip2)}
\arguments{
  \item{tip1}{A character string representing the label for tip 1}
  \item{tip2}{A character string representing the label for tip 2}
}
\value{A phyDat object of length 2}
\author{Martin R. Smith}
\seealso{\code{\link{single.taxon.tree}}
\code{\link{bind.tree}}}
\examples{two.tip.tree('Homo', 'Pan')}
\keyword{ tree }

