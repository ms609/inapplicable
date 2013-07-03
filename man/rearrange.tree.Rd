\name{rearrange.tree}
\alias{rearrange.tree}
\title{Rearrange phylogenetic tree}
\description{Performs one rearrangement of the specified method on the tree provided.}
\usage{
rearrange.tree(tree, data, rearrange, return.single = TRUE, iter = "<unknown>", cores = 4, trace = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{A rooted bifurcating phylogenetic tree with the desired outgroup, and the attributes:
    pscore: tree's parsimony score
    hits: The number of times the best score has been hit in the calling function
  
  }
  \item{data}{A data matrix defined by \code{prepare.data}}
  \item{rearrange}{A rearrangement function: probably one of rooted.nni, rooted.spr or rooted.tbr}
  \item{return.single}{Returns all trees if FALSE or a randomly selected tree if TRUE}
  \item{iter}{Iteration number of calling function, for reporting to user only}
  \item{cores}{Number of separate rearrangements to evaluate in parallel}
  \item{trace}{Determines how much information to output to screen}
}
\value{The most parsimonious of the trees generated, with attributes 'hits' and 'pscore'}
\author{Martin R. Smith}
\seealso{
  \code{\link{rooted.nni}}
  \code{\link{rooted.spr}}
  \code{\link{rooted.tbr}}
}
\examples{
data('SigSut')
random.tree <- rtree(34, tip.label=names(SigSut.data), br=NULL)
rearrange.tree(random.tree, SigSut.preparedata, rooted.nni)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
