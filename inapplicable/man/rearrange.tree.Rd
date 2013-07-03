\name{rearrange.tree}
\alias{rearrange.tree}
\title{Rearrange phylogenetic tree}
\description{Performs one tree rearrangement of the specified type.}
\usage{
rearrange.tree(tree, data, rearrange, return.single = TRUE,
 iter = "<unknown>", cores = 4, trace = 0)
}
\arguments{
  \item{tree}{a rooted bifurcating phylogenetic tree with the desired outgroup, and the attributes:
  \tabular{rl}{
    \code{pscore}\tab the tree's parsimony score\cr
    \code{hits}\tab the number of times the best score has been hit in the calling function;
  }
  }
  \item{data}{a data matrix defined by \code{\link{prepare.data}};}
  \item{rearrange}{a rearrangement function: probably one of 
    \code{\link{rooted.nni}}, \code{\link{rooted.spr}} or \code{\link{rooted.tbr}}}
  \item{return.single}{returns all trees if FALSE or a randomly selected tree if TRUE}
  \item{iter}{iteration number of calling function, for reporting to user only;}
  \item{cores}{number of separate rearrangements to evaluate in parallel;}
  \item{trace}{determines how much information to output to screen;}
}
\value{The most parsimonious of the trees generated, with attributes \code{hits} and \code{pscore}.}
\author{Martin R. Smith}
\seealso{
  \itemize{
    \item \code{\link{rooted.nni}}
    \item \code{\link{rooted.spr}}
  \item \code{\link{rooted.tbr}}
  }
}
\examples{
data('SigSut')
random.tree <- rtree(34, tip.label=names(SigSut.data), br=NULL)
rearrange.tree(random.tree, SigSut.preparedata, rooted.nni)
}