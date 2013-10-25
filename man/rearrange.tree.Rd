\name{rearrange.tree}
\alias{rearrange.tree}
\title{Rearrange phylogenetic tree}
\description{Performs one tree rearrangement of the specified type.}
\usage{
rearrange.tree(tree, data, rearrange, min.score = NULL, concavity = NULL, return.single = TRUE,
 iter = "<unknown>", cluster = NULL, trace = 0)
}
\arguments{
  \item{tree}{a rooted bifurcating phylogenetic tree with the desired outgroup, and the attributes
    \code{pscore}, the tree's parsimony score, and 
    \code{hits}, the number of times the best score has been hit in the calling function;
  }
  \item{data}{a data matrix defined by \code{\link{prepare.data}};}
  \item{rearrange}{a rearrangement function: probably one of 
    \code{\link{rooted.nni}}, \code{\link{rooted.spr}} or \code{\link{rooted.tbr}};}
  \item{min.score}{trees longer than \code{min.score}, probably the score of the starting tree,
    will be discarded;}
  \item{concavity}{concavity constant for implied weighting (not currently implemented!); 
    see \code{\link{parsimony.inapp}};}
  \item{return.single}{returns all trees if \kbd{FALSE} or a randomly selected tree if \kbd{TRUE};}
  \item{iter}{iteration number of calling function, for reporting to user only;}
  \item{cluster}{a cluster, prepared with \code{\link{prepare.cluster}}, to accelerate 
    searches on multicore machines;}
  \item{trace}{determines how much information to output to screen.}
}
\value{This function returns the most parsimonious of the trees generated, with attributes \code{hits} and \code{pscore}
 as described for argument \code{tree}.}
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