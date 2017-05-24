\name{RearrangeTree}
\alias{RearrangeTree}
\title{Rearrange phylogenetic tree}
\description{Performs one tree rearrangement of the specified type.}
\usage{
RearrangeTree(tree, data, rearrange, min.score = NULL, concavity = NULL, return.single = TRUE,
 iter = "<unknown>", cluster = NULL, trace = 0)
}
\arguments{
  \item{tree}{a rooted bifurcating phylogenetic tree with the desired outgroup, and the attributes
    \code{pscore}, the tree's parsimony score, and 
    \code{hits}, the number of times the best score has been hit in the calling function;
  }
  \item{data}{a data matrix in \code{morphyDat} format, perhaps created with \code{\link{MorphyData}};}
  \item{rearrange}{a rearrangement function: probably one of 
    \code{\link{RootedNNI}}, \code{\link{RootedSPR}} or \code{\link{RootedTBR}};}
  \item{min.score}{trees longer than \code{min.score}, probably the score of the starting tree,
    will be discarded;}
  \item{concavity}{concavity constant for implied weighting (not currently implemented!); 
    see \code{\link{InapplicableParsimony}};}
  \item{return.single}{returns all trees if \kbd{FALSE} or a randomly selected tree if \kbd{TRUE};}
  \item{iter}{iteration number of calling function, for reporting to user only;}
  \item{cluster}{a cluster, prepared with \code{\link{PrepareCluster}}, to accelerate 
    searches on multicore machines;}
  \item{trace}{determines how much information to output to screen.}
}
\value{This function returns the most parsimonious of the trees generated, with attributes \code{hits} and \code{pscore}
 as described for argument \code{tree}.}
\author{Martin R. Smith}
\seealso{
  \itemize{
    \item \code{\link{RootedNNI}}
    \item \code{\link{RootedSPR}}
    \item \code{\link{RootedTBR}}
  }
}
\examples{
data('SigSut')
random.tree <- rtree(34, tip.label=names(SigSut.data), br=NULL)
RearrangeTree(random.tree, SigSut.preparedata, RootedNNI)
}