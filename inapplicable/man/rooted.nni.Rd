\name{rooted.nni}
\alias{rooted.nni}
\alias{rooted.spr}
\alias{rooted.tbr}
\title{Rearrange a rooted tree}
\description{This function performs a rearrangement iteration on a tree, retaining the position of the root.}
\usage{
rooted.nni(tree)
rooted.spr(tree)
rooted.tbr(tree)
}
\arguments{
  \item{tree}{An object of class \code{\link{phylo}}, with all nodes resolved (bifurcating).}
}
\details{
A single \acronym{NNI}, \acronym{SPR} or \acronym{TBR} rearrangement is performed, subject to the constraint that 
no taxon may be moved to the opposite side of the root node.
Branch lengths are not (yet) supported.
}
\value{
This function returns a tree, in \code{phylo} format.
}
\author{
Martin R. Smith

\code{rooted.nni} is abridged from the \pkg{phangorn} function \code{nnin}
}

\seealso{
\itemize{
\item \code{\link{set.outgroup}}, set the outgroup of the phylogenetic tree
\item \code{\link{nni}}, unrooted \acronym{NNI} and \acronym{SPR}
\item \code{\link{tbr}}, unrooted \acronym{TBR}
}
}
\examples{
  tree <- read.tree(text='(((a,b),c),(d,(e,f)));')
  tree <- set.outgroup(tree, c('e', 'f'))
  plot(tree)
  dev.new()
  plot(rooted.nni(tree))
  plot(rooted.spr(tree))
  plot(rooted.tbr(tree))
}
