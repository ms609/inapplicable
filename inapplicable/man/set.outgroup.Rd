\name{set.outgroup}
\alias{set.outgroup}
\alias{root.robust}
\title{Root a phylogenetic tree}
\description{Sets the root of a phylogenetic tree, such that one child of the root node is \code{outgroup}.}
\usage{set.outgroup(tree, outgroup)}
\arguments{
  \item{tree}{a tree, in \code{\link{phyDat}} format, with all nodes resolved;}
  \item{outgroup}{a vector of mode numeric or character specifying the new outgroup.}
}
\value{This function returns a rooted tree with the new outgroup.}
\author{
Martin R. Smith
}

\seealso{
\code{\link{root}} ought to produce the same result as this function,
 but does not always do so in practice.
}
\examples{
  tree <- read.tree(text='(((a,b),c),(d,e));')
  plot(tree)
  plot(root.robust(tree, 'e'))
  plot(root.robust(tree, 4:5))
}