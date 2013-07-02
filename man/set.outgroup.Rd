\name{set.outgroup}
\alias{set.outgroup}
\alias{root.robust}

\title{
Root a phylogenetic tree
}
\description{
Sets the root of a phylogenetic tree, such that one child of the root node is \code{outgroup}.
}
\usage{
set.outgroup(tree, outgroup)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{A tree, in phyDat format, with all nodes resolved.}
  \item{outgroup}{A vector of mode numeric or character specifying the new outgroup.}
}
\details{
}
\value{A rooted tree with the new outgroup}
\author{
Martin R. Smith
}

\seealso{
\code{\link{ape,root}} \textit{should} produce the same result as this function, but does not always do so in practice.
}
\examples{
tree <- read.tree(text='(((a,b),c),(d,e));')
plot(tree)
plot(root.robust(tree, 'e'))
plot(root.robust(tree, 4:5))