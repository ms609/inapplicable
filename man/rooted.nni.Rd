\name{rooted.nni}
\alias{rooted.nni}
\alias{rooted.spr}
\alias{rooted.tbr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tree rearrangement
}
\description{
This function performs a rearrangement iteration on a tree, retaining the position of the root.
}
\usage{
rooted.nni(tree)
rooted.spr(tree)
rooted.tbr(tree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
A phyDat object, with all nodes resolved.
}
}
\details{
A single NNI, SPR or TBR rearrangement is performed, subject to the constrain that no taxon may be moved to the opposite side of the root node.
Branch lengths are not (yet) supported.
}
\value{
A tree, in phyDat format.
}
\author{
Martin R. Smith
rooted.nni is abridged from phangorn:::nnin
}

\seealso{
\code{\link{set.outgroup}}, set the outgroup of the phylogenetic tree
\code{\link{nni}}, unrooted NNI and SPR
\code{\link{tbr}}, unrooted TBR
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
