\name{visualise.character}
\alias{visualise.character}
\alias{visualise.char}
\alias{visualize.character}
\alias{visualize.char}

\title{Visualize optimization for a character}
\description{Determine and depict the possible states for a character on a tree under the most parsimonious conditions}
\usage{
visualise.character(tree, data, char.no, plot.fun = \link{plot})
}
\arguments{
  \item{tree}{a fully-resolved tree in \code{\link{phylo}} format, with the desired outgroup; edge lengths are not supported and will be deleted;}
  \item{data}{a matrix, ideally in the format produced by \code{\link{prepare.data}}, 
  alternatively as a \code{\link{phyDat}} object.  May contain inapplicable data;}
  \item{char.no}{number of the character to be displayed;}
  \item{plot.fun}{a function that plots a tree.}
}

\value{
The function plots a cladogram, annotating each node with the values that could occur there at the lowest parsimony score.
  Nodes highlighed in red potentially (but do not necessarily) contribute to tree length.
}
\author{
Martin R. Smith
}
\examples{
data('SigSut')
outgroup <- c('Lingula', 'Mickwitzia', 'Neocrania')
njtree <- root(nj(dist.hamming(SigSut.phy)), outgroup, resolve.root=TRUE)
njtree$edge.length <- NULL; njtree<-set.outgroup(njtree, outgroup)
visualise.character(njtree, SigSut.phy, 5, function(x) {plot(x); edgelabels();})
}
\keyword{ tree }
