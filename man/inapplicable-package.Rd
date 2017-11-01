\name{inapplicable-package}
\alias{inapplicable-package}
\alias{inapplicable}
\docType{package}
\title{Inapplicable data in parsimony analysis}
\description{
This package correctly calculates tree parsimony scores when datasets contain inapplicable data, and
provides functions to rearrange trees whilst preserving the outgroup. The package also contains heuristic
search methods to locate the most parsimonious tree.
}
\details{
This package calculates the parsimony score on phylogenetic trees.  It can also be used to find
the most parsimonious tree.  The \emph{Examples} section below provides a step-by-step guide to 
using this package on your own data.

Character data are read from a \emph{restricted} \acronym{NEXUS} format using the R function
\code{\link{read.nexus.data}}; see the latter function's \link[=read.nexus.data]{documentation} for details.

\acronym{NEXUS} files can be edited in any standard text editor,
 for example \href{http://notepad-plus-plus.org/}{Notepad++}.
 
A notable annoyance is that the parser cannot interpret curly braces, for example {01}.  
Ambiguous tokens of this nature should be replaced with a separate character, for example 
'A' to denote {01}, 'B' to denote {12}, perhaps using a search-and-replace operation in your 
favourite text editor.

Be aware that the algorithm distinguishes between character-states that are inapplicable; 
character-states that are ambiguous (but definitely not inapplicable);
and character-states that are ambiguous (and may or may not be applicable).  
These three cases should be denoted with different tokens: perhaps '-', '+' and '?'.
}

\author{
Martin R. Smith
}
\references{
BRAZEAU, M. D., GUILLERME, T. and SMITH, M. R. 2017. [Morphological phylogenetic analysis with inapplicable data](https://www.biorxiv.org/content/early/2017/10/26/209775). BioRxiv. doi:10.1101/209775
}
\keyword{ package }
\keyword{ tree }
\seealso{
\itemize{
\item \code{\link[ape:ape-package]{ape}}
\item \code{\link[phangorn:phangorn-package]{phangorn}}
}}
\examples{
\dontrun{

## Walkthrough of package functions

## NOT UP TO DATE. WILL BE UPDATED IN EARLY NOVEMBER 2017...

## To use this script on your own data, launch R, and type (or copy-paste) the following text
## into the R console.  Lines starting with '#' are comments and do not need to be copied.


## To install the package for the first time, type
install.packages('inapplicable')

## Once the package has been installed, load it using
library(inapplicable)

## Data can be read from a nexus file (note the restrictions detailed above):
my.data <- read.nexus.data('C:/path/to/filename.nex')

## Alternatively you can use a built-in dataset:
data(inapplicable.datasets); my.data <- inapplicable.datasets[[1]]

## Data should be converted to the phyDat format, with the inapplicable token listed as a level:
my.phyDat <- phyDat(my.data, type='USER', levels=c(0:9, '-'))

## Load a bifurcating tree,
tree <- read.nexus(treename.nex)
## or generate a random starting tree, 
tree <- rtree(length(my.phyDat), rooted=TRUE, tip.label=names(my.phyDat), br=NULL)
## or use neighbour joining to generate a starting tree
tree <- ape::nj(phangorn::dist.hamming(my.phyDat))
tree$edge.length <- NULL;

## Root the tree on an appropriate outgroup, if desired
my.outgroup <- names(my.phyDat[1])
tree <- ape::root(tree, my.outgroup, resolve.root=TRUE)

## View the starting tree by typing
par(mar=rep(0.5, 4)) # Shrinks the margins of the plot so there's more space for the tree
plot(tree)

## Calculate the tree's parsimony score
InapplicableFitch(tree, my.phyDat)

## Search for a better tree
better.tree <- TreeSearch(tree, my.phyDat, verbosity=5)
## Use Nearest Neighbour interchange to find a local optimum
better.tree <- TreeSearch(better.tree, my.phyDat, Rearrange = TreeSearch::RootedNNI,
                          maxIter=40000, maxHits=75, verbosity=2)
## Try the parsimony ratchet (Nixon, 1999) for a more exhaustive search of treespace
better.tree <- Ratchet(better.tree, my.phyDat, maxIt=25, maxIter=250, maxHits=20, k=5)
## The default parameters may not be enough to find the most parsimonious tree; type 
##    ?Ratchet or ?SectorialSearch to view all search parameters.

## View the results
plot(better.tree)

## Once you have reached the most parsimonious tree, 
## retain multiple trees to determine the consensus: 
best.trees <- RatchetConsensus(better.tree, my.phyDat, nSearch=10, maxIt=25, maxIter=250, maxHits=20, k=5)

## Calculate and display the consensus tree
plot(my.consensus <- consensus(best.trees))
}
}
