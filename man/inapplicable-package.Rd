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


## To use this script on your own data, launch R, and type (or copy-paste) the following text
## into the R console.  Lines starting with '#' are comments and do not need to be copied.


## To install the package for the first time, type
install.packages('inapplicable')

## Once the package has been installed, load it using
library(inapplicable)

## Data can be read from a nexus file (note the restrictions detailed above):
my.data <- read.nexus.data('C:/path/to/filename.nex')

## Alternatively you can use a built-in dataset:
data(SigSut); my.data <- SigSut.data

## A contrast matrix translates the tokens used in your dataset to the character states to 
##    which they correspond: for example decoding 'A' to {01}.
##    For more details, see the 'phangorn-specials' vignette in the phangorn package, accesible 
##    by typing '?phangorn' in the R prompt and navigating to index > package vignettes.

contrast.matrix <- matrix(data=c(
# 0 1 -  # Each column corresponds to a character-state
  1,0,0, # Each row corresponds to a token, here 0, denoting the character-state set {0} 
  0,1,0, # 1 | {1}
  0,0,1, # - | {-}
  1,1,0, # A | {01}
  1,1,0, # + | {01}
  1,1,1  # ? | {01-}
), ncol=3, byrow=TRUE); # ncol should correspond to the number of columns in the matrix
dimnames(contrast.matrix) <- list(
  # A list of the tokens corresponding to each rows in the contrast matrix
  c(0, 1, '-', 'A', '+', '?'),
  # A list of the character-states corresponding to the columns in the contrast matrix
  c(0, 1, '-') 
)

## To see the annotated contrast matrix, type
contrast.matrix

## Apply the contrast matrix, using the phyDat format...
my.phyDat <- phyDat(my.data, type='USER', contrast=contrast.matrix)

## ... and optimize the data for analysis 
my.prepdata <- PrepareData(my.phyDat)

## Specify the names of the outgroup taxa
my.outgroup <- c('taxon1', 'taxon2')

## Load a bifurcating tree,
tree <- read.nexus(treename.nex)
## or generate a random starting tree, 
tree <- rtree(length(my.phyDat), rooted=TRUE, tip.label=names(my.phyDat), br=NULL)
## or use neighbour joining to generate a starting tree
tree <- root(nj(dist.hamming(my.phyDat)), my.outgroup, resolve.root=TRUE)
tree$edge.length <- NULL;

## View the starting tree by typing
plot(tree)

## Calculate the tree's parsimony score
InapplicableParsimony(tree, my.prepdata)

## Search for a better tree
better.tree <- TreeSearch(tree, my.prepdata, my.outgroup, method='SPR')
## Try a sectorial search (Goloboff, 1999)
better.tree <- SectorialSearch(better.tree, my.prepdata, my.outgroup)
## Try the parsimony ratchet (Nixon, 1999)
better.tree <- InapplicablePratchet(better.tree, my.prepdata, my.outgroup, maxhits=50, k=20)
## The default parameters may not be enough to find the most parsimonious tree; type 
##    ?InapplicablePratchet or ?SectorialSearch to view all search parameters.


## View the results
plot(better.tree)

## Once you have reached the most parsimonious tree, 
## retain multiple trees to determine the consensus: 
best.trees <- TreeSearch(better.tree, my.prepdata, my.outgroup,
      method='SPR', forest.size=100, maxiter=10000, maxhits=100)

## Calculate and display the consensus tree
plot(my.consensus <- consensus(best.trees))
}
}
