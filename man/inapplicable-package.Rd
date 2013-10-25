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
the most parsimonious tree.  For usage guidance, see the \emph{Examples} section below.
}

\author{
Martin R. Smith
}
\references{
Smith, M. R., & Brazeau, M. D. (201X). \cite{Inapplicable data in parsimony analysis.} Systematic Biology XXX, doi: \href{http://dx.doi.org/10.1093/sysbio/syXXXX}{10.1093/sysbio/syXXXX}.
}
\keyword{ package }
\keyword{ tree }
\seealso{
\itemize{
\item \code{\link[ape:ape-package]{ape}}
\item \code{\link[phangorn:phangorn-package]{phangorn}}
}}
\examples{
\dontrun{ Walkthrough of package functions


## If necessary, install the package 
install.packages('inapplicable')
## Load the package 
library('inapplicable')

## Load your data: perhaps using the example datasets:
data(SigSut); my.data <- SigSut.data
## perhaps from a nexus file:
my.data <- read.nexus.data('path/to/filename.nex')
## Cells that are ambiguous (but definitely not inapplicable) should be denoted
##    with a different token (perhaps +) from cells that are ambiguous 
##    (and may or may not be applicable), which are usually denoted ?.

## Define the tokens used in your dataset using a contrast matrix: 
##    see the 'phangorn-specials' vignette  in the phangorn package,
##    accesible from help('phangorn')>index>package vignettes, for details.
my.contrast <- matrix(data=c(
1,0,0,
0,1,0,
0,0,1,
1,1,0,
1,1,1), ncol=3, byrow=TRUE); 
dimnames(my.contrast) <- list(
  c(0, 1, '-', '+', '?'), # rows in the contrast matrix correspond to tokens in the data matrix
  c(0, 1, '-') # columns correspond to the possible states of transformation series
 )

## Convert the data to phyDat format 
my.phyDat <- phyDat(my.data, type='USER', contrast=my.contrast)

## Prepare the data for analysis 
my.prepdata <- prepare.data(my.phyDat)

## Define the outgroup or root taxon 
my.outgroup <- c('taxon1', 'taxon2')

## Load a tree,
tree <- read.nexus(treename.nex)
## or generate a random starting tree, 
tree <- rtree(length(my.prepdata), names(my.prepdata), br=NULL)
## or use neighbour joining to generate a starting tree
tree <- root(nj(dist.hamming(my.phyDat)), outgroup, resolve.root=TRUE); tree$edge.length <- NULL;

## Calculate the tree's parsimony score
parsimony.inapp(tree, my.prepdata)

## Search for a better tree
better.tree <- tree.search(tree, my.prepdata, my.outgroup, method='SPR')
## Try a sectorial search (Goloboff, 1999)
better.tree <- sectorial.search(better.tree, my.prepdata, my.outgroup)
## Try the parsimony ratchet (Nixon, 1999)
better.tree <- pratchet.inapp(better.tree, my.prepdata, my.outgroup)

## View the results
plot(better.tree)

## Once you have reached the most parsimonious tree, retain multiple trees to determine the consensus: 
best.trees <- tree.search(better.tree, my.prepdata, my.outgroup,
      method='TBR', forest.size=100, maxiter=10000, maxhits=100)

## Calculate and display the consensus tree
plot(my.consensus <- consensus(best.trees))
}
}
