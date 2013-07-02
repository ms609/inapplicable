\name{SigSut}
\alias{SigSut}
\alias{SigSut.data}
\alias{SigSut.phy}
\alias{SigSut.preparedata}
\alias{SigSut.parstree}
\alias{SigSut.besttree}
\docType{data}
\title{Phylogenetic data from Sigwart and Sutton 2007}
\description{
Phylogenetic matrix from Sigwart and Sutton (2007), modified per Smith (submitted)
}
\usage{data('SigSut.RData')}
\format{
  The format is:
 chr "SigSut"
}
\details{
The original Sigwart & Sutton dataset contained ordered transformation series (which are not natively supported in R).
In this dataset I manually recoded the original ordered characters as a serial unordered transformation series following the formula:
a.  Original transformation series state: (0), greater than zero; (1), otherwise
b.  Original transformation series state: (0), greater than one; (1), otherwise
c.  Original transformation series state: (0), greater than two; (1), otherwise
The parsimony scores generated using this method are identical to those produced using a single ordered transformation series.

I have also replaced the ambiguous token ('?') with an 'ambiguous but applicable' token ('+') where appropriate.

The dataset contains:
  \item{SigSut.data}{The dataset, in the format produced by read.nexus.data('SigSut-fold.nex')}
  \item{SigSut.phy}{The dataset, in phyDat format}
  \item{SigSut.preparedata}{The dataset, in the format produced by prepare.data(SigSut.data)}
  \item{SigSut.parstree}{The phylogenetic tree (in phylo format) recovered by equal-weights analysis
    of the dataset under the standard parsimony method}
  \item{SigSut.besttree}{The phylogenetic tree (in phylo format) recovered by equal-weights analysis
    of the dataset under the revised method of Smith (submitted)}
}
\source{
Sigwart and Sutton (2007); modified per Smith (submitted)
}
\references{
Sigwart, J. D., & Sutton, M. D. (2007). Deep molluscan phylogeny: synthesis of palaeontological and neontological data. Proceedings of the Royal Society B: Biological Sciences, 274, 2413–2419. doi:10.1098/rspb.2007.0701
M. R. Smith (submitted). Incorporating inapplicable data in parsimony analyses. Systematic Biology XXX, doi: 10.1093/sysbio/syXXXX
}
\examples{
data('SigSut.RData')
contrast3 <- matrix(data=c(1,0,0, 0,1,0, 0,0,1, 1,1,0, 1,1,1), ncol=3, byrow=TRUE)
dimnames(contrast3) <- list(c(0:1, "x", '+', '?'), c(0:1, '-'))
as.phyDat(SigSut.data, type="USER", contrast = contrast3)
plot(SigSut.parstree)
dev.new()
plot(SigSut.besttree)
}
\keyword{datasets}