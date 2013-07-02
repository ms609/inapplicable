\name{VVRB}
\alias{VVRB.data}
\alias{VVRB.besttree}
\docType{data}
\title{Phylogenetic data from Vinther et al. 2008}
\description{
Phylogenetic matrix from Vinther, Van Roy and Briggs, modified per Smith (submitted)
}
\usage{data('VVRB.RData')}
\format{
  The format is:
 chr "VVRB"
}
\details{
The ambiguous token ('?') was replaced with an 'ambiguous but applicable' token ('+') where appropriate.

The dataset contains two items:
  \itemize{
    \item{VVRB.data}{The dataset, in the format produced by read.nexus.data('VVRB.nex')}
    \item{VVRB.besttree}{The phylogenetic tree (in phylo format) recovered by equal-weights analysis
      of the dataset under the revised method of Smith (submitted)}
  }
}
\source{
Vinther, Van Roy & Briggs (2008) modified per Smith (submitted)
}
\references{
Vinther, J., Van Roy, P., & Briggs, D. E. G. (2008). Machaeridians are Palaeozoic armoured annelids. Nature, 451(7175), 185-188. doi:10.1038/nature06474
M. R. Smith (submitted). Incorporating inapplicable data in parsimony analyses. Systematic Biology XXX, doi: 10.1093/sysbio/syXXXX
}
\examples{
data('VVRB')
contrast5 <- matrix(data=c(
 1,0,0,0,0,
 0,1,0,0,0,
 0,0,1,0,0, 
 0,0,0,1,0,
 0,0,0,0,1,
 1,1,0,0,0,
 0,1,1,0,0,
 1,1,1,1,0,
 1,1,1,1,1), ncol=5, byrow=TRUE);
 dimnames(contrast5) <- list(c(0:3, '-', "*", '&', '+', '?'), c(0:3, '-'))
vvrb.phy <- phyDat(vvrb.dat, type="USER", contrast=contrast5);
prepare.data(vvrb.phy)
outgroup <- c('Nemertean')
plot(VVRB.besttree)
}
\keyword{datasets}