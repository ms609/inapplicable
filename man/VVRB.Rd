\name{VVRB}
\alias{VVRB}
\alias{VVRB.data}
\alias{VVRB.phy}
\alias{VVRB.preparedata}
\alias{VVRB.besttree}
\docType{data}
\title{Data from Vinther et al. 2008}
\description{
Phylogenetic matrix from Vinther, Van Roy and Briggs, modified per Smith & Brazeau (201X),
and the strict consensus of the most parsimonious trees generated under traditional parsimony
and under the method proposed by Smith & Brazeau (201X).
}
\usage{data('VVRB')
VVRB.data
VVRB.phy
VVRB.preparedata
VVRB.besttree}
\format{
  The format is:
  \tabular{ll}{
    \code{VVRB.data       }  \tab \code{list}, as produced by \code{\link{read.nexus.data}}\cr
    \code{VVRB.phy        }  \tab \code{\link{phyDat}} \cr
    \code{VVRB.preparedata}  \tab \code{\link[=PrepareData]{*phyDat}} \cr
    \code{VVRB.besttree   }  \tab \code{\link{phylo}} \cr
  }
}
\details{
The ambiguous token ('?') was replaced with an 'ambiguous but applicable' token ('+') where appropriate.

The dataset contains:
\tabular{ll}{\code{VVRB.data}  \tab      The dataset, in the format produced by \code{read.nexus.data('VVRB.nex')}\cr
  \code{VVRB.phy}   \tab      The dataset, in \code{phyDat} format\cr
  \code{VVRB.preparedata}\tab The dataset, in the format produced by \code{PrepareData(VVRB.data)}\cr
  \code{VVRB.besttree}\tab     The phylogenetic tree recovered by equal-weights analysis\cr
  \tab of the dataset under the revised method of Smith & Brazeau (201X)\cr
}
}
\source{
Vinther, Van Roy & Briggs (2008) modified per Smith & Brazeau (201X)
}
\references{
\itemize{
\item Vinther, J., Van Roy, P., & Briggs, D. E. G. (2008). \cite{Machaeridians are Palaeozoic armoured annelids.} Nature, 451(7175), 185-188. doi:\href{http://dx.doi.org/10.1038/nature06474}{10.1038/nature06474}
\item Smith, M. R., & Brazeau, M. D. (201X). \cite{Incorporating inapplicable data in parsimony analyses.} Systematic Biology XXX, doi:\href{http://dx.doi.org/10.1093/sysbio/syXXXX}{10.1093/sysbio/syXXXX}
}}
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
identical(VVRB.phy, phyDat(VVRB.data, type="USER", contrast=contrast5))
VVRB.preparedata # <- PrepareData(VVRB.phy)
outgroup <- c('Nemertean')
plot(VVRB.besttree)
}
\keyword{datasets}