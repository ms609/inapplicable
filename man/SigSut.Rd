\name{SigSut}
\alias{SigSut}
\alias{SigSut.data}
\alias{SigSut.phy}
\alias{SigSut.preparedata}
\alias{SigSut.parstree}
\alias{SigSut.besttree}
\docType{data}
\title{Data from Sigwart and Sutton 2007}
\description{Contains the phylogenetic matrix from Sigwart and Sutton (2007), modified per Smith & Brazeau (201X),
and the strict consensus of the most parsimonious trees generated under traditional parsimony
and under the method proposed by Smith & Brazeau (201X).}
\usage{data('SigSut')
SigSut.data
SigSut.phy
SigSut.preparedata
SigSut.parstree
SigSut.besttree
}
\format{
  The format is:
  \tabular{ll}{
    \code{SigSut.data       }  \tab \code{list}, as produced by \code{\link{read.nexus.data}}\cr
    \code{SigSut.phy        }  \tab \code{\link{phyDat}} \cr
    \code{SigSut.preparedata}  \tab \code{\link[=PrepareData]{*phyDat}} \cr
    \code{SigSut.parstree   }  \tab \code{\link{phylo}} \cr
    \code{SigSut.besttree   }  \tab \code{\link{phylo}} \cr
    \code{SigSut.outgroup   }  \tab \code{\link{vector}} \cr
  }
}
\details{
Following Sigwart & Sutton (2007), I have used the dataset that follows the Brachiopod Fold Hypothesis.
The original dataset contained ordered transformation series (which are not natively supported in R).
In this dataset I manually recoded the original ordered characters as a serial unordered transformation series following the formula:
\tabular{llll}{a.\tab Original transformation series state:\tab (0), greater than zero;\tab (1), otherwise\cr
b.\tab Original transformation series state:\tab (0), greater than one; \tab (1), otherwise \cr
c.\tab Original transformation series state:\tab (0), greater than two; \tab (1), otherwise\cr
}
The parsimony scores generated using this method are identical to those produced using a single ordered transformation series.

I have also replaced the ambiguous token ('?') with an 'ambiguous but applicable' token ('+') where appropriate.

The dataset contains:
\tabular{ll}{\code{SigSut.data}  \tab      The dataset, in the format produced by \code{read.nexus.data('SigSut-fold.nex')}\cr
  \code{SigSut.phy}   \tab      The dataset, in phyDat format\cr
  \code{SigSut.preparedata}\tab The dataset, in the format produced by \code{PrepareData(SigSut.data)}\cr
  \code{SigSut.parstree}\tab    The phylogenetic tree (in phylo format) recovered by equal-weights analysis\cr
  \tab  of the dataset under the standard parsimony method\cr
  \code{SigSut.besttree}\tab     The phylogenetic tree (in phylo format) recovered by equal-weights analysis\cr
  \tab of the dataset under the revised method of Smith & Brazeau (201X)\cr
  \code{SigSut.outgroup}\tab     The outgroup, \code{c('Lingula', 'Neocrania', 'Mickwitzia')}\cr
}
}
\source{
Sigwart and Sutton (2007); modified per Smith & Brazeau (201X)
}
\references{
\itemize{
\item Sigwart, J. D., & Sutton, M. D. (2007). \cite{Deep molluscan phylogeny: synthesis of palaeontological and neontological data.} Proceedings of the Royal Society B: Biological Sciences, 274, 2413-2419. doi:\href{http://dx.doi.org/10.1098/rspb.2007.0701}{10.1098/rspb.2007.0701}
\item Smith, M. R., & Brazeau, M. D. (201X). \cite{Incorporating inapplicable data in parsimony analyses.} Systematic Biology XXX, doi: \href{http://dx.doi.org/10.1093/sysbio/syXXXX}{10.1093/sysbio/syXXXX}
}
}
\examples{
data('SigSut')
contrast3 <- matrix(data=c(1,0,0, 0,1,0, 0,0,1, 1,1,0, 1,1,1), ncol=3, byrow=TRUE)
dimnames(contrast3) <- list(c(0:1, "x", '+', '?'), c(0:1, '-'))
phyDat(SigSut.data, type="USER", contrast = contrast3)
plot(SigSut.parstree)
dev.new()
plot(SigSut.besttree)
}
\keyword{datasets}