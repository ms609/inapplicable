\name{SuttonEtAl}
\alias{SuttonEtAl}
\alias{SuttonEtAl.data}
\alias{SuttonEtAl.phy}
\alias{SuttonEtAl.preparedata}
\alias{SuttonEtAl.parstree}
\alias{SuttonEtAl.besttree}
\docType{data}
\title{Data from Sutton et al. 2012}
\description{
\code{#TODO!  } Re-generate ALL data to avoid typo of PolysacCos.  Re-generate preparedate with new function.

Phylogenetic matrix from Sutton et al. (2012), modified per Smith & Brazeau (submitted),
and the strict consensus of the most parsimonious trees generated under traditional parsimony
and under the method proposed by Smith & Brazeau (submitted).}
\usage{data('SuttonEtAl')
SuttonElAl.data
SuttonElAl.phy
SuttonElAl.preparedata
SuttonElAl.parstree
SuttonElAl.besttree}
\format{
  The format is:
  \tabular{ll}{
    \code{SuttonEtAl.data       }  \tab \code{list}, as produced by \code{\link{read.nexus.data}}\cr
    \code{SuttonEtAl.phy        }  \tab \code{\link{phyDat}} \cr
    \code{SuttonEtAl.preparedata}  \tab \code{\link[=prepare.data]{*phyDat}} \cr
    \code{SuttonEtAl.parstree   }  \tab \code{\link{phylo}} \cr
    \code{SuttonEtAl.besttree   }  \tab \code{\link{phylo}} \cr
  }
}
\details{
The original Sutton et al. dataset contained ordered transformation series (which are not natively supported in R).
In this dataset I manually recoded the original ordered characters as a serial unordered transformation series following the formula:
\tabular{llll}{a.\tab Original transformation series state:\tab (0), greater than zero;\tab (1), otherwise\cr
b.\tab Original transformation series state:\tab (0), greater than one; \tab (1), otherwise \cr
c.\tab Original transformation series state:\tab (0), greater than two; \tab (1), otherwise\cr
}
The parsimony scores generated using this method are identical to those produced using a single ordered transformation series.

I have also replaced the ambiguous token ('?') with an 'ambiguous but applicable' token ('+') where appropriate.

Some coding errors in the matrix were corrected:
  \itemize{
    \item \emph{Babinka, Fordilla, Pojetaia}: Character 12, 0->x (inapplicable)
    \item \emph{Babinka}: Character 24, 0->x (inapplicable)
    \item \emph{Chaetoderma}: Character 47, 1->x (inapplicable), to be consistent with character 22
    \item \emph{Pedanochiton}: Character 42, ?->1, to be consistent with character 48
  }
The dataset contains:
\tabular{ll}{\code{SuttonEtAl.data}  \tab      The dataset, in the format produced by \code{read.nexus.data('SuttonEtAl.nex')}\cr
  \code{SuttonEtAl.phy}   \tab      The dataset, in phyDat format\cr
  \code{SuttonEtAl.preparedata}\tab The dataset, in the format produced by \code{prepare.data(SuttonEtAl.data)}\cr
  \code{SuttonEtAl.parstree}\tab    The phylogenetic tree recovered by equal-weights analysis\cr
  \tab  of the dataset under the standard parsimony method\cr
  \code{SuttonEtAl.besttree}\tab     The phylogenetic tree recovered by equal-weights analysis\cr
  \tab of the dataset under the revised method of Smith & Brazeau (submitted)\cr
}
}
\source{
Original source, Sigwart and Sutton (2007); updated by Sutton et al. (2012); modified per Smith & Brazeau (submitted)
}
\references{
\itemize{
\item Sigwart, J. D., & Sutton, M. D. (2007). \cite{Deep molluscan phylogeny: synthesis of palaeontological and neontological data.} Proceedings of the Royal Society B: Biological Sciences, 274, 2413-2419. doi:\href{http://dx.doi.org/10.1098/rspb.2007.0701}{10.1098/rspb.2007.0701}
\item Sutton, M. D., Briggs, D. E. G., Siveter, D. J., Siveter, D. J., & Sigwart, J. D. (2012). \cite{A Silurian armoured aplacophoran and implications for molluscan phylogeny.} Nature, 490(7418), 94-97. doi:\href{http://dx.doi.org/10.1038/nature11328}{10.1038/nature11328}
\item Smith, M. R., & Brazeau, M. D. (submitted). \cite{Incorporating inapplicable data in parsimony analyses.} Systematic Biology XXX, doi:\href{http://dx.doi.org/10.1093/sysbio/syXXXX}{10.1093/sysbio/syXXXX}
}
}
\examples{
data('SuttonEtAl')
SuttonEtAl.phy
plot(SuttonEtAl.parstree)
dev.new()
plot(SuttonEtAl.besttree)
}
\keyword{datasets}