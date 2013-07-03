\name{SuttonEtAl}
\alias{SuttonEtAl.data}
\alias{SuttonEtAl.phy}
\alias{SuttonEtAl.preparedata}
\alias{SuttonEtAl.parstree}
\alias{SuttonEtAl.besttree}
\docType{data}
\title{Phylogenetic data from Sutton et al. 2012}
\description{Phylogenetic matrix from Sutton et al. (2012), modified per Smith (submitted).}
\usage{data('SuttonEtAl')}
\format{
  The format is:
 chr "SuttonEtAl"
}
\details{
The original Sutton et al. dataset contained ordered transformation series (which are not natively supported in R).
In this dataset I manually recoded the original ordered characters as a serial unordered transformation series following the formula:
a.  Original transformation series state: (0), greater than zero; (1), otherwise
b.  Original transformation series state: (0), greater than one; (1), otherwise
c.  Original transformation series state: (0), greater than two; (1), otherwise
The parsimony scores generated using this method are identical to those produced using a single ordered transformation series.

I have also replaced the ambiguous token ('?') with an 'ambiguous but applicable' token ('+') where appropriate.

Some coding errors in the matrix were corrected:
  \itemize{
    \item{Babinka, Fordilla, Pojetaia}{Character 12, 0->x (inapplicable)}
    \item{Babinka}{Character 24, 0->x (inapplicable)}
    \item{Chaetoderma}{Character 47, 1->x (inapplicable), to be consistent with character 22}
    \item{Pedanochiton}{Character 42, ?>1, to be consistent with character 48}
  }
The dataset contains:
  \itemize{
    \item{SuttonEtAl.data}{The dataset, in the format produced by read.nexus.data('SigSut-fold.nex')}
    \item{SuttonEtAl.phy}{The dataset, in phyDat format}
    \item{SuttonEtAl.preparedata}{The dataset, in the format produced by prepare.data(SigSut.data)}
    \item{SuttonEtAl.parstree}{The phylogenetic tree (in phylo format) recovered by equal-weights analysis
      of the dataset under the standard parsimony method}
    \item{SuttonEtAl.besttree}{The phylogenetic tree (in phylo format) recovered by equal-weights analysis
      of the dataset under the revised method of Smith (submitted)}
  }
}
\source{
Original source, Sigwart and Sutton (2007); updated by Sutton et al. (2012); modified per Smith (submitted)
}
\references{
Sigwart, J. D., & Sutton, M. D. (2007). \cite{Deep molluscan phylogeny: synthesis of palaeontological and neontological data.} Proceedings of the Royal Society B: Biological Sciences, 274, 2413-2419. doi:10.1098/rspb.2007.0701
Sutton, M. D., Briggs, D. E. G., Siveter, D. J., Siveter, D. J., & Sigwart, J. D. (2012). \cite{A Silurian armoured aplacophoran and implications for molluscan phylogeny.} Nature, 490(7418), 94-97. doi:10.1038/nature11328
M. R. Smith (submitted). \cite{Incorporating inapplicable data in parsimony analyses.} Systematic Biology XXX, doi: 10.1093/sysbio/syXXXX
}
\examples{
data('SuttonEtAl')
SuttonEtAl.phy
plot(SuttonEtAl.parstree)
dev.new()
plot(SuttonEtAl.besttree)
}
\keyword{datasets}