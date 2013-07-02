\name{fitch.combine}
\alias{fitch.combine}
\alias{fitch.combine.single}

\title{Modified Fitch downpass}
\description{
Generates the addition to parsimony score implied at a node where the child nodes have character 
states given by \code{'a'} and \code{'b'}.  fitch.combine.single calculates a single 
transformation series; fitch.combine calculates a matrix of multiple transformation series.
}
\usage{
fitch.combine(a, b, inapplicable.level)
}

\arguments{
  \item{a}{
In fitch.combine.single, a logical vector describing the possible character states present at the left child node.
In fitch.combine, a matrix where each row lists the possible character states present at the left child node.
}
  \item{b}{
As \code{a}, for the right child node
}
  \item{inapplicable.level}{
The level that represents the inapplicable token.
}
}
\details{
Calculates each character's contribution to the parsimony score of a tree at a specified node 
using the modified Fitch algorithm proposed in Smith (submitted).
The possible states of the each node are represented in a vector.  If the transformation series in 
a dataset can take the values \{0, 1, 2, -\}, then each transformation series will be represented at 
each node by a vector such as:
  node token is 0 ... \code{c(TRUE, FALSE, FALSE, FALSE)}
  node token is 1 ... \code{c(FALSE, TRUE, FALSE, FALSE)}
  node token is 0 or 2 ... \code{c(TRUE, FALSE, TRUE, FALSE)}
  node token is - ... code{c(FALSE, FALSE, FALSE, TRUE)}

}
\value{
  \item{comp1 }{Possible tokens at specified node}
  \item{comp2 }{Addition to parsimony score associated with each transformation series}
}
\references{
M. R. Smith (submitted). Incorporating inapplicable data in parsimony analyses. Systematic Biology XXX, doi: 10.1093/sysbio/syXXXX
}
\author{
Martin R. Smith
}
\examples{
## Child nodes have tokens '-' and '0'
fitch.combine.single(c(FALSE, FALSE, TRUE), c(TRUE, FALSE, FALSE), 3)
## Child nodes have tokens '1' and '+' (ambiuguous applicable); fourth token represents inapplicable state
fitch.combine.single(c(FALSE, TRUE, FALSE, FALSE), c(TRUE, TRUE, TRUE, FALSE), 4)

## Reconstruct the possible token states in a node whose children are the terminals 'Wiwaxia' and 'Halkieria'
data('SigSut')
fitch.combine(SigSut.preparedata$Wiwaxia, SigSut.preparedata$Halkieria, 3)

## The function is currently defined as
function (a, b, inapplicable.level) 
{
    shared <- a & b
    unions <- as.logical(colSums(shared))
    al <- a[-inapplicable.level, ]
    bl <- b[-inapplicable.level, ]
    one.child.inapp.only <- !colSums(al) | !colSums(bl)
    intersects <- !unions
    ret <- matrix(FALSE, nrow(a), nCol <- ncol(a))
    ret[, unions] <- shared[, unions]
    ret[, intersects] <- a[, intersects] | b[, intersects]
    ret[inapplicable.level, !one.child.inapp.only] <- FALSE
    col.score <- rep(TRUE, nCol)
    col.score[unions] <- FALSE
    col.score[one.child.inapp.only] <- FALSE
    return(list(ret, col.score))
  }

## Fitch.combine.single is currently defined as  
function (a, b, inapplicable.level) 
{
  shared <- a & b
  if (any(shared)) return (list(shared, 0))
  if (!any(a[-inapplicable.level]) || !any(b[-inapplicable.level])) return (list(a|b, 0))
  ret <- a|b
  ret[inapplicable.level] <- FALSE
  return (list(ret, 1))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ tre }
