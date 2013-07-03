\name{fitch.combine}
\alias{fitch.combine}
\alias{fitch.combine.single}

\title{Modified Fitch downpass}
\description{
Given two nodes, this function generates the tokens implied for the parent node by the Fitch 
algorithm as modified in Smith (submitted).  It also determines the addition to parsimony score
implied at this node.

\code{fitch.combine.single} should be provided with the values of a single transformation series;
\code{fitch.combine} handles multiple transformation series simeltaneously.
}
\usage{
fitch.combine(a, b, inapplicable.token.index)
fitch.combine.single(a, b, inapplicable.token.index)
}

\arguments{
  \item{a}{
\code{fitch.combine.single}: a logical vector describing the permissable tokens for the left child node;

\code{fitch.combine}: a matrix describing the permissable tokens for the left child node (see Details);
}
  \item{b}{As \code{a}, for the right child node;}
  \item{inapplicable.token.index}{Which token is the inapplicable token?}
}
\details{
Calculates each character's contribution to the parsimony score of a tree at a specified node 
using the modified Fitch algorithm proposed in Smith (submitted).

In \code{fitch.combine.single}, the possible tokens present at each node are represented by a 
logical vector.  Say that a transformation series can take a token from \code{{0, 1, 2, -}}.
This transformation series will be represented at each node by a vector such as:
\tabular{ll}{node token is 0       \tab \code{c(TRUE, FALSE, FALSE, FALSE)}\cr
node token is 1       \tab \code{c(FALSE, TRUE, FALSE, FALSE)}\cr
node token is 0 or 2  \tab \code{c(TRUE, FALSE, TRUE, FALSE)} \cr
node token is -       \tab \code{c(FALSE, FALSE, FALSE, TRUE)}
}

\code{fitch.combine} can handle multiple transformation series at once, and should be fed a 
matrix where each column represents a separate transformation series.

}
\value{
This function returns a list with the values:
  \item{ 1 }{Possible tokens at specified node;}
  \item{ 2 }{Addition to parsimony score associated with each transformation series.}
}
\references{
M. R. Smith (submitted). \cite{Incorporating inapplicable data in parsimony analyses.} Systematic Biology XXX, doi:\href{http://dx.doi.org/10.1093/sysbio/syXXXX}{10.1093/sysbio/syXXXX}
}
\author{Martin R. Smith}
\examples{
## Child nodes have tokens '-' and '0'; transformation series can take values {0, 1, -}
fitch.combine.single(c(FALSE, FALSE, TRUE), c(TRUE, FALSE, FALSE), 3)
## Child nodes have tokens 'B' and '+' (ambiuguous but applicable)
## Transformation series can take values {A, B, C, -}
fitch.combine.single(c(FALSE, TRUE, FALSE, FALSE), c(TRUE, TRUE, TRUE, FALSE), 4)

## Reconstruct the possible token states in a node whose children 
## are the terminals 'Wiwaxia' and 'Halkieria'
data('SigSut')
fitch.combine(SigSut.preparedata$Wiwaxia, SigSut.preparedata$Halkieria, 3)

## Fitch.combine.single is currently defined as  
function (a, b, inapplicable.token.index) 
{
  shared <- a & b
  if (any(shared)) return (list(shared, 0))
  if (!any(a[-inapplicable.token.index]) || !any(b[-inapplicable.token.index]))
    return (list(a|b, 0))
  ret <- a|b
  ret[inapplicable.token.index] <- FALSE
  return (list(ret, 1))
}
}
