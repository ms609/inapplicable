\name{as.binary}
\alias{as.binary}
\title{Convert a number to binary}
\description{Provides a (reversed) binary representation of a decimal integer}
\usage{
as.binary(x)
}
\arguments{
  \item{x}{Decimal integer to be converted to binary bits}
}
\details{
Provides an array corresponding to binary digits 1, 2, 4, 8, 16, ...

Binary number 0100 (= decimal 4) will be represented as 0 0 1.
}
\value{
An array corresponding to binary digits 1, 2, 4, 8, 16, ...

'Leading zeros' are not included.
}
\author{Martin R. Smith, adapted from code posted to R mailing list by Spencer Graves}

\examples{
  as.binary(4)  # 0 0 1
  as.binary(10) # 0 1 0 1
}