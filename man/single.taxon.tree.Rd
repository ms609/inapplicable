\name{single.taxon.tree}
\alias{single.taxon.tree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a phylogenetic 'tree' that comprises a single taxon
}
\description{
This function returns a phyDat object containing a single tip.
}
\usage{
single.taxon.tree(label)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{label}{A character vector specifying the label of the tip}
}
\value{
A phyDat object containing a single tip with the specified label
}
\seealso{
\code{\link{two.tip.tree}}
}
\examples{
single.taxon.tree('Homo_sapiens')
}
\keyword{ tree }
