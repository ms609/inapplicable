\name{prepare.cluster}
\alias{prepare.cluster}
\title{Prepare a cluster for use in tree search functions}
\description{...}
\usage{
prepare.cluster(cores)
}
\arguments{
  \item{cores}{Number of cores to include in the cluster.}
}
\details{
Creates a cluster of multiple cores and prepares it to analyse phylogenetic trees.
}
\value{
Returns a reference to a cluster.
}
\author{Martin R. Smith}

\examples{
  \dontrun{
  prepare.cluster(4)
  }
}