% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/morphyex.R
\name{mpl_get_symbols}
\alias{mpl_get_symbols}
\title{Retrieves the current list of symbols.}
\usage{
mpl_get_symbols(morphyobj)
}
\arguments{
\item{morphyobj}{An instance of the Morphy object.}
}
\value{
A C-style (null-terminated) string of the character state symbols
being used. NULL if failure.
}
\description{
Returns a pointer to the string of character state symbols
currently being used by Morphy (i.e. either the list of symbols extracted
from the matrix, or the caller-specified values).
}
\author{
Martin Brazeau
}
\keyword{internal}
