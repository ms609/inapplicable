#' @title PhyDat to MorphyDat
#'
#' @description Converts a PhyDat object to allow processing by MorphyDat
#'
#' @param phydat A \code{phyDat} object
#' 
#' @examples
#' morphy <- MorphyDat(SigSut.phy)
#' 
#' @value This function returns a matrix of class \code{*phyDat}.
#'    \itemize{
#' \item     Each row of the matrix corresponds to an operational taxonomic unit;
#' \item     Each column of the matrix represents a transformation series, with all transformation
#'           series containing three or more inapplicable tokens (and thus requiring the use of the 
#'           modified Fitch algorithm) listed first;
#' \item     Each cell in the matrix contains an integer that denotes the possible tokens.
#'           This integer is the decimal representation of the binary figure represented by a row
#'           in the contrast matrix (in reverse order), after
#'           phangorn:::prepareFitchData.          
#'   }
#'   In addition to the standard \link{phyDat} attributes, the returned object bears an attribute
#'   \code{inapp.level}, specifying the level that corresponds to the inapplicable token, \code{
#'   inapp.chars}, specifying the number of transformation series in which three or more inapplicable
#'   tokens are present (thus requiring the use of the modified Fitch algorithm), and [#TODO]
#'   \code{min.steps}, specifying the minimum number of steps necessary to obtain the observed
#'   distribution of tokens (necessary when calculating implied weights).
#' 
#' @seealso \code{\link{phyDat}}
#' 
#' @author Martin Smith
#' @export
MorphyDat <- function (phydat) {
  # Written with reference to phangorn:::prepareDataFitch
  at <- attributes(phydat) 
  nam <- at$names
  nLevel <- length(at$level)
  nChar <- at$nr
  if (nChar < 1) stop ("No (valid) characters sent to MorphyDat.")
  contrast <- attr(phydat, "contrast")
  nTip <- length(phydat)
  at$names <- NULL
  powersOf2 <- 2L^c(0L:(nLevel - 1L))
  tmp <- contrast %*% powersOf2
  tmp <- as.integer(tmp)
  if (nChar > 1) {
    ret <- t(vapply(phydat, function (x) as.integer(tmp[x]), integer(nChar))) 
  } else {
    ret <- as.integer(phydat)
    ret <- matrix(as.integer(tmp[ret]), ncol=1, dimnames=list(nam, NULL))
  }
  # as.integer so ready to send to C
  inapp.level <- which(at$levels == "-")
  needMorphyTreatment <- colSums(ret == tmp[inapp.level]) > 2
  newOrder <- c(which(needMorphyTreatment), which(!needMorphyTreatment))
  ret <- ret[, newOrder, drop=FALSE]
  attributes(ret) <- at
  class(ret) <- 'morphyDat'
  attr(ret, 'dim') <- c(nTip, nChar)
  dimnames(ret) <- list(nam, NULL)
  attr(ret, 'weight') <- at$weight[newOrder]
  attr(ret, 'index') <- newOrder[at$index]
  #attr(ret, 'nr') <- nChar
  #attr(ret, 'levels') <- at$levels
  #attr(ret, 'allLevels') <- at$allLevels
  #attr(ret, 'type') <- at$type
  #attr(ret, 'contrast') <- contrast
  attr(ret, 'inapp.level') <- 2 ^ (inapp.level - 1)
  attr(ret, 'inapp.chars') <- sum(needMorphyTreatment)
  ret
}

PrepareData <- function (data) {
# Written with reference to phangorn:::prepareDataFitch
  at <- attributes(data)
  nam <- at$names
  nLevel <- length(at$level)
  nChar <- at$nr
  cont <- attr(data, "contrast")
  nTip <- length(data)
  at$names <- NULL
  powers.of.2 <- 2L^c(0L:(nLevel - 1L))
  tmp <- cont %*% powers.of.2
  tmp <- as.integer(tmp)
  data <- unlist(data, FALSE, FALSE)
  ret <- tmp[data] 
  ret <- as.integer(ret)
  attributes(ret) <- at
  inapp.level <- which(at$levels == "-")
  if (!any(inapp.level)) warning ("    No inapplicable tokens detected in data.
    Use a hyphen ('-', not a dash) to denote taxa in which a transformation series is inapplicable.
    Does typing attributes(mydataname, 'contrast') return a column labelled '-'?")
  attr(ret, 'inapp.level') <- 2^(inapp.level - 1)
  attr(ret, 'dim') <- c(nChar, nTip)  
  attr(ret, 'unique.tokens') <- apply(ret, 1, function(x) quick.min(x, inapp.level))
  applicable.tokens <- setdiff(powers.of.2, 2^(inapp.level - 1))
  attr(ret, 'split.sizes') <- apply(ret, 1, function(x) vapply(applicable.tokens, function (y) sum(x==y), integer(1)))
  dimnames(ret) <- list(NULL, nam)
  class(ret) <- '*phyDat'
  ret
}

# No longer needed here - move to 'information' package
quick.min <- function (x, inapp.level) {
  if (length(inapp.level)) return(sum(2^(c(0:(inapp.level-2), inapp.level:12)) %in% unique(x)))
  return (sum(2^(0:12) %in% unique(x)))
}

#' @name AsBinary
#' @alias AsBinary
#' @title Convert a number to binary
#' @description Provides a (reversed) binary representation of a decimal integer
#' @usage AsBinary(x)
#'
#' @param x Decimal integer to be converted to binary bits
#' 
#' @details 
#' Provides an array corresponding to binary digits 1, 2, 4, 8, 16, ...
#' 
#' Binary number 0100 (= decimal 4) will be represented as 0 0 1.
#' 
#' @value 
#' An array corresponding to binary digits 1, 2, 4, 8, 16, ...
#' 
#' 'Leading zeros' are not included.
#' 
#' @author Martin R. Smith, adapted from code posted to R mailing list by Spencer Graves
#' 
#' @examples
#'   AsBinary(4)  # 0 0 1
#'   AsBinary(10) # 0 1 0 1
#' 
#' @export
AsBinary <- function(x) {
	N <- length(x)
	xMax <- max(x)	
	ndigits <- (floor(logb(xMax, base=2))+1)
	Base.b <- array(NA, dim=c(N, ndigits))
	for (i in 1:ndigits){#i <- 1
		Base.b[, i] <- (x %% 2)
		x <- (x %/% 2)
	}
	if(N == 1) Base.b[1, ] else Base.b
}