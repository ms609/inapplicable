#' @title PhyDat to MorphyDat
#'
#' @description Converts a PhyDat object to allow processing by MorphyDat
#'
#' @param phydat A \code{phyDat} object
#' 
#' @examples
#' morphy <- MorphyDat(SigSut.phy)
#' 
#' @return This function returns a matrix of class \code{*phyDat}.
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

#' @title String to MorphyDat
#'
#' @description Converts a PhyDat object to allow processing by MorphyDat
#'
#' @param string a string of tokens, optionally containing newlines, with no terminating semi-colon.  Polytomies not (yet) supported; each character must correspond to a unique state, ?, or the inapplicable token (-)
#' @param tips, a character vector corresponding to the names (in order) of each taxon in the matrix
#' @param byTaxon = TRUE, string is one TAXON's coding at a time; FALSE: one CHARACTER's coding at a time
#' 
#' @examples
#' morphy <- StringToMorphy("-?01231230?-", c('Lion', 'Gazelle'), byTaxon=TRUE)
#' # encodes the following matrix:
#' # Lion     -?0123
#' # Gazelle  1230?-
#' 
#' @return This function returns a matrix of class \code{morphyDat}; see \code{\link{MorphyData}}
#' @seealso \code{\link{phyDat}}
#' 
#' @author Martin Smith
#' @export
StringToMorphy <- function (x, tips, byTaxon = TRUE) {
  x <- strsplit(x, '')[[1]]
  x <- matrix(x[x != '\n'], nrow=length(tips), byrow=byTaxon)
  rownames(x) <- tips
  phy <- phyDat(x, levels=c(0:3, '-'), type='USER')
  MorphyDat(phy)
}

#' Extract character data from a phyDat object as a string
#' 
#' 
#' @param phy an object of class \code{\link{phyDat}}
#' @param ps text, perhaps ';', to append to the end of the string
#' 
#' @author Martin Smith
#' @importFrom phangorn phyDat
#' @export
PhyToString <- function (phy, ps='') {
  at <- attributes(phy)
  phyLevels <- at$allLevels
  phyChars <- at$nr
  phyContrast <- at$contrast == 1
  outLevels <- seq_len(ncol(phyContrast)) - 1
  if (any(inappLevel <- phyLevels == '-')) outLevels[which(phyContrast[inappLevel])] <- '-'
  levelTranslation <- apply(phyContrast, 1, function (x)  ifelse(sum(x) == 1, as.character(outLevels[x]), paste0(c('{', outLevels[x], '}'), collapse='')))
  if (any(ambigToken <- apply(phyContrast, 1, all))) levelTranslation[ambigToken] <- '?'
  ret <- paste0(c(t(vapply(phy, function (x) levelTranslation[x], character(phyChars))), ps), collapse='')
  return (ret)
}

#' Initialize a Morphy Object from a phyDat object
#' 
#' Creates a new Morphy object with the same size and characters as the phyDat object 
#' @param phy an object of class \code{\link{phyDat}}
#' @return a pointer to a Morphy object, with the attribute "weight" corres
#' 
#' @author Martin Smith
#' @importFrom phangorn phyDat
#' @export
LoadMorphy <- function (phy) {
  morphyObj <- mpl_new_Morphy()
  nTax <- length(phy)
  nChar <- length(phy[[1]])
  if(mpl_init_Morphy(nTax, nChar, morphyObj) -> error) stop("Error ", mpl_translate_error(error), " in mpl_init_Morphy")
  if(mpl_attach_rawdata(PhyToString(phy, ';'), morphyObj) -> error) stop("Error ", mpl_translate_error(error), " in mpl_attach_rawdata")
  if(mpl_set_num_internal_nodes(nTax + 1L, morphyObj) -> error) stop("Error ", mpl_translate_error(error), " in mpl_set_num_internal_nodes")
  if(mpl_apply_tipdata(morphyObj) -> error) stop("Error ", mpl_translate_error(error), " in mpl_apply_tipdata")
  weight <- attr(phy, 'weight')
  if (any(vapply(seq_along(weight), function (x) mpl_set_charac_weight(x, weight[x], morphyObj),
      integer(1)) -> error)) stop("Error ", mpl_translate_error(min(error)), "in mpl_set_charac_weight")
  
  mpl_get_charac_weight(3, morphyObj)
  
  return(morphyObj)
}

#' @name AsBinary
#' @aliases AsBinary
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
#' @return 
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