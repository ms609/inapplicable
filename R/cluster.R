#'  Prepare a cluster for use in tree search functions
#' \code{prepare.cluster}
#' @description ...
#' \usage{
#' prepare.cluster(cores)
#' }
#' \arguments{
#'   \item{cores}{Number of cores to include in the cluster.}
#' }
#' \details{
#' Creates a cluster of multiple cores and prepares it to analyse phylogenetic trees.
#' }
#' @return{
#' Returns a reference to a cluster.
#' }
#' @author Martin R. Smith
#' 
#' @examples{
#'   \dontrun{
#'   prepare.cluster(4)
#'   }
#' }
PrepareCluster <- function (cores) {
  cl <- makeCluster(getOption("cl.cores", cores))
  clusterEvalQ(cl, {library(inapplicable); NULL})
  setDefaultCluster(cl)
  attr(cl, 'cores') <- cores;
  return(cl)
}