#'  Prepare a cluster for use in tree search functions
#' 
#' @details \code{prepare.cluster} creates a cluster of multiple cores and prepares it to analyse phylogenetic trees.
#' @usage prepare.cluster(cores)
#' 
#' @param cores Number of cores to include in the cluster.
#' 
#' 
#' @return Returns a reference to a cluster.
#' 
#' @author Martin R. Smith
#' 
#' @examples
#' \dontrun{
#'   prepare.cluster(4)
#' }
#' 
#' @importFrom parallel setDefaultCluster clusterEvalQ makeCluster
#' @export
PrepareCluster <- function (cores) {
  cl <- makeCluster(getOption("cl.cores", cores))
  clusterEvalQ(cl, {library(inapplicable); NULL})
  setDefaultCluster(cl)
  attr(cl, 'cores') <- cores;
  return(cl)
}