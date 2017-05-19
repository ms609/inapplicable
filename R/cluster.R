PrepareCluster <- function (cores) {
  cl <- makeCluster(getOption("cl.cores", cores))
  clusterEvalQ(cl, {library(inapplicable); NULL})
  setDefaultCluster(cl)
  attr(cl, 'cores') <- cores;
  return(cl)
}