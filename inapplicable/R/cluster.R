prepare.cluster <- function (cores) {
  cl <- makeCluster(getOption("cl.cores", cores))
  clusterEvalQ(cl, {library(inapplicable);     NULL})
  attr(cl, 'cores') <- cores;
  return(cl)
}