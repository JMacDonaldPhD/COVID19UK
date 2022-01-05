# Population Structure Generator

#'@param M number of metapopulations
#'@param N population size of each metapopulation
#'@param within population mixing dynnamics within each metapopulation
#'@param between population mixing dynamics between each metapopulation
metaPopStruct <- function(N_M, within = "H", between = 0.5){
  M <- length(N_M)
  mixingMat <- matrix(between, nrow = M, ncol = M)
  if(within == "H"){
    diag(mixingMat) <- 1
  } else{
    diag(mixingMat) <- within
  }

  return(list(N_M = N_M, mixingMat = mixingMat))
}
