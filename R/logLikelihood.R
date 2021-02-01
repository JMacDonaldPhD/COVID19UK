# COVID-19 Cumulative Cases Observation Model (England UA areas)


#' @func
#' @param y cumulative case data up to current date.
#' @param N Populations of UK regions.
#' @param K kernel describing transmission between populations.
#' @param simulator Function which simulates progession of epidemic.
#' @noSims Number of simulations carried out.


llh = function(y = currentUACaseData(), simulator, noSims) {

  simulations = function(stateX0, epiParam){
    realisations = replicate(noSims, simulator(epiParam), simplify = FALSE)
    #' @param param a vector of model parameters c(beta, gamma, I0W, phi)
    #' @param visualise if TRUE then print out model parameters and a graph of
    # modelled case detections in Wuhan.  Useful for following progress of optimisers
    #' @return the log likelihood of the data conditional on the model and parameters
    logp_fn = function(obsParam, visualise=FALSE) {
      param = c(epiParam, obsParam)
      if(isTRUE(visualise)) {
        # pparam = c(beta=exp(param[1]), gamma=exp(param[2]),
        #            I0W=exp(param[3]), phi=invlogit(param[4]))
        print(param)
      }

      logp = sapply(X = realisations, function(X){
        p_detect = rep(param[4], length(N))

        # Observe increments in R in China
        cumCases = tail.matrix(X$R, n = 1)

        # China observation model
        assertthat::assert_that(all(dim(y) == dim(cumCases)))
        llik = dbinom(y, prob = param[4], size = cumCases, log=TRUE)
        llik = sum(llik)

        llik
      })
      log(mean(exp(logp)))
    }
    logp_fn
  }
  simulations
}


#' @func
#' @param y
#' @param simulator
#' @param noSims

llh.discrete = function(y, modelSim, noSims, parallel = TRUE){

  simulations = function(stateX0, epiParam, kernelParam){
    if(parallel){
      "%dopar%" = foreach::`%dopar%`
      noCores = parallel::detectCores() - 1
      cl = snow::makeCluster(noCores)
      snow::clusterExport(cl, list("modelSim", "transit_spat", "N", "n"), envir= environment() )
      doSNOW::registerDoSNOW(cl)
      realisations =
      foreach::foreach(i = 1:noSims) %dopar% {
        modelSim(stateX0, epiParam, kernelParam)
      }
      snow::stopCluster(cl)
    } else{
      realisations = replicate(noSims, modelSim(stateX0, epiParam, kernelParam), simplify = FALSE)
    }
    #' @param param a vector of model parameters c(beta, gamma, I0W, phi)
    #' @param visualise if TRUE then print out model parameters and a graph of
    # modelled case detections in Wuhan.  Useful for following progress of optimisers
    #' @return the log likelihood of the data conditional on the model and parameters
    logp_fn = function(obsParam, visualise=FALSE) {
      param = c(epiParam, obsParam)
      if(isTRUE(visualise)) {
        print(param)
      }

      logp = sapply(X = realisations, function(X){
        p_detect = rep(obsParam, length(N))

        # Observe increments in R in China
        cumCases = rowSums(X)

        # China observation model
        assertthat::assert_that(all(dim(y) == dim(cumCases)))
        llik = dbinom(y, prob = obsParam, size = cumCases, log=TRUE)
        llik = sum(llik)
        llik
      })
      log(mean(exp(logp)))
    }
    logp_fn
  }
  simulations
}
