# PM_MCMC for COVID Data allowing for different likelihoods.


PM_MCMCCOVIDLancsLlh.ADAPTIVE = function(totalHospAdmissions, fixedParam, alpha0, alphaPrior, kernelParam, simulator, likelihood, noSims,
                                      noIts, burnIn, lambda0, delta, runParallel = TRUE){

  start = as.numeric(Sys.time())


  # Adaptive Step
  alphaCurr = alpha0

  lambda = lambda0
  # Likelihood Function
  logPostCurr = -Inf

  while(logPostCurr == -Inf){

    if(runParallel){
      cl <- parallel::makePSOCKcluster(parallel::detectCores() - 1)
      doParallel::registerDoParallel(cl)
      sims <-
        foreach(1:noSims) %dopar% {
          simulator(param = fixedParam, kernelParam)$dailyNoCases
        }
      parallel::stopCluster(cl)

    } else{
      sims <- lapply(1:noSims, function(X) simulator(param = fixedParam, kernelParam)$dailyNoCases)
    }

    #sims = replicate(noSims, simulator(P0, epiParam = c(fixedParam[1:4], alphaCurr, fixedParam[5])), simplify = F)
    logPCurr = sapply(X = sims, function(X) totalHospLikelihood(colSums(X), totalHospAdmissions,
                                                                alphaCurr))
    logPCurr = log(mean(exp(logPCurr)))
    logPostCurr = logPCurr + alphaPrior(alphaCurr)
  }

  draws = matrix(nrow = noIts + 1, ncol = 2)
  draws[1, ] = c(alphaCurr, logPostCurr)
  noAccept = 0

  print("Sampling Progress")
  pb <- progress::progress_bar$new(total = noIts)
  for(i in 1:noIts){
    pb$tick()
    # Proposal

    u.delta = runif(1, 0, 1)
    if(u.delta > delta){
      alphaProp = RWMProposalUnit(alphaCurr, lambda)
    } else{
      alphaProp = RWMProposalUnit(alphaCurr, lambda0)
    }

    if(runParallel){
      cl <- parallel::makePSOCKcluster(parallel::detectCores() - 1)
      doParallel::registerDoParallel(cl)
      sims <-
        foreach(1:noSims) %dopar% {
          simulator(param = c(fixedParam, alphaProp), kernelParam)$dailyNoCases
        }
      parallel::stopCluster(cl)
    } else{
      sims <- lapply(1:noSims, function(X) simulator(param = c(fixedParam, alphaProp), kernelParam)$dailyNoCases)
    }


    #sims = replicate(noSims, simulator(P0, epiParam = c(fixedParam[1:4], alphaProp, fixedParam[5])), simplify = F)
    logPProp = sapply(X = sims, function(X) totalHospLikelihood(colSums(X), totalHospAdmissions,
                                                                alphaProp))
    logPProp = log(mean(exp(logPProp)))
    logPostProp= logPProp + alphaPrior(alphaProp)
    acceptProb = logPostProp - logPostCurr
    #print(acceptProb)

    u1 = runif(1)
    if(log(u1) < acceptProb){
      noAccept = noAccept + 1
      logPostCurr = logPostProp
      alphaCurr = alphaProp
      if(u.delta > delta){
        lambda = lambda + 0.93*(lambda/sqrt(i))
      }
    } else{
      if(u.delta > delta){
        lambda = lambda - 0.07*(lambda/sqrt((i)))
      }
    }
    draws[i+1, ] = c(alphaCurr, logPostCurr)
  }
  timeTaken = as.numeric(Sys.time()) - start

  # Calculate Effective Sample Sizes (and Per Second) and Acceptance Rates
  ESS <- coda::effectiveSize(draws[, 1])
  ESS.sec <- ESS/timeTaken
  acceptRate <-  noAccept/noIts
  print(c("Accept Rate:", acceptRate))

  return(list(draws = draws, lambda = lambda, noSims = noSims, ESS.sec = ESS.sec, acceptRate = acceptRate))
}


PM_MCMCCOVIDLancsLlh = function(totalHospAdmissions, fixedParam, alpha0, alphaPrior, kernelParam, simulator, lambda, noSims, noIts,
                             burnIn, runParallel = TRUE){

  start = as.numeric(Sys.time())
  # Adaptive Step

  alphaCurr = alpha0
  # Likelihood Function
  logPostCurr = -Inf

  while(logPostCurr == -Inf){

    if(runParallel){
      cl <- parallel::makePSOCKcluster(parallel::detectCores() - 1)
      doParallel::registerDoParallel(cl)
      sims <-
        foreach(1:noSims) %dopar% {
          simulator(param = c(fixedParam, alphaCurr), kernelParam)$dailyNoCases
        }
      parallel::stopCluster(cl)

    } else{
      sims = replicate(noSims, simulator(param = c(fixedParam, alphaCurr), kernelParam)$dailyNoCases, simplify = F)
    }


    logPCurr = sapply(X = sims, function(X) totalHospLikelihood(colSums(X), totalHospAdmissions,
                                                                alphaCurr))
    logPCurr = log(mean(exp(logPCurr)))
    logPostCurr = logPCurr + alphaPrior(alphaCurr)
  }

  draws = matrix(nrow = noIts + 1, ncol = 2)
  draws[1, ] = c(alphaCurr, logPostCurr)
  noAccept = 0

  print("Sampling Progress")
  pb <- progress::progress_bar$new(total = noIts)
  for(i in 1:noIts){
    pb$tick()
    # Proposal
    alphaProp = RWMProposalUnit(alphaCurr, lambda)

    if(runParallel){
      cl <- parallel::makePSOCKcluster(parallel::detectCores() - 1)
      doParallel::registerDoParallel(cl)
      sims <-
        foreach(1:noSims) %dopar% {
          simulator(param = c(fixedParam, alphaProp), kernelParam)$dailyNoCases
        }
      parallel::stopCluster(cl)
    } else{
      sims <- replicate(noSims, simulator(param = c(fixedParam, alphaProp), kernelParam)$dailyNoCases, simplify = F)
    }


    logPProp <- sapply(X = sims, function(X) totalHospLikelihood(colSums(X), totalHospAdmissions,
                                                                 alphaProp))
    logPProp <- log(mean(exp(logPProp)))
    logPostProp <- logPProp + alphaPrior(alphaProp)
    acceptProb <- logPostProp - logPostCurr
    #print(acceptProb)
    u1 <- runif(1)
    if(log(u1) < acceptProb){
      noAccept <- noAccept + 1
      logPostCurr <- logPostProp
      alphaCurr <- alphaProp
    }

    draws[i+1, ] <- c(alphaCurr, logPostCurr)
  }
  timeTaken <- as.numeric(Sys.time()) - start
  # Calculate Effective Sample Sizes (and Per Second) and Acceptance Rates
  ESS <- coda::effectiveSize(draws[, 1])
  ESS.sec <- ESS/timeTaken
  acceptRate <-  noAccept/noIts
  print(c("Accept Rate:", acceptRate))

  return(draws)
}



