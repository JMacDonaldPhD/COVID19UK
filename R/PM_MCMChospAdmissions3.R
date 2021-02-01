

PM_MCMCCOVIDLancs3.ADAPTIVE = function(totalHospAdmissions, P0, theta0, priors, simulator, noSims,
                                       noIts, burnIn, lambda0, delta){

  start = as.numeric(Sys.time())
  # Adaptive Step

  n = length(theta0)
  thetaCurr = theta0
  type = c(2, 2, 2, 3, 3, 3)
  V0 = diag(1, n)

  lambda = lambda0
  # Likelihood Function
  logPostCurr = -Inf

  while(logPostCurr == -Inf){
    sims = replicate(noSims, simulator(P0,
                                       epiParam = thetaCurr),
                     simplify = F)
    logPCurr = sapply(X = sims, function(X) totalHospLikelihood(colSums(X$dailyNoCases), totalHospAdmissions,
                                                                thetaCurr[5]))
    logPCurr = log(mean(exp(logPCurr)))
    logPostCurr = logPCurr + sum(evalPrior(thetaCurr, priors))
  }

  draws = matrix(nrow = noIts + 1, ncol = n + 1)
  draws[1, ] = c(thetaCurr, logPostCurr)
  noAccept = 0

  print("Sampling Progress")
  pb <- progress::progress_bar$new(total = noIts)
  for(i in 1:noIts){
    pb$tick()
    # Proposal

    u.delta = runif(1, 0, 1)
    if(u.delta > delta & noAccept >= 10){
      V = var(draws[,1:n], na.rm = T)
      thetaProp = RWMProposalMixed(thetaCurr, lambda, V, type)
    } else{
      thetaProp = RWMProposalMixed(thetaCurr, lambda0, V0, type)
    }
    #print(thetaProp)
    sims = replicate(noSims, simulator(P0,
                                       epiParam = thetaCurr),
                     simplify = F)
    logPProp = sapply(X = sims, function(X) totalHospLikelihood(colSums(X$dailyNoCases), totalHospAdmissions,
                                                                thetaProp[5]))
    logPProp = log(mean(exp(logPProp)))
    logPostProp = logPProp + sum(evalPrior(thetaProp, priors))
    acceptProb = logPostProp - logPostCurr
    #print(acceptProb)

    u1 = runif(1)
    if(log(u1) < acceptProb){
      noAccept = noAccept + 1
      logPostCurr = logPostProp
      thetaCurr = thetaProp
      if(u.delta > delta){
        lambda = lambda + 0.93*(lambda/sqrt(i))
      }
    } else{
      if(u.delta > delta){
        lambda = lambda - 0.07*(lambda/sqrt((i)))
      }
    }
    #print(thetaCurr)
    draws[i+1, ] = c(thetaCurr, logPostCurr)
  }
  timeTaken = as.numeric(Sys.time()) - start
  # Calculate Effective Sample Sizes (and Per Second) and Acceptance Rates
  ESS <- min(coda::effectiveSize(draws[, 1:n]))
  ESS.sec <- ESS/timeTaken
  acceptRate <-  noAccept/noIts


  print(c("Accept Rate:", acceptRate))

  return(list(draws = draws, lambda = lambda, V = V, noSims = noSims, ESS.sec = ESS.sec))
}

PM_MCMCCOVIDLancs3 = function(totalHospAdmissions, P0, theta0, priors, simulator, lambda, V, noSims, noIts,
                              burnIn){

  start = as.numeric(Sys.time())
  # Adaptive Step
  n = length(theta0)
  thetaCurr = theta0
  type = c(2, 2, 2, 3, 3, 3)
  # Likelihood Function
  logPostCurr = -Inf

  while(logPostCurr == -Inf){
    sims = replicate(noSims, simulator(P0,
                                       epiParam = thetaCurr),
                     simplify = F)
    logPCurr = sapply(X = sims, function(X) totalHospLikelihood(colSums(X$dailyNoCases), totalHospAdmissions,
                                                                thetaCurr[5]))
    logPCurr = log(mean(exp(logPCurr)))
    logPostCurr = logPCurr + sum(evalPrior(thetaCurr, priors))
  }

  draws = matrix(nrow = noIts + 1, ncol = n + 1)
  draws[1, ] = c(thetaCurr, logPostCurr)
  noAccept = 0

  print("Sampling Progress")
  pb <- progress::progress_bar$new(total = noIts)
  for(i in 1:noIts){
    pb$tick()
    # Proposal
    thetaProp = RWMProposalMixed(thetaCurr, lambda, V, type)

    sims = replicate(noSims, simulator(P0,
                                       epiParam = thetaProp),
                     simplify = F)
    logPProp = sapply(X = sims, function(X) totalHospLikelihood(colSums(X$dailyNoCases), totalHospAdmissions,
                                                                thetaProp[5]))
    logPProp = log(mean(exp(logPProp)))
    logPostProp = logPProp + sum(evalPrior(thetaProp, priors))
    acceptProb = logPostProp - logPostCurr
    #print(acceptProb)

    u1 = runif(1)
    if(log(u1) < acceptProb){
      noAccept = noAccept + 1
      logPostCurr = logPostProp
      thetaCurr = thetaProp
    }

    draws[i+1, ] = c(thetaCurr, logPostCurr)
  }
  timeTaken = as.numeric(Sys.time()) - start
  # Calculate Effective Sample Sizes (and Per Second) and Acceptance Rates
  ESS <- min(coda::effectiveSize(draws[, 1:n]))
  ESS.sec <- ESS/timeTaken
  acceptRate <-  noAccept/noIts


  print(c("Accept Rate:", acceptRate))

  return(list(draws = draws, ESS.sec = ESS.sec, acceptRate = acceptRate))
}

