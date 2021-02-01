# Modular MCMC


MCMC <- function(sampleData, theta0, priors, proposalType, simulator, likelihood, lambda, V, simParam, obsParam, noSims, noIts, burnIn,
                 runParallel = TRUE){

  start = as.numeric(Sys.time())
  # Adaptive Step


  k = length(theta0)
  s = length(simParam)
  o = length(obsParam)
  if(k != s + o){
    stop("Length of Parameter vector does not match number of utility arguments given.")
  } else if(!identical(unique(c(simParam, obsParam)),
                       sort(c(simParam, obsParam)))){
    stop("More than one purpose assigned to a parameter.")
  }


  thetaCurr = theta0

  # Likelihood Function
  logPostCurr = -Inf

  while(logPostCurr == -Inf){

    if(runParallel){
      cl <- parallel::makePSOCKcluster(parallel::detectCores() - 1)
      doParallel::registerDoParallel(cl)
      sims <-
        foreach(1:noSims) %dopar% {
          simulator(thetaCurr[simParam])
        }
      parallel::stopCluster(cl)

    } else{
      sims = replicate(noSims, simulator(thetaCurr[simParam]), simplify = F)
    }

    logPCurr = sapply(X = sims, function(X) likelihood(X, sampleData, thetaCurr[obsParam]))
    logPCurr = log(mean(exp(logPCurr)))
    logPostCurr = logPCurr + sum(evalPrior(thetaCurr, priors))
  }

  draws = matrix(nrow = noIts + 1, ncol = length(thetaCurr) + 1)
  draws[1, ] = c(thetaCurr, logPostCurr)
  noAccept = 0

  print("Sampling Progress")
  pb <- progress::progress_bar$new(total = noIts)
  for(i in 1:noIts){
    pb$tick()
    # Proposal
    thetaProp = RWMProposalMixed(thetaCurr, lambda, V, proposalType)

    if(runParallel){
      cl <- parallel::makePSOCKcluster(parallel::detectCores() - 1)
      doParallel::registerDoParallel(cl)
      sims <-
        foreach(1:noSims) %dopar% {
          simulator(param = thetaProp[simParam])
        }
      parallel::stopCluster(cl)
    } else{
      sims <- replicate(noSims, simulator(param = thetaProp[simParam]), simplify = F)
    }


    logPProp <- sapply(X = sims, function(X) likelihood(X, sampleData, thetaProp[obsParam]))
    logPProp <- log(mean(exp(logPProp)))
    logPostProp <- logPProp + sum(evalPrior(thetaProp, priors))
    acceptProb <- logPostProp - logPostCurr
    #print(acceptProb)
    u1 <- runif(1)
    if(log(u1) < acceptProb){
      noAccept <- noAccept + 1
      logPostCurr <- logPostProp
      thetaCurr <- thetaProp
    }
    draws[i+1, ] <- c(thetaCurr, logPostCurr)
  }
  timeTaken <- as.numeric(Sys.time()) - start
  # Calculate Effective Sample Sizes (and Per Second) and Acceptance Rates
  ESS <- coda::effectiveSize(draws[, 1])
  ESS.sec <- ESS/timeTaken
  acceptRate <-  noAccept/noIts
  print(c("Accept Rate:", acceptRate))

  return(draws)

}


adaptiveMCMC <- function(sampleData, theta0, priors, proposalType, simulator, likelihood, lambda0, V0, simParam,
                         obsParam, noSims, noIts, burnIn,
                         runParallel = TRUE, delta = 0.05){

  start <- as.numeric(Sys.time())
  # Adaptive Step

  #print("Hey, I got inside the function!")
  k <- length(theta0)
  s <- length(simParam)
  o <- length(obsParam)
  if(k != s + o){
    stop("Length of Parameter vector does not match number of utility arguments given.")
  } else if(!identical(unique(c(simParam, obsParam)),
                       sort(c(simParam, obsParam)))){
    stop("More than one purpose assigned to a parameter.")
  }

  varParam <- which(!apply(X = V0, 2, function(X) all(X == 0)))
  thetaCurr = theta0
  n <- length(thetaCurr)

  # Likelihood Function
  logPostCurr <- -Inf

  while(logPostCurr == -Inf){
    #print("I'm Stuck!")
    if(runParallel){
      cl <- parallel::makePSOCKcluster(parallel::detectCores() - 1)
      doParallel::registerDoParallel(cl)
      sims <-
        foreach(1:noSims) %dopar% {
          simulator(thetaCurr[simParam])
        }
      parallel::stopCluster(cl)

    } else{
      sims <- replicate(noSims, simulator(thetaCurr[simParam]), simplify = F)
    }

    logPCurr <- sapply(X = sims, function(X) likelihood(X, sampleData, thetaCurr[obsParam]))
    logPCurr <- log(mean(exp(logPCurr)))
    logPostCurr <- logPCurr + sum(evalPrior(thetaCurr, priors))
  }


  lambda <- lambda0
  V <- V0
  draws <- matrix(nrow = noIts + 1, ncol = length(thetaCurr) + 2 + length(as.vector(V[varParam, varParam])))
  draws[1, ] <- c(thetaCurr, logPostCurr, lambda, as.vector(V[varParam, varParam]))
  noAccept <- 0

  print("Sampling Progress")
  pb <- progress::progress_bar$new(total = noIts)
  for(i in 1:noIts){
    pb$tick()
    # Proposal

    u.delta <- runif(1)
    if(u.delta > delta & noAccept >= 10){
      V <- var(draws[,1:n], na.rm = T)
      thetaProp <- RWMProposalMixed(thetaCurr, lambda, V, proposalType)
    } else{
      thetaProp <- RWMProposalMixed(thetaCurr, lambda0, V0, proposalType)
    }


    if(runParallel){
      cl <- parallel::makePSOCKcluster(parallel::detectCores() - 1)
      doParallel::registerDoParallel(cl)
      sims <-
        foreach(1:noSims) %dopar% {
          simulator(param = thetaProp[simParam])
        }
      parallel::stopCluster(cl)
    } else{
      sims <- replicate(noSims, simulator(param = thetaProp[simParam]), simplify = F)
    }


    logPProp <- sapply(X = sims, function(X) likelihood(X, sampleData, thetaProp[obsParam]))
    logPProp <- log(mean(exp(logPProp)))
    logPostProp <- logPProp + sum(evalPrior(thetaProp, priors))
    acceptProb <- logPostProp - logPostCurr
    #print(acceptProb)

    u1 <- runif(1)
    if(log(u1) < acceptProb){
      noAccept <- noAccept + 1
      logPostCurr <- logPostProp
      thetaCurr <- thetaProp
      if(u.delta > delta){
        lambda <- lambda + 0.93*(lambda/((i)^(1/3)))
      }
    } else{
      if(u.delta > delta){
        lambda <- lambda - 0.07*(lambda/((i)^(1/3)))
      }
    }
    draws[i+1, ] <- c(thetaCurr, logPostCurr, lambda, as.vector(V[varParam, varParam]))
  }
  timeTaken <- as.numeric(Sys.time()) - start
  # Calculate Effective Sample Sizes (and Per Second) and Acceptance Rates
  ESS <- coda::effectiveSize(draws[, 1])
  ESS.sec <- ESS/timeTaken
  acceptRate <-  noAccept/noIts
  print(c("Accept Rate:", acceptRate))

  if(!exists("V")){
    V <- var(draws[,1:n], na.rm = T)
  }

  return(list(draws = draws, V = V, lambda = lambda))

}



