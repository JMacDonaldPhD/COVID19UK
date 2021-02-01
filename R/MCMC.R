#' MCMC Sampler

#'@param caseData, number of daily ascertained cases of 2019-nCoV by city

# Markov Discrete Time Approx parameter prior
logGammaPrior = function(rate, shape){
  dgammaPrior = function(x){
    dgamma(x, rate = rate, shape = shape, log = T)
  }
}

#' Geometric Event parameter prior
logUnifPrior = function(){
  duinfPrior = function(x){
    dunif(x, log = TRUE)
  }
}

# Proposal for Geometric Event Parameters
FoldedNormProposal = function(paramCurr, lambda, V, unit){
  p = length(paramCurr)
  if(missing(V)){
    V = diag(1, p)
  }
  paramProp = paramCurr + lambda*mvtnorm::rmvnorm(1, rep(0, p), V)

  for(i in (1:p)[-unit]){
    paramProp[i] = abs(paramProp[i])
  }
  for(i in unit){
    while(paramProp[i] < 0 | paramProp[i] > 1){
      if(paramProp[i] < 0){
        paramProp[i] = abs(paramProp[i])
      }
      if(paramProp[i] > 1){
        paramProp[i] = 2 - paramProp[i]
      }
    }
  }
  paramProp
}

# ==== Bayesian Inference ====

BayesianModel = function(y, simulations, posiPrior, unitPrior, param0, lambda01, V0, lambda02, delta = 0.05,
                         stateX0, spark, kernelParam){

  Start = as.numeric(Sys.time())
  # Adaptive Step
  noIts = 10000

  paramCurr = param0[1:5]
  phiCurr = param0[6]
  # Likelihood Function
  logPostCurr = -Inf

  unit = 3:5
  posi = 1:2

  while(logPostCurr == -Inf){
    likelihoodCurr = simulations(stateX0, c(paramCurr, spark), kernelParam)
    logPCurr = likelihoodCurr(phiCurr)
    logPostCurr = logPCurr + sum(unitPrior(paramCurr[unit])) + sum(posiPrior(paramCurr[posi])) + unitPrior(phiCurr)
  }

  # Create Storage Matrix
  draws = matrix(NA, nrow = noIts + 1, ncol = length(param0) + 1)
  draws[1,] = c(paramCurr, phiCurr, logPostCurr)

  lambda1  = lambda01
  lambda2  = lambda02
  accParam = 0
  accPhi = 0

  print("Sampling Progress")
  pb <- progress::progress_bar$new(total = noIts)

  for(i in 1:noIts){
    pb$tick()

    #print(paramCurr)
    # Proposal Acceptance Counter
    udelta1 = runif(1, 0, 1)
    if(udelta1 > delta & accParam > 10){
      Vi = var(draws[,1:5], na.rm = T)
      paramProp = FoldedNormProposal(paramCurr, lambda1, Vi, unit)
    } else{
      paramProp = FoldedNormProposal(paramCurr, lambda01, V0, unit)
    }
    #print(paramProp)
    likelihoodProp = simulations(stateX0, c(paramProp, spark), kernelParam)
    logLProp = likelihoodProp(phiCurr)
    logPostProp = logLProp + sum(posiPrior(paramProp[posi])) + sum(unitPrior(paramProp[unit])) + unitPrior(phiCurr)

    accProb = (logPostProp) - (logPostCurr)

    u1 = log(runif(1))
    reject = !(u1 < accProb)
    if(u1 < accProb){
      logPostCurr = logPostProp
      paramCurr = paramProp
      likelihoodCurr = likelihoodProp
      accParam = accParam + 1
      if(udelta1 > delta & accParam > 10){
        lambda1 = lambda1 + 0.93*(lambda1/sqrt(i))
      }
    } else{
      if(udelta1 > delta & accParam > 10){
        lambda1 = lambda1 - 0.07*(lambda1/sqrt((i)))
      }
    }

    # Multiplicative RWM Proposal for phi
    udelta2 = runif(1)
    if(udelta2 > delta){
      phiProp = FoldedNormProposal(phiCurr, lambda2, unit = 1)
    } else{
      phiProp = FoldedNormProposal(phiCurr, lambda02, unit = 1)
    }

    logLProp = likelihoodCurr(phiProp)
    logPostProp = logLProp + sum(posiPrior(paramCurr[posi])) + sum(unitPrior(paramCurr[unit])) + unitPrior(phiProp)

    accProb = (logPostProp) - (logPostCurr)

    u1 = log(runif(1))
    if(u1 < accProb){
      logPostCurr = logPostProp
      phiCurr = phiProp
      accPhi = accPhi + 1
      if(udelta2 > delta){
        lambda2 = lambda2 + 0.75*(lambda2/sqrt(i))
      }
    } else{
      if(udelta2 > delta){
        lambda2 = lambda2 - 0.25*(lambda2/sqrt((i)))
      }
    }
    draws[i + 1, ] = c(paramCurr, phiCurr, logPostCurr)
    print(accParam)
  }

  par(mfrow = c(4, 2))
  plot(draws[, 1], type = 'l')
  acf(draws[,1])

  plot(draws[, 2], type = 'l')
  acf(draws[,2])

  plot(draws[, 3], type = 'l')
  acf(draws[,3])

  plot(draws[, 4], type = 'l')
  acf(draws[,4])

  V = var(draws[, 1:2], na.rm = T)


  # Posterior (MCMC Sampler)
  MCMC = function(param0, noIts){
    Start = as.numeric(Sys.time())
    paramCurr = param0

    # Likelihood Function
    logPCurr = -Inf

    while(logPostCurr == -Inf){
      likelihoodCurr = simulations(paramCurr[1:3])
      logLCurr = likelihoodCurr(paramCurr[4])
      logPostCurr = logLCurr + sum(epiPrior1(paramCurr[1:2])) + epiPrior2(paramCurr[3]) +
        phiPrior(paramCurr[4])
    }

    draws = matrix(NA, nrow = noIts + 1, ncol = length(paramCurr) + 1)
    draws[1,] = c(paramCurr, logPostCurr)
    accParam = 0
    accPhi = 0

    print("Sampling Progress")
    pb <- progress::progress_bar$new(total = noIts)

    for(i in 1:noIts){
      pb$tick()
      # Folded Normal Proposal for beta, gamma, I0W
      paramProp = c(abs(paramCurr[1:2] + mvtnorm::rmvnorm(1, mean = rep(0, 2), sigma = lambda1*V)),
                    paramCurr[3:4])

      paramProp[3] = paramCurr[3] + sample(c(1, -1), size = 1)

      likelihoodProp = simulations(paramProp[1:3])
      logLProp = likelihoodProp(paramProp[4])
      logPostProp = logLProp + sum(epiPrior1(paramProp[1:2])) + epiPrior2(paramProp[3]) +
        phiPrior(paramProp[4])
      accProb = (logPostProp) - (logPostCurr)

      u1 = log(runif(1))
      if(u1 < accProb){
        logPostCurr = logPostProp
        paramCurr = paramProp
        likelihoodCurr = likelihoodProp
        accParam = accParam + 1
      }

      # Multiplicative RWM Proposal for phi
      paramProp = c(paramCurr[1:3], unitFoldedNormProposal(paramCurr[4], lambda2))
      logLProp = likelihoodCurr(paramProp[4])
      logPostProp = logLProp + sum(epiPrior1(paramProp[1:2])) + epiPrior2(paramProp[3]) +
        phiPrior(paramProp[4])

      accProb = (logPostProp) - (logPostCurr)

      u1 = log(runif(1))
      if(u1 < accProb){
        logPostCurr = logPostProp
        paramCurr = paramProp
        accPhi = accPhi + 1
      }
      draws[i + 1, ] = c(paramCurr, logPostCurr)
    }
    list(draws = draws, accRateParam = accParam/noIts, accRatePhi = accPhi/noIts, time = as.numeric(Sys.time()) - Start)
  }
  list(MCMCSampler = MCMC, draws = draws, lambda1, lambda2, V,
       time = as.numeric(Sys.time()) - Start)
}


