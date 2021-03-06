# Household SIR Experiment

householdSIRexperiment <- function(experimentName, simulator, obsModel,
                                   simParam, obsParam, dir, varSymbols, seed = NULL){

  # Create Directory for Experiments
  # If a directory already exists with name "dir", then the directory will not be created
  # but the wd will be set to dir
  if(missing(dir)){
    dir = "Experiments/"
  }

  if(!dir.exists(dir)){
    dir.create(dir)
  }
  setwd(dir)

  # Create Directory for this experiment
  if(!dir.exists(experimentName)){
    dir.create(experimentName)
  }
  setwd(experimentName)




  # ==== Simulate Positive test data ====

  if(!is.null(seed)){
    set.seed(seed)
  }

  res <- simulator(param = simParam)

  jpeg(paste0(c(experimentName, ".jpeg"), collapse = ""), width = 480*2)
  par(mfrow = c(1, 2))
  plot(res$SIRsummary[,1], col = "blue", type = 'l', ylim = c(0, N), xlab = "Time (Days)", ylab = "# Individuals")
  lines(res$SIRsummary[,2], col = "red", type = 'l', ylim = c(0, N))
  lines(res$SIRsummary[,3], col = "grey", type = 'l', ylim = c(0, N))

  finalHHsize <- factor(res$householdFinalSize, levels = 0:4, labels = as.character(0:4))
  plot(finalHHsize, xlab = "Household Epidemic Size", ylab = "Frequency")
  dev.off()
  print(c("FINAL SIZE:", tail(res$SIRsummary[,2] + res$SIRsummary[,3], n  = 1)))
  sampleData <- obsModel$NGF(res, obsParam)
  #totalPositiveCases <- colSums(dailyPositiveCases)
  likelihood <- obsModel$llh

  # ==== Epidemic Summary Output

  allParam <- c(simParam, obsParam)
  p <- length(allParam)
  simParamIndex <- 1:length(simParam)
  obsParamIndex <- (1:length(obsParam)) + length(simParam)


  likelihoodCalc_function <- function(simParam.part, obsParam, noSims = 1, calc.sd = F){
    llh.est <- function(llhParam){
      simParam <- c(llhParam, simParam.part)
      simData <- replicate(noSims, simulator(simParam), simplify = F)
      llh <- simplify2array(lapply(X = simData, function(X) likelihood(X, sampleData, obsParam)))

      llh.mean <- log(mean(exp(llh)))

      if(calc.sd){
        llh.sd <- log(sd(exp(llh)))
        return(list(mean = llh.mean, sd = sd.mean))
      }
      return(llh.mean)
    }
  }


  # ==== MCMC Functions ====
  MCMC_functions <- function(priors, proposalType, varParam, runParallel = TRUE){

    # ==== Convert Parameter Names into Indices
    paramNames <- c(names(simParam), names(obsParam))

    for(i in 1:length(varParam)){
      varParam[i] <- which(paramNames == varParam[i])
    }
    varParam <- as.numeric(varParam)
    # ==== Assign fixed priors to fixed parameters

    finalPriors <- rep(list(NA), p)
    finalPriors[varParam] <- priors

    fixedParam <- (1:p)[-varParam]
    for(i in fixedParam){
      finalPriors[[i]] <- fixedPrior(allParam[i])
      print(finalPriors[[i]](allParam[i])) # This fixes problem with assign prior closures
    }

    finalProposalType <- rep(NA, p)
    finalProposalType[varParam] <- proposalType
    finalProposalType[fixedParam] <- 1


    graphOutput <- function(MCMCdraws, outputName, density = TRUE){
      #varParam <- (1:length(theta0))[-fixedParam]


      jpeg(outputName)

      if(density){
        ncol.graphs <- 3
      } else{
        ncol.graphs <- 2
      }

      par(mfrow = c(length(varParam), ncol.graphs))
      if(density){
        for(i in varParam){
          plot(MCMCdraws[, i], type = 'l', ylab = varSymbols[i])
          acf(MCMCdraws[, i], main = '')
          plot(density(MCMCdraws[, i]), col = 'blue', xlab = varSymbols[i])
          abline(h = 0, v = allParam[i], col = "red", lty = 2)
        }

      } else{
        for(i in varParam){
          plot(MCMCdraws[, i], type = 'l', ylab = varSymbols[i])
          acf(MCMCdraws[, i], main = '')

        }
      }
      dev.off()
    }

    # profileMCMC <- function(lambda0 = 1e-3, V0 = proposal.VarCovMat0(length(theta0)),
    #                         noSims, runParallel = F){
    #   # Estimate time
    #   tmp <- "TestRunProfile.out"
    #   Rprof(tmp, interval = 0.01)
    #   noIts <- 100
    #   run <-
    #     adaptiveMCMC(sampleData, theta0, priors, proposalType, simulator, likelihood, lambda0, V0, simParamIndex, obsParamIndex,
    #                  noSims, noIts = noIts, burnIn = 0, runParallel, delta = 0.05)
    #
    #   Rprof(NULL)
    #
    #   capture.output(summaryRprof(tmp), file = "TestRunProfile.txt")
    #
    #   print("==== Estimated Time for 10000 iterations ====")
    #   print(paste0(10000*summaryRprof(tmp)$sampling.time/noIts/60, " minutes"))
    # }

    adaptMCMC <- function(theta0, lambda0 = 1e-3, V0 = proposal.VarCovMat0(p, varParam), noSims,
                          noIts, burnIn){

      run <-
        adaptiveMCMC(sampleData, theta0, priors = finalPriors, proposalType = finalProposalType, simulator,
                     likelihood, lambda0, V0,
                     simParamIndex, obsParamIndex,
                     noSims, noIts, burnIn, runParallel, delta = 0.05)


      save(run, file = "adaptRun.Rdata")

      # MCMC Diagnostic
      graphOutput(run$draws[-(1:burnIn), ], outputName = "AdaptMCMCOutput.jpeg", density = TRUE)
      return(list(lambda = run$lambda, V = run$V, draws = run$draws))
    }

    inferenceMCMC <- function(theta, lambda, V, noSims, noIts, burnIn){

      run <-
        MCMC(sampleData, theta, priors = finalPriors, proposalType = finalProposalType, simulator, likelihood,
             lambda, V, simParamIndex, obsParamIndex,
             noSims, noIts, burnIn, runParallel)
      save(run, file = "inferenceRun.Rdata")

      # MCMC Diagnostic
      graphOutput(run[-(1:burnIn), ], outputName = "inferenceMCMCOutput.jpeg", density = TRUE)
      return(run)
    }

    return(list(adaptMCMC = adaptMCMC,
                inferenceMCMC = inferenceMCMC))
  }
  print("Experiment set up. Do not change directory! Files will not save where you may expect them to")

  return(list(MCMC_functions = MCMC_functions, likelihoodCalc_function = likelihoodCalc_function))
}

