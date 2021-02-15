#' Epidemic Model

#' Takes an underlying epidemic model, observation model and a set of parameters. Generates noise based on these parameters,
#' creates a closure function which then attempts to make inference on the noise through a Psuedo-Marginal MCMC.
#' This closure will take arguements related to the PM-MCMC, including which parameters are fixed and which are
#' to be made inference on.


#' @param simulator function which simulates the underlying epidemic process, given a set of parameters
#' @param obsModel A list of two functions.One named `NGF` (Noise Generating Function), which generates data,
#'                 a realisation of the underlying epidemic process (i.e output from `simulator`) and a set of
#'                 observation parameters. The other, `llh`, calculates the likelihood of a given underlying
#'                 epidemic process resulting in a given sample, given a set of observational parameters.
#' @param simParam "True" values of the underlying epidemic process
#' @param obsParam "True" values of the observation model
#' @param varSymbols A vector of strings/expressions allocating names for the model parameters.
#' @param seed Gives a seed to set, which will determine the underlying epidemic process generated, and the
#'             sample generated, given a set of model parameters. If NULL, no seed is set and random number
#'             generation will be based on current R session.
#' @return A likelihood estimation function and MCMC set up function
epiModel <- function(simulator, obsModel,
                     simParam, obsParam, varNames, seed = NULL, conditional = TRUE){
  #print(ls(environment(fun = simulator)))
  # Create Directory for Experiments
  # If a directory already exists with name "dir", then the directory will not be created
  # but the wd will be set to dir
  # if(missing(dir)){
  #   dir = "Experiments/"
  # }
  #
  # if(!dir.exists(dir)){
  #   dir.create(dir)
  # }
  # setwd(dir)
  #
  # # Create Directory for this experiment
  # if(!dir.exists(experimentName)){
  #   dir.create(experimentName)
  # }
  # setwd(experimentName)
  #
  #


  # ==== Simulate Positive test data ====

  if(!is.null(seed)){
    set.seed(seed)
  }

  res <- simulator(param = simParam)

  #print(c("FINAL SIZE:", sum(res)))
  sampleData <- obsModel$NGF(res, obsParam)
  #totalPositiveCases <- colSums(dailyPositiveCases)
  likelihood <- obsModel$llh


  # Conditional Simulator
  if(conditional){
    simulator <- conditionalTestingHouseholdSIR(sampleData, pop, endTime = environment(fun = simulator)$endTime)
  }



  # ==== Epidemic Summary Output ====


  # ==== MCMC Functions ====

  allParam <- c(simParam, obsParam)
  p <- length(allParam)
  simParamIndex <- 1:length(simParam)
  obsParamIndex <- (1:length(obsParam)) + length(simParam)


  #' Creates a function which calculates a Monte Carlo estimate for a given subset of parameters for the underlying epidemic process.

  #' @param simParam_names Character vector describing which simulation parameters are to be profiled (according to the given value of)
  #' @param obsParam A set of observation parameters to condition the likilhood calculation on
  #' @param noSims Number of particles to simulate for Monte Carlo estimation of the likelihood
  #' @param calc.sd If TRUE, the standard deviation of the estimate will also be calculated.
  likelihoodCalc_function <- function(simParam_names, obsParam, noSims = 1, calc.sd = F){
    simParam.subset <- which(varNames %in% simParam_names)
    simParam.part <- simParam[-simParam.subset]
    llh.est <- function(llhParam){
      simParam[simParam.subset] <- llhParam
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


  #'

  #' @param priors
  #' @param proposalType
  #' @param varParam
  #' @param runParallel
  #' @return An adaptive MCMC and MCMC sampler

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









