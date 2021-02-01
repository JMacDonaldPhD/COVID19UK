#
# SMC - Covid 19
#

#
#' Discrete time daily model.
#


discreteModel = function(N, kernel, mx, daySim = "unconditional", y = NULL, noDays = 1, outputStates = F){
  n = length(N)
  dailyCases = matrix(nrow = n, ncol = noDays)
  if(outputStates){
    tot.states = sum(mx) + 2
    simOut = rep(list(matrix(nrow = n, ncol = noDays)), tot.states)

    ma=c(1,mx[1]+1,sum(mx[1:2])+1,sum(mx)+1)

    names = c("S", rep(NA, sum(mx)), "R")

    stateNames = c("E", "P", "I")
    for(i in 1:length(mx)){
      stageIndex = (ma[i] + 1):ma[i+1]
      for(j in 1:length(stageIndex)){
        names[stageIndex[j]] = paste(c(stateNames[i], j), sep = "", collapse = "")
      }
      names(simOut) = names
    }
    if(daySim == "unconditional"){
      simModel = function(StateX0, lambda, theta, kernelParam){
        StateX = StateX0
        K = kernel(kernelParam)
        day = 1
        while(day <= noDays){
          simDay = transit_spat_unc(StateX, mx, lambda, theta, K, output = output)
          StateX = simDay$StateX

          for(i in 1:tot.states){
            simOut[[i]][,day] = StateX[,i]
          }

          dailyCases[,day] = simDay$noCases
          day = day + 1
        }
        return(list(noCases = dailyCases, StateX = simOut))
      }
    } else if(daySim == "conditional"){

      if(is.null(y)) stop("Provide case data")

      simModel = function(StateX0, lambda, theta, kernelParam){
        StateX = StateX0
        K = kernel(kernelParam)
        day = 1
        while(day <= noDays){
          simDay = transit_spat_cond(y[,day], StateX, mx, lambda, theta, K, output = output)
          StateX = simDay$StateX

          for(i in 1:tot.states){
            simOut[[i]][,day] = StateX[,i]
          }

          dailyCases[,day] = simDay$noCases
          day = day + 1
        }
        return(list(noCases = dailyCases, StateX = simOut))
      }

    } else{
      stop("Choose valid simulation option; conditional/unconditional")
    }
  } else{
    if(daySim == "unconditional"){
      simModel = function(StateX0, lambda, theta, kernelParam){
        StateX = StateX0
        K = kernel(kernelParam)
        day = 1
        while(day <= noDays){
          simDay = transit_spat_unc(StateX, mx, lambda, theta, K, output = output)
          StateX = simDay$StateX
          dailyCases[,day] = simDay$noCases
          day = day + 1
        }
        return(dailyCases)
      }
    } else if(daySim == "conditional"){

      if(is.null(y)) stop("Provide case data")

      simModel = function(StateX0, lambda, theta, kernelParam){
        StateX = StateX0
        K = kernel(kernelParam)
        day = 1
        while(day <= noDays){
          print(day)
          #simDay = transit_spat_cond1(y[,day], StateX, mx, lambda, theta, K, output = output)
          simDay = transit_spat_cond2(y[,c(day, day + 1)], StateX, mx, lambda, theta, K, output = output)
          StateX = simDay$StateX
          dailyCases[,day] = simDay$noCases
          day = day + 1
        }
        return(dailyCases)
      }

    } else{
      stop("Choose valid simulation option; conditional/unconditional")
    }
  }


  simModel
}


#' Noise Generating Function

NoiseGeneratingFunction = function(N, StateX0, lambda, theta, kernelParam,
                                   obsParam, simulator=NA){
  realisation = simulator(StateX0, lambda, theta, kernelParam)
  n = length(N)

  # # Simulated Cumulative Cases (X)
  # X = rowSums(realisation)

  # Noise (y)
  y = apply(X = realisation, MARGIN = 2, FUN = function(X) rbinom(n, size = X, prob = obsParam))
  rownames(y) = names(N)
  y
}






# transit_spat - simulates forward number of cases for a day allowing for N and K being vectors and matrices of dependence
# obs - states how many people becoming infected are observed
#
#
#
# NN=c(2,2,2)*10000000
# KK=matrix(c(1,0.2,0.2,0.2,1,0.2,0.2,0.2,1),ncol=3)
# SSX=matrix(sample(1:10,21,replace=T),nrow=3,ncol=7)
# SSX[,1]=NN
#
# for(i in 1:10)
# {
# SSX=transit_spat(NN,KK,SSX, kappa = 0.4, gamma = 1, delta = 0.7, lambda_P = 2, lambda_I = 1, spark = 0.5, phi = 1)
# print(SSX)
# SSX=SSX$StateX
#
# points()
#
#
# }

