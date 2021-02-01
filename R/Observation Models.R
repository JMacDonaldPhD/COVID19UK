# Observation Models


# More abstracted function for creating the observation model.
# Takes complete daily cases and simulates observation data.
# Parameters of the observation process
# ascertainment rate of positive cases, alpha
# What is the structure of aplha? is it group based or homogeneous between groups
# period of accumulation of positive cases.

# A test is done at a certain time. If it is positive, it will be logged in the observation data
# Because of the nature of the inference (psusedo-marginal), it is impossible to calculate likelihood
# based on exact test times (i/e the time that the positive case was ascertained). Instead, we make
# inference of the cumulation of positive test over a chosen time window. The size/frequency of these
# time windows will dictate how much informatino we have about hte evolution about the epidemic, as
# well as the ascertainment rate. Psuedo-Marginal methods work extremely well when there is less data,
# however, if there is too little data, we cannot identify the model parameters. The question is, therefire,
# how much information do we need to make adequate inference on the model parameters. In addition, is it
# possible to do Psuedo-Marginal inference to make adequate inference or what techniques, i.e particle
# filtering, do we need to employ to make this possible.


# alphaStruct will specify how alpha is assumed to be assigned to the different groups
# If NULL, it is assumed that one alpha parameter pertains to the positive case ascertianment
# in all groups

# obsWindows a list of vectors (size two) that will specify the time periods which positive test data is accumulated. List?
observationModel <- function(alphaStruct = NULL, obsWindows = "Total"){
  OWT <- identical(obsWindows, "Total")

  # Determine whether alpha is homogeneous across all groups
  if(is.list(obsWindows)){
    lapply(X = obsWindows, FUN = function(X) if(length(X) != 2) stop("Observation windows must be defined by vector
                                                                     of length 2"))
    n_obs <- length(obsWindows)
  } else if(!OWT){
    stop("Invalid argument given to ObsWindows")
  }


  # Noise Generating Function
  if(is.null(alphaStruct)){
    # ==== Total cases, no alpha struct ====
    if(OWT){
      NGF <- function(completeData, alpha){
        n_alpha <- length(alpha)

        if(n_alpha != 1){
          stop("Number of alpha parameters is > 1, but no group assignment structure is given!")
        }

        cumulativeCases <- colSums(completeData)

        sampleData <- rbinom(length(cumulativeCases), cumulativeCases, prob = alpha)

        return(sampleData)
      }
    } else{
      NGF <- function(completeData, alpha){
        n_alpha <- length(alpha)

        if(n_alpha != 1){
          stop("Number of alpha parameters is > 1, but no group assignment structure is given!")
        }

        cumulativeCases <-  t(sapply(X = obsWindows, function(X) colSums(completeData[(X[1] + 1):X[2], ])))

        sampleData <- matrix(t(rbinom(length(cumulativeCases), t(cumulativeCases), prob = alpha)),
                             nrow = n_obs, byrow = TRUE)

        return(sampleData)
      }
    }
  } else{

    if(OWT){
      NGF <- function(completeData, alpha){
        n_alpha <- length(alpha)
        k_alpha <- length(unique(alphaStruct))
        n <- ncol(completeData)

        if(n_alpha != k_alpha){
          stop("Number of alpha parameters, does not match the number of unique group assignments!")
        } else{
          alphaValues <- matrix(rep(alpha[alphaStruct], n), ncol = n, byrow = TRUE)
          #alphaValues <- alpha[alphaStruct]
        }

        cumulativeCases <- colSums(completeData)



        #sampleData <- t(apply(X = cumulativeCases, MARGIN = 1, FUN = rbinom, n = n, prob = alphaValues))
        sampleData <- rbinom(length(cumulativeCases), cumulativeCases, probs = alphaValues)

        return(sampleData)
      }
    } else{
      NGF <- function(completeData, alpha){
        n_alpha <- length(alpha)
        k_alpha <- length(unique(alphaStruct))
        n <- ncol(completeData)

        if(n_alpha != k_alpha){
          stop("Number of alpha parameters, does not match the number of unique group assignments!")
        } else{
          alphaValues <- matrix(rep(alpha[alphaStruct], n), ncol = n, byrow = TRUE)
          #alphaValues <- alpha[alphaStruct]
        }

        cumulativeCases <- t(sapply(X = obsWindows, function(X) colSums(completeData[(X[1] + 1):X[2], ])))

        #sampleData <- t(apply(X = cumulativeCases, MARGIN = 1, FUN = rbinom, n = n, prob = alphaValues))
        sampleData <- matrix(t(rbinom(length(cumulativeCases), t(cumulativeCases), probs = t(alphaValues))),
                             nrow = n_obs, byrow = TRUE)

        return(sampleData)
      }
    }
  }

  # Likelihood Calculation
  if(is.null(alphaStruct)){
    if(OWT){
      likelihoodFunction <- function(simulatedData, sampleData, alpha, log = TRUE){
        n_alpha <- length(alpha)

        if(n_alpha != 1){
          stop("Number of alpha parameters is > 1, but no group assignment structure is given!")
        }

        cumulativeCases <- colSums(simulatedData)

        llh <- sum(dbinom(sampleData, cumulativeCases, alpha, log = log))

        return(llh)
      }
    } else{
      likelihoodFunction <- function(simulatedData, sampleData, alpha, log = TRUE){
        n_alpha <- length(alpha)

        if(n_alpha != 1){
          stop("Number of alpha parameters is > 1, but no group assignment structure is given!")
        }

        cumulativeCases <- t(sapply(X = obsWindows, function(X) colSums(simulatedData[(X[1] + 1):X[2], ])))

        llh <- sum(dbinom(sampleData, cumulativeCases, alpha, log = log))

        return(llh)
      }
    }

  } else{
    if(OWT){
      likelihoodFunction <- function(simulatedData, sampleData, alpha, log = TRUE){

        n_alpha <- length(alpha)
        k_alpha <- length(unique(alphaStruct))
        n <- ncol(completeData)

        if(n_alpha != k_alpha){
          stop("Number of alpha parameters, does not match the number of unique group assignments!")
        } else{
          alphaValues <- matrix(rep(alpha[alphaStruct], n), ncol = n, byrow = TRUE)
          #alphaValues <- alpha[alphaStruct]
        }

        cumulativeCases <- colSums(simulatedData)

        llh <- sum(dbinom(sampleData, cumulativeCases, alphaValues, log = log))
        return(llh)

      }
    } else{
      likelihoodFunction <- function(simulatedData, sampleData, alpha, log = TRUE){

        n_alpha <- length(alpha)
        k_alpha <- length(unique(alphaStruct))
        n <- ncol(completeData)

        if(n_alpha != k_alpha){
          stop("Number of alpha parameters, does not match the number of unique group assignments!")
        } else{
          alphaValues <- matrix(rep(alpha[alphaStruct], n), ncol = n, byrow = TRUE)
          #alphaValues <- alpha[alphaStruct]
        }

        cumulativeCases <- t(sapply(X = obsWindows, function(X) colSums(simulatedData[(X[1] + 1):X[2], ])))

        llh <- sum(dbinom(sampleData, cumulativeCases, alphaValues, log = log))
        return(llh)

      }
    }

  }
  return(list(NGF = NGF, llh = likelihoodFunction))
}


# Positive case rate is constant between groups
commonAlphaObsModel <- function(completeData, alpha, daily = TRUE){
  if(is.matrix(completeData)){
    n = ncol(noCases)
    noPositiveTest <- t(apply(X = noCases, MARGIN = 1, FUN = rbinom, n = n, prob = alpha))
  } else{
    noPositiveTest <- rbinom(length(completeData), completeData, prob = alpha)
  }
  return(noPositiveTest)
}


groupAlphaObsModel <- function(completeData, alpha, alphaGroups){
  if(is.matrix(completeData)){
    n = ncol(completeData)
    noPositiveTest <- t(apply(X = noCases, MARGIN = 1, FUN = rbinom, n = n, prob = alpha[alphaGroups]))
  } else{
    noPositiveTest <- rbinom(length(completeData), completeData, prob = alpha[alphaGroups])
  }

  return(noPositiveTest)
}
