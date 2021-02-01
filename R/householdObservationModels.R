
dHyperGeom = function(x, Y, k, log = TRUE){
  if(!(is.list(x) & is.list(Y))){
    return(extraDistr::dmvhyper(x, Y, k = sum(x), log = log))
  } else if(is.list(x) & is.list(Y)){
    return(sum(mapply(extraDistr::dmvhyper, x, Y, MoreArgs = list(k, log))))
  }
}
# Hypergeometric Sampling of households
#' @param completeData The size of epidemic in all of the households on the days of observation
#' @param sampleData The size of epidemic in sampled households
HG_likelihood <- function(completeData, sampleData, obsParam){
  CD <- table(factor(completeData, levels = 0:4))
  SD <- table(factor(sampleData, levels = 0:4))
  k <- sum(SD)
  llh <- dHyperGeom(x = as.numeric(SD), Y = as.numeric(CD), k = obsParam, log = T)
  return(llh)
}

householdHyperNGF <- function(completeData, sampleSize){
  return(completeData[sample(1:length(completeData), size = sampleSize, replace = F)])
}



Biased_likelihood <- function(completeData, sampleData, obsParam){

  # = CHECK 1 =
  # In observed households, is the number of infected individuals the same in each dataset
  check1 <- sum((completeData$Hstate$I[-1, ] - sampleData$observedInfectious) * sampleData$observed_HH) == 0
  if(!check1){
    return(-Inf)
  }

  # = CHECK 2 =
  # Is the number of cases at least the number of observed cases on any given day? The following calculation
  # will return -Inf if this is not the case
  llh <- sum(dbinom(sampleData$observedCases, size = completeData$Infections, prob = obsParam, log = TRUE))

  return(llh)
}

householdBiasedNGF <- function(completeData, alpha){
  n_trials <- length(completeData$Infections)
  observedCases <- matrix(rbinom(n_trials, size = completeData$Infections, prob = alpha),
                          nrow = nrow(completeData$Infections), ncol = ncol(completeData$Infections))

  observed_HH <- observedCases > 0

  observedInfectious <- observed_HH * completeData$Hstate$I[-1, ]

  sampleData <- list(observedCases = observedCases, observed_HH = observed_HH,
                     observedInfectious = observedInfectious)

  return(sampleData)

  # Table which records which houses were observed on which day, how many cases where
  # ascertained and how many individuals are infected
  #print(sum(observedCases != 0))
  caseTable <- matrix(nrow = sum(observedCases != 0), ncol = 4)
  colnames(caseTable) <- c("Day", "Household", "Observed Cases", "Infected Residents")
  caseNo <- 1
  for(day in 1:nrow(observedCases)){
    obs_hh_index <- which(observedCases[day, ] > 0)
    if(!identical(obs_hh_index, integer(0))){
      for(j in 1:length(obs_hh_index)){
        hh_index <- obs_hh_index[j]
        hh_no_inf <- completeData$Hstate$I[day + 1, hh_index]
        caseTable[caseNo, ] <- c(day, hh_index, observedCases[day, hh_index], hh_no_inf)
        caseNo <- caseNo + 1
        #print(caseNo)
      }
    }
  }
  return(list(caseTable = caseTable, observedCases = observedCases))
}




observationModel <- function(householdSize){

  # Hypergeometric
  # Observe a proportion of households
  # For households of size 4, there are 5 possible states regarding size of epidemic in the household
  # at any given time.
  # Any cross-sectional sample can therefore be seen as a hypergeometric sample from these 5 different
  # states

  possibleStates <- diag(TRUE, householdSize + 1)



  biased_likelihood <- function(sampleData, epiParam, obsParam){
    # Like the COVID model, the biased model will only follow-up households who aquire infection, but with
    # a probability
    # The likelihood of follow-up will be proportional to the number of people who become infected in the
    # household.
    # This means the observation process is informed by the epidemic process itself. If infections spreads
    # between households, more households will be followed up. If the within household infection is higher,
    # the likelihood that households are followed up increases




  }

}
