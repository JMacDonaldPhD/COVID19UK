#' Hospital Admission Likelihood



totalHospLikelihood_commonAlpha = function(completeData, sampleData, alpha, log = TRUE){
  if(!identical(dim(as.matrix(completeData)), dim(as.matrix(sampleData)))){
    stop("Likelihood Error: Complete dataset and sample dataset not compatible")
  }

  if(log){
    return(sum(dbinom(sampleData, completeData, prob = alpha, log = log)))
  } else{
    return(prod(dbinom(sampleData, completeData, prob = alpha, log = log)))
  }
}

groupAlphaLlh <- function(alphaGroup){
  llh_fn <- function(completeData, sampleData, alpha, log = TRUE){
    if(!identical(dim(as.matrix(completeData)), dim(as.matrix(sampleData)))){
      stop("Likelihood Error: Complete dataset and sample dataset not compatible")
    }
    if(log){
      return(sum(dbinom(sampleData, completeData, prob = alpha[alphaGroup], log = log)))
    } else{
      return(prod(dbinom(sampleData, completeData, prob = alpha[alphaGroup], log = log)))
    }
  }
}


dailyHospLikelihood = function(dailyCases, dailyHospAdmissions, alpha, log = TRUE){
  if(!identical(dim(dailyCases), dim(dailyHospAdmissions))){
    stop("Likelihood Error: Complete dataset and sample dataset not compatible")
  }

  if(log){
    return(sum(dbinom(dailyHospAdmissions, dailyCases, prob = alpha, log = log)))
  } else{
    return(prod(dbinom(dailyHospAdmissions, dailyCases, prob = alpha, log = log)))
  }
}
