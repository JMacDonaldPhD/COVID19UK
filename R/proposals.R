#' Proposals

unitTransform = function(x){
  while(x < 0 | x > 1){
    if(x < 0){
      x = abs(x)
    }
    if(x > 1){
      x = 2 - x
    }
  }
  return(x)
}


RWMProposalUnit = function(alphaCurr, lambda){
  alphaProp = alphaCurr + lambda*rnorm(1, 0, 1)
  while(alphaProp < 0 | alphaProp > 1){
    if(alphaProp < 0){
      alphaProp = abs(alphaProp)
    }
    if(alphaProp > 1){
      alphaProp = 2 - alphaProp
    }
  }
  return(alphaProp)
}



RWMProposalptve = function(thetaCurr, lambda){
  thetaProp = abs(thetaCurr + lambda*rnorm(1, 0, 1))
  return(thetaProp)
}


# type 1: Regular RWM proposal
# type 2: Positive folded RWM proposal
# type 3: Unit folded RWM proposal
RWMProposalMixed = function(thetaCurr, lambda, V = diag(1, length(thetaCurr)), type = rep(1, length(thetaCurr))){

  type2 = (type == 2)
  type3 = (type == 3)

  thetaProp = thetaCurr + lambda*mvtnorm::rmvnorm(1, mean = rep(0, length(thetaCurr)), sigma = V)

  thetaProp[type2] = abs(thetaProp[type2])

  thetaProp[type3] = sapply(thetaProp[type3], FUN = unitTransform, simplify = "array")

  return(simplify2array(thetaProp))
}


# Create covariance matrix which zeros parameters that are considered fixed

proposal.VarCovMat0 <- function(noParam, whichVar){

  V <- diag(1, noParam)

  V[-whichVar, ] <- 0
  V[, -whichVar] <- 0

  return(V)
}





