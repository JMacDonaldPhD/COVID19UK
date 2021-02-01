#' Priors


betaPrior = function(a, b, log = TRUE){
  bPrior = function(x){
    dbeta(x, a, b, log = log)
  }
}

unifPrior <- function(a = 0, b = 1, log = TRUE){
  uPrior = function(x){
    dunif(x, a, b, log)
  }
}

gammaPrior <- function(shape, rate, log = TRUE){
  gPrior = function(x){
    dgamma(x, shape, rate, log = log)
  }
}

fixedPrior <- function(trueValue, log = TRUE){

  if(log){
    fPrior <- function(x){
      if(x == trueValue){
        0
      } else{
        -Inf
      }
    }
  } else{
    fPrior <- function(x){
      if(x == trueValue){
        1
      } else{
        0
      }
    }
  }
}

evalPrior = function(xs, priors){
  eval = mapply(FUN = function(prior, x) prior(x), prior = priors, x = xs)

  return(eval)
}
