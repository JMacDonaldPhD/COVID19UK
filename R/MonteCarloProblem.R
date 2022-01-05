#' Monte Carlo Problem

#' Define a Monte Carlo problem, in turn creating function which
#' enables the user to carry out crude Monte Carlo estimation or
#' define an importance density and carry out Importance Sampling.
#' An "optimal" importance density for a given family of distributions
#' can be found and passed to the Importance Sampling function.

#' @param H function of underlying random variables which is an unbiased estimate of the target value
#' @param f The nominal probability distribution which the random variables are beleived to originate from.
#'        List of 2, the first of which is sampling function, and the second which is a density function.
#'        f should already have any parameters it relies on defined inside it, unless a non-centered representation is
#'        chosen. In this case, the non-centering transformation should be provided, along with the relevant
#'        parameters.
#' @param non_centered function which defines a non-centering parameterisation. NULL if non-centered parameterisation
#'        is not used. Should take two arguments. One particle and the parameters.
#' @param par The parameters to be used if a non_centered parameterisation is defined.

MonteCarloProblem <- function(H, f, centre_fn = identity, par){

  if(!identical(centre_fn, identity)){
    centre_fn_part <- function(particle){
      return(centre_fn(particle, par))
    }
  }
  #' @param N Number of particles to be drawn from f
  #' @return  Returns a crude Monte Carlo estimate
  #'          with estimate variance and relative error.
  crudeSampling <- function(N){
    particles <- replicate(N, f$sample(), simplify = F)

    centred_particles <- lapply(X = particles, FUN = centre_fn_part)
    H_calc <- sapply(centred_particles, FUN = H)
    MCest <- mean(H_calc)

    MCest.var <- (1/(N-1))*sum((MCest - H_calc)^2)
    MCest.RE <- sqrt(MCest.var/N)/MCest
    return(list(MCest = MCest, MCest.var = MCest.var, MCest.RE = MCest.RE))
  }

  #' @param N Number of particles to be drawn from g
  #' @param g Distribution comprising of sampling and
  #'          density calculation functions. Should
  #'          generate adn calculate densities for
  #'          particles of the same form associated
  #'          with f.
  #' @return  Returns an Importance Sampling Monte Carlo
  #'          estimate based on the choice of g, along
  #'          with estimate variance and relative error.
  importanceSampling <- function(N, g){

    particles <- replicate(N, g$sample(), simplify = F)

    centred_particles <- lapply(X = particles, FUN = centre_fn_part)


    H_calc <- sapply(centred_particles, FUN = H)

    g_calc <- sapply(particles, FUN = g$density)

    f_calc <- sapply(particles, FUN = f$density)

    W <- f_calc/g_calc
    Y <- H_calc*W


    MCest <- mean(Y)

    MCest.var <- (1/(N-1))*sum((MCest - Y)^2)
    MCest.RE <- sqrt(MCest.var/N)/MCest
    return(list(MCest = MCest, MCest.var = MCest.var, MCest.RE = MCest.RE))
  }

  #' @param N Number of particles drawn from f
  #' @param g_family The family of distributions to be optimised over
  #' @param StartPar Starting parameters for optim. Should be in the
  #'        correct form for use with g_family (i.e correct length and
  #'        valid values)
  #' @return Returns a distribution comprising of sampling and density
  #'         calculation functions for use in Importance Sampling.
  ISoptim <- function(N, g_family, startPar){

    particles <- replicate(N, f$sample(), simplify = F)
    centred_particles <- lapply(particles, FUN = centre_fn_part)

    f_calc <- sapply(particles, FUN  = f$density)

    H_calc <- sapply(centred_particles, FUN = H)


    approxVariance <- function(eta){
      g_calc <- sapply(particles, FUN = g_family(eta)$density)
      W <- f_calc/g_calc
      V_hat <- (1/N)*sum((H_calc^2)*W)
      return(V_hat)
    }
    par <- optim(par = startPar, fn = approxVariance)$par

    return(g_family(par))
  }

  return(list(IS = importanceSampling, crudeMC = crudeSampling,
              ISoptim = ISoptim))
}
