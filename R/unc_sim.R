#' Unconditional one-day forward simulation of SEPIR meta-popn epidemic


#' @func
#' @param N population of UA
#' @param K matrix describing mixing of between populations
#' @param mx, number of stages for each theta state state i.e mx[1] will be the number of exposure states an individual must go
#'            through before reaching the pre-symptomatic state
#' @param stateX length(N) by 8 matrix representing current state of epidemics
#' @param kappa Exit rate from Exposure state 1 (E_1) and Exposure state 2 (E_2)
#' @param delta Exit rate from pre-symptomatic state 1 (P_1) and pre-symptomatic state 2 (P_2)
#' @param gamma Exit rate from Infectious state 1 (I_1) and Infectious state 2 (I_2)
#' @param lambda_P human-human infection rate between an (syptomatic) infectious and susceptible indiviudal
#' @param lambda_I base human-human infection rate between an (pre-syptomatic) infectious and susceptible indiviudal
#' @param spark Background infectious pressure
#' @param phi Ascertainment rate

# lambda_P,lambda_I, kappa, delta, gamma, spark

# lambdaParam
# thetaParam

transit_spat_unc = function(StateX, mx, lambda, theta, K, output = "cases"){

  # Total theta states
  tot.mx = sum(mx)

  # check StateX matches upto total states

  # column indices for the last stage of each state (excluding susceptible).
  ma=c(2,mx[1]+1,sum(mx[1:2])+1,tot.mx+1)

  # Total Infectious Pressure from pre-symp.
  lambda_P=lambda[1]*K%*%(rowSums(StateX[,(ma[2] + 1):ma[3], drop = F]))


  # Total Infectious pressure from infectious individuals on susceptible people
  lambda_I=lambda[2]*rowSums(StateX[,(ma[3] + 1):ma[4], drop = F])

  # Spark (background transmisson)
  lambda_S = lambda[3]*N/mean(N)

  # Total Infectious Pressure
  Lam=(lambda_S + lambda_I + lambda_P)/N

  trans=ink=matrix(0,ncol= ma[4] + 1,nrow=n)

  # Binomial draw for number of infectious in each population
  trans[,1]=rbinom(n, StateX[,1], 1 - exp(-Lam))

  # Binomial Draw for number of numbers exits from E_1, E_2, P_1, P_2, I_1, I_2
  # etr=matrix(rbinom(n*6, StateX[,2:7],
  #                   c(rep(1 - exp(kappa) ,n*2),rep(1 - exp(gamma), n*2), rep(1 - exp(delta),n*2))), ncol=6, byrow=F)
  etr = matrix(rbinom(n*tot.mx, StateX[,ma[1]:ma[4]],
                      c(rep(theta[1] ,n*mx[1]),rep(theta[2], n*mx[2]), rep(theta[3],n*mx[3]))), ncol= tot.mx, byrow=F)

  # Leaving States
  trans[,ma[1]:ma[4]] = etr
  # Entering New State
  ink[,ma[1]:ma[4]] = trans[,ma[1]:ma[4] - 1]

  # New epidemic state
  StateX = StateX + ink - trans

  # New cases are defined as those which transition from the pre-symptomatic state to the symptomatic stage.
  return(list(noCases = trans[,ma[3]], StateX = StateX))

}
