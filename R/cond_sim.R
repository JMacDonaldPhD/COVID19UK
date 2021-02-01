#' One-day forward simulation of SEPIR meta-popn epidemic conditional on the number of detected cases
#' in each region (I -> R event is a case)


#' @func
#' @param y Number of ascertained cases on that day
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

# Conditions on number of ascertained cases in that day, y[day]
# Fails if there is not enough individuals in the pre-symptomatic state on a particular day.
transit_spat_cond1 = function(y, StateX, mx, lambda, theta, K, output = "cases"){

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
  Lam= (lambda_S + lambda_I + lambda_P)/N

  trans = ink = matrix(0, ncol= ma[4] + 1, nrow = n)


  # Remove ascertained cases so number of unobserved is simulated
  StateX[,ma[3]] = StateX[, ma[3]] - y

  if(any(StateX[, ma[3]] < 0)){
    stop("Insufficent Number of Pre-symptomatics to support observed data")
  }


  # Binomial draw for number of infectious in each population
  trans[,1] = rbinom(n, StateX[,1], 1 - exp(-Lam))

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
  StateX[,ma[3]+1]=StateX[,ma[3]+1] + y
  return(list(noCases = trans[,ma[3]] + y, StateX = StateX))

}

# Conditions on the number of new cases using the number of ascertained cases for that day, y[day]

# Makes sure there are enough pre-symptomatic people for the number of ascertained cases in the following day

transit_spat_cond2 = function(y, StateX, mx, lambda, theta, K, output = "cases"){

  # Total theta states
  tot.mx = sum(mx)

  # check StateX matches upto total states

  #
  today = 1
  tomorrow = 2

  # column indices for the first stage of each state
  mf = c(2, mx[1] + 2,  sum(mx[1:2]) + 2)

  # column indices for the last stage of each state (excluding susceptible).
  ma = c(mx[1]+1, sum(mx[1:2])+1, tot.mx+1)

  # Total Infectious Pressure from pre-symp.
  lambda_P = lambda[1]*K%*%(rowSums(StateX[,mf[2]:ma[2], drop = F]))

  # Total Infectious pressure from infectious individuals on susceptible people
  lambda_I=lambda[2]*rowSums(StateX[,mf[3]:ma[3], drop = F])

  # Spark (background transmisson)
  lambda_S = lambda[3]*N/mean(N)

  # Total Infectious Pressure
  Lam= (lambda_S + lambda_I + lambda_P)/N

  trans=ink=matrix(0,ncol= ma[3] + 1,nrow=n)

  # Remove ascertained cases so number of unobserved is simulated (last P -> First I)
  StateX[,ma[2]] = StateX[, ma[2]] - y[today]

  if(any(StateX[,ma[2]] < 0)){
    stop("Not enough Pre-symptomatics to support observed cases on this day")
  }

  # How many unascertained cases leave P
  unascertained = rbinom(n, StateX[,ma[2]], rep(theta[2], n))

  # Do some populations currentlly have a lack of Pre-symptomatics for the next day?
  lackOfP = StateX[, ma[2]] - unascertained - y[tomorrow] < 0

  # Lowest amount of movement from state previous of last P state that gives sufficient pre-symptomatics
  # for next day
  bound = lackOfP*abs(StateX[, ma[2]] - unascertained - y[tomorrow])

  # Is the previous stage an exposed or pre-symptomatic state
  prevStageE = (ma[2] - 1) %in% mf[1]:ma[1]

  # Condition on the sufficient pre-symptomatics we need for the next day
  if(any(lackOfP > 0)){
    if(all(StateX[, ma[2] - 1] >= bound)){
      StateX[, ma[2] - 1] = StateX[, ma[2] - 1] - bound
    } else{
      stop("Not enough Stage 2 exposed individuals to progress to next day")
    }
  }

  # Binomial draw for number of infectious in each population
  trans[,1] = rbinom(n, StateX[,1], 1 - exp(-Lam))

  # Binomial Draw for number of numbers exits from E_1, E_2, P_1, P_2, I_1, I_2

  Edraws = matrix(rbinom(n*mx[1], StateX[, mf[1]:ma[1]], rep(theta[1], n*mx[1])), ncol = mx[1], nrow = n, byrow = F)

  Pdraws = matrix(rbinom(n*(mx[2] - 1), StateX[, mf[2]:(ma[2] - 1)], rep(theta[2], n*(mx[2] - 1))),
                  ncol = mx[2] - 1, nrow = n, byrow = F)
  Idraws = matrix(rbinom(n*mx[3], StateX[, mf[3]:ma[3]], rep(theta[3], n*mx[3])), ncol = mx[3], nrow = n, byrow = F)


  # etr = matrix(rbinom(n*(tot.mx - 1), StateX[,c(mf[1]:ma[1], (mf[2]:ma[2])[-mx[2]], mf[3]:ma[3])],
  #                     c(rep(theta[1] ,n*mx[1]),rep(theta[2], n*(mx[2] - 1)), rep(theta[3],n*mx[3]))), ncol= tot.mx - 1, byrow=F)


  # Leaving States
  trans[, mf[1]:ma[3]] = cbind(Edraws, Pdraws, unascertained, Idraws)

  # Entering New State
  ink[, mf[1]:ma[3]] = trans[, mf[1]:ma[3] - 1]

  # New epidemic state
  StateX = StateX + ink - trans
  StateX[, ma[2] + 1] = StateX[, ma[2] + 1] + y[today]
  StateX[, ma[2]] = StateX[, ma[2]] + bound
  return(list(noCases = trans[,ma[2]] + y[today] , StateX = StateX))
}
