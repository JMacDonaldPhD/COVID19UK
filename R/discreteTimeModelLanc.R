

# Simulates one day in the future.

discreteModelLanc = function(N, kernel, stoch = matrix(c(-1, 1, 0, 0, 0, -1, 1, 0, 0, 0, -1, 1), byrow = T, ncol = 4, nrow = 3),
                             oneObs = TRUE, lockdownDate = lubridate::dmy('23/03/2020'), startDate = lubridate::dmy('01/03/2020'),
                             endDate = lubridate::today(), kernelPreLD, kernelPostLD){

  n <- length(N)
  noDays = as.numeric(difftime(endDate, startDate, unit = 'days'))
  lockdownDay = as.numeric(difftime(lockdownDate, startDate, unit = 'days'))
  # projects epidemic one day using binomial draws
  # Either takes P0 as an argument for the intial number of presymptomatic individuals
  # OR StateX if simulating from a timepoint which is mid-epidemic.


  # epiParam = c(beta.base, beta70plus, betaI, rho, alpha, gamma)
  # beta.base, infection rate between those not in over 70 category
  # beta70plus, infection rate for people in 70 plus category
  # betaI, reduction in infection rate due to quarantine or hospitalisation
  # rho, probability of moving from presymtomatic to symptomatic in one day
  # alpha, probability of being hospitalised given you have become symptomatic
  # gamma, probability of dying or becoming immune to COVID19 in one day

  daily_household_prog <- function(StateX = NA, P0, beta_k){


  }

  daily_pop_prog <- function(StateX = NA, P0, epiParam, K_P, K_I){

    # if(is.na(StateX)){
    #   if(missing(P0)){
    #     stop("Must provide StateX or P0")
    #   } else{
    #     S <- N - P0
    #     I <- R <- rep(0, n)
    #     stateX <- cbind(S, I, P0, R)
    #   }
    # }

    # beta = rep(epiParam[1], n)
    # beta[n] = epiParam[2]
    # beta.I = epiParam[3]
    # K = kernel(epiParam[1], beta[n])
    rho = rep(epiParam[1], n)
    alpha = c(epiParam[2], rep(epiParam[4], n - 2), epiParam[3])
    gamma = rep(epiParam[5], n)



    # Total Infectious Pressure from pre-symptomatic people
    Lambda_P= K_P%*%(StateX[,2])


    # Total Infectious pressure from infectious individuals on susceptible people
    Lambda_I= K_I%*%(StateX[,3])


    # Total Infectious Pressure
    Lam=(Lambda_I+Lambda_P)/N

    trans <- ink <- matrix(0,ncol=4,nrow=n)

    # Binomial draw for number of infectious in each population
    trans[,1] = rbinom(n,StateX[,1],1-exp(-Lam))

    # Binomial Draw for number of numbers exits from E_1, E_2, P_1, P_2, I_1, I_2
    # etr=matrix(rbinom(n*6, StateX[,2:7],
    #                   c(rep(1 - exp(kappa) ,n*2),rep(1 - exp(gamma), n*2), rep(1 - exp(delta),n*2))), ncol=6, byrow=F)

    # Pre-symptomatic -> Symptomatic
    PtoI = rbinom(n, StateX[,2], rho)

    # Hospital Admissions (Positive Testing)
    hospAdmission = rbinom(n, PtoI, alpha)

    # Leaving States
    trans[,2] = PtoI
    trans[,3] = rbinom(n, StateX[,3], gamma)
    # Entering New State
    ink[,2:4] = trans[,1:3]

    # New epidemic state
    StateX = StateX + ink - trans

    # # How many people becoming symptomatic are asertained
    # obs = rbinom(n,trans[,5], phi)

    # Return state and observation
    return(list(StateX = StateX, noCases = PtoI, hospAdmission = hospAdmission))
  }

  # Return S, P, I, R and Hospital Admissions
  # S, P, I, R is a simulation of the underlying epidemic process
  # Hospital Admissions is a draw from the observation process given S, P, I, R, the total number of people who were admitted
  # to hospital.
  sim = function(P0, epiParam, LDparam){

    # Pre-Lockdown infectious tranmission matrices
    K_P = kernelPreLD(epiParam[1])
    K_I = kernelPreLD(epiParam[2])

    S <- P <- I <- R <- matrix(nrow = noDays + 1, ncol = n)

    S0 = N - P0
    I0 <- R0 <- rep(0, n)

    S[1,] <- S0
    P[1,] <- P0
    I[1,] <- I0
    R[1,] <- R0

    stateX <- cbind(S0, P0, I0, R0)
    dailyHospAdmissions <- dailyNoCases <- matrix(nrow = noDays, ncol = n)

    for(day in 1:noDays){
      if(day == lockdownDay + 1){
        # Pre-Lockdown infectious tranmission matrices
        K_P = kernelPostLD(epiParam[1], LDparam)
        K_I = kernelPostLD(epiParam[2], LDparam)
      }

      nextDay <- transit_spat(StateX = stateX, epiParam = epiParam)
      stateX <- nextDay$StateX

      dailyHospAdmissions[day, ] <- nextDay$hospAdmission
      dailyNoCases[day, ] <- nextDay$noCases

      S[day + 1,] <- stateX[,1]
      P[day + 1,] <- stateX[,2]
      I[day + 1,] <- stateX[,3]
      R[day + 1,] <- stateX[,4]
    }

    return(list(S = S, P = P, I = I, R = R, dailyNoCases = dailyNoCases,
                dailyHospAdmissions = dailyHospAdmissions))
  }
}



# projects epidemic one day using binomial draws
transit_spat=function(StateX, epiParam){

  beta = rep(epiParam[1], n)
  beta[n] = epiParam[2]
  beta.I = epiParam[3]
  K = kernel(epiParam[1], beta[n])
  rho = epiParam[4]
  alpha = epiParam[5]
  gamma = epiParam[6]


  # Total Infectious Pressure from pre-symptomatic people
  Lambda_P= K%*%(StateX[,2])


  # Total Infectious pressure from infectious individuals on susceptible people
  Lambda_I= beta.I*K%*%(StateX[,3])


  # Total Infectious Pressure
  Lam=(Lambda_I+Lambda_P)/N

  trans=ink=matrix(0,ncol=3,nrow=n)

  # Binomial draw for number of infectious in each population
  trans[,1] = rbinom(n,StateX[,1],1-exp(-Lam))

  # Binomial Draw for number of numbers exits from E_1, E_2, P_1, P_2, I_1, I_2
  # etr=matrix(rbinom(n*6, StateX[,2:7],
  #                   c(rep(1 - exp(kappa) ,n*2),rep(1 - exp(gamma), n*2), rep(1 - exp(delta),n*2))), ncol=6, byrow=F)

  # Pre-symptomatic -> Symptomatic
  PtoI = rbinom(n, State[,2], rep(rho, n))

  # Hospital Admissions
  hospAdmission = rbinom(n, PtoI, rep(alpha, n))

  # Leaving States
  trans[,2] = PtoI
  trans[,3] = rbinom(n, StateX[,3], rep(gamma, n))
  # Entering New State
  ink[,2:3] = trans[,1:2]

  # New epidemic state
  StateX = StateX + ink - trans

  # # How many people becoming symptomatic are asertained
  # obs = rbinom(n,trans[,5], phi)

  # Return state and observation
  return(list(StateX = StateX, noCases = PtoI, hospAdmission = hospAdmission))
}
