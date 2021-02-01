# discrete time COVID19 model with local and global mixing steps.


# pop will be a dataframe containing information about all individuals in the population of interest
# In particular, it will include household ID, Age group, Key Worker status, Current Infectious status.


discreteModel_HG = function(pop, stoch = matrix(c(0, 0, 0, 0,
                                                  -1, 1, 0, 0,
                                                  0, -1, 1, 0,
                                                  0, 0, -1, 1), byrow = T, ncol = 4, nrow = 4),
                             oneObs = TRUE, lockdownDate = lubridate::dmy('23/03/2020'), startDate = lubridate::dmy('01/03/2020'),
                             endDate = lubridate::today(), kernelPreLD, kernelPostLD){


  noDays <- as.numeric(difftime(endDate, startDate, unit = 'days'))
  lockdownDay <- as.numeric(difftime(lockdownDate, startDate, unit = 'days'))
  # projects epidemic one day using binomial draws
  # Either takes P0 as an argument for the intial number of presymptomatic individuals
  # OR StateX if simulating from a timepoint which is mid-epidemic.

  # ==== Set up (Only needs to be done once per population) ====

  hh_sizes <- as.numeric(table(pop$hh_id))

  n_hh = length(unique(pop$hh_id))
  N = nrow(pop)

  group <- as.numeric(pop$ageBin)
  n <- length(unique(group))


  hh_block_mat <- lapply(X = hh_sizes, FUN = function(X){
    mat <- matrix(1, ncol = X, nrow = X)
    diag(mat) <- 0
    mat
  })
  hh_mat <- Matrix::bdiag(hh_block_mat)


  # Individual Group Status
  indGroup = matrix(F, nrow = N, ncol = n)
  for(i in 1:N){
    indGroup[i, group[i]] = T
  }



  susc = c(T, F, F, F)
  pre = c(F, T, F, F)
  inf = c(F, F, T, F)
  rem = c(F, F, F, T)
  states = matrix(c(susc, pre, inf, rem), nrow = 4, ncol = 4, byrow = F)
  noStates = ncol(states)

  indState = states[pop$state, ]

  # Construct Population State Matrix using infectious status and groups
  StateX = matrix(nrow = ncol(indGroup), ncol = nrow(states))

  for(i in 1:ncol(indGroup)){
    for(j in 1:nrow(states)){
      StateX[i, j] = sum(indState[, j] & indGroup[, i])
    }
  }




  # epiParam = c(beta.base, beta70plus, betaI, rho, alpha, gamma)
  # beta.base, infection rate between those not in over 70 category
  # beta70plus, infection rate for people in 70 plus category
  # betaI, reduction in infection rate due to quarantine or hospitalisation
  # rho, probability of moving from presymtomatic to symptomatic in one day
  # alpha, probability of being hospitalised given you have become symptomatic
  # gamma, probability of dying or becoming immune to COVID19 in one day

  dailyHouseholdProg <- function(StateX = NA, indState, group, beta_H){
    #print("==== Household Contacts ====")
    # Calculates the number of infected individuals exerting pressure on each susceptible individual
    infPressure <- beta_H*textTinyR::sparse_Sums(hh_mat[indState[,1], indState[,2] | indState[,3], drop = F], rowSums = T)

    # Draws whether susceptibles become infected
    infections <- rbinom(sum(indState[,1]), size = 1, prob = 1 - exp(-infPressure)) + 1


    # Add infections to population level state information.
    whoInfected = infections == 2
    if(sum(whoInfected) > 0){
      group[indState[,1]][infections == 2]
      x = table(group[indState[,1]][infections == 2])
      StateX[as.numeric(names(x)), ] <- StateX[as.numeric(names(x)), ] + t(sapply(x, FUN = `*`, e2 = stoch[2,]))

      # Change States at individual level
      if(length(whoInfected) == 1){
        indState[indState[,1], ] <- matrix(states[2, ], nrow = sum(whoInfected), ncol = noStates, byrow = T)
      } else{
        indState[indState[,1], ][whoInfected, ] <- matrix(states[2, ], nrow = sum(whoInfected), ncol = noStates, byrow = T)
      }
    }

    # StateCheck = matrix(nrow = ncol(indGroup), ncol = nrow(states))

    # for(i in 1:ncol(indGroup)){
    #   for(j in 1:nrow(states)){
    #     StateCheck[i, j] = sum(indState[, j] & indGroup[, i])
    #   }
    # }
    # if(!sum(StateCheck == StateX) == noStates * n){
    #   print("Individual States and Population States no longer match")
    # }


    return(list(StateX = StateX, indState = indState))

  }

  dailyGlobalProg <- function(StateX = NA, indState, epiParam, K_P, K_I){
    ##print("==== Global Contacts ====")
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
    #K = kernel(epiParam[1], beta[n])

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
    postiveTests = rbinom(n, PtoI, alpha)

    # Leaving States
    trans[,2] = PtoI
    trans[,3] = rbinom(n, StateX[,3], gamma)


    #print(trans)
    #print(StateX)
    symptoms = c()
    infected = c()
    recover = c()

    fn <- function(){

      for(i in 1:n){
        x1 = (1:N)[indGroup[,i] & indState[,2]]
        l1 = length(x1)
        if(l1 > 0){
          #print(i)
          #print(x1)
          symptoms <<- c(symptoms, sample(x = c(0, x1), size = trans[i, 2], prob = c(0, rep(1, length(x1))), replace = FALSE))
        }

        x2 = (1:N)[indGroup[,i] & indState[,1]]
        l2 = length(x2)
        if(l2 > 0){
          infected <<- c(infected, sample(x = c(0, x2), size = trans[i, 1], prob = c(0, rep(1, length(x2))), replace = FALSE))
        }

        x3 = (1:N)[indGroup[,i] & indState[,3]]
        l3 = length(x3)
        if(l3 > 0){
          recover <<- c(recover, sample(x = c(0, x3), size = trans[i, 3], prob = c(0, rep(1, length(x3))), replace = FALSE))
        }
      }
    }

    fn()

    # Entering New State
    ink[,2:4] = trans[,1:3]

    # New epidemic state
    StateX = StateX + ink - trans

    indState[symptoms, ] <- matrix(states[3, ], nrow = length(symptoms), ncol = noStates, byrow = T)

    indState[infected, ] <- matrix(states[2, ], nrow = length(infected), ncol = noStates, byrow = T)
    indState[recover, ] <- matrix(states[4, ], nrow = length(recover), ncol = noStates, byrow = T)



    # StateCheck = matrix(nrow = ncol(indGroup), ncol = nrow(states))
    #
    # for(i in 1:ncol(indGroup)){
    #   for(j in 1:nrow(states)){
    #     StateCheck[i, j] = sum(indState[, j] & indGroup[, i])
    #   }
    # }
    # if(!sum(StateCheck == StateX) == noStates * n){
    #   print("Individual States and Population States no longer match")
    # }
    # # How many people becoming symptomatic are asertained
    # obs = rbinom(n,trans[,5], phi)

    # Return state and observation
    return(list(StateX = StateX, indState = indState, noCases = PtoI, postiveTests = postiveTests))
  }

  # Return S, P, I, R and Hospital Admissions
  # S, P, I, R is a simulation of the underlying epidemic process
  # Hospital Admissions is a draw from the observation process given S, P, I, R, the total number of people who were admitted
  # to hospital.
  sim = function(beta_H, epiParam, beta_G, l){

    # Pre-Lockdown infectious tranmission matrices
    K_P = kernelPreLD(beta_G)
    K_I = kernelPreLD(beta_G*0.1)

    S <- P <- I <- R <- matrix(nrow = noDays + 1, ncol = n)

    S0 <- StateX[,1]
    P0 <- StateX[,2]
    I0 <- StateX[,3]
    R0 <- StateX[,4]

    S[1,] <- S0
    P[1,] <- P0
    I[1,] <- I0
    R[1,] <- R0

    dailyPositiveTests <- dailyNoCases <- matrix(nrow = noDays, ncol = n)

    for(day in 1:noDays){
      #print(c("======== day", day))
      if(day == lockdownDay + 1){
        # Pre-Lockdown infectious tranmission matrices
        K_P = kernelPostLD(beta_G, l)
        K_I = kernelPostLD(beta_G*0.1, l)
      }


      # ==== Household contact ====

      # if(day == 18){
      #   debug(dailyHouseholdProg)
      # }
      #debug(dailyHouseholdProg)
      householdContact <- dailyHouseholdProg(StateX, indState, group, beta_H)
      StateX <- householdContact$StateX
      indState <- householdContact$indState

      # ==== Global Contact ====

      #debug(dailyGlobalProg)
      globalContact <- dailyGlobalProg(StateX = StateX, indState = indState, epiParam = epiParam, K_P, K_I)
      StateX <- globalContact$StateX
      indState <- globalContact$indState

      dailyPositiveTests[day, ] <- globalContact$postiveTests
      dailyNoCases[day, ] <- globalContact$noCases

      S[day + 1,] <- StateX[,1]
      P[day + 1,] <- StateX[,2]
      I[day + 1,] <- StateX[,3]
      R[day + 1,] <- StateX[,4]
    }

    return(list(S = S, P = P, I = I, R = R, dailyNoCases = dailyNoCases,
                dailyPositiveTests = dailyPositiveTests))
  }
  return(sim)
}



# projects epidemic one day using binomial draws
# #transit_spat=function(StateX, epiParam){
#
#   beta = rep(epiParam[1], n)
#   beta[n] = epiParam[2]
#   beta.I = epiParam[3]
#   K = kernel(epiParam[1], beta[n])
#   rho = epiParam[4]
#   alpha = epiParam[5]
#   gamma = epiParam[6]
#
#
#   # Total Infectious Pressure from pre-symptomatic people
#   Lambda_P= K%*%(StateX[,2])
#
#
#   # Total Infectious pressure from infectious individuals on susceptible people
#   Lambda_I= beta.I*K%*%(StateX[,3])
#
#
#   # Total Infectious Pressure
#   Lam=(Lambda_I+Lambda_P)/N
#
#   trans=ink=matrix(0,ncol=3,nrow=n)
#
#   # Binomial draw for number of infectious in each population
#   trans[,1] = rbinom(n,StateX[,1],1-exp(-Lam))
#
#   # Binomial Draw for number of numbers exits from E_1, E_2, P_1, P_2, I_1, I_2
#   # etr=matrix(rbinom(n*6, StateX[,2:7],
#   #                   c(rep(1 - exp(kappa) ,n*2),rep(1 - exp(gamma), n*2), rep(1 - exp(delta),n*2))), ncol=6, byrow=F)
#
#   # Pre-symptomatic -> Symptomatic
#   PtoI = rbinom(n, State[,2], rep(rho, n))
#
#   # Hospital Admissions
#   hospAdmission = rbinom(n, PtoI, rep(alpha, n))
#
#   # Leaving States
#   trans[,2] = PtoI
#   trans[,3] = rbinom(n, StateX[,3], rep(gamma, n))
#   # Entering New State
#   ink[,2:3] = trans[,1:2]
#
#   # New epidemic state
#   StateX = StateX + ink - trans
#
#   # # How many people becoming symptomatic are asertained
#   # obs = rbinom(n,trans[,5], phi)
#
#   # Return state and observation
#   return(list(StateX = StateX, noCases = PtoI, hospAdmission = hospAdmission))
# }
