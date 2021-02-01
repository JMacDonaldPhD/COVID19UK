# Seperate Epidemic Model and Observation Model

COVID19Model <- function(pop, startDate = lubridate::dmy('01/03/2020'),
                         endDate = lubridate::today(), lockdownDate = lubridate::dmy('23/03/2020'),
                         kernel, kernelParam, output = "daily"){


  #'
  #' S -> E -> P -> I  -> R
  #'            |
  #'             -> A  -> R
  #'

  stoch <- matrix(c(-1, 1, 0, 0, 0, 0,
                    0, -1, 1, 0, 0, 0,
                    0, 0, -1, 1, 0, 0,
                    0, 0, -1, 0, 1, 0,
                    0, 0, 0, -1, 0, 1,
                    0, 0, 0, 0, -1, 1), nrow = 6, ncol = 6, byrow = T)
  noDays <- as.numeric(difftime(endDate, startDate, unit = 'days'))
  lockdownDay <- as.numeric(difftime(lockdownDate, startDate, unit = 'days'))
  # projects epidemic one day using binomial draws
  # Either takes P0 as an argument for the intial number of presymptomatic individuals
  # OR StateX if simulating from a timepoint which is mid-epidemic.

  # ==== Set up (Only needs to be done once per population) ====

  hh_sizes <- as.numeric(table(pop$hh_id))

  n_hh <- length(unique(pop$hh_id))
  N <- nrow(pop)
  n <- max(pop$uniqueGroup)
  group <- factor(as.integer(pop$uniqueGroup), levels = 1:14,
                  labels = as.character(1:14))


  hh_block_mat <- lapply(X = hh_sizes, FUN = function(X){
    mat <- matrix(1, ncol = X, nrow = X)
    diag(mat) <- 0
    mat
  })
  hh_mat <- Matrix::bdiag(hh_block_mat)


  # Individual Group Status
  indGroup <- as.integer(pop$uniqueGroup)


  # individual infectious status
  susc <- c(T, F, F, F, F, F)
  exp <- c(F, T, F, F, F, F)
  pre <- c(F, F, T, F, F, F)
  inf <- c(F, F, F, T, F, F)
  asym <- c(F, F, F, F, T, F)
  rem <- c(F, F, F, F, F, T)

  states <- matrix(c(susc, exp, pre, inf, asym, rem), nrow = length(susc),
                   ncol = length(susc), byrow = F)
  noStates <- ncol(states)

  indState <- states[pop$state, ]

  # Construct Population State Matrix using infectious status and groups
  StateX <- matrix(nrow = n, ncol = noStates)

  for(i in 1:n){
    for(j in 1:noStates){
      StateX[i, j] <- sum(group[indState[, j]] == i)
    }
  }

  fast.tabulate <- function(draws, groups){
    a <- data.table::data.table(V1 = groups, V2 = draws)
    b <- a[, .N, by = list(V1, V2)]
    c <- tapply(b$N, list(b$V1, b$V2), function(x) sum(x))
    c[] <- sapply(c, function(x) replace(x, is.na(x), 0))
    return(c)
  }

  # epiParam = c(beta.base, beta70plus, betaI, rho, gamma)
  # beta.base, infection rate between those not in over 70 category
  # beta70plus, infection rate for people in 70 plus category
  # betaI, reduction in infection rate due to quarantine or hospitalisation
  # rho, probability of moving from presymtomatic to symptomatic in one day
  # gamma, probability of dying or becoming immune to COVID19 in one day

  dailyProg <- function(StateX, indState, beta_H, beta_G, epiParam, K_P, K_I){
    # Calculates the number of infected individuals exerting pressure on each susceptible individual
    infPressure_H <- beta_H*textTinyR::sparse_Sums(hh_mat[indState[,1], (indState[,3] | indState[,4]) | indState[,5],
                                                          drop = F], rowSums = T)

    # Total Infectious Pressure from pre-symptomatic and asymptomatic people
    Lambda_P <- K_P%*%(StateX[,3] + StateX[,5])


    # Total Infectious pressure from infectious individuals on susceptible people
    Lambda_I <- K_I%*%(StateX[,4])

    groupInfPressure <- beta_G*(Lambda_P + Lambda_I)/N

    # Should be vector of size noS (i.e number of susceptible individuals)
    infPressure_G <- groupInfPressure[indGroup][indState[,1]]

    noS <- sum(StateX[,1])
    # Decided whether each individual gets individual is infected (equivalent to a binomial draw with prob
    # equal to the probability of infection over one day, given the state of the population at the start
    # of the day.)
    infections <- runif(noS) < (1 - exp(-infPressure_G - infPressure_H))
    noInfection <- sum(infections)

    if(noInfection > 0){
      infectionSum <- fast.tabulate(infections,  group[indState[,1]])
      if(dim(infectionSum)[2] == 2){
        infectionSum <- infectionSum[,2]
      }
    } else{
      infectionSum <- rep(0, n)
    }

    eta <- epiParam[1] # Exposed -> Pre-Symptomatic
    rho <- epiParam[2] # Pre-symptomatic -> Symptomatic
    phi <- epiParam[3] # Pre-symptomatic -> Asymptomatic
    gamma <- epiParam[4] # A/symtomatic -> Removed



    # Exposed -> Pre-symptomatic
    infectiousRate <- eta
    noE <- sum(StateX[,2])
    infectious <- rbinom(noE, size = 1, infectiousRate)

    noInfectious <- sum(infectious == 1)

    if(noInfectious > 0){
      infectiousSum <- fast.tabulate(infectious,  group[indState[,2]])
      if(dim(infectiousSum)[2] == 2){
        infectiousSum <- infectiousSum[,2]
      }
    } else{
      infectiousSum <- rep(0, n)
    }

    # Pre-symptomatic -> Symptomatic/Asymptomatic
    symptomRate <- rho
    noP <- sum(StateX[,3])
    draws <- runif(noP)
    symptoms <- draws < (symptomRate)
    noSymptoms <- sum(symptoms)
    if(noSymptoms > 0){
      symptomSum <- fast.tabulate(symptoms,  group[indState[,3]])
      if(dim(symptomSum)[2] == 2){
        symptomSum <- symptomSum[,2]
      }
    } else{
      symptomSum <- rep(0, n)
    }

    noCases <- symptomSum

    asymptomRate <- phi
    asymptoms <- !symptoms & (draws < symptomRate + asymptomRate)
    noAsymptoms <- sum(asymptoms)
    if(noAsymptoms > 0){
      asymptomSum <- fast.tabulate(asymptoms,  group[indState[,3]])
      if(dim(asymptomSum)[2] == 2){
        asymptomSum <- asymptomSum[,2]
      }
    } else{
      asymptomSum <- rep(0, n)
    }

    # symptoms/asymptoms -> removal
    removalRate <- gamma

    noI <- sum(StateX[,4])
    noA <- sum(StateX[,5])
    finalState <- runif(noI)
    removalI <- finalState < gamma
    noRemovalI <- sum(removalI)
    if(noRemovalI > 0){
      removalISum <- fast.tabulate(removalI,
                                   group[indState[,4]])

      if(dim(removalISum)[2] == 2){
        removalISum <- removalISum[,2]
      }
    } else{
      removalISum <- rep(0, n)
    }


    noA <- sum(StateX[,5])
    finalState <- runif(noA)

    removalA <- finalState < gamma

    noRemovalA <- sum(removalA)
    if(noRemovalA > 0){
      removalASum <- fast.tabulate(removalA,
                                   group[indState[,5]])

      if(dim(removalASum)[2] == 2){
        removalASum <- removalASum[,2]
      }
    } else{
      removalASum <- rep(0, n)
    }

    # Updating individual states

    savedIndState <- indState

    if(noInfection > 0){
      if(noS == 1){
        indState[savedIndState[,1], ] <- states[2, ]
      } else{
        indState[savedIndState[,1], ][infections == 1, ] <- matrix(states[2, ],
                                                                   nrow = noInfection,
                                                                   ncol = noStates,
                                                                   byrow = T)
      }
    }

    if(noInfectious > 0){
      if(noE == 1){
        indState[savedIndState[,2], ] <- states[3, ]
      } else{
        indState[savedIndState[,2], ][infectious == 1, ] <- matrix(states[3, ],
                                                                   nrow = noInfectious,
                                                                   ncol = noStates,
                                                                   byrow = T)
      }
    }

    if(noSymptoms > 0){
      if(noP == 1){
        indState[savedIndState[,3], ] <- states[4, ]
      } else{
        indState[savedIndState[,3], ][symptoms, ] <- matrix(states[4, ],
                                                            nrow = noSymptoms,
                                                            ncol = noStates,
                                                            byrow = T)
      }
    }
    if(noAsymptoms > 0){
      if(noP == 1){
        indState[savedIndState[,3], ] <- states[5, ]
      } else{
        indState[savedIndState[,3], ][asymptoms, ] <- matrix(states[5, ],
                                                             nrow = noAsymptoms,
                                                             ncol = noStates,
                                                             byrow = T)
      }
    }

    if(noRemovalI > 0){
      if(noI == 1){
        indState[savedIndState[,4], ] <- states[6, ]
      } else{
        indState[savedIndState[,4], ][removalI, ] <- matrix(states[6, ],
                                                            nrow = noRemovalI,
                                                            ncol = noStates,
                                                            byrow = T)
      }
    }
    if(noRemovalA > 0){
      if(noA == 1){
        indState[savedIndState[,5], ] <- states[6, ]
      } else{
        indState[savedIndState[,5], ][removalA, ] <- matrix(states[6, ],
                                                            nrow = noRemovalA,
                                                            ncol = noStates,
                                                            byrow = T)
      }
    }



    StateX <- StateX + t(sapply(infectionSum, FUN = `*`, e2 = stoch[1,])) +
      t(sapply(infectiousSum, FUN = `*`, e2 = stoch[2,])) +
      t(sapply(symptomSum, FUN = `*`, e2 = stoch[3,])) +
      t(sapply(asymptomSum, FUN = `*`, e2 = stoch[4,])) +
      t(sapply(removalISum, FUN = `*`, e2 = stoch[5,])) +
      t(sapply(removalASum, FUN = `*`, e2 = stoch[6,]))


    return(list(indState = indState,
                noCases = noCases, StateX = StateX))
  }

  dailyProg <- compiler::cmpfun(dailyProg)
  # Return S, P, I, R and Hospital Admissions
  # S, P, I, R is a simulation of the underlying epidemic process
  # Hospital Admissions is a draw from the observation process given S, P, I, R, the total number of people who were admitted
  # to hospital.
  sim = function(param){
    beta_G <- param[1]
    beta_H <- param[2]
    epiParam <- param[3:6]
    a <- kernelParam[1]
    b <- kernelParam[2]
    # Pre-Lockdown infectious tranmission matrices
    K_P = kernel(a = 1, b = 1)
    K_I = 0.1*kernel(a = 1, b = 1)

    S <- E <- P <- I <- A <- R  <- matrix(nrow = noDays + 1, ncol = n)

    S0 <- StateX[,1]
    E0 <- StateX[,2]
    P0 <- StateX[,3]
    I0 <- StateX[,4]
    A0 <- StateX[,5]
    R0 <- StateX[,6]
    #D0 <- StateX[,6]

    S[1,] <- S0
    E[1,] <- E0
    P[1,] <- P0
    I[1,] <- I0
    A[1,] <- A0
    R[1,] <- R0
    #D[1,] <- D0

    dailyNoCases <- matrix(nrow = noDays, ncol = n)

    for(day in 1:noDays){
      if(day == lockdownDay + 1){
        # Pre-Lockdown infectious tranmission matrices
        K_P = 0.1*kernel(a, b)
        K_I = 0.1*kernel(a, b)
      }


      # ==== Daily Contact ====

      #debug(dailyProg)
      #print(c("==== Day", day))
      # if(day == 2){
      #   debug(dailyProg)
      # }
      dayProgression <- dailyProg(StateX, indState, beta_H, beta_G,
                                  epiParam, K_P, K_I)
      StateX <- dayProgression$StateX
      indState <- dayProgression$indState

      dailyNoCases[day, ] <- dayProgression$noCases

      S[day + 1,] <- StateX[,1]
      E[day + 1,] <- StateX[,2]
      P[day + 1,] <- StateX[,3]
      I[day + 1,] <- StateX[,4]
      A[day + 1,] <- StateX[,5]
      R[day + 1,] <- StateX[,6]
      #D[day + 1,] <- StateX[,6]
    }

    if(output == "Daily Cases"){
      return(dailyNoCases)
    } else if(output == "Epidemic"){
      return(list(S = S, E = E, P = P, I = I, A = A, R = R))
    }
  }
  return(sim)
}
