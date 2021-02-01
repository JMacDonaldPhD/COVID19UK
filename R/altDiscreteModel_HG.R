# Alt Discrete Model HG



altDiscreteModel_HG = function(pop, stoch = matrix(c(-1, 1, 0, 0, 0, 0,
                                                  0, -1, 1, 0, 0, 0,
                                                  0, 0, -1, 1, 0, 0,
                                                  0, 0, 0, -1, 1, 0,
                                                  0, 0, 0, -1, 0, 1), byrow = T, ncol = 6, nrow = 5),
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

  group <- factor(as.integer(pop$ageBin), levels = 1:7,
                  labels = as.character(1:7))
  n <- length(unique(group))


  hh_block_mat <- lapply(X = hh_sizes, FUN = function(X){
    mat <- matrix(1, ncol = X, nrow = X)
    diag(mat) <- 0
    mat
  })
  hh_mat <- Matrix::bdiag(hh_block_mat)


  # Individual Group Status
  indGroup <- as.integer(pop$ageBin)


  # individual infectious status
  susc <- c(T, F, F, F, F, F)
  exp <- c(F, T, F, F, F, F)
  pre <- c(F, F, T, F, F, F)
  inf <- c(F, F, F, T, F, F)
  rem <- c(F, F, F, F, T, F)
  death <- c(F, F, F, F, F, T)
  states <- matrix(c(susc, exp, pre, inf, rem, death), nrow = length(susc),
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

  # epiParam = c(beta.base, beta70plus, betaI, rho, alpha, gamma)
  # beta.base, infection rate between those not in over 70 category
  # beta70plus, infection rate for people in 70 plus category
  # betaI, reduction in infection rate due to quarantine or hospitalisation
  # rho, probability of moving from presymtomatic to symptomatic in one day
  # alpha, probability of being hospitalised given you have become symptomatic
  # gamma, probability of dying or becoming immune to COVID19 in one day

  dailyProg <- function(StateX, indState, beta_H, beta_G, epiParam, K_P, K_I){
    #print("==== Hous ehold Contacts ====")
    # Calculates the number of infected individuals exerting pressure on each susceptible individual
    infPressure_H <- beta_H*textTinyR::sparse_Sums(hh_mat[indState[,1], indState[,2] | indState[,3], drop = F], rowSums = T)

    # Total Infectious Pressure from pre-symptomatic people
    Lambda_P <- K_P%*%(StateX[,2])


    # Total Infectious pressure from infectious individuals on susceptible people
    Lambda_I <- K_I%*%(StateX[,3])

    #infPressure_G <- c(0, N)
    groupInfPressure <- beta_G*(Lambda_P + Lambda_I)/N

    infPressure_G <- groupInfPressure[indGroup][indState[,1]]

    #indInfPressure <- -infPressure_G + infPressure_H
    noS <- sum(StateX[,1])
    infections <- rbinom(noS, size = 1,
                         prob = 1 - exp(-infPressure_G - infPressure_H))
      #runif(noS) < (1 - exp(-infPressure_G - infPressure_H))
      #rbinom(noS, size = 1, prob = 1 - exp(-infPressure_G - infPressure_H))

    noInfection <- sum(infections == 1)
    #noInfection <- sum(infections)

    if(noInfection > 0){

      #infectionSum <- table(group[indState[,1]], infections)

      infectionSum <- fast.tabulate(infections,  group[indState[,1]])
      # if(!identical(infectionSum, infectionSum2)){
      #   stop("Methods not equivalent")
      # }


      if(dim(infectionSum)[2] == 2){
        infectionSum <- infectionSum[,2]
      }
    } else{
      infectionSum <- rep(0, n)
    }


    # Draws whether susceptibles become infected
    #infections <- rbinom(sum(indState[,1]), size = 1, prob = 1 - exp(-infPressure)) + 1

    eta <- epiParam[1]
    rho <- epiParam[2]
    gamma <- epiParam[3]
    phi <- epiParam[4]
    alpha <- epiParam[5]


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

    # Pre-symptomatic -> Symptomatic
    symptomRate <- rho
    noP <- sum(StateX[,3])
    symptoms <- runif(noP) < (symptomRate)
      #rbinom(noP, size = 1, symptomRate)

    #noSymptoms <- sum(symptoms == 1)
    noSymptoms <- sum(symptoms)
    if(noSymptoms > 0){
      symptomSum <- fast.tabulate(symptoms,  group[indState[,3]])
      #symptomSum <- table(group[indState[,2]], symptoms)
      if(dim(symptomSum)[2] == 2){
        symptomSum <- symptomSum[,2]
      }
    } else{
      symptomSum <- rep(0, n)
    }

    positiveTests <- rbinom(n, size = symptomSum, prob = alpha)

    noCases <- symptomSum

    recoveryRate <- gamma
    deathRate <- phi
    noI <- sum(StateX[,4])
    finalState <- runif(noI)
    deaths <- finalState < phi
    recovery <- (finalState > phi) & (finalState < phi + gamma)
    #recovery = rbinom(noI, size = 1, recoveryRate)
      #runif(noI) < (recoveryRate)
      #rbinom(noI, size = 1, recoveryRate)

    #noRecovery <- sum(recovery == 1)
    noRecovery <- sum(recovery)
    noDeath <- sum(deaths)
    if(noRecovery > 0){
      #recoverySum <- table(group[indState[,3]], recovery)
      recoverySum <- fast.tabulate(recovery,  group[indState[,4]])
      # if(!identical(recoverySum, recoverySum2)){
      #   stop("Methods not equivalent")
      # }

      if(dim(recoverySum)[2] == 2){
        recoverySum <- recoverySum[,2]
      }
    } else{
      recoverySum <- rep(0, n)
    }
    if(noDeath > 0){
      #recoverySum <- table(group[indState[,3]], recovery)
      deathSum <- fast.tabulate(deaths,  group[indState[,4]])
      # if(!identical(recoverySum, recoverySum2)){
      #   stop("Methods not equivalent")
      # }

      if(dim(deathSum)[2] == 2){
        deathSum <- deathSum[,2]
      }
    } else{
      deathSum <- rep(0, n)
    }



    # Hospital Admissions (Positive Testing)

    # recovery[indState[,3]]


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

    if(noRecovery > 0){
      if(noI == 1){
        indState[savedIndState[,4], ] <- states[5, ]
      } else{
        indState[savedIndState[,4], ][recovery, ] <- matrix(states[5, ],
                                                                 nrow = noRecovery,
                                                                 ncol = noStates,
                                                                 byrow = T)
      }
    }


    if(noDeath > 0){
      if(noI == 1){
        indState[savedIndState[,4], ] <- states[6, ]
      } else{
        indState[savedIndState[,4], ][deaths, ] <- matrix(states[6, ],
                                                                 nrow = noRecovery,
                                                                 ncol = noStates,
                                                                 byrow = T)
      }
    }

    StateX <- StateX + t(sapply(infectionSum, FUN = `*`, e2 = stoch[1,])) +
      t(sapply(infectiousSum, FUN = `*`, e2 = stoch[2,])) +
      t(sapply(symptomSum, FUN = `*`, e2 = stoch[3,])) +
      t(sapply(recoverySum, FUN = `*`, e2 = stoch[4,])) +
      t(sapply(deathSum, FUN = `*`, e2 = stoch[5,]))


    # Construct Population State Matrix using infectious status and groups
    StateCheck = matrix(nrow = n, ncol = noStates)

    for(i in 1:n){
      for(j in 1:noStates){
        StateCheck[i, j] = sum(group[indState[, j]] == i)
      }
    }

    if(sum(StateCheck != StateX) > 0){
      stop("Individual States and Population State differ")
    }
    if(sum(rowSums(indState) > 1)){
      stop("invalid Individual State(s)")
    }

    return(list(positiveTests = positiveTests, indState = indState,
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
    epiParam <- param[3:7]
    # Pre-Lockdown infectious tranmission matrices
    K_P = kernelPreLD(beta_G)
    K_I = kernelPreLD(beta_G*0.1)

    S <- E <- P <- I <- R <- D <- matrix(nrow = noDays + 1, ncol = n)

    S0 <- StateX[,1]
    E0 <- StateX[,2]
    P0 <- StateX[,3]
    I0 <- StateX[,4]
    R0 <- StateX[,5]
    D0 <- StateX[,6]

    S[1,] <- S0
    E[1,] <- E0
    P[1,] <- P0
    I[1,] <- I0
    R[1,] <- R0
    D[1,] <- D0

    dailyPositiveTests <- dailyNoCases <- matrix(nrow = noDays, ncol = n)

    for(day in 1:noDays){
      #print(c("======== day", day))
      if(day == lockdownDay + 1){
        # Pre-Lockdown infectious tranmission matrices
        K_P = kernelPostLD(beta_G*0.01)
        K_I = kernelPostLD(beta_G*0.1*0.01)
      }


      # ==== Daily Contact ====

      #debug(dailyProg)
      # print(c("==== Day", day))
      # if(day == 2){
      #   debug(dailyProg)
      # }
        dayProgression <- dailyProg(StateX, indState, beta_H, beta_G,
                                  epiParam, K_P, K_I)
      StateX <- dayProgression$StateX
      indState <- dayProgression$indState

      dailyPositiveTests[day, ] <- dayProgression$positiveTests
      dailyNoCases[day, ] <- dayProgression$noCases

      S[day + 1,] <- StateX[,1]
      E[day + 1,] <- StateX[,2]
      P[day + 1,] <- StateX[,3]
      I[day + 1,] <- StateX[,4]
      R[day + 1,] <- StateX[,5]
      D[day + 1,] <- StateX[,6]
    }
    return(list(dailyNoCases = dailyNoCases,
                dailyPositiveTests = dailyPositiveTests))
  }
  return(sim)
}

