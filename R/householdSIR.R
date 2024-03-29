#' Household SIR


HouseholdSIR <- function(pop, startTime = 0, endTime = Inf, PRINT = FALSE){



  # S -> I  -> R


  stoch <- matrix(c(-1, 1, 0,
                    0, -1, 1), nrow = 2, ncol = 3, byrow = T)
  #noDays <- endTime

  if(startTime > 0){
    cutData <- 1:startTime
  }

  # projects epidemic one day using binomial draws
  # Either takes P0 as an argument for the intial number of presymptomatic individuals
  # OR StateX if simulating from a timepoint which is mid-epidemic.



  # ==== Set up (Only needs to be done once per population) ====

  #hh_sizes <- as.numeric(table(pop$householdID))
  hh_size <- pop$N_M
  #print(c("Household Sizes", hh_sizes[1:10]))
  M <- length(pop$N_M)
  #n_hh <- length(unique(pop$householdID))

  #N <- nrow(pop)
  N <- sum(pop$N_M)

  #hhIndexMat <- sapply(1:n_hh, function(X) pop$householdID == X)

  #Hstate <- t(apply(hhIndexMat, MARGIN = 2, function(X) with(pop, fast.tabulate(state[X], householdID[X]))))

  #startingSize <- rowSums(Hstate[, 2:3])


  # Construct Population State Matrix using infectious status and groups
  #StateX <- colSums(Hstate)

  #print(StateX)


  # epiParam = c(beta.base, beta70plus, betaI, rho, gamma)
  # beta.base, infection rate between those not in over 70 category
  # beta70plus, infection rate for people in 70 plus category
  # betaI, reduction in infection rate due to quarantine or hospitalisation
  # rho, probability of moving from presymtomatic to symptomatic in one day
  # gamma, probability of dying or becoming immune to COVID19 in one day

  dailyProg <- function(StateX, Hstate, beta_G, beta_H, gamma){
    # Calculates the number of infected individuals exerting pressure on each susceptible individual

    globalInfPressure <- beta_G*StateX[2]/N

    householdInfPressure <- beta_H*Hstate[, 2]

    infProb <- 1 - exp(-(globalInfPressure + householdInfPressure))

    # S -> I
    noInf <- rbinom(M, size = Hstate[, 1], prob = infProb)
    logProb <- dbinom(noInf, size = Hstate[,1], prob = infProb)
    # I -> R
    removalRate <- gamma
    noRem <- rbinom(M, size = Hstate[, 2], prob = gamma)
    logProb <- logProb + dbinom(noRem, size = Hstate[,2], prob = gamma)

    # Update Household States
    Hstate <- Hstate +
      noInf*matrix(stoch[1, ], nrow = M, ncol = 3, byrow = T) +
      noRem*matrix(stoch[2, ], nrow = M, ncol = 3, byrow = T)

    # Update population state
    totInf <- sum(noInf)
    totRem <- sum(noRem)
    StateX <- StateX + c(-totInf, totInf - totRem, totRem)

    # Population State and Household State Check
    if(!identical(colSums(Hstate), StateX)){
      warning("Population State and Household States do not line up")
    }

    return(list(Hstate = Hstate, StateX = StateX, Infections = noInf, logProb = logProb))
  }

  dailyProg <- compiler::cmpfun(dailyProg)
  # Return S, P, I, R and Hospital Admissions
  # S, P, I, R is a simulation of the underlying epidemic process
  # Hospital Admissions is a draw from the observation process given S, P, I, R, the total number of people who were admitted
  # to hospital.
  #debug(sim)
  sim <- function(I0, param){

    Hstate <- cbind(pop$N_M - I0, I0, rep(0, M))
    Hstate <- unname(Hstate)
    StateX <- c(N - sum(I0), sum(I0), 0)
    beta_G <- param[1]
    beta_H <- param[2]
    gamma <- param[3]

    maxSimDays <- 1e3
    # Stores summry of each epidemic state at the end of each day
    if(is.infinite(endTime)){
      SIRsummary <- matrix(nrow = maxSimDays + 1, ncol = 3)
      A <- matrix(nrow = maxSimDays + 1, ncol = M)
    } else{
      SIRsummary <- matrix(nrow = endTime + 1, ncol = 3)
      A <- matrix(nrow = endTime + 1, ncol = M)
    }

    time <- 0
    SIRsummary[time + 1, ] <- StateX
    epidemic <- list(S = A, I = A, R = A)
    epidemic$S[time + 1, ] <- Hstate[,1]
    epidemic$I[time + 1, ] <- Hstate[,2]
    epidemic$R[time + 1, ] <- Hstate[,3]

    Infections <- A[-1, , drop = F]

    logProb <- 0
    while(time < endTime  & time < maxSimDays){

      # ==== Daily Contact ====

      #debug(dailyProg)
      #print(c("==== Day", day))
      # if(day == 2){
      #   debug(dailyProg)
      # }
      #debug(dailyProg)
      dayProgression <- dailyProg(StateX, Hstate, beta_G, beta_H,
                                  gamma)
      StateX <- dayProgression$StateX
      Hstate <- dayProgression$Hstate
      time <- time + 1

      epidemic$S[time + 1, ] <- Hstate[,1]
      epidemic$I[time + 1, ] <- Hstate[,2]
      epidemic$R[time + 1, ] <- Hstate[,3]
      SIRsummary[time + 1, ] <- dayProgression$StateX
      Infections[time, ] <- dayProgression$Infections
      logProb <- logProb + dayProgression$logProb
      if(StateX[2] == 0 & time < endTime){
        for(i in (time + 1):(nrow(epidemic$S) - 1)){
          epidemic$S[i + 1, ] <- Hstate[,1]
          epidemic$I[i + 1, ] <- Hstate[,2]
          epidemic$R[i + 1, ] <- Hstate[,3]
          SIRsummary[i + 1, ] <- dayProgression$StateX
          Infections[i, ] <- rep(0, nrow(Hstate))
        }
        time <- endTime + 1
      }
    }
    if(PRINT){
      print("Day by day household States:")
      print(epidemic)
      print("Daily Infections by Household:")
      print(Infections)
    }

    # Cut Data

    if(startTime > 0){
      epidemic$S <- epidemic$S[-cutData, ]
      epidemic$I <- epidemic$I[-cutData, ]
      epidemic$R <- epidemic$R[-cutData, ]
      SIRsummary <- SIRsummary[-cutData, ]
      Infections[-cutData[-length(cutData)], ]

    }
    #Hstate2 <- Hstate
    return(list(Hstate = epidemic, SIRsummary = SIRsummary,
                householdFinalSize = rowSums(Hstate[, 2:3]),
                Infections = Infections, logProb = logProb))
  }
  return(sim)
}
