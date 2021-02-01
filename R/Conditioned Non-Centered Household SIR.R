#' Conditioned Non-Centered Household SIR simulation

#' Non-centered Household SIR Simulation (for data augmentation approach)



Conditioned_NonCentered_HouseholdSIR <- function(pop, sampleData, startTime = 0, endTime){



  # S -> I  -> R
  stoch <- matrix(c(-1, 1, 0,
                    0, -1, 1), nrow = 2, ncol = 3, byrow = T)
  noDays <- abs(endTime - startTime)


  # projects epidemic one day using binomial draws
  # Either takes P0 as an argument for the intial number of presymptomatic individuals
  # OR StateX if simulating from a timepoint which is mid-epidemic.



  # ==== Set up (Only needs to be done once per population) ====

  hh_sizes <- as.numeric(table(pop$householdID))
  #print(c("Household Sizes", hh_sizes[1:10]))
  n_hh <- length(unique(pop$householdID))

  N <- nrow(pop)

  hhIndexMat <- sapply(1:n_hh, function(X) pop$householdID == X)

  Hstate <- t(apply(hhIndexMat, MARGIN = 2, function(X) with(pop, fast.tabulate(state[X], householdID[X]))))

  startingSize <- rowSums(Hstate[, 2:3])


  # Construct Population State Matrix using infectious status and groups
  StateX <- colSums(Hstate)

  # What is easiet form to have Sample Data in?
  # Split into elements of list by household?
  obs_index <- which(colSums(sampleData$observed_HH) > 0)

  # beta_G, Global Rate of Infection
  # beta_H, Within Household of Infection
  # gamma, probability of dying or becoming immune to COVID19 in one day

  dailyProg <- function(day, StateX, Hstate, beta_G, beta_H, gamma, U, V, P_bias){
    # Calculates the number of infected individuals exerting pressure on each susceptible individual

    globalInfPressure <- beta_G*StateX[2]/N

    householdInfPressure <- beta_H*Hstate[, 2]




    infProb <- 1 - exp(-(globalInfPressure + householdInfPressure))

    # S -> I
    # noInf <- c()
    # for(i in 1:n_hh){
    #   c <- pbinom(0:Hstate[i,1], size = Hstate[i,1], prob = infProb[i])
    #   noInf[i] <- min(which(U[i] < c)) - 1
    # }
    noInf <-sapply(X = 1:n_hh, FUN = function(X){
      c <- pbinom(0:Hstate[X,1], size = Hstate[X,1], prob = infProb[X])
      return(which(U[X] < c)[1] - 1)
    })

    for(i in obs_index){
      # Introduce Bias to make sure simulation lines up with observed data
      # x_it
      no_obs <- sampleData$observedCases[day, i] # Need at least no_obs infections

      # If there are too many infections on this day, there will not be enough susceptibles
      # for infections on later days.
      # How many infections are there on later days?

      # s > t
      futureDays <- 1:noDays > day

      # \sum_{s > t} x_is
      futureKnownInfections <- sum(sampleData$observedCases[futureDays, i])


      # Conditions on number of infections
      # x_it < y_it < S_i(t-1) - \sum_{s > t} x_is

      possibleNoInfections <- no_obs:(Hstate[i, 1] - futureKnownInfections)

      # Conditions on number of recoveries
      # z_it <= S_i(t-1) + I_i(t-1) - I_is for all s > t
      # -> z_it <= S_i(t-1) + I_i(t-1) - max_{s>t}(I_is)
      possibleNoRecoveries <- 0:(Hstate[i, 1] + Hstate[i, 2] - max(sampleData$observedInfectious[futureDays, i]))

      possibleComb <- expand.grid(possibleNoInfections, possibleNoRecoveries)

      allowedComb <- apply(possibleComb, MARGIN = 1, FUN = `-`) == sampleData$observedInfectious[day, i] - Hstate[i, 2]

      c <- pbinom((no_obs - 1):Hstate[i, 1], size = Hstate[i, 1], prob = infProb[i])
      bounds <- c(min(c), max(c))
      P_contribution <- max(c) - min(c) # Probability that the bias occured
      P_bias <- P_bias * P_contribution # Add the bias on to past bias

      # Simulate within the condition(s)
      U_cond <- runif(1, a  = bounds[1], b = bounds[2])
      noInf <- ((no_obs - 1):Hstate[i, 1])[which(U[X] < c)[1]]

    }

    if(!identical(noInf, noInf2)){
      stop("Two Methods don't give same result!")
    }
    # I -> R
    removalRate <- gamma
    # noRem <- c()
    # for(i in 1:n_hh){
    #   c <- pbinom(0:Hstate[i, 2], size = Hstate[i, 2], prob = gamma)
    #   noRem[i] <- min(which(V[i] < c)) - 1
    # }
    noRem <-sapply(X = 1:n_hh, FUN = function(X){
      c <- pbinom(0:Hstate[X,2], size = Hstate[X,2], prob = gamma)
      return(min(which(V[X] < c)) - 1)
    })

    # Update Household States
    Hstate <- Hstate +
      noInf*matrix(stoch[1, ], nrow = n_hh, ncol = 3, byrow = T) +
      noRem*matrix(stoch[2, ], nrow = n_hh, ncol = 3, byrow = T)

    # Update population state
    totInf <- sum(noInf)
    totRem <- sum(noRem)
    StateX <- StateX + c(-totInf, totInf - totRem, totRem)

    # Population State and Household State Check
    if(!identical(colSums(Hstate), StateX)){
      warning("Population State and Household States do not line up")
    }

    return(list(Hstate = Hstate, StateX = StateX, Infections = noInf))
  }

  dailyProg <- compiler::cmpfun(dailyProg)
  # Return S, P, I, R and Hospital Admissions
  # S, P, I, R is a simulation of the underlying epidemic process
  # Hospital Admissions is a draw from the observation process given S, P, I, R, the total number of people who were admitted
  # to hospital.
  # The non-centered simulator will take 2 matrices (T by n_hh) of random unit uniform variables which determine how many
  # individuals transition from each household on each day, based on the CDF of the binomial draw which would
  # determine this usually.
  NCsim = function(param, U, V){
    beta_G <- param[1]
    beta_H <- param[2]
    gamma <- param[3]

    # Stores summry of each epidemic state at the end of each day
    SIRsummary <- matrix(nrow = noDays + 1, ncol = 3)
    SIRsummary[1, ] <- StateX

    A <- matrix(nrow = noDays + 1, ncol = n_hh)
    epidemic <- list(S = A, I = A, R = A)
    epidemic$S[1, ] <- Hstate[,1]
    epidemic$I[1, ] <- Hstate[,2]
    epidemic$R[1, ] <- Hstate[,3]

    Infections <- A[-1, ]

    P_bias <- 0
    for(day in 1:noDays){

      # ==== Daily Contact ====

      #debug(dailyProg)
      #print(c("==== Day", day))
      # if(day == 2){
      #   debug(dailyProg)
      # }
      debug(dailyProg)
      dayProgression <- dailyProg(day, StateX, Hstate, beta_G, beta_H,
                                  gamma, U[day, ], V[day, ], P_bias)
      StateX <- dayProgression$StateX
      Hstate <- dayProgression$Hstate
      epidemic$S[day + 1, ] <- Hstate[,1]
      epidemic$I[day + 1, ] <- Hstate[,2]
      epidemic$R[day + 1, ] <- Hstate[,3]
      SIRsummary[day + 1, ] <- dayProgression$StateX
      Infections[day, ] <- dayProgression$Infections
    }
    return(list(Hstate = epidemic, SIRsummary = SIRsummary,
                householdFinalSize = rowSums(Hstate[, 2:3]) - startingSize,
                Infections = Infections))
  }
  return(NCsim)
}
