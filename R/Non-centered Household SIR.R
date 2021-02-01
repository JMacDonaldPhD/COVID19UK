#' Non-centered Household SIR Simulation (for data augmentation approach)



NonCentered_HouseholdSIR <- function(pop, startTime = 0, endTime){



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

  #print(StateX)

  # epiParam = c(beta.base, beta70plus, betaI, rho, gamma)
  # beta.base, infection rate between those not in over 70 category
  # beta70plus, infection rate for people in 70 plus category
  # betaI, reduction in infection rate due to quarantine or hospitalisation
  # rho, probability of moving from presymtomatic to symptomatic in one day
  # gamma, probability of dying or becoming immune to COVID19 in one day

  dailyProg2 <- function(StateX, Hstate, beta_G, beta_H, gamma){
    # Calculates the number of infected individuals exerting pressure on each susceptible individual

    globalInfPressure <- beta_G*StateX[2]/N

    householdInfPressure <- beta_H*Hstate[, 2]




    infProb <- 1 - exp(-(globalInfPressure + householdInfPressure))

    # S -> I
    U <- runif(n_hh)
    noInf <-sapply(X = 1:n_hh, FUN = function(X){
      c <- pbinom(0:Hstate[X,1], size = Hstate[X,1], prob = infProb[X])
      return(min(which(U[X] < c)) - 1)
    })

    # I -> R
    removalRate <- gamma
    V <- runif(n_hh)
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

    return(list(Hstate = Hstate, StateX = StateX, Infections = noInf,
                U = U, V = V))
  }



  dailyProg1 <- function(StateX, Hstate, beta_G, beta_H, gamma, U, V){
    # Calculates the number of infected individuals exerting pressure on each susceptible individual

    globalInfPressure <- beta_G*StateX[2]/N

    householdInfPressure <- beta_H*Hstate[, 2]




    infProb <- 1 - exp(-(globalInfPressure + householdInfPressure))

    # S -> I

    noInf <-sapply(X = 1:n_hh, FUN = function(X){
      c <- pbinom(0:Hstate[X,1], size = Hstate[X,1], prob = infProb[X])
      return(min(which(U[X] < c)) - 1)
    })

    # I -> R
    removalRate <- gamma

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


  dailyProg1 <- compiler::cmpfun(dailyProg1)
  dailyProg2 <- compiler::cmpfun(dailyProg2)

  # Return S, P, I, R and Hospital Admissions
  # S, P, I, R is a simulation of the underlying epidemic process
  # Hospital Admissions is a draw from the observation process given S, P, I, R, the total number of people who were admitted
  # to hospital.
  # The non-centered simulator will take 2 matrices (T by n_hh) of random unit uniform variables which determine how many
  # individuals transition from each household on each day, based on the CDF of the binomial draw which would
  # determine this usually.
  NCsim = function(param, U = NULL, V = NULL){
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

    if(is.null(U) & is.null(V)){
      U <- V <- matrix(nrow = noDays, ncol = n_hh)


      for(day in 1:noDays){

        # ==== Daily Contact ====

        #debug(dailyProg)
        #print(c("==== Day", day))
        # if(day == 2){
        #   debug(dailyProg)
        # }
        #debug(dailyProg)
        dayProgression <- dailyProg2(StateX, Hstate, beta_G, beta_H,
                                    gamma)
        U[day, ] <- dayProgression$U
        V[day, ] <- dayProgression$V
        StateX <- dayProgression$StateX
        Hstate <- dayProgression$Hstate
        epidemic$S[day + 1, ] <- Hstate[,1]
        epidemic$I[day + 1, ] <- Hstate[,2]
        epidemic$R[day + 1, ] <- Hstate[,3]
        SIRsummary[day + 1, ] <- dayProgression$StateX
        Infections[day, ] <- dayProgression$Infections
      }
    } else{
      for(day in 1:noDays){

        # ==== Daily Contact ====

        #debug(dailyProg)
        #print(c("==== Day", day))
        # if(day == 2){
        #   debug(dailyProg)
        # }
        #debug(dailyProg)
        dayProgression <- dailyProg1(StateX, Hstate, beta_G, beta_H,
                                    gamma, U[day, ], V[day, ])
        StateX <- dayProgression$StateX
        Hstate <- dayProgression$Hstate
        epidemic$S[day + 1, ] <- Hstate[,1]
        epidemic$I[day + 1, ] <- Hstate[,2]
        epidemic$R[day + 1, ] <- Hstate[,3]
        SIRsummary[day + 1, ] <- dayProgression$StateX
        Infections[day, ] <- dayProgression$Infections
      }
    }

    return(list(Hstate = epidemic, SIRsummary = SIRsummary,
                householdFinalSize = rowSums(Hstate[, 2:3]) - startingSize,
                Infections = Infections,
                U = U, V = V))
  }
  return(NCsim)
}
