# Continuous-time epidemic simulation model thing


CT_Homo_SIR <- function(N, startTime = 0, endTime = Inf){

  # S -> I -> R
  stoch <- matrix(c(-1, 1, 0,
                    0, -1, 1), nrow = 2, ncol = 3, byrow = T)

  simEvent <- function(StateX, beta, gamma){
    # Calculates the number of infected individuals exerting pressure on each susceptible individual

    InfPressure <- beta*StateX[1]*StateX[2]
    removalRate <- StateX[2]*gamma
    totalEventRate <- InfPressure + removalRate

    # Time to next event
    #E <- rexp(1, rate = 1)
    timeToEvent <- rexp(1, rate = totalEventRate)

    # Event nature
    U <- runif(1, 0, 1)
    event <- (U < removalRate/totalEventRate) + 1

    # Update population state
    StateX <- StateX + stoch[event, ]

    logProb <- dexp(timeToEvent, totalEventRate, T) + dunif(U, 0, 1, T)

    return(list(StateX = StateX, timeToEvent = timeToEvent, logProb = logProb))
  }

  #dailyProg <- compiler::cmpfun(dailyProg)
  # Return S, P, I, R and Hospital Admissions
  # S, P, I, R is a simulation of the underlying epidemic process
  # Hospital Admissions is a draw from the observation process given S, P, I, R, the total number of people who were admitted
  # to hospital.
  #debug(sim)
  sim <- function(I0, param){

    # Define Initial State
    StateX <- c(N - I0, I0, 0)

    beta <- param[1]
    gamma <- param[2]

    # Store Epidemic State
    SIRsummary <- matrix(nrow = 1000, ncol = 4)

    # Counters
    time <- startTime
    eventNo <- 1
    SIRsummary[eventNo, ] <- c(StateX, time)

    logProb <- 0
    while(time < endTime  & StateX[2] > 0){
      # if(eventNo == 136){
      #   debug(simEvent)
      # }
      event <- simEvent(StateX, beta, gamma)
      StateX <- event$StateX
      s <- event$timeToEvent
      logProb <- logProb + event$logProb
      time <- time + s
      if(StateX[1] < 0){
        print(eventNo)
      }
      SIRsummary[eventNo + 1, ] <- c(StateX, time)
      eventNo <- eventNo + 1
      if(eventNo + 1 > 1000){
        SIRsummary <- rbind(SIRsummary, matrix(nrow = 1000, ncol = 4))
      }

    }

    # Remove NA rows
    rowsNA <- apply(SIRsummary, MARGIN = 1, FUN = function(X) prod(is.na(X)))
    SIRsummary <- SIRsummary[!rowsNA, ]
    finalSize <- sum(SIRsummary[nrow(SIRsummary),2:3]) - I0

    return(list(SIRsummary = SIRsummary, logProb = logProb,
                finalSize = finalSize))
  }
  return(sim)
}


#' Example
#' N <- 100
#' simulator <- CT_Homo_SIR(N)
#
#' set.seed(1)
#' simulator(10, param = c(1, 1))
