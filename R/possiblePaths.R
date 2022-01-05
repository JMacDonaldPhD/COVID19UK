# Counting paths


# How many different events per step
# Number of events is random to though
# How many states can you get to from the current state?

# How many of these states line up with the data?

# What is the likelihood of getting to these states from the given current
# state?

# Hypothesis: In general, the likelihood of generating a valid path will decrease
#             the more events which are assumed to occur.

# SIR with fully exchangable population

# X_0 to X_t, is it possible to count the paths?



paths <- function(X0, Xt){


  noInfections <- X0[1] - Xt[1]

  noRemovals <- Xt[3] - X0[3]
  # For every event, there is 2 choices, infection or removal
  # Until the total number of one event has been reached, then the rest of the
  # events are the other one.

  # Multisets of infection and removals
  # Create Env for the possible sequence of infections and removal
  I <- iterpc::iterpc(c(noInfections, noRemovals), labels = c(1, 2), ordered = T)
  possiblePaths <- iterpc::getall(I)

  return(possiblePaths)
}

probPath <- function(path, beta, gamma, X0, t){
  stoch <- matrix(c(-1, 1, 0,
                    0, -1, 1),
                  nrow = 2, ncol = 3, byrow = T)
  infRate <- function(X, beta){
    return(X[1]*X[2]*beta)
  }
  remRate <- function(X, gamma){
    return(X[2]*gamma)
  }

  # Efficiently getting the states
  # Naively do a for loop

  N <- length(path)
  states <- matrix(nrow = N + 1, ncol = 3)


  # Does this give right shape matrix (want N+1 by 3)
  states[1:(N + 1),] <- matrix(X0, byrow = T, nrow = N + 1, ncol = 3)
  states[2:(N + 1), ] <- states[2:(N + 1), ] + t(sapply(1:length(path), FUN = function(X){
    return(colSums(stoch[path[1:X], , drop = F]))
  }))


  # Does this give right shape matrix? (want N + 1 by 2)
  rates <- t(apply(X = states[1:N, ], MARGIN = 1, FUN = function(X, beta, gamma){
    return(c(infRate(X, beta), remRate(X, gamma)))
  }, beta = beta, gamma = gamma))


  lambda <- rowSums(rates)


  pEvent <- prod(rates[, path]/lambda)
  pTime <- pPathTime(path_lambda = lambda, t)
  #print(c(event = pEvent, time = pTime))
  p <- prod(c(pEvent, pTime))
  return(p)

}

pSim <- function(X0, Xt, beta, gamma, t){

  possiblePaths <- paths(X0, Xt)

  probPaths <- apply(X = possiblePaths, MARGIN = 1, FUN = probPath, beta = beta,
                     gamma = gamma, X0 = X0, t = t)
  p <- sum(probPaths)

  return(p)

}

# Each event must have a time which is less than the time the next state is
# at.
# Last event will probably not reach exactly t. The condition is that the next
# event will pass t.

# What is the probability that the sum of all event times is less than t

# each event time T_i ~ Exp() dependent on the current state.
# So if we have a sequence of events up until time t, we know the distribution
# of every event time. We then condition that the sum of this event times is
# less than t. WORKED OUT


pPathTime <- function(path_lambda, pathTime){
  N <- length(path_lambda)
  p <- prod((path_lambda[-N] - path_lambda[N])/path_lambda[-N])*(1 - exp(-path_lambda[N]*pathTime))
  return(p)
}













