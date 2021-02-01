#' MetaPopulation Markov Epidemic Model for COVID-19 in UK
#' Structure based on 'Wuhan' pkg by Chris Jewell

#' @func  Markov Model
#' @param N populations of UK unified authorities
#' @param K contact matrix describing the movement between populations
#' @param stoch Stochiometry matrix describing individual disease progression within an individual (SEIR by default)
#' @param alpha parameter desribing reciprocal average exposure period
#' @param max_t time which simulator will simulate upto
#' @return A simulator which, when given epidemic parameters,
#'         will simulate the spread of COVID-19 in the UK according to the
#'         assumed model.

MarkovModel = function(N, K, stoch = matrix(c(-1, 1, 0, 0,
                                              0, -1, 1, 0,
                                              0, 0, -1, 1), byrow = T, ncol = 4, nrow = 3),
                                init_loc = "York", alpha = 1/4,
                                max_t = lubridate::as.period(lubridate::today() -
                                                               lubridate::as_date("2020-01-31"))@day){

  # Number of Cities
  n = length(N)

  # Days
  t = 0:max_t
  x =
  logDensity = 0
  simulate = function(param){
    # Initial infectives in Wuhan
    I0W = round(param[3])
    I0 = rep(0, length(N))
    I0[rownames(N)==init_loc] = I0W
    S0 = N - I0
    E0 = rep(0, length(N))
    R0 = E0

    y = c(S0, E0, I0, R0)
    ymat = matrix(y, nrow = n, byrow = FALSE)
    S = matrix(ncol = n, nrow = length(t))
    E = S
    I = S
    R = S

    S[1, ] = S0
    E[1, ] = E0
    I[1, ] = I0
    R[1, ] = R0

    day = t[1]
    currentTime = day
    while(sum(ymat[,2] > 0 | ymat[,3] > 0) >= 1 & day < max_t){
      # Calculate rate of each possible event
      exposureRates = ymat[,1]*param[1]*(ymat[,3]/N + (K%*%(ymat[,3]/N)/N))
      infectiousRates = ymat[,2]*alpha
      removalRates = ymat[,3]*param[2]
      rates = c(exposureRates, infectiousRates, removalRates)
      # Total rate will dictate the waiting time
      # distribution.
      cs = cumsum(c(exposureRates, infectiousRates, removalRates))
      totalRate = cs[length(cs)]

      waitingTime = rexp(1, rate = totalRate)
      logDensity = logDensity + dexp(waitingTime, rate = totalRate, log = T)
      currentTime = currentTime + waitingTime


      daysPassed = (currentTime > 0:max_t) & (day < 0:max_t)
      dP = sum(daysPassed)
      if(dP > 0){
        day = day + dP
        S[daysPassed,] = matrix(ymat[,1], nrow = dP, ncol = n, byrow = T)
        E[daysPassed,] = matrix(ymat[,2], nrow = dP, ncol = n, byrow = T)
        I[daysPassed,] = matrix(ymat[,3], nrow = dP, ncol = n, byrow = T)
        R[daysPassed,] = matrix(ymat[,4], nrow = dP, ncol = n, byrow = T)
      }

      # Which event occured and in what population?
      u1 = runif(1, 0, totalRate)
      index = sum(u1 > cs) + 1
      event = ceiling(index/n)
      population = index - n*(event - 1)
      logDensity = logDensity + punif(cs[index], 0, totalRate, log = T)

      # Alter state of population
      ymat[population, ] = ymat[population, ] + stoch[event, ]
    }
    if(day < max_t){
      daysNotPassed = 0:max_t > day
      dNP = sum(daysNotPassed)
      S[daysNotPassed,] = matrix(ymat[,1], nrow = dNP, ncol = n, byrow = T)
      E[daysNotPassed,] = matrix(ymat[,2], nrow = dNP, ncol = n, byrow = T)
      I[daysNotPassed,] = matrix(ymat[,3], nrow = dNP, ncol = n, byrow = T)
      R[daysNotPassed,] = matrix(ymat[,4], nrow = dNP, ncol = n, byrow = T)
    }
    list(t=t, S=S, E=E, I=I, R=R, logDensity = logDensity)
  }
  simulate
}



NoiseGeneratingFunction = function(param, N, K, max_t, simulator=NA) {
  realisation = simulator(param)
  n = length(N)
  p_detect = param[4]

  # Simulated Cumulative Cases (X)
  X = tail(realisation$R, n = 1)

  # Noise (y)
  y = rbinom(n, size = cumCases, prob = p_detect)

  names(y) = names(N)

  y
}



