


MarkovModelLanc = function(N, P0, kernel, stoch = matrix(c(-1, 1, 0, 0,
                                                       0, -1, 1, 0,
                                                       0, 0, -1, 1), byrow = T, ncol = 4, nrow = 3),
                       max_t){

  # Number of Metapopulations
  n = length(N)

  # Days
  t = 0:max_t
  logDensity = 0

  # param = c(beta.P, beta.70plus, beta.I, rho, alpha, gamma)
  simulate = function(param){

    K = kernel(param[1], param[2])

    S0 = N - P0
    I0 = rep(0, n)
    R0 = 0

    hospMat = matrix(0, nrow = 2, ncol = n)

    y = c(S0, P0, I0, R0)
    ymat = matrix(y, nrow = n, byrow = FALSE)
    S = matrix(ncol = n, nrow = length(t))
    P = S
    I = S
    R = S

    S[1, ] = S0
    P[1, ] = P0
    I[1, ] = I0
    R[1, ] = R0

    beta = rep(param[1], n)
    beta[n] = param[2]

    day = t[1]
    currentTime = day
    while(sum(ymat[,2] > 0 | ymat[,3] > 0) >= 1 & day < max_t){
      # Calculate rate of each possible event
      exposureRates = ymat[,1]*(beta*ymat[,2]/N + (K%*%(ymat[,2]/N)/N) + param[3]*(beta*ymat[,3]/N + (K%*%(ymat[,3]/N)/N)))
      infectiousRates = ymat[,2]*param[4]
      removalRates = ymat[,3]*param[6]
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
        P[daysPassed,] = matrix(ymat[,2], nrow = dP, ncol = n, byrow = T)
        I[daysPassed,] = matrix(ymat[,3], nrow = dP, ncol = n, byrow = T)
        R[daysPassed,] = matrix(ymat[,4], nrow = dP, ncol = n, byrow = T)
      }

      # Which event occured and in what population?
      u1 = runif(1, 0, totalRate)
      index = sum(u1 > cs) + 1
      event = ceiling(index/n)


      population = index - n*(event - 1)

      if(event == 2){
        hosp = rbinom(1, n = 1, prob = parma[5])

        if(hosp){
          hospMat[1, population] = hospMat[2, population] + 1
        } else{
          hospMat[2, population] = hospMat[2, population] + 1
        }
      }


      logDensity = logDensity + punif(cs[index], 0, totalRate, log = T)


      # Alter state of population
      ymat[population, ] = ymat[population, ] + stoch[event, ]
    }
    if(day < max_t){
      daysNotPassed = 0:max_t > day
      dNP = sum(daysNotPassed)
      S[daysNotPassed,] = matrix(ymat[,1], nrow = dNP, ncol = n, byrow = T)
      P[daysNotPassed,] = matrix(ymat[,2], nrow = dNP, ncol = n, byrow = T)
      I[daysNotPassed,] = matrix(ymat[,3], nrow = dNP, ncol = n, byrow = T)
      R[daysNotPassed,] = matrix(ymat[,4], nrow = dNP, ncol = n, byrow = T)
    }
    list(t=t, S=S, P=P, I=I, R=R, logDensity = logDensity)
  }
  simulate
}
