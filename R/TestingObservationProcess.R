# Testing

#alpha <- 0.05 # testing rate

#pi <- 1 # Sensitivity, True Positive Rate
#psi <- 1 # Specitivity, True Negative Rate

# On day of infection, an individual has probablity alpha of being tested
# If the test is positive, each individual in the household is tested.

# Bernoulli Draw for deciding to test

# Bernoulli Draw for ascertaining infection

# If infection ascertained, sequence of bernouilli draws which should aim
# to ascertain infection if individual infected or no infection is individual
# is not infected.

# So observed data will consist of true positive tests with further household tests (which may be true or false)
# AND false negative tests where a new infection arises and is tested but fails to be detected.

# If pi = 1, psi = 1, then we arrive at the original biased household sampling scheme.
testingNGF <- function(completeData, obsParam){
  alpha <- obsParam[1]
  pi <- obsParam[2]
  psi <- obsParam[3]
  n_trials <- length(completeData$Infections)

  NROWS <- nrow(completeData$Infections)
  # Testing; Which of the new cases which are tested are detected.
  NCOLS <- ncol(completeData$Infections)

  # Which new cases are tested
  firstStageTests <- matrix(rbinom(n_trials, size = completeData$Infections, prob = alpha),
                  nrow = NROWS, ncol = NCOLS)
  whichHHfirstStage <- firstStageTests > 0

  # Testing; Which of the new cases which are tested are detected.
  ptveFirstStage <- matrix(rbinom(n_trials, size = firstStageTests, prob = pi),
                          nrow = NROWS, ncol = NCOLS)
  whichHHptveFirstStage <- ptveFirstStage > 0

  # Household Testing; For those who tested tpositive, all other residents are tested

  # The test will carry different uncertainty on returning the currect result depending
  # on the individuals actual disease status.

  # Testing those who are infected
  secondStageInfectedTests <- whichHHptveFirstStage * (completeData$Hstate$I[-1, ] - firstStageTests)
  ptveSecondStageInfectedTests <- matrix(rbinom(n_trials, size = secondStageInfectedTests,
                                           prob = pi), nrow = NROWS, ncol = NCOLS)

  # Testing those who are not infected
  secondStageNoninfectedTests <- whichHHptveFirstStage * (completeData$Hstate$S[-1, ] +
                                        completeData$Hstate$R[-1, ])
  ntveSecondStageNoninfectedTests <- matrix(rbinom(n_trials, size = secondStageNoninfectedTests, prob = psi),
                                    nrow = NROWS, ncol = NCOLS)

  # When testing, the disease status is not known, so test results are collapsed into positive and negative results.
  ptveSecondStageTests <- ptveSecondStageInfectedTests + (secondStageNoninfectedTests - ntveSecondStageNoninfectedTests)
  ntveSecondStageTests <- ntveSecondStageNoninfectedTests + (secondStageInfectedTests - ptveSecondStageInfectedTests)
  #totalPositives <- positiveTests + positive_infected_tests + (noninfected_tests - negative_noninfected_tests)
  #totalNegatives <- (tests - PositiveTests) + negative_noninfected_tests + (infected_tests - positive_infected_tests)


  sampleData <- list(firstStageTests = firstStageTests, ptveFirstStage = ptveFirstStage,
                    ptveSecondStageTests = ptveSecondStageTests,
                    ntveSecondStageTests = ntveSecondStageTests)

  return(sampleData)
}

# Density of positive and negative tests given infectious status of household.

ptveTestDensity_hh <- function(ntveTests, ptveTests, S, I, R, x, pi, psi){
  if(ntveTests + ptveTests == 0){
    return(0)
  }
  # possible True Positives, lpp
  poss_tp <- (0:(I - x))[0:(I - x) <= ptveTests]
  # True Negatives, lnn
  poss_tn <- (0:(S + R))[0:(S + R) <= ntveTests]

  combinations <- expand.grid(poss_tp, poss_tn)

  combinations$poss <- ptveTests == combinations[,1] + (S + R - combinations[,2])

  llh <- log(sum(apply(combinations[combinations$poss, 1:2], MARGIN = 1,
                      FUN = function(X) prod(dbinom(X, size = c(I-x, S + R), prob = c(pi, psi))))
                 )
             )
  return(llh)
}


testingLlh <- function(completeData, sampleData, obsParam){
  alpha <- obsParam[1]
  pi <- obsParam[2]
  psi <- obsParam[3]
  NROWS <- nrow(completeData$Infections)
  NCOLS <- ncol(completeData$Infections)
  llh_alpha <- sum(dbinom(sampleData$firstStageTests, size = completeData$Infections, prob = alpha, log = TRUE))
  #print(llh_alpha)

  # First stage testing probability, only reliant on sample dataset
  llh_firstStageTest <- sum(dbinom(sampleData$ptveFirstStage, size = sampleData$firstStageTests, prob = pi, log = TRUE))
  #print(llh_firstStageTest)

  if(llh_alpha == -Inf | llh_firstStageTest == -Inf){
    return(-Inf)
  }

  func <- function(X, i){
    result <- ptveTestDensity_hh(sampleData$ntveSecondStageTests[i,X],
                       sampleData$ptveSecondStageTests[i,X],
                       completeData$Hstate$S[i + 1, X],
                       completeData$Hstate$I[i + 1, X],
                       completeData$Hstate$R[i + 1, X],
                       sampleData$firstStageTests[i, X],
                       pi, psi)
    if(is.nan(result)){
      print(c(i, X))
    }
    return(result)
  }
  llh_secondStageTest <- matrix(nrow = NROWS, ncol = NCOLS)

  for(i in 1:NROWS){
    for(X in 1:NCOLS){
      llh_secondStageTest[i, X] <- func(X, i)
    }
    #p <- sum(sapply(X = 1:NCOLS, FUN = func, i = i))
    #llh_secondStageTest <- p + llh_secondStageTest

    #print(c(i, llh_secondStageTest))
  }
  return(sum(llh_secondStageTest) + llh_alpha + llh_firstStageTest)
}









