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
testingNGF <- function(completeData, obsParam, obsWindow = c(1, Inf), PRINT = FALSE){

  # 1st entry in these matrices is the intial state before day 1 of pandemic
  completeData$Hstate$S <- completeData$Hstate$S[-1, ]
  completeData$Hstate$I <- completeData$Hstate$I[-1, ]
  completeData$Hstate$R <- completeData$Hstate$R[-1, ]



  # Cut down complete data to the days which epidemic is observed
  if(is.infinite(obsWindow[2])){
    obsWindow[2] <- nrow(completeData$Infections)
  }

  daysObserved <- obsWindow[1]:obsWindow[2]
  completeData$Hstate$S <- completeData$Hstate$S[daysObserved, ]
  completeData$Hstate$I <- completeData$Hstate$I[daysObserved, ]
  completeData$Hstate$R <- completeData$Hstate$R[daysObserved, ]

  completeData$Infections <- completeData$Infections[daysObserved, ]

  alpha <- obsParam[1]
  pi <- obsParam[2]
  psi <- obsParam[3]
  n_trials <- length(completeData$Infections)

  NROWS <- nrow(completeData$Infections)
  # Testing; Which of the new cases which are tested are detected.
  NCOLS <- ncol(completeData$Infections)

  # Which new cases are tested
  ascertainedIndividuals <- matrix(rbinom(n_trials, size = completeData$Infections, prob = alpha),
                  nrow = NROWS, ncol = NCOLS)



  whichHHascertained <- ascertainedIndividuals > 0

  # Testing; Which of the new cases which are tested are detected.
  # ptveFirstStage <- matrix(rbinom(n_trials, size = firstStageTests, prob = pi),
  #                         nrow = NROWS, ncol = NCOLS)
  # whichHHptveFirstStage <- ptveFirstStage > 0

  # Household Testing; For those who tested tpositive, all other residents are tested

  # The test will carry different uncertainty on returning the currect result depending
  # on the individuals actual disease status.

  # Testing those who are infected
  secondStageInfectedTests <- whichHHascertained * (completeData$Hstate$I - ascertainedIndividuals)
  ptveSecondStageInfectedTests <- matrix(rbinom(n_trials, size = secondStageInfectedTests,
                                           prob = pi), nrow = NROWS, ncol = NCOLS)

  # Testing those who are not infected
  secondStageNoninfectedTests <- whichHHascertained * (completeData$Hstate$S+
                                        completeData$Hstate$R)
  ntveSecondStageNoninfectedTests <- matrix(rbinom(n_trials, size = secondStageNoninfectedTests, prob = psi),
                                    nrow = NROWS, ncol = NCOLS)

  # When testing, the disease status is not known, so test results are collapsed into positive and negative results.
  ptveSecondStageTests <- ptveSecondStageInfectedTests + (secondStageNoninfectedTests - ntveSecondStageNoninfectedTests)
  ntveSecondStageTests <- ntveSecondStageNoninfectedTests + (secondStageInfectedTests - ptveSecondStageInfectedTests)
  #totalPositives <- positiveTests + positive_infected_tests + (noninfected_tests - negative_noninfected_tests)
  #totalNegatives <- (tests - PositiveTests) + negative_noninfected_tests + (infected_tests - positive_infected_tests)


  sampleData <- list(ascertainedIndividuals = ascertainedIndividuals,
                    ptveSecondStageTests = ptveSecondStageTests,
                    ntveSecondStageTests = ntveSecondStageTests)
  if(PRINT){
    print(sampleData)
  }
  return(sampleData)
}

# Density of positive and negative tests given infectious status of household.

ptveTestDensity_hh <- function(ntveTests, ptveTests, S, I, R, x, pi, psi){
  if(ntveTests + ptveTests == 0){
    return(0) # Second Stage Testing did not occur in this household
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


testingLlh <- function(completeData, sampleData, obsParam, PRINT = FALSE){
  alpha <- obsParam[1]
  pi <- obsParam[2]
  psi <- obsParam[3]
  NROWS <- nrow(completeData$Infections)
  NCOLS <- ncol(completeData$Infections)


  llh_alpha_breakdown <- dbinom(sampleData$ascertainedIndividuals, size = completeData$Infections, prob = alpha, log = TRUE)
  if(PRINT){
    print("Recruitment Process")
    print(llh_alpha_breakdown)
  }
  llh_alpha <- sum(llh_alpha_breakdown)

  # First stage testing probability, only reliant on sample dataset
  # llh_firstStageTest_breakdown <- dbinom(sampleData$ptveFirstStage, size = sampleData$firstStageTests, prob = pi, log = TRUE)
  # if(PRINT){
  #   print("First Stage Testing")
  #   print(llh_firstStageTest_breakdown)
  # }
  # llh_firstStageTest <- sum(llh_firstStageTest_breakdown)
  #print(llh_firstStageTest)

  if(llh_alpha == -Inf){
    if(PRINT){
      print(c("Testing Likelihood Not Calculated"))
      print(c("Total Log-likelihood:", -Inf))
    }

    return(-Inf)
  }

  func <- function(X, i){
    result <- ptveTestDensity_hh(sampleData$ntveSecondStageTests[i,X],
                       sampleData$ptveSecondStageTests[i,X],
                       completeData$Hstate$S[i + 1, X],
                       completeData$Hstate$I[i + 1, X],
                       completeData$Hstate$R[i + 1, X],
                       sampleData$ascertainedIndividuals[i, X],
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
  }
  if(PRINT){
    print("Second Stage Testing")
    print(llh_secondStageTest)

    print(c("Total Log-likelihood:", sum(llh_secondStageTest) + llh_alpha))
  }
  return(sum(llh_secondStageTest) + llh_alpha)
}









