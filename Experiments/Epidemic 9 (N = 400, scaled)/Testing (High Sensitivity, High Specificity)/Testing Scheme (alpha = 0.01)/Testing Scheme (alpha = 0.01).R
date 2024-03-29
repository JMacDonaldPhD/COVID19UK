rm(list = ls())
setwd("C:/Users/Work/Documents/git/COVID19UK/Experiments/Epidemic 9 (N = 400, scaled)/")

load("epidemic 9.RData")



setwd("C:/Users/Work/Documents/git/COVID19UK/Experiments/Epidemic 9 (N = 400, scaled)/Testing (High Sensitivity, High Specificity)/Testing Scheme (alpha = 0.01)/")

`%<<%` <- function(x, y){
  return(x > y[1] & x < y[2])
}



# ==== Simulating representative testing Data ====
obsParam <- c(alpha = .01, pi = .9, psi = .98)
obsWindow <- c(0, 5)
simulator <- COVID19UK::HouseholdSIR(pop, endTime = obsWindow[2], PRINT = FALSE) # Default start time is startTime = 0

samples <- replicate(n = 1e4, COVID19UK::testingNGF(simulatedEpidemic, obsParam = obsParam, obsWindow = obsWindow, PRINT = FALSE), simplify = F)

totalTests <- sapply(X = samples, FUN = function(X) sum(X$ntveSecondStageTests + X$ptveSecondStageTests))
totalPosTests <- sapply(X = samples, FUN = function(X) sum(X$ptveSecondStageTests))
totalNegTests <- sapply(X = samples, FUN = function(X) sum(X$ntveSecondStageTests))

#png(filename = "testing_summary_densities.png", width = 480*3)


#par(mfrow = c(1, 3))
#
# plot(density(totalTests, from = 0, to = N), xlab = "Total Number of Tests", main = "", xlim = c(0, max(totalTests) + 5))
# plot(density(totalPosTests, from = 0, to = N), xlab = "Total Number of Positive Tests", main = "", xlim = c(0, max(totalPosTests) + 5))
# plot(density(totalNegTests, from = 0, to = N), xlab = "Total Number of Negative Tests", main = "", xlim = c(0, max(totalNegTests) + 5))
#
#
# dev.off()

png(filename = "testing_summary_hists.png", width = 480*3)

par(mfrow = c(1, 3))

hist(totalTests, breaks = seq(0, max(totalTests) + 5, by = 2.5))
hist(totalPosTests, breaks = seq(0, max(totalPosTests) + 5, by = 2))
hist(totalNegTests, breaks = seq(0, max(totalNegTests)+ 5, by = 2))

dev.off()

totalNegTestsBounds <- rep(mean(totalNegTests), 2) + c(-5, 5)
totalPosTestsBounds <- rep(mean(totalPosTests), 2) + c(-5, 5)

totalTestsBounds <- rep(mean(totalTests)) + c(-5, 5)

simulatedSamples <- list()


for(i in 1:6){

  totalNegTestsSample <- -Inf
  totalPosTestsSample <- -Inf
  totalTestsSample <- -Inf

  while(!(totalNegTestsSample %<<% totalNegTestsBounds &
          totalPosTestsSample %<<% totalPosTestsBounds &
          totalTestsSample %<<% totalTestsBounds)){
    simulatedSample <- COVID19UK::testingNGF(simulatedEpidemic, obsParam = obsParam, obsWindow = obsWindow, PRINT = FALSE)
    totalTestsSample <- sum(simulatedSample$ntveSecondStageTests + simulatedSample$ptveSecondStageTests)
    totalPosTestsSample <- sum(simulatedSample$ptveSecondStageTests)
    totalNegTestsSample <- sum(simulatedSample$ntveSecondStageTests)
  }

  dailyNegTests <- rowSums(simulatedSample$ntveSecondStageTests)
  dailyPosTests <- rowSums(simulatedSample$ptveSecondStageTests)

  dailyTests <- dailyPosTests + dailyNegTests




  cumNegTests <- c(0, cumsum(dailyNegTests))
  cumPosTests <- c(0, cumsum(dailyPosTests))
  cumNoTests <- c(0, cumsum(dailyTests))

  simulatedSamples[[i]] <- list(dailyTests = dailyTests, dailyNegTests = dailyNegTests, dailyPosTests = dailyPosTests,
                                cumNoTests = cumNoTests, cumNegTests = cumNegTests, cumPosTests = cumPosTests,
                                sample = simulatedSample)
  rm(list = c("dailyTests", "dailyNegTests", "dailyPosTests",
              "cumNegTests", "cumPosTests", "cumNoTests", "simulatedSample"))

}


png(filename = "sampled_data_tests.png", width = 2*480, height = 6*480)

par(mfrow = c(6, 2))
days <- (obsWindow[1]:obsWindow[2])[-1]
timepoints <- obsWindow[1]:obsWindow[2]
for(i in 1:6){

  attach(simulatedSamples[[i]])
  plot(days, dailyTests, col = 'blue', ylim = c(0, max(cumNoTests) + 1), ylab = "Cumulative # Tests",
       xlab = "Day")
  points(days, dailyNegTests, col = 'red', pch = 2)
  points(days, dailyPosTests, col = 'green', pch = 3)



  plot(timepoints, cumNoTests, col = 'blue', type = 'l', ylim = c(0, max(cumNoTests) + 1), ylab = "Cumulative # Tests",
       xlab = "Time (Days)")

  lines(timepoints, cumNegTests, col = 'red', lty = 2)
  lines(timepoints, cumPosTests, col = 'green', lty = 3)

  legend("topleft", legend = c("Total", "Positive", "Negative"), col = c("blue", "green", "red"), lty = c(1, 3, 2))
  detach(simulatedSamples[[i]])
}




dev.off()





# ==== Set up experiment? ====

varNames <- c("beta_G", "beta_H", "gamma", "alpha", "pi", "psi")

for(i in 1:6){
  attach(simulatedSamples[[i]])
  householdModel <- COVID19UK::epiModel(simulator = simulator, obsModel = list(NGF = testingNGF, llh = testingLlh), simParam, obsParam,
                                        varNames = varNames, seed = NULL, conditional = F , simulatedSample = sample, simulatedEpidemic = simulatedEpidemic)

  llh_calc_beta <- householdModel$likelihoodCalc_function(simParam_names = c("beta_G", "beta_H"), obsParam = obsParam, noSims = 1)

  simulatedSamples[[i]]$MCsamples <- replicate(n = 1e5, llh_calc_beta(simParam[1:2]) , simplify = T)

  rm(list = c("householdModel", "llh_calc_beta"))

  detach(simulatedSamples[[i]])
}

# ==== Distribution of Likelihood Estimate for true parameters ====


png("observation_density_over_underlying_epidemic.png")
par(mfrow = c(2,3))

for(i in 1:6){
  attach(simulatedSamples[[i]])
  hist.default(MCsamples, main = "", xlab = "Log probability of sample being generated from a simulation")

  legend("topright", legend = c("prop. impossible epidemics:", sum(is.infinite(MCsamples))/length(MCsamples) ))

  detach(simulatedSamples[[i]])
}

dev.off()

fileConnection <- file("parameters.txt")

writeLines(c(paste0(c("Population size",
                      N,
                      "",
                      "Household Size",
                      N_h,
                      "",
                      "Observation Period",
                      obsWindow,
                      "",
                      "beta[G] (Global Infection Rate)",
                      simParam[1],
                      "",
                      "beta[H] (Within-household Infection Rate)",
                      simParam[2],
                      "",
                      "gamma (Removal rate)",
                      simParam[3],
                      "",
                      "alpha (new case ascertainment probability)",
                      obsParam[1],
                      "",
                      "pi (Test sensitivity)",
                      obsParam[2],
                      "",
                      "psi (Test specitivity)",
                      obsParam[3]))),
           fileConnection)

file.show("parameters.txt")

close.connection(fileConnection)


rm(list = ls()[!ls() %in% c("simParam", "simulatedEpidemic", "simulatedDSample", "N", "N_h", "obsParam", "householdModel"
                            )])

save.image(file = "Testing Scheme (alpha = 0.01).RData")

rm(list = ls())

