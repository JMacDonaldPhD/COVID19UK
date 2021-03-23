rm(list = ls())
setwd(paste0(c(proj_wd, "/Experiments/Epidemic 8 (N = 400, scaled)/"), collapse = ""))

load("epidemic 8.RData")


setwd(paste0(c(proj_wd, "/Experiments/Epidemic 8 (N = 400, scaled)/Testing (Straight sims vs. Permutations)"), collapse = ""))
`%<<%` <- function(x, y){
  return(x > y[1] & x < y[2])
}





# ==== Simulating representative testing Data ====

obsParam <- c(alpha = 0.1, pi = .9, psi = .98)
obsWindow <- c(0, 1)
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

png(filename = "sampled_data_tests.png", width = 480)

par(mfrow = c(1, 2))

days <- (obsWindow[1]:obsWindow[2])[-1]
timepoints <- obsWindow[1]:obsWindow[2]


plot(days, dailyTests, col = 'blue', ylim = c(0, max(cumNoTests) + 1), ylab = "Cumulative # Tests",
     xlab = "Day")
points(days, dailyNegTests, col = 'red', pch = 2)
points(days, dailyPosTests, col = 'green', pch = 3)



plot(timepoints, cumNoTests, col = 'blue', type = 'l', ylim = c(0, max(cumNoTests) + 1), ylab = "Cumulative # Tests",
     xlab = "Time (Days)")

lines(timepoints, cumNegTests, col = 'red', lty = 2)
lines(timepoints, cumPosTests, col = 'green', lty = 3)

legend("topleft", legend = c("Total", "Positive", "Negative"), col = c("blue", "green", "red"), lty = c(1, 3, 2))

dev.off()


# ==== benchmark one-time llh calc. vs. household permutation ====

simulations <- replicate(n = 100, simulator(simParam), simplify = F)

perm_testingLlh <- testingLlh_exchangeability(initialState, maxPermutations = 20)

perm_testingLlh1 <- testingLlh_exchangeability(initialState, maxPermutations = 10)

perm_testingLlh2 <- testingLlh_exchangeability(initialState, maxPermutations = 5)

perm_testingLlh3 <- testingLlh_exchangeability(initialState, maxPermutations = 1)

perm_testingLlh4 <- testingLlh_exchangeability(initialState, maxPermutations = 0)

perm_testingLlh5 <- testingLlh_exchangeability(initialState, maxPermutations = 1000)


perm_testingLlh2(simulations[[1]], sampleData = simulatedSample, obsParam)




microbenchmark::microbenchmark(lapply(simulations, FUN = testingLlh, sampleData = simulatedSample, obsParam = obsParam),
                               lapply(simulations, FUN = perm_testingLlh, sampleData = simulatedSample, obsParam = obsParam),
                               lapply(simulations, FUN = perm_testingLlh1, sampleData = simulatedSample, obsParam = obsParam),
                               lapply(simulations, FUN = perm_testingLlh2, sampleData = simulatedSample, obsParam = obsParam),
                               lapply(simulations, FUN = perm_testingLlh3, sampleData = simulatedSample, obsParam = obsParam),
                               lapply(simulations, FUN = perm_testingLlh4, sampleData = simulatedSample, obsParam = obsParam),
                               lapply(simulations, FUN = perm_testingLlh5, sampleData = simulatedSample, obsParam = obsParam),
                               times = 10L)

# With a maximum of 10 permutations, llh_calc takes on average 20 times longer.  It does a maximum of 10 extra calls to testingLlh, which should account for about 20 seconds


# ==== Straight Simulations for likelihood calculation ====

varNames <- c("beta_G", "beta_H", "gamma", "alpha", "pi", "psi")


llh_func <- testingLlh

householdModel <- COVID19UK::epiModel(simulator = simulator, obsModel = list(NGF = testingNGF, llh = llh_func), simParam, obsParam,
                                      varNames = varNames, seed = NULL, conditional = F , simulatedSample = simulatedSample, simulatedEpidemic = simulatedEpidemic)


# ==== Distribution of Likelihood Estimate for true parameters (straight simulations) ====
llh_calc_beta <- householdModel$likelihoodCalc_function(simParam_names = c("beta_G", "beta_H"), obsParam = obsParam, noSims = 1)

samples <- replicate(n = 1e5, llh_calc_beta(simParam[1:2]) , simplify = T)


# ==== Permutated Simulations for likelihood calculation ====

varNames <- c("beta_G", "beta_H", "gamma", "alpha", "pi", "psi")


llh_func <- testingLlh_exchangeability(initialState, maxPermutations = 1)

householdModel <- COVID19UK::epiModel(simulator = simulator, obsModel = list(NGF = testingNGF, llh = llh_func), simParam, obsParam,
                                      varNames = varNames, seed = NULL, conditional = F , simulatedSample = simulatedSample, simulatedEpidemic = simulatedEpidemic)


# ==== Distribution of Likelihood Estimate for true parameters (straight simulations) ====
llh_calc_beta <- householdModel$likelihoodCalc_function(simParam_names = c("beta_G", "beta_H"), obsParam = obsParam, noSims = 1)

samples2 <- replicate(n = 1e5, llh_calc_beta(simParam[1:2]) , simplify = T)









png("observation_density_over_underlying_epidemic.png")
par(mfrow = c(1,2))
hist.default(samples, main = "", xlab = "Log probability of sample being generated from a simulation")

legend("topright", legend = c("prop. impossible epidemics:", sum(is.infinite(samples))/length(samples) ))


hist.default(samples2, main = "", xlab = "Log probability of sample being generated from a simulation")

legend("topright", legend = c("prop. impossible epidemics:", sum(is.infinite(samples))/length(samples) ))


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

save.image(file = "Testing Scheme 1.RData")

rm(list = ls())

