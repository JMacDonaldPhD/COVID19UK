
setwd("C:/Users/Josh/git/COVID19UK/Experiments/Epidemic 1 (N = 100)/")

load("epidemic 1.RData")



setwd("C:/Users/Josh/git/COVID19UK/Experiments/Epidemic 1 (N = 100)/Testing Scheme 1/")

`%<<%` <- function(x, y){
  return(x > y[1] & x < y[2])
}



# ==== Simulating representative testing Data ====



set.seed(1)
obsParam <- c(alpha = .1, pi = .9, psi = .6)
obsWindow <- c(1, 2)
simulator <- COVID19UK::HouseholdSIR(pop, endTime = 2, PRINT = FALSE) # Default start time is startTime = 0


samples <- replicate(n = 1e4, COVID19UK::testingNGF(simulatedEpidemic, obsParam = obsParam, obsWindow = obsWindow, PRINT = FALSE), simplify = F)

totalTests <- sapply(X = samples, FUN = function(X) sum(X$ntveSecondStageTests + X$ptveSecondStageTests))
totalPosTests <- sapply(X = samples, FUN = function(X) sum(X$ptveSecondStageTests))
totalNegTests <- sapply(X = samples, FUN = function(X) sum(X$ntveSecondStageTests))

png(filename = "testing_summary_densities.png", width = 480*3)


par(mfrow = c(1, 3))

plot(density(totalTests, from = 0, to = N), xlab = "Total Number of Tests", main = "", xlim = c(0, max(totalTests) + 5))
plot(density(totalPosTests, from = 0, to = N), xlab = "Total Number of Positive Tests", main = "", xlim = c(0, max(totalPosTests) + 5))
plot(density(totalNegTests, from = 0, to = N), xlab = "Total Number of Negative Tests", main = "", xlim = c(0, max(totalNegTests) + 5))


dev.off()

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

cumNoTests <- cumsum(dailyPosTests + dailyNegTests)

png(filename = "sampled_data_tests.png", width = 480)

par(mfrow = c(1, 1))
plot(cumNoTests, col = 'blue', type = 'l', ylim = c(0, max(cumNoTests) + 1), ylab = "Cumulative # Tests",
     xlab = "Time (Days)")

lines(cumsum(dailyNegTests), col = 'red', lty = 2)
lines(cumsum(dailyPosTests), col = 'green', lty = 3)

legend("topleft", legend = c("Total", "Positive", "Negative"), col = c("blue", "green", "red"), lty = c(1, 3, 2))

dev.off()
# ==== Set up experiment? ====

varNames <- c("beta_G", "beta_H", "gamma", "alpha", "pi", "psi")

householdModel <- COVID19UK::epiModel(simulator = simulator, obsModel = list(NGF = testingNGF, llh = testingLlh), simParam, obsParam,
                                      varNames = varNames, seed = NULL, conditional = F , simulatedSample = simulatedSample, simulatedEpidemic = simulatedEpidemic)




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

