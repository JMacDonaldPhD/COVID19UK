

# Create 'true' epidemic given a set of 'true' epidemic parameters


setwd("C:/Users/Josh/git/COVID19UK/Experiments/Epidemic 1 (N = 100)")


`%<<%` <- function(x, y){
  return(x > y[1] & x < y[2])
}


# ==== Defining Population Structure ====

# Size of each household
N_h <- 4

# Number of households
numberHouseholds <- 25

# Size of the population
N <- N_h*numberHouseholds

# Define initial state of individuals
indState0 <- rep(1, N)
indState0[(0:9)*4 + 1] <- 2 # Setting one person out of the first 10% of households to be infected

# dataframe structure for population
pop <- data.frame(ID = 1:N,
                  householdID = as.vector(sapply(X = 1:(N/N_h), function(X) return(rep(X, N_h)))),
                  state = indState0)

# Should probably bake this into HouseholdSIR function. Why do I need to do it? Some matrix sum breaks if I don't do this.
pop$state <- factor(pop$state, levels = 1:3, labels = c("S", "I", "R"))


# ==== Creating Simulator ====

# Ensures same 'true' epidemic and sample data is drawn each time. Not true as epidemic size is based on endTime
# So endTime = 5 will result in a different 'true' epidemic to endTime = 10
set.seed(1)
simulator <- COVID19UK::HouseholdSIR(pop, endTime = 100, PRINT = FALSE) # Default start time is startTime = 0

# ==== Simulating a representative Epidemic ====

simParam <- c(beta_G = .075, beta_H = .1, gamma = .1)
simulations <- replicate(n = 1e4, simulator(simParam), simplify = F)

finalSizes <- sapply(X = simulations, FUN = function(X) sum(X$SIRsummary[nrow(X$SIRsummary), 2:3]))
averageHouseholdFinalSizes <- sapply(X = simulations, FUN = function(X) mean(X$householdFinalSize))
propHousholdsFullyInfected <- sapply(X = simulations, FUN = function(X) sum(X$householdFinalSize == N_h)/numberHouseholds)


png(filename = "Summary Plots.png", width = 480*3)
par(mfrow = c(1, 3))

plot(density(finalSizes, from = 0, to = N), xlab = "Size of epidemic", main = "")
plot(density(averageHouseholdFinalSizes, from = 0, to = N_h), xlab = "Average size of household epidemic", main = "")

hist.default(propHousholdsFullyInfected, breaks = (0:10)/10, main = "", xlab = "Prop. of Households Fully Infected")
#plot(density(propHousholdsFullyInfected, from = 0, to = 1), main = "Prop. of Households Fully Infected (Day 10)")

dev.off()
epidemicSizeBounds <- rep(mean(finalSizes),2) + c(-5, 5)

avgHouseholdSizeBounds <- rep(mean(averageHouseholdFinalSizes), 2) + c(-.5, .5)


epidemicSize <- 0
avgHouseholdSize <- 0
while(!((epidemicSize %<<% epidemicSizeBounds) &
        (avgHouseholdSize %<<% avgHouseholdSizeBounds))){
  simulatedEpidemic <- simulator(simParam)
  epidemicSize <- sum(simulatedEpidemic$SIRsummary[nrow(simulatedEpidemic$SIRsummary), 2:3])
  avgHouseholdSize <- mean(simulatedEpidemic$householdFinalSize)
}


# Simulated Epidemic that sample data will be derived from

png(filename = "True_Epidemic_Curve.png")

par(mfrow = c(1, 1))
time <- 0:(nrow(simulatedEpidemic$SIRsummary) - 1)
plot(time, simulatedEpidemic$SIRsummary[, 1], type = 'l', col = 'blue', ylim = c(0, N))
lines(time, simulatedEpidemic$SIRsummary[, 2], col = 'red')
lines(time, simulatedEpidemic$SIRsummary[, 3], col = 'grey')

legend("topright", legend = c("S", "I", "R"), col = c("blue", "red", "grey"),
       lty = rep(1, 3))

dev.off()


fileConnection <- file("parameters.txt")

writeLines(c(paste0(c("Population size",
                      N,
                      "",
                      "Household Size",
                      N_h,
                      "",
                      "beta[G] (Global Infection Rate)",
                      simParam[1],
                      "",
                      "beta[H] (Within-household Infection Rate)",
                      simParam[2],
                      "",
                      "gamma (Removal rate)",
                      simParam[3],
                      ""))),
           fileConnection)

#file.show("parameters.txt")

close.connection(fileConnection)


rm(list = ls()[!ls() %in% c("simParam", "simulator", "simulatedEpidemic", "N", "N_h", "pop")])

save.image(file = "epidemic 1.RData")

rm(list = ls())

