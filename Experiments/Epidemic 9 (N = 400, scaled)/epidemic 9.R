

# Create 'true' epidemic given a set of 'true' epidemic parameters

rm(list = ls())
# Test this on other machines
setwd(paste0(c(proj_wd, "/Experiments/Epidemic 9 (N = 400, scaled)"), collapse = ""))


`%<<%` <- function(x, y){
  return(x > y[1] & x < y[2])
}

# ==== Defining Population Structure ====

# Size of each household
N_h <- 4

# Number of households
numberHouseholds <- 100

# Size of the population
N <- N_h*numberHouseholds

# Define initial state of individuals
indState0 <- rep(1, N)
indState0[(0:(N/10 - 1))*4 + 1] <- 2 # Setting one person out of the first 10% of households to be infected

# dataframe structure for population
pop <- data.frame(ID = 1:N,
                  householdID = as.vector(sapply(X = 1:(N/N_h), function(X) return(rep(X, N_h)))),
                  state = indState0)

# Should probably bake this into HouseholdSIR function. Why do I need to do it? Some matrix sum breaks if I don't do this.
pop$state <- factor(pop$state, levels = 1:3, labels = c("S", "I", "R"))


states <- matrix(c(1, 0, 0,
                   0, 1, 0,
                   0, 0, 1), nrow = 3, ncol = 3, byrow = F)

initialState <- matrix(0, nrow = numberHouseholds, ncol = ncol(states))


#initialState[as.numeric(X[2]), ] <<- initialState[as.numeric(X[2]), ] + states[as.numeric(X[3]), ]
pop_copy <- pop
attr(pop_copy$state, "levels") <- NULL
attr(pop_copy$state, "class") <- NULL
apply(X = pop_copy, MARGIN = 1, function(X) initialState[X[2], ] <<- initialState[X[2], ] + states[X[3], ])

# ==== Creating Simulator ====

# Ensures same 'true' epidemic and sample data is drawn each time. Not true as epidemic size is based on endTime
# So endTime = 5 will result in a different 'true' epidemic to endTime = 10
set.seed(10)
simulator <- COVID19UK::HouseholdSIR(pop, startTime = 0, endTime = 100, PRINT = FALSE) # Default start time is startTime = 0

# ==== Simulating a representative Epidemic ====

simParam <- c(beta_G = 0.1*100/N, beta_H = 0.05, gamma = 1/14)
simulations <- replicate(n = 1e4, simulator(simParam), simplify = F)

finalSizes <- sapply(X = simulations, FUN = function(X) sum(X$SIRsummary[nrow(X$SIRsummary), 2:3]))
averageHouseholdFinalSizes <- sapply(X = simulations, FUN = function(X) mean(X$householdFinalSize))
propHousholdsFullyInfected <- sapply(X = simulations, FUN = function(X) sum(X$householdFinalSize == N_h)/numberHouseholds)
noInfectionsInFirstDay <- sapply(X = simulations, FUN = function(X) diff(c(X$SIRsummary[2, 1], X$SIRsummary[1, 1])))
noInfectionsInFirstDay2 <- sapply(X = simulations, FUN = function(X) sum(X$Infections[1, ]))

png(filename = "Summary Plots.png", width = 480*2, height =  480*2)
par(mfrow = c(2, 2))

plot(density(finalSizes, from = 0, to = N), xlab = "Size of epidemic", main = "")
plot(density(averageHouseholdFinalSizes, from = 0, to = N_h), xlab = "Average size of household epidemic", main = "")

hist.default(propHousholdsFullyInfected, breaks = (0:10)/10, main = "", xlab = "Prop. of Households Fully Infected")
hist.default(noInfectionsInFirstDay, breaks = 0:max(noInfectionsInFirstDay), main = "Number of Infections on Day 1",
             right = F)
dev.off()



epidemicSizeBounds <- rep(mean(finalSizes),2) + c(-5, 5)

avgHouseholdSizeBounds <- rep(mean(averageHouseholdFinalSizes), 2) + c(-.5, .5)
noInfectionsInFirstDayBounds <- rep(mean(noInfectionsInFirstDay), 2) + c(-2, +2)



epidemicSize <- -Inf
avgHouseholdSize <- -Inf
noInfectionsInFirstDay <- -Inf
while(!((epidemicSize %<<% epidemicSizeBounds) &
        (avgHouseholdSize %<<% avgHouseholdSizeBounds) &
        (noInfectionsInFirstDay %<<% noInfectionsInFirstDayBounds))){
  simulatedEpidemic <- simulator(simParam)
  epidemicSize <- sum(simulatedEpidemic$SIRsummary[nrow(simulatedEpidemic$SIRsummary), 2:3])
  avgHouseholdSize <- mean(simulatedEpidemic$householdFinalSize)
  noInfectionsInFirstDay <- sum(simulatedEpidemic$Infections[1, ])
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


rm(list = ls()[!ls() %in% c("simParam", "simulator", "simulatedEpidemic", "N", "N_h", "pop", "initialState")])

save.image(file = "epidemic 9.RData")

rm(list = ls())

