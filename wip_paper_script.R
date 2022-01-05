# WIP Paper Simulation Script


pop <- metaPopStruct(N = rep(4, 100), within = "H", between = 0.5)

sim <- HouseholdSIR(pop, startTime = 0, endTime = Inf, PRINT = FALSE)
I0 <- c(rep(1, 5), rep(0, 95))
param <- c(0.5, 0.75, 1)




# ====  Multiple Simulations ====

setwd("C:/Users/Work/Lancaster University/Grp-Josh PhD - General/Papers/figures")
jpeg("simulations.jpg", width = 2*480, height = 480)
par(mfrow = c(1,2))
noSims <- 1000
simulations <- replicate(n = noSims, sim(I0, param), simplify = F)

avgInfections <- c()
for(i in 1:20){
  avgInfections[i] <- mean(sapply(X = simulations, function(X) return(X$SIRsummary[i, 2])))
}
avgInfections

CImatrix <- matrix(nrow = 2, ncol = 20)
for(i in 1:20){
  CImatrix[, i] <- avgInfections[i] + c(1, -1)*1.96*sd(sapply(X = simulations, function(X) return(X$SIRsummary[i, 2])))
}

plot.new()
plot.window(xlim = c(0, 20), ylim = c(0, 50))

Axis(x = c(0,20), at = c(0, 5, 10, 15, 20), side = 1)
mtext("Day", side=1, line=3, cex.lab=1, adj = NA, las=1, col="black")

Axis(x = c(0,50), at = c(0, 10, 20, 30, 40, 50), side = 2)
mtext("Number Infectives", side=2, line=3, cex.lab=1, adj = 0.5, las=3, col="black")

lines(avgInfections, col = 'red')
lines(CImatrix[1, ], col = 'gray', lty = 2)
lines(CImatrix[2,], col = 'gray', lty = 2)

plot.new()
plot.window(xlim = c(0, 20), ylim = c(0, 50))

Axis(x = c(0,20), at = c(0, 5, 10, 15, 20), side = 1)
mtext("Day", side=1, line=3, cex.lab=1, adj = NA, las=1, col="black")

Axis(x = c(0, 50), at = c(0, 10, 20, 30, 40, 50), side = 2)
mtext("Number Infectives", side=2, line=3, cex.lab=1, adj = 0.5, las=3, col="black")

rand_sample <- sample(noSims, size = 10, replace = F)
lapply(X = rand_sample, function(X) lines(simulations[[X]]$SIRsummary[1:20, 2], col = "red"))
dev.off()
# ==== One Simulation ====

set.seed(1)
simulation <- sim(I0, param)

# = Summary Figures =

setwd("C:/Users/Work/Lancaster University/Grp-Josh PhD - General/Papers/figures")

jpeg("simulation.jpg", width = 2*480, height = 480)

par(mfrow = c(1, 2))

# Overall Epidemic Time Series
plot(simulation$SIRsummary[1:20,2], type = 'l', col = 'red', ylim = c(0, 50),
     xlab = "Day", ylab = "# Individuals")
#lines(simulation$SIRsummary[1:20, 1], col = 'blue')
#lines(simulation$SIRsummary[1:20, 3], col = 'grey')
#legend(x = "topleft", legend = c("S", "I", "R"), col = c("blue", "red", "grey"), lty = 1)



# Cumulative Infections
plot(cumsum(simulation$SIRsummary[1:20,2]), type = 'l', col = 'red', ylim = c(0, 400),
     xlab = "Day", ylab = "Total Infections")

dev.off()

# average household infection time series? Gives same curve of no infectives
#plot(rowMeans(simulation$Hstate$I), type = 'l', col = 'red', ylim = c(0,4), xlim = c(1,10))

