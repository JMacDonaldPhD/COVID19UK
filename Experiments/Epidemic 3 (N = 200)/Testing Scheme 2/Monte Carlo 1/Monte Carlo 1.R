
# Monte Carlo Estimation of likelihood
rm(list = ls())
setwd("C:/Users/Josh/git/COVID19UK/Experiments/Epidemic 2 (N = 100)/Testing Scheme 1/")
load("Testing Scheme 1.RData")

setwd("C:/Users/Josh/git/COVID19UK/Experiments/Epidemic 2 (N = 100)/Testing Scheme 1/Monte Carlo 1")

 # ==== Estimate and plot Likelihood ====
noSims <- 1
hh_llhCalc <- householdModel$likelihoodCalc_function(simParam_names = c("beta_G", "beta_H"), obsParam = obsParam, noSims = noSims)

beta_G.seq <- seq(0, .5, length = 500)
beta_H.seq <- seq(0, .5, length = 500)
length_beta_G.seq <- length(beta_G.seq)
betaGrid <- expand.grid(beta_G.seq, beta_H.seq)
start <- as.numeric(Sys.time())
testCalc <- apply(betaGrid[(1:length_beta_G.seq)*length_beta_G.seq,], MARGIN = 1, FUN = hh_llhCalc)
timeTaken <- as.numeric(Sys.time()) - start

print(paste0("Time started: ", lubridate::now(tzone = "GMT"), collapse = ""))
print(c("Estimated Time to calculate likelihood values (mins.): ", timeTaken*(nrow(betaGrid)/length_beta_G.seq/60)))

start <- as.numeric(Sys.time())
betaGrid$llh <- apply(betaGrid[, 1:2], MARGIN = 1, FUN = hh_llhCalc)
timeTaken <- as.numeric(Sys.time()) - start

print(c("Actual time taken (mins.)", timeTaken/60))


z <- matrix(betaGrid$llh, nrow = length(beta_G.seq), ncol = length(beta_H.seq))

beta_G.profile <- log(colMeans(exp(z)))
beta_H.profile <- log(rowMeans(exp(z)))
# Small beta values not return small likelihood values? why?

png(filename = "Log_Profile_Likelihood_Estimates (0-1.5).png", width = 480*2)

par(mfrow = c(1, 2))

plot(beta_G.seq, beta_G.profile, type = 'l')
abline(v = simParam[1], lty = 2, col = "red")
legend("topleft", legend =  "True value", col = "red", lty = 2)


plot(beta_H.seq, beta_H.profile, type = 'l')
abline(v = simParam[2], lty = 2, col = "red")


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
                      "",
                      "alpha (new case ascertainment probability)",
                      obsParam[1],
                      "",
                      "pi (Test sensitivity)",
                      obsParam[2],
                      "",
                      "psi (Test specitivity)",
                      obsParam[3],
                      "",
                      "Monte Carlo Samples:",
                      noSims))),
           fileConnection)

file.show("parameters.txt")

close.connection(fileConnection)

rm(list = ls()[!ls() %in% c("betaGrid")])

save.image("betaGrid.RData")
load("betaGrid.RData")

rm(list = ls())

