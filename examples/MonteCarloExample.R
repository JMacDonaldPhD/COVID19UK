#' Define likelihood estimatinon in terms of Monto Carlo Problem
setwd(paste0(c(proj_wd, "/Experiments/Epidemic 9 (N = 400, scaled)"), collapse = ""))
load("Epidemic 9.Rdata")

obsParam <- c(0.05, 0.9, 0.98)
obsWindow <- c(0, 1)
simulatedSample <- testingNGF(simulatedEpidemic, obsParam = obsParam)


model <- epiModel(simulator, obsModel = list(testingNGF, testingLlh), simParam = simParam, obsParam = obsParam,
                  varNames = c("bG", "bH", "gamma", "alpha", "pi", "psi"),
                  simulatedSample = simulatedSample, simulatedEpidemic = simulatedEpidemic)

llh_calc <- model$likelihoodCalc_function(simParam_names = c("bG", "bH"), obsParam = obsParam, noSims = 1)
simulator(simParam)
