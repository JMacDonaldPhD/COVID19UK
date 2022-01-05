# Particle Filter


# 1. Importance Sampling step.
# 2. Starting from different initial conditions.
# 3. Calculating sequential likelihood.


# ==== Importance Sampling ====

IS <- function(particles, weights){
  resample <- sample(1:length(particles), size = length(particles))
  particles <- particles[resample]
  weights <- weights[resample]
  return(particles = particles, weights = weights)
}


# setwd(paste0(c(proj_wd, "/Experiments/Epidemic 9 (N = 400, scaled)"), collapse = ""))
# load("Epidemic 9.Rdata")
#
# # Simulate sample data
# obsParam <- c(0.05, 0.9, 0.98)
# obsWindow <- c(0, 1)
# simulatedSample <- testingNGF(simulatedEpidemic, obsParam = obsParam)
#
#
# model <- epiModel(simulator, obsModel = list(testingNGF, testingLlh), simParam = simParam, obsParam = obsParam,
#                   varNames = c("bG", "bH", "gamma", "alpha", "pi", "psi"),
#                   simulatedSample = simulatedSample, simulatedEpidemic = simulatedEpidemic)
# model$likelihoodCalc_function()
