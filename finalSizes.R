# Epidemic final size distribution with respect to R_t

N <- 500
model <- CT_Homo_SIR(N)

R <- c(0.5, 0.75, 1, 1.25, 1.5, 5)
finalSizes <- rep(list(NA), 3)
# Simulate the model many times and make a log of the final size and the log
# prob

K <- 1e3

I0 <- 10
gamma <- 1
beta <- gamma*R/(N - I0)

for(i in 1:length(R)){
  finalSizes[[i]] <- replicate(K, model(I0, c(beta, gamma))$finalSize)
}

save.image(file = "finalSizes.RData")




