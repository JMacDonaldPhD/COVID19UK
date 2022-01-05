model <- CT_Homo_SIR(20, endTime = 5)
beta <- 1
gamma <- 10
sim <- model(2, c(beta, gamma))

X0 <- sim$SIRsummary[1, 1:3]
Xt <- sim$SIRsummary[11, 1:3]
t <- sim$SIRsummary[11, 4]

nrow(paths(X0, Xt))

debug(probPath)
pSim(X0, Xt, beta, gamma, t)


# z is the sum of N independent but not necessarily identically distributed
# exponential random variables.
# Lambda is a vector of size N of the rate parameter for each exponential
# random variable

# DONT THINK THE EXPRESSION I HAVE IS NORMALISED!
dz <- function(z, lambda){

  N <- length(lambda)
  equal_lambdaN <- lambda[-N] == lambda[N]
  num_equal_lambdaN <- sum(equal_lambdaN)

  if(num_equal_lambdaN == N - 1){
    return(dgamma(z, rate = lambda[1], shape = N))
  } else if(num_equal_lambdaN == 0){
    lambdaN_minus_lambda <- lambda[N] - lambda[-N]
    p1 <- prod(lambda)
    p3 <- prod((1/(lambdaN_minus_lambda))*(exp(z*(lambdaN_minus_lambda)) - 1))
    return(p1*p3)
  } else{

    p1 <- prod(lambda)
    p2 <- z^(num_equal_lambdaN)
    lambdaN_minus_lambda <- lambda[N] - lambda[c(!equal_lambdaN, F)]
    if(length(lambdaN_minus_lambda) + num_equal_lambdaN != N){
      stop("Something Went Wrong!")
    }
    p3 <- prod((1/(lambdaN_minus_lambda))*(exp(z*(lambdaN_minus_lambda)) - 1))
    return(p1*p2*p3)
  }

}

dz(100, lambda = c(1, 2))
z.seq <- 0
plot()
