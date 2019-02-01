#What is the significance of output?
#Why are we getting the D values as 1 every time?
#Binomial, Geometric, Exponential distribution are not continuous distributions, so how can we use ks.test on it?

z <<- NULL
general.gof.solver <- function(y, string.distribution, n00 = null) {
  nboot <- 100
  statvec <- NULL
  z <<- y
  
  if (string.distribution == "Binomial") {
    estimatedValueBinomial = mean(z)/n00
    stat0 <- ks.test(z, pbinom, n00, estimatedValueBinomial)$statistic
    for (i in 0:nboot) {
      zboot <- rbinom(length(z), n00, estimatedValueBinomial)
      phatstar <- mean(zboot)/n00
      statboot <- ks.test(zboot, pbinom, n00, phatstar)$statistic
      statvec <- c(statvec, statboot)
    }
    output <- sum(statvec>stat0)/nboot
    list(output)
  }
  else if (string.distribution == "Geometric") {
    pHat <- 1/mean(z)
    stat0 <- ks.test(z, pgeom, pHat)$statistic
    for (i in 0:nboot) {
      zboot <- rgeom(length(z), pHat)
      pHatStar <- 1/mean(zboot)
      statboot <- ks.test(zboot, pgeom, pHatStar)$statistic
      statvec <- c(statvec, statboot)
    }
    output <- sum(statvec>stat0)/nboot
    list(output)
  }
  else if (string.distribution == "Poisson") {
    pHat <- mean(z)
    stat0 <- ks.test(z, ppois, pHat)$statistic
    for (i in 0:nboot) {
      zboot <- rpois(length(z), pHat)
      pHatStar <- mean(zboot)
      statboot <- ks.test(zboot, ppois, pHatStar)$statistic
      statvec <- c(statvec, statboot)
    }
    output <- sum(statvec>stat0)/nboot
    list(output)
  }
  else if (string.distribution == "Uniform") {
    pHat <- max(z)
    stat0 <- ks.test(z, punif, 0, pHat)$statistic
    for (i in 0:nboot) {
      zboot <- runif(length(z), 0, pHat)
      phatstar <- max(zboot)
      statboot <- ks.test(zboot, punif, 0, phatstar)$statistic
      statvec <- c(statvec, statboot)
    }
    output <- sum(statvec>stat0)/nboot
    list(output)
  }
  else if (string.distribution == "Normal") {
    muHat <- mean(z)
    sigma_square_hat <- mean(z^2)
    stat0 <- ks.test(z, pnorm, muHat, sigma_square_hat)$statistic
    for (i in 0:nboot) {
      zboot <- rnorm(length(z), muHat, sigma_square_hat)
      muHatStar <- mean(zboot)
      sigma_square_hat_star <- mean(zboot^2)
      statboot <- ks.test(zboot, pnorm, muHatStar, sigma_square_hat_star)$statistic
      statvec <- c(statvec, statboot)
    }
    output <- sum(statvec>stat0)/nboot
    list(output)
  }
  else if (string.distribution == "Exponential") {
    pHat <- mean(z)
    stat0 <- ks.test(z, pexp, pHat)$statistic
    for (i in 0:nboot) {
      zboot <- rexp(length(z), pHat)
      pHatStar <- mean(zboot)
      statboot <- ks.test(zboot, pexp, pHatStar)$statistic
      statvec <- c(statvec, statboot)
    }
    output <- sum(statvec>stat0)/nboot
    list(output)
  }
  else if (string.distribution == "Gamma") {
    mle <- Beta_Gamma_MLE(z, "Gamma")
    stat0 <- ks.test(z, pgamma, mle[1], mle[2])$statistic
    for (i in 0:nboot) {
      zboot <- rgamma(length(z), mle[1], mle[2])
      mleStar <- Beta_Gamma_MLE(zboot, "Gamma")
      statboot <- ks.test(zboot, pgamma, mleStar[1], mleStar[2])$statistic
      statvec <- c(statvec, statboot)
    }
    output <- sum(statvec>stat0)/nboot
    list(output)
  }
  else if (string.distribution == "Beta") {
    mle <- Beta_Gamma_MLE(z, "Beta")
    stat0 <- ks.test(z, pbeta, mle[1], mle[2])$statistic
    for(i in 0:nboot) {
      zboot <- rbeta(length(z), mle[1], mle[2])
      mleStar <- Beta_Gamma_MLE(zboot, "Beta")
      statboot <- ks.test(zboot, pbeta, mleStar[1], mleStar[2])$statistic
      statvec <- c(statvec, statboot)
    }
    output <- sum(statvec>stat0)/nboot
    list(output)
  }
}











#MLE for beta and gamma
x <<- NULL
Beta_Gamma_MLE <- function(z, distribution) {
  if (distribution == "Beta") {
    x <<- z
    optim(theta <- c(2,2), beta_log, hessian=TRUE)$par
  }
  else if (distribution == "Gamma") {
    x <<- z
    optim(theta <- c(2,3), gamma_log)$par
  }
}
# alpha - theta 1
# beta - theta 2
beta_log <- function(theta) {
  value1 <- (length(x)*log(gamma(sum(theta[1], theta[2]))))
  value2 <- (length(x)*log(gamma(theta[1])))
  value3 <- (length(x)*log(gamma(theta[2])))
  logx <- log(x)
  log1x <- log(1-x)
  value4 <- (theta[1] - 1)*(sum(logx))
  value5 <- (theta[2] - 1)*(sum(log1x))
  sum(-(value1 - value2 - value3
        + value4
        + value5))
}

gamma_log <- function(theta) {
  if (theta[1] < 0 || theta[2] < 0 || theta[1] > 50) {
    return(NA)
  }
  value1 <- -(theta[1]*length(x)*log(theta[2]))
  value2 <- -(length(x)*log(gamma(theta[1])))
  value3 <- -(sum(x)/theta[2])
  value4 <- (theta[1]-1)*(sum(log(x)))
  sum(-(value1 + value2 + value3 + value4))
}