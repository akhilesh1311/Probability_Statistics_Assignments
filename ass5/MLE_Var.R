#n00 is important only for binomial distribution, neglect it for the rest of the distributions
MLE_Variance <- function(z, distribution, n00) {
  n <- length(z)
  p <- MLE(z, distribution)
  if (distribution == "Binomial") {
    fisher = n00/((1-p)*p)
    fisher_n = n*fisher
    var = 1/fisher_n
    list(var = var)
  }
  else if (distribution == "Geometric") {
    fisher = 1/((p^2)*(1-p))
    fisher_n = n*fisher
    var = 1/fisher_n
    list(var = var)
  }
  else if (distribution == "Poisson") {
    fisher = 1/p
    fisher_n = n*fisher
    var = 1/fisher_n
    list(var = var)
  }
  else if (distribution == "Normal") {
    fisher_n <- matrix(c(-n/p[2], 0, 0, -2*n/p[2]), nrow = 2, ncol = 2)
    var = solve(-1*fisher_n)
    list(var = var)
  }
  else if (distribution == "Exponential") {
    fisher = 1/(p[1]^2)
    fisher_n = n*fisher
    var = 1/fisher_n
    list(var = var)
  }
  else if (distribution == "Gamma") {
    beta <- 1/p[2]
    value_repeat <- -n/beta
    fisher_n <- matrix(c(-n*trigamma(p[1]), value_repeat, value_repeat, -(p[1]*n/(beta^2))), nrow = 2, ncol = 2)
    var = solve(-1*fisher_n)
    list(var = var)
  }
  else if (distribution == "Beta") {
    trigamma_sum_alpha_beta <- trigamma(p[1] + p[2])
    trigamma_alpha <- trigamma(p[1])
    trigamma_beta <- trigamma(p[2])
    fisher_n <- matrix(c(n*trigamma_sum_alpha_beta - n*trigamma_alpha, n*trigamma_sum_alpha_beta, n*trigamma_sum_alpha_beta, 
                         n*trigamma_sum_alpha_beta - n*trigamma_beta), nrow = 2, ncol = 2)
    var <- solve(-1*fisher_n)
    list(var = var)
  }
}


x <<- NULL
MLE <- function(z, distribution) {
  if (distribution == "Binomial") {
    mu_1 <- mean(z)
    mu_2 <- mean((z-mu_1)^2)
    p <- 1-((mu_2)/(mu_1))
    n <- mu_1/p
    phat <- (1/length(z))*(n*p)
    return(p)
  }
  #done
  else if (distribution == "Geometric") {
    mu_1 <- mean(z)
    p <- 1/mu_1
    phat <- 1/mu_1
    return(p)
  }
  #done
  else if (distribution == "Poisson") {
    mu_1 <- mean(z)
    lambda <- mu_1
    lambdahat <- lambda
    return(lambdahat)
  }
  else if (distribution == "Uniform") {
    mu_1 <- mean(z)
    mu_2 <- mean(z^2)
    a <- mu_1 - sqrt(3*(mu_2 - (mu_1)^2))
    b <- mu_1 + sqrt(3*(mu_2 - (mu_1)^2))
    thetahat <- max(z)
    return(thetahat)
  }
  else if (distribution == "Normal") {
    mu <- mean(z)
    n<-length(z)
    sigma_square <- var(z)*(n/(n-1))
    muhat <- (sum(z))/n
    sigma_squarehat <- var(z)*(n/(n-1))
    return(c(muhat, sigma_squarehat))
  }
  #done
  else if (distribution == "Exponential") {
    mu_1 <- mean(z)
    beta <- mu_1
    thetahat <- beta
    return(thetahat)
  }
  else if (distribution == "Beta") {
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