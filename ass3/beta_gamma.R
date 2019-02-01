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