SD_by_n <- function(v1, nboot, alpha) {
  n1 <- length(v1)
  jackvec <- NULL
  sd_hat <- sqrt(sd(v1)^2*((n1-1)/n1))
  for (i in 1:n1) {
    sd_jack <- sqrt(sd(v1[-i])^2*((n1-1)/n1))
    jackvec <- c(jackvec, n1*sd_hat - (n1-1)*sd_jack)
  }
  
  sample_variance <- var(jackvec)*(n1-1)/n1
  confidence_interval <- vector(length=2)
  confidence_interval[1] <- mean(jackvec) - qnorm(1-alpha/2)*sqrt(sample_variance/n1)
  confidence_interval[2] <- mean(jackvec) + qnorm(1-alpha/2)*sqrt(sample_variance/n1)
  
  for(i in 1:nboot) {
    x1 <- sample(v1, n1, replace = TRUE, prob = NULL)
    sd_boot <- c(sd_boot, sqrt(sd(x1)^2*((n1-1)/n1)))
  }
}