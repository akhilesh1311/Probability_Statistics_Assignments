Bayesian_Confidence_Interval <- function(z, string.distribution, pos_dist_param, alpha) {
  n <- length(z)
  # pass Beta function for Bernoulli
  if (string.distribution == "Bernoulli") {
    alphaa <- sum(z) + pos_dist_param[1]
    betaa <- n - sum(z) + pos_dist_param[2]
    sd <- sqrt(alphaa*betaa/(((alphaa+betaa)^2)*(alphaa+betaa+1)*n))
    mu <- alphaa/(alphaa + betaa)
    C1 <- mu - 1.96*sd
    C2 <- mu + 1.96*sd
    list(C1, C2)
  }
  # pass Gamma function for Exponential
  else if (string.distribution == "Exponential") {
    alphaa <- n + pos_dist_param[1]
    betaa <- sum(z) + pos_dist_param[2]
    sd <- sqrt(alphaa*(betaa^2)/n)
    mu <- (alphaa*betaa)
    C1 <- mu - 1.96*sd
    C2 <- mu + 1.96*sd
    list(C1, C2)
  }
  # pass Gamma function for Poisson
  else if (string.distribution == "Poisson") {
    alphaa <- pos_dist_param[1] + sum(z)
    betaa <- pos_dist_param[2] + n
    sd <- sqrt(alphaa*(betaa^2)/n)
    mu <- (alphaa*betaa)
    C1 <- mu - 1.96*sd
    C2 <- mu + 1.96*sd
    list(C1, C2)
  }
  
  #pos_dist_param[1] <- mu
  #pos_dist_param[2] <- tau
  #pos_dist_param[3] <- alpha
  #pos_dist_param[4] <- beta
  else if (string.distribution == "Normal") {
    #confidence interval of precision
    alphaa <- pos_dist_param[3] + (n/2)
    betaa <- pos_dist_param[4] + (1/2)(sum((z-mean(z))^2)) + (pos_dist_param[2]*n*(mean(z)-pos_dist_param[1])^2)/(2*(pos_dist_param[2]+n))
    sd <- sqrt(alphaa*(betaa^2)/n)
    mu <- (alphaa*betaa)
    C1 <- mu - 1.96*sd
    C2 <- mu + 1.96*sd
    list(C1, C2)

    #confidence interval for mean of Normal
    
  }
}