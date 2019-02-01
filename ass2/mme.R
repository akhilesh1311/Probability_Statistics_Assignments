MethodsOfMomentsEstimator <- function(z, distribution) {
  if (distribution == "Normal") {
    mu <- mean(z)
    n <- length(z)
    sigma_square <- var(z)*(n/(n-1))
    list(mu = mu, sigma_square = sigma_square)
  }
  else if (distribution == "Bernoulli") {
    mu <- mean(z)
    list(mu = mu)
  }
  else if (distribution == "Point") {
    a <- mean(z)
    list(a=a)
  }
  else if (distribution == "Binomial") {
    alpha_one <- mean(z)
    alpha_two <- mean(z^2)
    p <- (alpha_one - alpha_two + alpha_one^2)/alpha_one
    n <- (alpha_one^2)/(alpha_one - alpha_two + alpha_one^2)
    list(p = p, n = n)
  }
  else if (distribution == "Geometric") {
    alpha_one <- mean(z)
    p <- 1/alpha_one
    list(p = p)
  }
  else if (distribution == "Poisson") {
    alpha_one <- mean(z)
    lambda <- alpha_one
    list(lambda = lambda)
  }
  else if (distribution == "Uniform") {
    alpha_one <- mean(z)
    alpha_two <- mean(z^2)
    a <- alpha_one - sqrt(3*(alpha_two - alpha_one^2))
    b <- alpha_one + sqrt(3*(alpha_two - alpha_one^2))
    list(a=a,b=b)
  }
  else if (distribution == "Exponential") {
    alpha_one <- mean(z)
    beta <- alpha_one
    list(beta = beta)
  }
  else if (distribution == "Gamma") {
    alpha_one <- mean(z)
    alpha_two <- mean(z^2)
    alpha <- (alpha_one^2)/(alpha_two - alpha_one^2)
    beta <- (alpha_two - alpha_one^2)/(alpha_one)
    list(alpha = alpha, beta = beta)
  }
  else if (distribution == "Beta") {
    alpha.start <- 1
    beta.start <- 1
    out <- beta_solver1(z,alpha.start, beta.start, 1000)
  }
  else if (distribution == "T") {
    alpha_two <- mean(z^2)
    v <- 2*alpha_two/(alpha_two - 1)
    list(v = v)
  }
  else if (distribution == "ChiSquare") {
    alpha_one <- mean(z)
    p <- alpha_one
    list(p = p)
  }
  else if (distribution == "Multivariate Normal") {
    alpha_one <- mean(z)
    alpha_two <- mean(z^2)
    mu <- alpha_one
    summation <- alpha_two - alpha_one^2
    list(mu = my, summation = summation)
  }
}

beta_solver <- function(x,x2, alpha0, beta0) {
  res1 <- x-(alpha0/(alpha0 + beta0))
  res2 <- x2-((alpha0+1) * (alpha0)/((alpha0+beta0)*(alpha0+beta0+1)))
  c(res1, res2)
}

beta_solver1 <- function(z, alpha0, beta0, ntry = 100) {
  x <- mean(z)
  y <- mean(z^2)
  ss1 <- sum(beta_solver(x,y,alpha0, beta0)^2)
  ss1 <- ss1
  for (i in 1:ntry) {
    alpha1 <- (1/2 + runif(1))*alpha0
    beta1 <- (1/2 + runif(1))*beta0
    ss2 <- sum(beta_solver(x,y,alpha1,beta1)^2)
    if (ss2 < ss1) {
      print(c(ss1, ss2, i))
      alpha0 <- alpha1
      beta0 <- beta1
      ss1 <- ss2
    }
  }
}