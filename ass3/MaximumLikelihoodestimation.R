#MAXIMUM LIKELIHOOD ESTIMATOR

MaximumLikelihood <- function(z,distribution){
  if (distribution == "Binomial") {
    mu_1 <- mean(z)
    mu_2 <- mean((z-mu_1)^2)
    p <- 1-((mu_2)/(mu_1))
    n <- mu_1/p
    phat <- (1/length(z))*(n*p)
    list(p = p, n = n, mle = phat)
  }
  else if (distribution == "Geometric") {
    mu_1 <- mean(z)
    p <- 1/mu_1
    phat <- 1/mu_1
    list(p = p, mle = phat)
  }
  else if (distribution == "Poisson") {
    mu_1 <- mean(z)
    lambda <- mu_1
    lambdahat <- lambda
    list(lambda = lambda, mle = lambdahat)
  }
  else if (distribution == "Uniform") {
    mu_1 <- mean(z)
    mu_2 <- mean(z^2)
    a <- mu_1 - sqrt(3*(mu_2 - (mu_1)^2))
    b <- mu_1 + sqrt(3*(mu_2 - (mu_1)^2))
    thetahat <- max(z)
    list(a=a,b=b, mle=thetahat)
  }
  else if (distribution == "Normal") {
    mu <- mean(z)
    n<-length(z)
    sigma_square <- var(z)*(n/(n-1))
    muhat <- (sum(z))/n
    sigma_squarehat <- var(z)*(n/(n-1))
    list(mu = mu, sigma_square = sigma_square,mle1=muhat,mle2=sigma_squarehat)
  }
  else if (distribution == "Exponential") {
    mu_1 <- mean(z)
    beta <- mu_1
    thetahat <- beta
    list(beta = beta, mle = beta)
  }
  else if (distribution == "Gamma") {
    mu_1 <- mean(z)
    mu_2 <- mean(z^2)
    alpha <- (mu_1^2)/(mu_2 - (mu_1)^2)
    beta <- (mu_2 - (mu_1)^2)/(mu_1)
    list(alpha = alpha, beta = beta)
  }
  else if (distribution == "Beta") {
    alpha.start <- 1
    beta.start <- 1
    out <- beta_solver1(z,alpha.start, beta.start, 1000)
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

# Could not find maximum likelihood estimate for beta and gamma distributions. The beta distribution maximum likelihood estimation was going into double derivates and
# I did not know whether the maximum likelihood approach in case of beta and gamma are right.

OUTPUT

MaximumLikelihood(z,"Binomial")
$p
[1] -6.219015

$n
[1] -4.092962

$mle
[1] 0.2545419

MaximumLikelihood(z,"Geometric")
$p
[1] 0.03928626

$mle
[1] 0.03928626

MaximumLikelihood(z,"Poisson")
$lambda
[1] 25.45419

$mle
[1] 25.45419

MaximumLikelihood(z,"Uniform")
$a
[1] 1.975211

$b
[1] 48.93318

$mle
[1] 49.76183

MaximumLikelihood(z,"Normal")
$mu
[1] 25.45419

$sigma_square
[1] 187.4852

$mle1
[1] 25.45419

$mle2
[1] 187.4852

MaximumLikelihood(z,"Exponential")
$beta
[1] 25.45419

$mle
[1] 25.45419