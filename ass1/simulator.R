Simulator <- function(n, nboot, alpha) {
	v1 <- rlnorm(n, meanlog = 0, sdlog = 1)
	
	#central limit theorem based confidence interval
	confidence_interval <- vector(length=2)
	confidence_interval[1] <- mean(v1) - qnorm(1-alpha/2)*sd(v1)
	confidence_interval[2] <- mean(v1) + qnorm(1-alpha/2)*sd(v1)
	
	answer_bootstrap <- Bootstrap(v1, nboot, mean, alpha)
	answer_jackknife <- Jackknife(v1, nboot, mean, alpha)
	
	list(answer_bootstrap=answer_bootstrap, answer_jackknife=answer_jackknife)
}
