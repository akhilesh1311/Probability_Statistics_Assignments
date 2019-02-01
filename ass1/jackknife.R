library(e1071)

Jackknife <- function(v1, nboot, statfunc, alpha) {
	n1 <- length(v1)
	jackvec <- NULL
	mu0 <- statfunc(v1)
	for(i in 1:n1) {
		mua <- statfunc(v1[-i])
		jackvec <- c(jackvec, n1*mu0 - (n1-1)*mua)
		}
	jackbias <- mean(jackvec) - mu0
	jacksd <- sd(jackvec)
	
	sample_variance <- var(jackvec)
	confidence_interval <- vector(length=2)
	confidence_interval[1] <- mean(jackvec) - qnorm(1-alpha/2)*sqrt(sample_variance/n1)
	confidence_interval[2] <- mean(jackvec) + qnorm(1-alpha/2)*sqrt(sample_variance/n1)
	
#	list(mu0=mu0, jackbias=jackbias, jacksd=jacksd, confidence_interval)
#	return(answer)
#	list(answer)
	return(jacksd)
}
