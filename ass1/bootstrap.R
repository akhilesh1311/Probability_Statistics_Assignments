library(e1071)

Bootstrap <- function(v1, nboot, statfunc, alpha) {
	n1 <- length(v1)
	statfunc_hat <- statfunc(v1)
	statfunc_boot <- NULL
	answer <- NULL

	for(i in 1:nboot) {
		x1 <- sample(v1, n1, replace = TRUE, prob = NULL)
		statfunc_boot <- c(statfunc_boot, statfunc(x1))
	}
	bias <- mean(statfunc_boot) - statfunc_hat
	boot_sd <- sd(statfunc_boot)
	se <- sqrt(var(statfunc_boot))

	lq <- quantile(statfunc_boot, alpha/2)
	uq <- quantile(statfunc_boot, 1-alpha/2)

	#The first index in any of the following vector denotes the lower limit, and the second
	#denotes the upper limit
	bootstrap_confidence <- vector(length = 2)
	bootstrap_confidence[1] <- statfunc_hat - qnorm(1-alpha/2)*se
	bootstrap_confidence[2] <- statfunc_hat + qnorm(1-alpha/2)*se

	normal_confidence <- vector(length = 2)
	normal_confidence[1] <- statfunc_hat - qnorm(1-alpha/2)*(sd(v1)/sqrt(n1))
	normal_confidence[2] <- statfunc_hat + qnorm(1-alpha/2)*(sd(v1)/sqrt(n1))

	percentile_confidence <- vector(length = 2)
	percentile_confidence[1] <- lq
	percentile_confidence[2] <- uq

	pivotal_confidence <- vector(length = 2)
	pivotal_confidence[1] <- 2*statfunc_hat - uq
	pivotal_confidence[2] <- 2*statfunc_hat - lq
	
	answer <- list(bias=bias,
	  bootstrap_confidence=c(bootstrap_confidence[1],bootstrap_confidence[2])
	, normal_confidence=c(normal_confidence[1], normal_confidence[2])
	, percentile_confidence=c(percentile_confidence[1], percentile_confidence[2])
	, pivotal_confidence=c(pivotal_confidence[1], pivotal_confidence[2]))
	return(answer)
#	list(answer)
}
