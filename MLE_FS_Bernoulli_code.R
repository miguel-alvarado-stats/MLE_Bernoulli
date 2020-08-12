# MAXIMUM LIKELIHOOD ESTIMATION: FISHER-SCORING
# BERNOULLI DISTRIBUTION

MLE_FS_Bernoulli <- function(y, maxiter = 100, epsilon = 0.000000001, stop_criteria = 10000){

  n <- length(y)
  sumy <- sum(y)

  U <- function(n, sumy, theta){
    (1/(theta*(1 - theta)))*sumy - (n/(1 - theta))
  }

  J <- function(n, theta){
    n/(theta*(1 - theta))
  }

  # first guess:
  theta <- 0.5
  Estimator <- matrix(NA, ncol = 1)
  Estimator[1, 1] <- theta

  iter <- 1
  while( (stop_criteria > epsilon) & (iter <= maxiter) ){

    num <- U(n, sumy, Estimator[iter, 1])
    den <- J(n, Estimator[iter, 1])

    UPD <- (num/den)
    Estimator_iter <- as.matrix(Estimator[iter, 1]) + UPD
    Estimator <- rbind(Estimator, Estimator_iter)
    stop_criteria <- UPD^2
    iter <- iter + 1
  }

  Likelihood <- matrix(NA, ncol = 1)
  LogLikelihood <- matrix(NA, ncol = 1)

  for(i in 1:length(Estimator)){
    Likelihood[i] <- prod( (Estimator[i,1]^y)*((1 - Estimator[i,1])^(1 - y)) )
  }

  for(i in 1:length(Estimator)){
    LogLikelihood[i] <- sum( y*log(Estimator[i,1]) + (1 - y)*log(1 - Estimator[i,1]))
  }

  results <- cbind(format(Estimator, nsmall = 6), Likelihood, LogLikelihood)
  colnames(results) <- c("ML Estimator", "Likelihood", "Log-Likelihood")
  return(results)

}
