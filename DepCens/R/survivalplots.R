plot.weibull.t <- function(object){
  X_T <- object$X_T
  bi <- object$bi
  risco_a_T <- object$risco_a_T
  tempo <- object$time
  p <- object$p

  beta_T <- c()
  for (i in 1:p){
    beta_T[i] <- cbind(object[[1]][c(2+i)])
  }

  w_i <- rowMeans(bi)
  sobrev_T <- exp(-risco_a_T*exp(X_T%*%beta_T + w_i))
  sobrev_T <- as.data.frame(sobrev_T)
  names(sobrev_T)[1] <- "sobrev_T"

  teste_frame <- as.data.frame(cbind(tempo, X_T, risco_a_T, w_i, sobrev_T))
  o <- order(teste_frame$tempo)
  teste_frame <- teste_frame[o,]

  return(plot(teste_frame$tempo,teste_frame$sobrev_T, xlab = "Follow up time", ylab = "Survival probability of T", main = "Fit with Weibull distribution") + lines(lowess(teste_frame$sobrev_T~teste_frame$tempo), col='red', lwd=2))
}

plot.weibull.c <- function(object){
  X_C <- object$X_C
  bi <- object$bi
  risco_a_C <- object$risco_a_C
  tempo <- object$time
  p <- object$p
  q <- object$q

  beta_C <- c()
  for (j in 1:q){
    beta_C[j] <- cbind(object[[1]][c(4+p+j)])
  }

  w_i <- rowMeans(bi)
  alpha <- object[[1]][(5+p+q)]

  sobrev_C <- exp(-risco_a_C*exp(X_C%*%beta_C + alpha*w_i))

  sobrev_C <- as.data.frame(sobrev_C)
  names(sobrev_C)[1] <- "sobrev_C"

  teste_frame <- as.data.frame(cbind(tempo, X_C, risco_a_C, w_i, sobrev_C))
  o <- order(teste_frame$tempo)
  teste_frame <- teste_frame[o,]
  tempo <- teste_frame$tempo

  return(plot(teste_frame$tempo,teste_frame$sobrev_C, xlab = "Follow up time", ylab = "Survival probability of C", main = "Fit with Weibull distribution") + lines(lowess(teste_frame$sobrev_C~teste_frame$tempo), col='red', lwd=2))
}

plot.mep.t <- function(object){
  X_T <- object$X_T
  bi <- object$bi
  risco_a_T <- object$risco_a_T
  tempo <- object$time
  p <- object$p

  beta_T <- c()
  for (i in 1:p){
    beta_T[i] <- cbind(object[[1]][c(i)])
  }

  w_i <- rowMeans(bi)
  sobrev_T <- exp(-risco_a_T*exp(X_T%*%beta_T + w_i))
  sobrev_T <- as.data.frame(sobrev_T)
  names(sobrev_T)[1] <- "sobrev_T"

  teste_frame <- as.data.frame(cbind(tempo, X_T, risco_a_T, w_i, sobrev_T))
  o <- order(teste_frame$tempo)
  teste_frame <- teste_frame[o,]

  return(plot(teste_frame$tempo,teste_frame$sobrev_T, xlab = "Follow up time", ylab = "Survival probability of T", main = "Fit with piecewise exponential distribution") + lines(lowess(teste_frame$sobrev_T~teste_frame$tempo), col='red', lwd=2))
}

plot.mep.c <- function(object){
  X_C <- object$X_C
  bi <- object$bi
  risco_a_C <- object$risco_a_C
  tempo <- object$time
  p <- object$p
  q <- object$q
  bmax <- object$bmax

  beta_C <- c()
  for (j in 1:q){
    beta_C[j] <- cbind(object[[1]][c(p+bmax+j)])
  }

  w_i <- rowMeans(bi)
  alpha <- object[[1]][(p+bmax+q+1)]

  sobrev_C <- exp(-risco_a_C*exp(X_C%*%beta_C + alpha*w_i))

  sobrev_C <- as.data.frame(sobrev_C)
  names(sobrev_C)[1] <- "sobrev_C"

  teste_frame <- as.data.frame(cbind(tempo, X_C, risco_a_C, w_i, sobrev_C))
  o <- order(teste_frame$tempo)
  teste_frame <- teste_frame[o,]
  tempo <- teste_frame$tempo

  return(plot(teste_frame$tempo,teste_frame$sobrev_C, xlab = "Follow up time", ylab = "Survival probability of C", main = "Fit with piecewise exponential distribution") + lines(lowess(teste_frame$sobrev_C~teste_frame$tempo), col='red', lwd=2))
}
