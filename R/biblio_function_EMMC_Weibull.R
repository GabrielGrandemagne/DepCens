
###--------------------------------------------------------------------------------------------------------------------------------------------------------
# Function to generate frailties of the full conditional log-likelihood function
###--------------------------------------------------------------------------------------------------------------------------------------------------------

log_cond_wi <- function(w_k_aux,Alpha,delta_T,delta_C,X_T,X_C,Betas_T,Betas_C,risco_a_T,risco_a_C,Sigma2){

  w_k <- log(w_k_aux/(1-w_k_aux))

  if (ncol(t(t(X_T))) == 1){
    pred_linear_T <- exp((t(t(X_T))*Betas_T)+w_k)
  }
  if (ncol(t(t(X_T))) != 1){
    pred_linear_T <- exp((X_T[,]%*%Betas_T)+w_k)
  }
  if (ncol(t(t(X_C))) == 1){
    pred_linear_C <- exp((t(t(X_C))*Betas_C)+Alpha*w_k)
  }
  if (ncol(t(t(X_C))) != 1){
    pred_linear_C <- exp((X_C[,]%*%Betas_C)+Alpha*w_k)
  }

  log_vero_w <- sum(delta_T*(w_k) - risco_a_T*pred_linear_T + delta_C*Alpha*(w_k) - risco_a_C*pred_linear_C) -((w_k^2)/(2*Sigma2))- log(w_k_aux*(1-w_k_aux))

  return(log_vero_w)
}

support_wi <-  function(w_k_aux,Alpha,delta_T,delta_C,X_T,X_C,Betas_T,Betas_C,risco_a_T,risco_a_C,Sigma2){(w_k_aux>0)*(w_k_aux<1)}

# _    _      _ _           _ _
#| |  | |    (_) |         | | |
#| |  | | ___ _| |__  _   _| | |
#| |/\| |/ _ \ | '_ \| | | | | |
#\  /\  /  __/ | |_) | |_| | | |
# \/  \/ \___|_|_.__/ \__,_|_|_|
###--------------------------------------------------------------------------------------------------------------------------------------------------------
# Function to estimate the parameters of failure times of the Weibull model, score functions of each parameter, a function used in multiroot
###--------------------------------------------------------------------------------------------------------------------------------------------------------

modelo_T_Weibull <-  function(param_t, t, X_T, delta.t, bi){
  p <- ncol(X_T)
  n <- nrow(X_T)

  risco_a_T <- (t^exp(param_t[1]))*exp(param_t[2])

  if(ncol(t(t(X_T))) == 1){
    pred_T <- exp(t(t(X_T))*param_t[3])*rowMeans(exp(bi[,]))
    w_kl_beta_T <- risco_a_T*pred_T
  }

  else{
    pred_T <- exp(X_T[,]%*%param_t[3:(p+2)])*rowMeans(exp(bi[,]))
    w_kl_beta_T <- risco_a_T*pred_T
  }

  w_kl_beta_T_num <- cbind(matrix(w_kl_beta_T, nrow = n, ncol = p))*X_T

  U_T_1 <- colSums(X_T*delta.t - w_kl_beta_T_num)
  w_kl_alpha_T_num <- w_kl_beta_T*log(t)*exp(param_t[1])
  U_alphaT <-  colSums(delta.t*(1 + exp(param_t[1])*log(t)) - w_kl_alpha_T_num)

  w_kl_lambda_T_num  <- w_kl_beta_T
  U_lambdaT <-  colSums(delta.t*(1) - w_kl_lambda_T_num)

  c(U_T_1 = U_T_1,U_alphaT=U_alphaT,U_lambdaT=U_lambdaT)
}

###--------------------------------------------------------------------------------------------------------------------------------------------------------
# Function to estimate the parameters of dependent censoring times of the Weibull model, score functions of each parameter, a function used in multiroot
###--------------------------------------------------------------------------------------------------------------------------------------------------------

modelo_C_Weibull <-  function(param_c, t, X_C, delta.c, bi){
  q <- ncol(X_C)
  n <- nrow(X_C)
  risco_a_C <- (t^exp(param_c[1]))*exp(param_c[2])

  if(ncol(t(t(X_C))) == 1){
    pred_C <- exp(t(t(X_C))*param_c[4])*rowMeans(exp(param_c[3]*bi[,]))
    w_kl_beta_C <- risco_a_C*pred_C
    w_kl_alpha_num <- risco_a_C*exp(t(t(X_C))*param_c[4])*rowMeans(bi*exp(param_c[3]*bi[,]))
  }

  else{
    pred_C <- exp(X_C[,]%*%param_c[4:(q+3)])*rowMeans(exp(param_c[3]*bi[,]))
    w_kl_beta_C <- risco_a_C*pred_C
    w_kl_alpha_num <- risco_a_C*exp(X_C[,]%*%param_c[4:(q+3)])*rowMeans(bi*exp(param_c[3]*bi[,]))
  }

  w_kl_beta_C_num <- cbind(matrix(w_kl_beta_C, nrow = n, ncol = q))*X_C

  U_betas <- colSums(X_C*delta.c - w_kl_beta_C_num)
  U_alpha <- colSums(delta.c*rowMeans(bi[,]) - w_kl_alpha_num)
  w_kl_alpha_C_num <- w_kl_beta_C*log(t)*exp(param_c[1])
  U_alphaC <-  colSums(delta.c*(1 + exp(param_c[1])*log(t)) - w_kl_alpha_C_num)

  w_kl_lambda_C_num   <- w_kl_beta_C
  U_lambdaC <-  colSums(delta.c*(1) - w_kl_lambda_C_num)

  c(U_alphaC=U_alphaC,U_lambdaC=U_lambdaC,U_alpha=U_alpha,U_betas = U_betas)
}

crit_weibull <- function(X_T,X_C,bi,n,alpha_T,lambda_T,beta_T,alpha_C,lambda_C,beta_C,alpha,delta_t,delta_c,time,m,ident) {
  w <- bi                                                                       # pega a matriz de fragilidades salva na ultima itera??o
  L <- ncol(w)
  num_param <- length(c(alpha_T,lambda_T,beta_T,alpha_C,lambda_C,beta_C,alpha)) # numero de parametros, parametros salvos em fit (fit <- c((out[s,]+out[s-1,]+out[s-2,])/3))
  log_vero <- matrix(NA,n,L)                                                    #log-verossimilhan?a, L numero de replicas monte carlo de w

  for ( l in 1:L){
    log_vero[,l] <-  delta_t*(log(alpha_T) + (alpha_T-1)*log(time) + log(lambda_T) + X_T%*%beta_T + w[,l]) - (time^alpha_T)*lambda_T*exp(X_T%*%beta_T + w[,l])
    + delta_c*(log(alpha_C) + (alpha_C-1)*log(time) + log(lambda_C) + X_C%*%beta_C + alpha*w[,l]) - (time^alpha_C)*lambda_C*exp(X_C%*%beta_C + alpha*w[,l])
  }

  vero <- exp(log_vero)
  vero_grupo <- matrix(NA,m,L)  #m eh o numero de grupos/cluster

  for (i in 1:m){
    vero_grupo[i,] <- matrixStats::colProds(vero, rows = which(ident==i))  # verossimilhan?a para cada grupo
  }

  log_lik <- sum(log(rowSums(vero_grupo)/L))

  AIC <- 2*(-log_lik + num_param)

  BIC <- 2*(-log_lik + 0.5*log(n)*num_param)

  HQ <- 2*(-log_lik + log(log(n))*num_param)

  return(cbind(AIC,BIC,HQ))
}

###--------------------------------------------------------------------------------------------------------------------------------------------------------
# Function to calculate the vector with the first order derivatives of the Weibull model
###--------------------------------------------------------------------------------------------------------------------------------------------------------

Esp_Deriv_Weibull <-  function(t,delta_T, delta_C, X_T, X_C, beta_T, beta_C,alpha,Sigma2, alpha_T, alpha_C,lambda_T, lambda_C, w_k_grupo, ident){
  wk = w_k_grupo[ident,]
  num_param <- length(beta_T)+2+length(beta_C)+2+length(alpha)+length(Sigma2)
  deriv1 <- matrix(NA,num_param,ncol(wk))  #vetor que com as derivadas de primeira ordem
  p <- ncol(X_T)
  q <- ncol(X_C)
  pred_T <- as.vector(exp(X_T%*%beta_T))
  pred_C <- as.vector(exp(X_C%*%beta_C))
  Delta_t <- delta_T%*%t(rep(1,ncol(wk)))
  Delta_c <- delta_C%*%t(rep(1,ncol(wk)))

  deriv1[1,] <- colSums(Delta_t*(alpha_T^(-1) + log(t)) - (t^(alpha_T))*log(t)*lambda_T*pred_T*exp(wk)) #derivada de alpha.t

  deriv1[2,] <- colSums(Delta_t*(lambda_T^(-1)) - (t^(alpha_T))*pred_T*exp(wk)) #derivada de lambda.t

  for (i in 1:p){
    deriv1[2+i,] <- colSums(Delta_t*(X_T[,i]) - (t^(alpha_T))*lambda_T*pred_T*X_T[,i]*exp(wk)) #derivada do vetor beta_T
  }

  deriv1[3+p,] <- colSums(Delta_c*(alpha_C^(-1) + log(t)) - (t^(alpha_C))*log(t)*lambda_C*pred_C*exp(alpha*wk))  #derivada de alpha.c

  deriv1[4+p,] <- colSums(Delta_c*(lambda_C^(-1)) - (t^(alpha_C))*pred_C*exp(alpha*wk))  #derivada de lambda.c

  for (i in 1:q){
    deriv1[4+p+i,] <- colSums(Delta_c*(X_C[,i]) - (t^(alpha_C))*lambda_C*pred_C*X_C[,i]*exp(alpha*wk))  #derivada do vetor beta_C
  }

  deriv1[4+p+q+1,] <- colSums(Delta_c*(wk) - (t^(alpha_C))*lambda_C*pred_C*exp(alpha*wk)*wk)  #derivada de alpha (parametro da dependencia)

  deriv1[4+p+q+2,] <- colSums(-0.5*Sigma2^(-1) + 0.5*(Sigma2^(-2))*w_k_grupo^(2)) #derivada de sigma2

  aux <- deriv1[,1]%*%t(deriv1[,1])

  for( i in 2:ncol(wk)){
    aux <- aux + deriv1[,i]%*%t(deriv1[,i])
  }

  return(aux/ncol(wk))
}

###--------------------------------------------------------------------------------------------------------------------------------------------------------
# Function to calculate the matrix with the second order derivatives of the Weibull model
###--------------------------------------------------------------------------------------------------------------------------------------------------------

Esp_DerivParc_Weibull <-  function(t,delta_T, delta_C, X_T, X_C, beta_T, beta_C, alpha,Sigma2, alpha_T, alpha_C, lambda_T, lambda_C,w_k_grupo, ident){

  num_param <- length(beta_T)+2+length(beta_C)+2+length(alpha)+length(Sigma2)
  deriv2 <- matrix(0,num_param,num_param)  #vetor que com as derivadas de segunda ordem,
  wk = w_k_grupo[ident,]
  p <- ncol(X_T)
  q <- ncol(X_C)
  pred_T <- as.vector(exp(X_T%*%beta_T))
  pred_C <- as.vector(exp(X_C%*%beta_C))
  Delta_t <- delta_T%*%t(rep(1,ncol(wk)))
  Delta_c <- delta_C%*%t(rep(1,ncol(wk)))

  deriv2[1,1] <- mean(colSums(Delta_t*(-alpha_T^(-2)) - (t^(alpha_T))*(log(t)^(2))*lambda_T*pred_T*exp(wk)))  #derivada de alpha.t da derivada de alpha.t

  deriv2[1,2] <- mean(colSums(- (t^(alpha_T))*(log(t))*pred_T*exp(wk)))   #derivada de lambda.t da derivada de alpha.t
  deriv2[2,1] <- deriv2[1,2]

  for (i in 1:p){
    deriv2[1,2+i] <- mean(colSums(- (t^(alpha_T))*(log(t))*lambda_T*pred_T*X_T[,i]*exp(wk)))  #derivada do vetor beta_T da derivada de alpha.t
    deriv2[2+i,1] <- deriv2[1,2+i]
  }

  deriv2[2,2] <- mean(colSums( Delta_t*(-lambda_T^(-2))))   #derivada de lambda.t da derivada de lambda.t

  for (i in 1:p){
    deriv2[2,2+i] <- mean(colSums( - (t^(alpha_T))*pred_T*X_T[,i]*exp(wk)))  #derivada do vetor beta_T da derivada de lambda.t
    deriv2[2+i,2] <- deriv2[2,2+i]
  }

  for (i in 1:p){
    deriv2[2+i,2+i] <- mean(colSums( - ((t^(alpha_T))*lambda_T*pred_T*(X_T[,i]^(2))*exp(wk))))   #derivada do vetor beta_T da derivada do vetor beta_T
  }

  for (j in 1:(p-1)){
    for (i in 1:(p-1)){
      if (ncol(t(t(X_T))) == 1){
        next
      }
      deriv2[2+j,3+i] <- mean(colSums(- (t^(alpha_T))*lambda_T*pred_T*X_T[,1+i]*X_T[,j]*exp(wk))) #derivada de beta_T_i da derivada de beta_T_(i-1)
      deriv2[3+i,2+j] <- deriv2[2+j,3+i]
    }
  }

  deriv2[3+p,3+p] <- mean(colSums(Delta_c*(-alpha_C^(-2)) - (t^(alpha_C))*(log(t)^(2))*lambda_C*pred_C*exp(alpha*wk)))    #derivada de alpha.c da derivada de alpha.c

  deriv2[3+p,4+p] <- mean(colSums( - (t^(alpha_C))*log(t)*pred_C*exp(alpha*wk)))  #derivada de lambda.c da derivada de alpha.c
  deriv2[4+p,3+p] <- deriv2[3+p,4+p]

  for (i in 1:q){
    deriv2[3+p,(4+p+i)] <- mean(colSums( - (t^(alpha_C))*log(t)*lambda_C*pred_C*X_C[,i]*exp(alpha*wk)))   #derivada do vetor beta_C da derivada de alpha.c
    deriv2[(4+p+i),3+p] <- deriv2[3+p,(4+p+i)]
  }

  deriv2[3+p,(4+p+q+1)] <- mean(colSums( - (t^(alpha_C))*log(t)*lambda_C*pred_C*exp(alpha*wk)*wk))    #derivada de alpha da derivada de alpha.c
  deriv2[(4+p+q+1),3+p] <- deriv2[3+p,(4+p+q+1)]

  deriv2[4+p,4+p] <- mean(colSums(Delta_c*(-lambda_C^(-2)) ))  #derivada de lambda.c da derivada de lambda.c

  for (i in 1:q){
    deriv2[4+p,(4+p+i)] <- mean(colSums( - (t^(alpha_C))*pred_C*X_C[,i]*exp(alpha*wk)))    #derivada de beta_C da derivada de lambda.c
    deriv2[(4+p+i),4+p] <- deriv2[4+p,(4+p+i)]
  }

  deriv2[4+p,(4+p+q+1)] <- mean(colSums( - (t^(alpha_C))*pred_C*exp(alpha*wk)*wk))  #derivada de alpha da derivada de lambda.c
  deriv2[(4+p+q+1),4+p] <- deriv2[4+p,(4+p+q+1)]

  for (i in 1:q){
    deriv2[4+p+i,4+p+i] <- mean(colSums(-((t^(alpha_C))*lambda_C*pred_C*(X_C[,i]^(2))*exp(alpha*wk))))  #derivada de beta_C da derivada de beta_C
  }

  for (j in 1:(q-1)){
    for (i in 1:(q-1)){
      if (ncol(t(t(X_C))) == 1){
        next
      }
      deriv2[4+p+j,4+p+1+i] <- mean(colSums( - (t^(alpha_C))*lambda_C*pred_C*X_C[,1+i]*X_C[,j]*exp(alpha*wk)))  #derivada de beta_C_i da derivada de beta_C_(i-1)
      deriv2[4+p+1+i,4+p+j] <- deriv2[4+p+j,4+p+1+i]
    }
  }

  for (i in 1:q){
    deriv2[4+p+i,(4+p+q+1)] <- mean(colSums(- t^(alpha_C)*lambda_C*pred_C*X_C[,i]*exp(alpha*wk)*wk))   #derivada de alpha da derivada de beta_C
    deriv2[(4+p+q+1),4+p+i] <- deriv2[4+p+i,(4+p+q+1)]
  }

  deriv2[(4+p+q+1),(4+p+q+1)] <- mean(colSums( - (t^(alpha_C))*lambda_C*pred_C*exp(alpha*wk)*wk^(2)))     #derivada de alpha da derivada de alpha

  deriv2[(4+p+q+2),(4+p+q+2)] <- mean(colSums( 0.5*Sigma2^(-2) - (Sigma2^(-3))*w_k_grupo^(2)))   #derivada de sigma2 da derivada de sigma2

  return((as.matrix(deriv2)))
}

###---------------------------------------------------------------------------------------------------
# Funcao para ajustar o modelo completo
###---------------------------------------------------------------------------------------------------

model_Weibull_dep <-  function(formula, data, delta_t, delta_c, ident){

  formula <- Formula::Formula(formula)
  mf <- stats::model.frame(formula=formula, data=data)
  Terms <- stats::terms(mf)
  Z <- stats::model.matrix(formula, data = mf, rhs = 1)
  X <- stats::model.matrix(formula, data = mf, rhs = 2)
  X_T <- Z[,-1]
  X_T <- as.matrix(X_T)
  X_C <- X[,-1]
  X_C <- as.matrix(X_C)
  Xlabels <- colnames(X_T)
  Zlabels <- colnames(X_C)
  time <- stats::model.response(mf)
  q <- ncol(X_C)
  p <- ncol(X_T)
  n <- nrow(data)
  m <- max(ident)

  # Common mistakes
  if(any(time < 0 | time == 0)) stop('time must be greater than 0')

  ###----------------------------------------------------------------------------------------------------
  # Initial value of the parameters betas_T, alpha_T, lambda_T, sigma2

  control = coxph.control()
  ajuste_coxph_T <- coxph(Surv(time, delta_t) ~  X_T, method="breslow")
  #risco_a_T <- basehaz(ajuste_coxph_T, centered=FALSE)$hazard  #cumulative hazard
  beta_T <- ajuste_coxph_T$coef

  ajuste_coxph_C <- coxph(Surv(time, delta_c) ~  X_C, method="breslow")
  #risco_a_C <- basehaz(ajuste_coxph_C, centered=FALSE)$hazard  #cumulative hazard
  beta_C <- ajuste_coxph_C$coef

  alpha <- 0
  sigma2 <- 1
  alpha_T <- 1
  lambda_T <- 1
  alpha_C <- 1
  lambda_C <- 1

  param <- c(alpha_T,lambda_T,beta_T,alpha_C,lambda_C,beta_C,alpha,sigma2)

  risco_a_T <- rep(0.05,n)
  risco_a_C <- rep(0.05,n)

  ###----------------------------------------------------------------------------------------------------
  # Especificacoes do algorimo EMMC

  maxit <- 100 #numero maximo de iteracoes
  eps1= rep(1e-3, length(param))
  eps2= rep(1e-4, length(param))
  n_intMCs = c(rep(10,40),rep(25,30), rep(50,20), rep(100,10))

  ## Iniciando os objetos utilizados
  out =  matrix(NA,maxit+1,length(param))
  dif =  matrix(NA,maxit+1,length(param))
  final = length(param)
  count = rep(0,maxit+1)
  s=1
  continue=TRUE

  ###-------------------------------------------------------------------------------------------------

  while (continue == TRUE) {

    count = rep(0,maxit+1)
    out[s,] =  c(alpha_T,lambda_T,beta_T,alpha_C,lambda_C,beta_C,alpha,sigma2)
    n_intMC = n_intMCs[s]
    w_chapeu_grupo <- matrix(NA, m, n_intMC)

    for ( k in 1:m){  #k eh grupo
      w_trans <- arms(0.5, myldens=log_cond_wi, indFunc=support_wi,n.sample=n_intMC,Alpha=alpha,delta_T=delta_t[ident==k],delta_C=delta_c[ident==k],X_T=X_T[ident==k,], X_C=X_C[ident==k,],Betas_T=beta_T,Betas_C=beta_C, risco_a_T=risco_a_T[ident==k],risco_a_C=risco_a_C[ident==k],Sigma2=sigma2)
      w_auxi <- log(w_trans)-log(1-w_trans)
      w_chapeu_grupo[k,] <- w_auxi
    }
    bi = w_chapeu_grupo[ident,]

    sigma2 <- mean(w_chapeu_grupo^2)  #variancia da fragilidade
    ###------------------------------------------------------------------------------------
    # estimando os parametros da dist. dos tempos de falha

    S_T <- multiroot(f = modelo_T_Weibull, start = rep(0.01, 2+p),t=time, X_T=X_T, delta.t=delta_t,bi=bi)
    param_T <- S_T$root
    alpha_T <- exp(param_T[1])
    lambda_T <- exp(param_T[2])
    beta_T <- param_T[3:(2+p)]
    risco_a_T <- (time^alpha_T)*(lambda_T)

    ###-----------------------------------------------------------------------------------
    #para censura

    S_C <- multiroot(f = modelo_C_Weibull, start = rep(0.01,3+q),t=time, X_C=X_C, delta.c=delta_c,bi=bi)
    param_C <- S_C$root
    alpha_C <- exp(param_C[1])
    lambda_C <- exp(param_C[2])
    alpha <- param_C[3]
    beta_C <- param_C[4:(3+q)]
    risco_a_C <- (time^alpha_C)*(lambda_C)


    ###-----------------------------------------------------------------------------------
    # criterio de parada
    out[s+1,]=  c(alpha_T,lambda_T,beta_T,alpha_C,lambda_C,beta_C,alpha,sigma2)
    #print(out[s+1,])
    dif[s,] <- (abs(out[s+1,]-out[s,]))/(abs(out[s,])-eps1)

    for (z in 1: length(param)){
      if (dif[s,z]<eps2[z]) {
        count[s] = count[s] + 1
      }
    }

    s=s+1
    if (s>3) {
      continue=((s <= maxit) & (count[s] < length(param))&(count[s-1] < length(param))&(count[s-2] < length(param))) #a diferenca precisa ser menor que um erro por 3 iteracaoes consecutivas
    }
  } ## fim do EMMC

  param_est <- c((out[s,]+out[s-1,]+out[s-2,])/3)


  Esp_deriv_ordem1 <- Esp_Deriv_Weibull(t=time,delta_T=delta_t, delta_C=delta_c,X_T=X_T, X_C=X_C,beta_T=beta_T, beta_C=beta_C,
                                        alpha=alpha, Sigma2=sigma2, alpha_T=alpha_T, alpha_C=alpha_C,
                                        lambda_T=lambda_T,lambda_C=lambda_C, w_k_grupo=w_chapeu_grupo, ident=ident)


  Esp_deriv_ordem2 <- Esp_DerivParc_Weibull(t=time,delta_T=delta_t, delta_C=delta_c, X_T=X_T, X_C=X_C,beta_T=beta_T, beta_C=beta_C,
                                            alpha=alpha, Sigma2=sigma2, alpha_T=alpha_T, alpha_C=alpha_C,lambda_T=lambda_T, lambda_C=lambda_C, w_k_grupo=w_chapeu_grupo, ident=ident)


  InfFisher <- (-Esp_deriv_ordem2) - Esp_deriv_ordem1
  Var <- solve(InfFisher)
  ErroPadrao <- sqrt(diag(Var))
  if (any(is.na(ErroPadrao))) warning("The algorithm did not converge. It might converge if you run the function again.", call. = FALSE)

  p_value_alpha <- 2*stats::pnorm(-abs(param_est[(5+p+q)]/ErroPadrao[(5+p+q)]))
  p_value_t <- vector()
  for (i in 1:p){
    p_value_t[i] <- 2*stats::pnorm(-abs((param_est[c(2+i)])/(ErroPadrao[c(2+i)])))
  }
  p_value_c <- vector()
  for (j in 1:q){
    p_value_c[j] <- 2*stats::pnorm(-abs((param_est[c(4+p+j)])/(ErroPadrao[c(4+p+j)])))
  }

  criterios <- crit_weibull(X_T,X_C,bi,n,alpha_T,lambda_T,beta_T,alpha_C,lambda_C,beta_C,alpha,delta_t,delta_c,time,m,ident)

  # Ajustando a classe/lista
  fit <- param_est
  fit <- list(fit=fit)
  fit$stde <- ErroPadrao
  fit$crit <- criterios
  fit$pvalue <- c(p_value_alpha, p_value_t, p_value_c)
  fit$n <- n
  fit$p <- p
  fit$q <- q
  fit$call <- match.call()
  fit$formula <- stats::formula(Terms)
  fit$terms <- stats::terms.formula(formula)
  fit$labels1 <- Zlabels
  fit$labels2 <- Xlabels
  fit$risco_a_T <- risco_a_T
  fit$risco_a_C <- risco_a_C
  fit$bi <- bi
  fit$X_T <- X_T
  fit$X_C <- X_C
  fit$time <- time
  class(fit) <- "dcensoring"
  return(fit)
}

#  ___  ___ ___________
#  |  \/  ||  ___| ___ \
#  | .  . || |__ | |_/ /
#  | |\/| ||  __||  __/
#  | |  | || |___| |
#  \_|  |_/\____/\_|
###---------------------------------------------------------------------------------------------------------------------------

modelo_T_MEP <-  function(beta_T, X_T, delta.t,risco_a_T,bi,n_intMC){
  p <- ncol(X_T)
  n <- nrow(X_T)

  if(ncol(t(t(X_T))) == 1){
    w_kl_beta_T <- risco_a_T*exp((X_T[,]*beta_T))*(rowSums(exp(bi[,]))/n_intMC)
  }
  else{
    w_kl_beta_T <- risco_a_T*exp((X_T[,]%*%beta_T))*(rowSums(exp(bi[,]))/n_intMC)
  }

  w_kl_beta_T_num <- cbind(matrix(w_kl_beta_T, nrow = n, ncol = p))*X_T
  U_T_1 <- colSums((X_T*delta.t - w_kl_beta_T_num))
  c(U_T_1 = U_T_1)
}

###---------------------------------------------------------------------------------------------------------------------------
# funcao para estimar os coeficientes de regress?o dos tempos de censura dependente

modelo_C_MEP <-  function(beta_C, X_C,delta.c,risco_a_C,bi,n_intMC){
  q <- ncol(X_C)
  n <- nrow(X_C)

  if(ncol(t(t(X_C))) == 1){
    w_kl_beta_C <- risco_a_C*exp((X_C[,]*beta_C[1:q]))*(rowSums(exp(beta_C[q+1]*bi[,]))/n_intMC)
    w_kl_beta_C_alpha_num <- risco_a_C*exp(X_C[,]*beta_C[1:q])*(rowSums(bi[,]*exp(beta_C[q+1]*bi[,]))/n_intMC)
  }

  else{
    w_kl_beta_C <- risco_a_C*exp((X_C[,]%*%beta_C[1:q]))*(rowSums(exp(beta_C[q+1]*bi[,]))/n_intMC)
    w_kl_beta_C_alpha_num <- risco_a_C*exp(X_C[,]%*%beta_C[1:q])*(rowSums(bi[,]*exp(beta_C[q+1]*bi[,]))/n_intMC)
  }

  w_kl_beta_C_num <- cbind(matrix(w_kl_beta_C, nrow = n, ncol = q))*X_C
  U_C_1 <- colSums((X_C*delta.c - w_kl_beta_C_num))

  U_C_alpha <- sum(delta.c*(rowSums(bi[,])/n_intMC) - w_kl_beta_C_alpha_num)

  c(U_C_1 = U_C_1,U_C_alpha=U_C_alpha)
}

###-------------------------------------------------------------------------------------------
# funcao que calcula a grade

time_grid <- function(time, event, n.int=NULL)
{
  o <- order(time)
  time <- time[o]
  event <- event[o]
  time.aux <- unique(time[event==1])
  if(is.null(n.int))
  {
    n.int <- length(time.aux)
  }

  m <- length(time.aux)
  if(n.int > m)
  {
    a <- c(0,unique(time[event==1]))
    a[length(a)] <- Inf
  }
  else
  {
    b <- min(m,n.int)
    k1 <- trunc(m/b)
    r <- m-b*k1
    k2 <- k1+1
    idf1 <- seq(k1,(b-r)*k1, k1)
    idf2 <- sort(seq(m,max(idf1),-k2))
    idf <- unique(c(idf1,idf2))
    a_inf <- c(0,time.aux[idf])
    a_inf[length(a_inf)] <- Inf
    a_s_inf  <- c(0,time.aux[idf])
  }
  saida <- list(a_inf,a_s_inf)

  return(saida)
}

###-----------------------------------------------------------------------------------
# funcao que calculo a variavel inidicadora conforme a grade de intervalos

ID_a <- function(a,U,t){

  for( i in 1:(length(U)) ){
    if( U[length(U)-i+1]== 1 ){
      a[(length(U)-i+2)] <- a[(length(U)-i+3)]
    }
  }
  a2 <- c(unique(a))
  id <- as.numeric(cut(t,a2))

  return(a2)
}

###-------------------------------------------------------------------------------------------------
# fun??o que calcula as diferen?as (t-a[j]) para cada individuo, que sera usado no calculo do risco acumulado

deltas <- function(a_inf,t,IDE,b,n){
  tij <- matrix(0,n,b)
  deltaij <- matrix(0,n,b)
  a <- a_inf
  for(i in 1:n){
    for(j in 1:b){
      deltaij[i,j] <- (min(t[i], a[j+1])-a[j])*((t[i]-a[j])>0)
    }
  }
  return(deltaij)
}

###------------------------------------------------------------------------------------------------------
#H0ti eh a fun??o risco acumulado para o MEP

H0ti <- function(a_inf,t,deltaij,lambda,b,n){
  Hoti <- matrix(0,n,b)
  H0ti <- NULL
  for(i in 1:n){
    for(j in 1:b){
      Hoti[i,j] <- (deltaij[i,j]*lambda[j] )
    }
  }
  for(i in 1:n){
    H0ti[i]<- sum(Hoti[i,])
  }
  return(H0ti)
}

crit_mep <- function(X_T,X_C,bi,n,beta_T,lambda_T_j,beta_C,alpha,lambda_C_j,delta_t,delta_c,risco_a_T,risco_a_C,id_T,id_C,m,ident) {
  w <- bi                                                           # pega a matriz de fragilidades salva na ultima itera??o
  L <- ncol(w)
  num_param <- length(c(beta_T,lambda_T_j,beta_C,alpha,lambda_C_j)) # numero de parametros, parametros salvos em fit (fit <- c((out[s,]+out[s-1,]+out[s-2,])/3))
  log_vero <- matrix(NA,n,L)                                        # log-verossimilhan?a, L numero de replicas monte carlo de w


  for (l in 1:L){
    if (ncol(t(t(X_T))) == 1 && ncol(t(t(X_C))) == 1){
      log_vero[,l] <-  delta_t*(t(t(log(lambda_T_j[id_T]))) + X_T*beta_T + w[,l]) - risco_a_T*exp(X_T*beta_T + w[,l])
      + delta_c*(t(t(log(lambda_C_j[id_C]))) + X_C*beta_C + alpha*w[,l]) - risco_a_C*exp(X_C*beta_C + alpha*w[,l])
    }

    if (ncol(t(t(X_T))) == 1 && ncol(t(t(X_C))) != 1){
      log_vero[,l] <-  delta_t*(t(t(log(lambda_T_j[id_T]))) + X_T*beta_T + w[,l]) - risco_a_T*exp(X_T*beta_T + w[,l])
      + delta_c*(t(t(log(lambda_C_j[id_C]))) + X_C%*%beta_C + alpha*w[,l]) - risco_a_C*exp(X_C%*%beta_C + alpha*w[,l])
    }

    if (ncol(t(t(X_T))) != 1 && ncol(t(t(X_C))) == 1){
      log_vero[,l] <-  delta_t*(t(t(log(lambda_T_j[id_T]))) + X_T%*%beta_T + w[,l]) - risco_a_T*exp(X_T%*%beta_T + w[,l])
      + delta_c*(t(t(log(lambda_C_j[id_C]))) + X_C*beta_C + alpha*w[,l]) - risco_a_C*exp(X_C*beta_C + alpha*w[,l])

    }
    if (ncol(t(t(X_T))) != 1 && ncol(t(t(X_C))) != 1){
      log_vero[,l] <-  delta_t*(t(t(log(lambda_T_j[id_T]))) + X_T%*%beta_T + w[,l]) - risco_a_T*exp(X_T%*%beta_T + w[,l])
      + delta_c*(t(t(log(lambda_C_j[id_C]))) + X_C%*%beta_C + alpha*w[,l]) - risco_a_C*exp(X_C%*%beta_C + alpha*w[,l])
    }
  }

  vero <- exp(log_vero)
  vero_grupo <- matrix(NA,m,L)  #m eh o numero de grupos/cluster

  for (i in 1:m){
    vero_grupo[i,] <- matrixStats::colProds(vero, rows = which(ident==i))  # verossimilhan?a para cada grupo
  }

  log_lik <- sum(log(rowSums(vero_grupo)/L))

  AIC <- 2*(-log_lik + num_param)
  BIC <- 2*(-log_lik + 0.5*log(n)*num_param)
  HQ <- 2*(-log_lik + log(log(n))*num_param)

  return(cbind(AIC,BIC,HQ))
}

###-------------------------------------------------------------------------------------------
# calculo das derivas de primeira ordem

Esp_Deriv_MEP <-  function(X_T,X_C,delta_T,delta_C,beta_T,beta_C,alpha,theta,nu_ind_j_T,nu_ind_j_C,lambda_T_j,lambda_C_j,risco_a_T,risco_a_C,deltaij_T,deltaij_C,w_k_grupo,ident){

  n <- nrow(X_T)
  p <- ncol(X_T)
  q <- ncol(X_C)
  wk <- w_k_grupo[ident,]
  deriv1 <- matrix(NA,(p+q+2+length(lambda_T_j)+length(lambda_C_j)),ncol(wk))

  pred_T <- as.vector(exp(X_T%*%beta_T))*exp(wk)
  pred_C <- as.vector(exp(X_C%*%beta_C))*exp(alpha*wk)
  Delta_t <- delta_T%*%t(rep(1,ncol(wk)))
  Delta_c <- delta_C%*%t(rep(1,ncol(wk)))
  lambda_t <- t(lambda_T_j%*%t(rep(1,n)))
  lambda_c <- t(lambda_C_j%*%t(rep(1,n)))

  for (i in 1:p){
    deriv1[i,] <- colSums(X_T[,i]*(Delta_t - risco_a_T*pred_T))  #derivadas de beta_T's
  }

  for ( j in 1:length(lambda_T_j)){
    deriv1[(j+p),] <- nu_ind_j_T[j]*(lambda_T_j[j]^(-1)) - colSums(deltaij_T[,j]*pred_T)  #derivada de lambda_T_j
  }

  for (i in 1:q){
    deriv1[(length(lambda_T_j)+p+i),] <- colSums(X_C[,i]*(Delta_c - risco_a_C*pred_C))  #derivadas de beta_C's
  }

  deriv1[(length(lambda_T_j)+p+q+1),] <- colSums(wk*(Delta_c - risco_a_C*pred_C))  #derivada de alpha

  for ( j in 1:length(lambda_C_j)){
    deriv1[(length(lambda_T_j)+p+q+1+j),] <- nu_ind_j_C[j]*(lambda_C_j[j]^(-1)) - colSums(deltaij_C[,j]*pred_C)   #derivada de lambda_C_j
  }

  deriv1[(length(lambda_C_j)+length(lambda_T_j)+p+q+2),] <- colSums(-0.5*theta^(-1) + 0.5*(theta^(-2))*w_k_grupo^(2))    #derivada de tau

  aux <- deriv1[,1]%*%t(deriv1[,1])

  for( i in 2:ncol(wk)){
    aux <- aux + deriv1[,i]%*%t(deriv1[,i])
  }
  return(aux/ncol(wk))
}

###---------------------------------------------------------------------------------------------------------------------------
#funcoes para calcular o vetor com as derivadas de segunda ordem

Esp_DerivParc_MEP <-  function( X_T, X_C,delta_T,delta_C, beta_T, beta_C, alpha, theta, nu_ind_j_T,nu_ind_j_C,lambda_T_j,lambda_C_j, risco_a_T, risco_a_C,deltaij_T,deltaij_C, w_k_grupo, ident){

  n <- nrow(X_T)
  p <- ncol(X_T)
  q <- ncol(X_C)
  deriv2 <- matrix(0,(p+q+2+length(lambda_T_j)+length(lambda_C_j)),(p+q+2+length(lambda_T_j)+length(lambda_C_j)))
  wk = w_k_grupo[ident,]

  pred_T <- as.vector(exp(X_T%*%beta_T))*exp(wk)
  pred_C <- as.vector(exp(X_C%*%beta_C))*exp(alpha*wk)
  Delta_t <- delta_T%*%t(rep(1,ncol(wk)))
  Delta_c <- delta_C%*%t(rep(1,ncol(wk)))
  lambda_t <- t(lambda_T_j%*%t(rep(1,n)))
  lambda_c <- t(lambda_C_j%*%t(rep(1,n)))

  for (i in 1:p){
    deriv2[i,i] <- - mean(colSums(X_T[,i]*X_T[,i]*risco_a_T*pred_T))  #derivada de beta_T da derivada de beta_T
  }

  for (j in 1:(p-1)){
    for (i in 1:(p-1)){
      if (ncol(t(t(X_T))) == 1){
        next
      }
      deriv2[j,1+i] <- - mean(colSums(X_T[,1+i]*X_T[,j]*risco_a_T*pred_T))   #derivada de beta_T_i da derivada de beta_T_(i-1)
      deriv2[1+i,j] <- deriv2[j,1+i]
    }
  }

  for (j in 1:length(lambda_T_j)){
    deriv2[(j+p),(j+p)] <- - nu_ind_j_T[j]*(lambda_T_j[j]^(-2))  #derivada de lambda_T_j da derivada de lambda_T_j
  }

  for (i in 1:p){
    for (j in 1:length(lambda_T_j)){
      deriv2[i,(j+p)] <- -mean(colSums(deltaij_T[,j]*pred_T*X_T[,i]))    #derivada de beta_T da derivada de lambda_T_j
      deriv2[(j+p),i] <- deriv2[i,(j+p)]
    }
  }

  for (i in 1:q){
    deriv2[(length(lambda_T_j)+p+i),(length(lambda_T_j)+p+i)] <- - mean(colSums(X_C[,i]*X_C[,i]*risco_a_C*pred_C))  #derivada de beta_C da derivada de beta_C
  }

  for (j in 1:(q-1)){
    for (i in 1:(q-1)){
      if (ncol(t(t(X_C))) == 1){
        next
      }
      deriv2[(length(lambda_T_j)+p+j),(length(lambda_T_j)+p+1+i)] <- - mean(colSums(X_C[,1+i]*X_C[,j]*risco_a_C*pred_C))
      deriv2[(length(lambda_T_j)+p+1+i),(length(lambda_T_j)+p+j)] <-  deriv2[(length(lambda_T_j)+p+j),(length(lambda_T_j)+p+1+i)]   #derivada de beta_C_i da derivada de beta_C_(i-1)
    }
  }

  deriv2[(length(lambda_T_j)+p+q+1),(length(lambda_T_j)+p+q+1)] <- - mean(colSums(wk*wk*risco_a_C*pred_C))  #derivada de alpha da derivada de alpha

  for (i in 1:q){
    deriv2[(length(lambda_T_j)+p+i),(length(lambda_T_j)+p+q+1)] <- - mean(colSums(X_C[,i]*wk*risco_a_C*pred_C))
    deriv2[(length(lambda_T_j)+p+q+1),(length(lambda_T_j)+p+i)] <- deriv2[(length(lambda_T_j)+p+i),(length(lambda_T_j)+p+q+1)]    #derivada de alpha da derivada de beta_C's
  }

  for (j in 1:length(lambda_C_j)){
    deriv2[(j+length(lambda_T_j)+p+q+1),(j+length(lambda_T_j)+p+q+1)] <- - nu_ind_j_C[j]*(lambda_C_j[j]^(-2))    #derivada de lambda_C_j da derivada de lambda_C_j
  }

  for (i in 1:q){
    for (j in 1:length(lambda_C_j)){
      deriv2[(length(lambda_T_j)+p+i),(j+length(lambda_T_j)+p+q+1)] <- -mean(colSums(deltaij_C[,j]*pred_C*X_C[,i]))
      deriv2[(j+length(lambda_T_j)+p+q+1),(length(lambda_T_j)+p+i)] <- deriv2[(length(lambda_T_j)+p+i),(j+length(lambda_T_j)+p+q+1)]    #derivada de beta_C's da derivada de lambda_C_j
    }
  }

  for (j in 1:length(lambda_C_j)){
    deriv2[(length(lambda_T_j)+p+q+1),(j+length(lambda_T_j)+p+q+1)] <- -mean(colSums(deltaij_C[,j]*pred_C*wk))
    deriv2[(j+length(lambda_T_j)+p+q+1),(length(lambda_T_j)+p+q+1)] <- deriv2[(length(lambda_T_j)+p+q+1),(j+length(lambda_T_j)+p+q+1)]   #derivada de alpha da derivada de lambda_C_j
  }

  deriv2[(length(lambda_T_j)+length(lambda_C_j)+p+q+2),(length(lambda_T_j)+length(lambda_C_j)+p+q+2)] <- mean(colSums( 0.5*theta^(-2) - (theta^(-3))*w_k_grupo^(2)))  #derivada de theta da derivada de theta

  return((as.matrix(deriv2)))
}

model_MEP_dep <-  function(formula, data, delta_t, delta_c, ident, Num_intervals){

  formula <- Formula::Formula(formula)
  mf <- stats::model.frame(formula=formula, data=data)
  Terms <- stats::terms(mf)
  Z <- stats::model.matrix(formula, data = mf, rhs = 1)
  X <- stats::model.matrix(formula, data = mf, rhs = 2)
  X_T <- Z[,-1]
  X_T <- as.matrix(X_T)
  X_C <- X[,-1]
  X_C <- as.matrix(X_C)
  Xlabels <- colnames(X_T)
  Zlabels <- colnames(X_C)
  time <- stats::model.response(mf)
  q <- ncol(X_C)
  p <- ncol(X_T)
  n <- nrow(data)
  m <- max(ident)

  ## Erros mais comuns
  if(any(time < 0 | time == 0)) stop('time must be greater than 0')

  ###----------------------------------------------------------------------------------------------------
  # chute inicial para os betas_T,betas_C, betas_R, alpha2,alpha3 e sigma2
  control = coxph.control()
  ajuste_coxph_T <- coxph(Surv(time, delta_t) ~  X_T, method="breslow")
  #risco_a_T <- basehaz(ajuste_coxph_T, centered=FALSE)$hazard  #cumulative hazard
  beta_T <- ajuste_coxph_T$coef
  # beta_T <- rep(0.1,p)

  ajuste_coxph_C <- coxph(Surv(time, delta_c) ~  X_C, method="breslow")
  #risco_a_C <- basehaz(ajuste_coxph_C, centered=FALSE)$hazard  #cumulative hazard
  beta_C <- ajuste_coxph_C$coef
  # beta_C <- rep(0.1,q)
  ###----------------------------------------------------------------------------------------------------

  alpha <- 0
  sigma2 <- 1
  bmax <- Num_intervals #numero de intervalos
  lambda_T_j <- rep(0.1,bmax)
  lambda_C_j <- rep(0.1,bmax)

  param <- c(beta_T,lambda_T_j,beta_C,alpha,lambda_C_j,sigma2)

  risco_a_T <- rep(0.05,n)
  risco_a_C <- rep(0.05,n)

  ###----------------------------------------------------------------------------------------------------
  # Especificacoes do algorimo EMMC

  maxit <- 100 #numero maximo de iteracoes
  eps1= rep(1e-3, length(param))
  eps2= rep(1e-4, length(param))
  n_intMCs = c(rep(10,40),rep(25,30), rep(50,20), rep(100,10))

  ## Iniciando os objetos utilizados
  out =  matrix(NA,maxit+1,length(param))
  dif =  matrix(NA,maxit+1,length(param))
  final = length(param)
  count = rep(0,maxit+1)
  s=1
  continue=TRUE


  while (continue == TRUE) {

    count = rep(0,maxit+1)
    out[s,] =  c(beta_T,lambda_T_j,beta_C,alpha,lambda_C_j,sigma2)
    n_intMC = n_intMCs[s]
    w_chapeu_grupo <- matrix(NA, m, n_intMC)

    for (k in 1:m){
      w_trans <- arms(0.5, myldens=log_cond_wi, indFunc=support_wi,n.sample=n_intMC,Alpha=alpha,delta_T=delta_t[ident==k],delta_C=delta_c[ident==k],X_T=X_T[ident==k,], X_C=X_C[ident==k,],Betas_T=beta_T,Betas_C=beta_C, risco_a_T=risco_a_T[ident==k],risco_a_C=risco_a_C[ident==k],Sigma2=sigma2)
      w_auxi <- log(w_trans)-log(1-w_trans)
      w_chapeu_grupo[k,] <- w_auxi
    }
    bi = w_chapeu_grupo[ident,]
    sigma2 <- mean(w_chapeu_grupo^2)

    ###------------------------------------------------------------------------------------
    # estimando a taxa de falha para tempo de falha
    a_T <- time_grid(time=time, event=delta_t, n.int=bmax)
    a_inf_T   <- a_T[[1]]
    a_s_inf_T <- a_T[[2]]
    num_int_T  <- (length(a_T[[1]])-1)
    id_T <- as.numeric(cut(time,a_inf_T))  # id: grade com inf na ultima casela
    U_T <- rep(0,length(unique(id_T))-1)
    nu_T <- tapply(delta_t,id_T ,sum)
    pred_linear_T <- exp(X_T%*%beta_T)*rowMeans(exp(bi[,]))
    A_T <- ID_a(a_s_inf_T,U_T,time)

    xi <- NULL
    for(j in 1:num_int_T){
      y <- time
      y[id_T<j] <- A_T[j]
      y[id_T>j] <- A_T[j+1]
      xi[j] <- sum((y-A_T[j])*pred_linear_T)
    }
    sumX_T <- xi
    lambda_T_j <- nu_T/sumX_T

    ###-----------------------------------------------------------------------------------
    # estimando a taxa de falha para tempo de censura
    a_C <- time_grid(time=time, event=delta_c, n.int=bmax)
    a_inf_C   <- a_C[[1]]
    a_s_inf_C <- a_C[[2]]
    num_int_C  <- (length(a_C[[1]])-1)
    id_C <- as.numeric(cut(time,a_inf_C))
    U_C <- rep(0,length(unique(id_C))-1)
    nu_C <- tapply(delta_c,id_C ,sum)
    pred_linear_C <- exp(X_C%*%beta_C)*rowMeans(exp(alpha*bi[,]))
    A_C <- ID_a(a_s_inf_C,U_C,time)

    xi <- NULL
    for(j in 1:num_int_C){
      y <- time
      y[id_C<j] <- A_C[j]
      y[id_C>j] <- A_C[j+1]
      xi[j] <- sum((y-A_C[j])*pred_linear_C)
    }
    sumX_C <- xi
    lambda_C_j <- nu_C/sumX_C

    ###----------------------------------------------------------------------
    # calculo do risco acumulado
    deltaij_T <- deltas(A_T,time,id_T,num_int_T,n)
    risco_a_T <- H0ti(a_inf_T,time,deltaij_T,lambda_T_j,num_int_T,n)

    deltaij_C <- deltas(A_C,time,id_C,num_int_C,n)
    risco_a_C <- H0ti(a_inf_C,time,deltaij_C,lambda_C_j,num_int_C,n)

    ###---------------------------------------------------------------------------
    ### estimacao dos coeficientes de regressao para os tempos de falha
    S_T <- multiroot(f = modelo_T_MEP, start = rep(0.1,p),  X_T=X_T, delta.t=delta_t,risco_a_T=risco_a_T,bi=bi,n_intMC=n_intMC)
    beta_T <- S_T$root

    ###----------------------------------------------------------------------------------------
    ### estimacao dos coeficientes de regressao para os tempos de censura dependente e para alpha
    S_C <- multiroot(f = modelo_C_MEP, start = rep(0.1,q+1), X_C=X_C,delta.c=delta_c,risco_a_C=risco_a_C,bi=bi,n_intMC=n_intMC)
    beta_C <- S_C$root[1:q]
    alpha <- S_C$root[q+1]

    ###--------------------------------------------------------------------------------
    # criterio de parada
    out[s+1,]= c(beta_T,lambda_T_j,beta_C,alpha,lambda_C_j,sigma2)
    dif[s,] <- (abs(out[s+1,]-out[s,]))/(abs(out[s,])-eps1)

    for (r in 1:length(param)){
      if (dif[s,r]<eps2[r]) {
        count[s] = count[s] + 1
      }
    }
    s=s+1
    if (s>3) {
      continue=((s <= maxit) & (count[s] < length(param))&(count[s-1] < length(param))&(count[s-2] < length(param)))
    }
  } ## fim do EMMC

  param_est <- c((out[s,]+out[s-1,]+out[s-2,])/3)

  ###--------------------------------------------------------------------------------
  # calculo das derivadas de primeira e segunda ordem, para calcular o erro padrao
  Esp_deriv_ordem1 <- Esp_Deriv_MEP( X_T=X_T, X_C=X_C,delta_T=delta_t,delta_C=delta_c, beta_T=beta_T, beta_C=beta_C, alpha=alpha,
                                     theta=sigma2,nu_ind_j_T=nu_T,nu_ind_j_C=nu_C, lambda_T_j=lambda_T_j,lambda_C_j=lambda_C_j, risco_a_T=risco_a_T,
                                     risco_a_C=risco_a_C,deltaij_T=deltaij_T,deltaij_C=deltaij_C, w_k_grupo=w_chapeu_grupo, ident=ident)

  Esp_deriv_ordem2 <- Esp_DerivParc_MEP( X_T=X_T, X_C=X_C,delta_T=delta_t,delta_C=delta_c, beta_T=beta_T, beta_C=beta_C, alpha=alpha,
                                         theta=sigma2,nu_ind_j_T=nu_T,nu_ind_j_C=nu_C, lambda_T_j=lambda_T_j,lambda_C_j=lambda_C_j, risco_a_T=risco_a_T,
                                         risco_a_C=risco_a_C,deltaij_T=deltaij_T,deltaij_C=deltaij_C, w_k_grupo=w_chapeu_grupo, ident=ident)


  InfFisher <- (-Esp_deriv_ordem2) - Esp_deriv_ordem1
  Var <- solve(InfFisher)
  ErroPadrao <- sqrt(diag(Var))
  if (any(is.na(ErroPadrao))) warning("The algorithm did not converge. It might converge if you run the function again.", call. = FALSE)

  # p-value
  p_value_alpha <- 2*stats::pnorm(-abs(param_est[(p+bmax+q+1)]/ErroPadrao[(p+bmax+q+1)]))
  p_value_t <- vector()
  for (i in 1:p){
    p_value_t[i] <- 2*stats::pnorm(-abs((param_est[c(i)])/(ErroPadrao[c(i)])))
  }
  p_value_c <- vector()
  for (j in 1:q){
    p_value_c[j] <- 2*stats::pnorm(-abs((param_est[c(p+bmax+j)])/(ErroPadrao[c(p+bmax+j)])))
  }

  # information criteria
  criterios <- crit_mep(X_T,X_C,bi,n,beta_T,lambda_T_j,beta_C,alpha,lambda_C_j,delta_t,delta_c,risco_a_T,risco_a_C,id_T,id_C,m,ident)

  # creating the list/object
  fit <- param_est
  fit <- list(fit=fit)
  fit$stde <- ErroPadrao
  fit$crit <- criterios
  fit$pvalue <- c(p_value_alpha, p_value_t, p_value_c)
  fit$n <- n
  fit$p <- p
  fit$q <- q
  fit$call <- match.call()
  fit$formula <- stats::formula(Terms)
  fit$terms <- stats::terms.formula(formula)
  fit$labels1 <- Zlabels
  fit$labels2 <- Xlabels
  fit$bmax <- bmax
  fit$risco_a_T <- risco_a_T
  fit$risco_a_C <- risco_a_C
  fit$bi <- bi
  fit$X_T <- X_T
  fit$X_C <- X_C
  fit$time <- time
  class(fit) <- "dcensoring"
  return(fit)
}
