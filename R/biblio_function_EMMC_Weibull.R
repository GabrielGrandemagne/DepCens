
###------------------------------------------------------------------------------------------------
# Funcao para gerar da log condicional completa da fragilidade, wk (sendo k o indice de grupo)
###------------------------------------------------------------------------------------------------

log_cond_wi <- function(w_k_aux,Alpha,delta_T,delta_C,X_T,X_C, Betas_T,Betas_C,risco_a_T,risco_a_C,Sigma2){


  w_k <- log(w_k_aux/(1-w_k_aux))

  if (ncol(t(t(X_T))) == 1)
  {
    pred_linear_T <- exp((t(t(X_T))*Betas_T)+w_k)
    pred_linear_C <- exp((t(t(X_C))*Betas_C)+Alpha*w_k)
  }

  else
  {
    pred_linear_T <- exp((X_T[,]%*%Betas_T)+w_k)
    pred_linear_C <- exp((X_C[,]%*%Betas_C)+Alpha*w_k)
  }

  log_vero_w <- sum(delta_T*(w_k) - risco_a_T*pred_linear_T + delta_C*Alpha*(w_k) - risco_a_C*pred_linear_C) -((w_k^2)/(2*Sigma2))- log(w_k_aux*(1-w_k_aux))

  return(log_vero_w)
}

# funcao utilizada no arms
support_wi <-  function(w_k_aux,Alpha,delta_T,delta_C,X_T,X_C,Betas_T,Betas_C,risco_a_T,risco_a_C,Sigma2){(w_k_aux>0)*(w_k_aux<1)}

###-----------------------------------------------------------------------------------
# Funcoes para estimar os parametros dos tempos de falha, funcao usando multiroot
###-----------------------------------------------------------------------------------

modelo_T_Weibull <-  function(param_t, t, X_T, delta.t, bi)
{
  risco_a_T <- (t^exp(param_t[1]))*exp(param_t[2])

  if(ncol(t(t(X_T))) == 1)
  {
    pred_T <- exp(t(t(X_T))*param_t[3])*rowMeans(exp(bi[,])) #MUDAR MULT VETOR
    w_kl_beta_T <- risco_a_T*pred_T
    w_kl_beta_T_num <- cbind(w_kl_beta_T)*X_T
  }

  if(ncol(X_T) == 2)
  {
    pred_T <- exp(X_T[,]%*%param_t[3:4])*rowMeans(exp(bi[,]))
    w_kl_beta_T <- risco_a_T*pred_T
    w_kl_beta_T_num <- cbind(w_kl_beta_T, w_kl_beta_T)*X_T
  }

  if(ncol(X_T) == 3)
  {
    pred_T <- exp(X_T[,]%*%param_t[3:5])*rowMeans(exp(bi[,]))
    w_kl_beta_T <- risco_a_T*pred_T
    w_kl_beta_T_num <- cbind(w_kl_beta_T, w_kl_beta_T, w_kl_beta_T)*X_T
  }

  if(ncol(X_T) == 4)
  {
    pred_T <- exp(X_T[,]%*%param_t[3:6])*rowMeans(exp(bi[,]))
    w_kl_beta_T <- risco_a_T*pred_T
    w_kl_beta_T_num <- cbind(w_kl_beta_T, w_kl_beta_T, w_kl_beta_T, w_kl_beta_T)*X_T
  }

  if(ncol(X_T) == 5)
  {
    pred_T <- exp(X_T[,]%*%param_t[3:7])*rowMeans(exp(bi[,]))
    w_kl_beta_T <- risco_a_T*pred_T
    w_kl_beta_T_num <- cbind(w_kl_beta_T, w_kl_beta_T, w_kl_beta_T, w_kl_beta_T, w_kl_beta_T)*X_T
  }

  if(ncol(X_T) == 6)
  {
    pred_T <- exp(X_T[,]%*%param_t[3:8])*rowMeans(exp(bi[,]))
    w_kl_beta_T <- risco_a_T*pred_T
    w_kl_beta_T_num <- cbind(w_kl_beta_T, w_kl_beta_T, w_kl_beta_T, w_kl_beta_T, w_kl_beta_T, w_kl_beta_T)*X_T
  }

  if(ncol(X_T) == 7)
  {
    pred_T <- exp(X_T[,]%*%param_t[3:9])*rowMeans(exp(bi[,]))
    w_kl_beta_T <- risco_a_T*pred_T
    w_kl_beta_T_num <- cbind(w_kl_beta_T, w_kl_beta_T, w_kl_beta_T, w_kl_beta_T, w_kl_beta_T, w_kl_beta_T, w_kl_beta_T)*X_T
  }

  if(ncol(X_T) == 8)
  {
    pred_T <- exp(X_T[,]%*%param_t[3:10])*rowMeans(exp(bi[,]))
    w_kl_beta_T <- risco_a_T*pred_T
    w_kl_beta_T_num <- cbind(w_kl_beta_T, w_kl_beta_T, w_kl_beta_T, w_kl_beta_T, w_kl_beta_T, w_kl_beta_T, w_kl_beta_T, w_kl_beta_T)*X_T
  }

  if(ncol(X_T) == 9)
  {
    pred_T <- exp(X_T[,]%*%param_t[3:11])*rowMeans(exp(bi[,]))
    w_kl_beta_T <- risco_a_T*pred_T
    w_kl_beta_T_num <- cbind(w_kl_beta_T, w_kl_beta_T, w_kl_beta_T, w_kl_beta_T, w_kl_beta_T, w_kl_beta_T, w_kl_beta_T, w_kl_beta_T, w_kl_beta_T)*X_T
  }

  if(ncol(X_T) == 10)
  {
    pred_T <- exp(X_T[,]%*%param_t[3:12])*rowMeans(exp(bi[,]))
    w_kl_beta_T <- risco_a_T*pred_T
    w_kl_beta_T_num <- cbind(w_kl_beta_T, w_kl_beta_T, w_kl_beta_T, w_kl_beta_T, w_kl_beta_T, w_kl_beta_T, w_kl_beta_T, w_kl_beta_T, w_kl_beta_T, w_kl_beta_T)*X_T
  }

  U_T_1 <- colSums(X_T*delta.t - w_kl_beta_T_num)
  w_kl_alpha_T_num <- w_kl_beta_T*log(t)*exp(param_t[1])
  U_alphaT <-  colSums(delta.t*(1 + exp(param_t[1])*log(t)) - w_kl_alpha_T_num)

  w_kl_lambda_T_num  <- w_kl_beta_T
  U_lambdaT <-  colSums(delta.t*(1) - w_kl_lambda_T_num)

  c(U_T_1 = U_T_1,U_alphaT=U_alphaT,U_lambdaT=U_lambdaT)
}

###-----------------------------------------------------------------------------------------------
# Funcoes para estimar os parametros dos tempos de censura, funcao usando multiroot
###-----------------------------------------------------------------------------------------------

modelo_C_Weibull <-  function(param_c, t, X_C, delta.c, bi)
{
  risco_a_C <- (t^exp(param_c[1]))*exp(param_c[2])

  if(ncol(t(t(X_T))) == 1)
  {
    pred_C <- exp(t(t(X_C))*param_c[4])*rowMeans(exp(param_c[3]*bi[,]))
    w_kl_beta_C <- risco_a_C*pred_C
    w_kl_beta_C_num <- cbind(w_kl_beta_C)*X_C
    w_kl_alpha_num <- risco_a_C*exp(t(t(X_C))*param_c[4])*rowMeans(bi*exp(param_c[3]*bi[,]))
  }

  if(ncol(X_C) == 2)
  {
    pred_C <- exp(X_C[,]%*%param_c[4:5])*rowMeans(exp(param_c[3]*bi[,]))
    w_kl_beta_C <- risco_a_C*pred_C
    w_kl_beta_C_num <- cbind(w_kl_beta_C, w_kl_beta_C)*X_C
    w_kl_alpha_num <- risco_a_C*exp(X_C[,]%*%param_c[4:5])*rowMeans(bi*exp(param_c[3]*bi[,]))
  }

  if(ncol(X_C) == 3)
  {
    pred_C <- exp(X_C[,]%*%param_c[4:6])*rowMeans(exp(param_c[3]*bi[,]))
    w_kl_beta_C <- risco_a_C*pred_C
    w_kl_beta_C_num <- cbind(w_kl_beta_C, w_kl_beta_C, w_kl_beta_C)*X_C
    w_kl_alpha_num <- risco_a_C*exp(X_C[,]%*%param_c[4:6])*rowMeans(bi*exp(param_c[3]*bi[,]))
  }

  if(ncol(X_C) == 4)
  {
    pred_C <- exp(X_C[,]%*%param_c[4:7])*rowMeans(exp(param_c[3]*bi[,]))
    w_kl_beta_C <- risco_a_C*pred_C
    w_kl_beta_C_num <- cbind(w_kl_beta_C, w_kl_beta_C, w_kl_beta_C, w_kl_beta_C)*X_C
    w_kl_alpha_num <- risco_a_C*exp(X_C[,]%*%param_c[4:7])*rowMeans(bi*exp(param_c[3]*bi[,]))
  }

  if(ncol(X_C) == 5)
  {
    pred_C <- exp(X_C[,]%*%param_c[4:8])*rowMeans(exp(param_c[3]*bi[,]))
    w_kl_beta_C <- risco_a_C*pred_C
    w_kl_beta_C_num <- cbind(w_kl_beta_C, w_kl_beta_C, w_kl_beta_C, w_kl_beta_C, w_kl_beta_C)*X_C
    w_kl_alpha_num <- risco_a_C*exp(X_C[,]%*%param_c[4:8])*rowMeans(bi*exp(param_c[3]*bi[,]))
  }

  if(ncol(X_C) == 6)
  {
    pred_C <- exp(X_C[,]%*%param_c[4:9])*rowMeans(exp(param_c[3]*bi[,]))
    w_kl_beta_C <- risco_a_C*pred_C
    w_kl_beta_C_num <- cbind(w_kl_beta_C, w_kl_beta_C, w_kl_beta_C, w_kl_beta_C, w_kl_beta_C, w_kl_beta_C)*X_C
    w_kl_alpha_num <- risco_a_C*exp(X_C[,]%*%param_c[4:9])*rowMeans(bi*exp(param_c[3]*bi[,]))
  }

  if(ncol(X_C) == 7)
  {
    pred_C <- exp(X_C[,]%*%param_c[4:10])*rowMeans(exp(param_c[3]*bi[,]))
    w_kl_beta_C <- risco_a_C*pred_C
    w_kl_beta_C_num <- cbind(w_kl_beta_C, w_kl_beta_C, w_kl_beta_C, w_kl_beta_C, w_kl_beta_C, w_kl_beta_C, w_kl_beta_C)*X_C
    w_kl_alpha_num <- risco_a_C*exp(X_C[,]%*%param_c[4:10])*rowMeans(bi*exp(param_c[3]*bi[,]))
  }

  if(ncol(X_C) == 8)
  {
    pred_C <- exp(X_C[,]%*%param_c[4:11])*rowMeans(exp(param_c[3]*bi[,]))
    w_kl_beta_C <- risco_a_C*pred_C
    w_kl_beta_C_num <- cbind(w_kl_beta_C, w_kl_beta_C, w_kl_beta_C, w_kl_beta_C, w_kl_beta_C, w_kl_beta_C, w_kl_beta_C, w_kl_beta_C)*X_C
    w_kl_alpha_num <- risco_a_C*exp(X_C[,]%*%param_c[4:11])*rowMeans(bi*exp(param_c[3]*bi[,]))
  }

  if(ncol(X_C) == 9)
  {
    pred_C <- exp(X_C[,]%*%param_c[4:12])*rowMeans(exp(param_c[3]*bi[,]))
    w_kl_beta_C <- risco_a_C*pred_C
    w_kl_beta_C_num <- cbind(w_kl_beta_C, w_kl_beta_C, w_kl_beta_C, w_kl_beta_C, w_kl_beta_C, w_kl_beta_C, w_kl_beta_C, w_kl_beta_C, w_kl_beta_C)*X_C
    w_kl_alpha_num <- risco_a_C*exp(X_C[,]%*%param_c[4:12])*rowMeans(bi*exp(param_c[3]*bi[,]))
  }

  if(ncol(X_C) == 10)
  {
    pred_C <- exp(X_C[,]%*%param_c[4:13])*rowMeans(exp(param_c[3]*bi[,]))
    w_kl_beta_C <- risco_a_C*pred_C
    w_kl_beta_C_num <- cbind(w_kl_beta_C, w_kl_beta_C, w_kl_beta_C, w_kl_beta_C, w_kl_beta_C, w_kl_beta_C, w_kl_beta_C, w_kl_beta_C, w_kl_beta_C, w_kl_beta_C)*X_C
    w_kl_alpha_num <- risco_a_C*exp(X_C[,]%*%param_c[4:13])*rowMeans(bi*exp(param_c[3]*bi[,]))
  }

  U_betas <- colSums(X_C*delta.c - w_kl_beta_C_num)
  U_alpha <- colSums(delta.c*rowMeans(bi[,]) - w_kl_alpha_num)
  w_kl_alpha_C_num <- w_kl_beta_C*log(t)*exp(param_c[1])
  U_alphaC <-  colSums(delta.c*(1 + exp(param_c[1])*log(t)) - w_kl_alpha_C_num)

  w_kl_lambda_C_num   <- w_kl_beta_C
  U_lambdaC <-  colSums(delta.c*(1) - w_kl_lambda_C_num)

  c(U_alphaC=U_alphaC,U_lambdaC=U_lambdaC,U_alpha=U_alpha,U_betas = U_betas)
}

###-----------------------------------------------------------------------
# Funcoes para calcular o vetor com as derivadas de primeira ordem
###-----------------------------------------------------------------------

Esp_DerivPrimOrdem <-  function(t,delta_T, delta_C, X_T, X_C, beta_T, beta_C,alpha,Sigma2, alpha_T, alpha_C,lambda_T, lambda_C, w_k_grupo, ident){
  wk = w_k_grupo[ident,]
  num_param <- length(beta_T)+2+length(beta_C)+2+length(alpha)+length(Sigma2)
  deriv1 <- matrix(NA,num_param,ncol(wk))  #vetor que com as derivadas de primeira ordem
  p <- ncol(X_T)
  q <- ncol(X_C)
  pred_T <- as.vector(exp(X_T%*%beta_T))
  pred_C <- as.vector(exp(X_C%*%beta_C))
  Delta_t <- delta_T%*%t(rep(1,ncol(wk)))
  Delta_c <- delta_C%*%t(rep(1,ncol(wk)))

  deriv1[1,] <- colSums(Delta_t*(alpha_T^(-1) + log(t)) - (t^(alpha_T))*log(t)*lambda_T*pred_T*exp(wk))

  deriv1[2,] <- colSums(Delta_t*(lambda_T^(-1)) - (t^(alpha_T))*pred_T*exp(wk))

  for (i in 1:p){
    deriv1[2+i,] <- colSums(Delta_t*(X_T[,i]) - (t^(alpha_T))*lambda_T*pred_T*X_T[,i]*exp(wk))
  }

  deriv1[3+p,] <- colSums(Delta_c*(alpha_C^(-1) + log(t)) - (t^(alpha_C))*log(t)*lambda_C*pred_C*exp(alpha*wk))

  deriv1[4+p,] <- colSums(Delta_c*(lambda_C^(-1)) - (t^(alpha_C))*pred_C*exp(alpha*wk))

  for (i in 1:q){
    deriv1[4+p+i,] <- colSums(Delta_c*(X_C[,i]) - (t^(alpha_C))*lambda_C*pred_C*X_C[,i]*exp(alpha*wk))
  }

  deriv1[4+p+q+1,] <- colSums(Delta_c*(wk) - (t^(alpha_C))*lambda_C*pred_C*exp(alpha*wk)*wk)

  deriv1[4+p+q+2,] <- colSums(-0.5*Sigma2^(-1) + 0.5*(Sigma2^(-2))*w_k_grupo^(2))

  aux <- deriv1[,1]%*%t(deriv1[,1])

  for( i in 2:ncol(wk)){
    aux <- aux + deriv1[,i]%*%t(deriv1[,i])
  }

  return(aux/ncol(wk))
}

###----------------------------------------------------------------------
# Funcoes para calcular o vetor com as derivadas de segunda ordem
###----------------------------------------------------------------------

Esp_DerivParciais <-  function(t,delta_T, delta_C, X_T, X_C, beta_T, beta_C, alpha,Sigma2=sigma2, alpha_T, alpha_C, lambda_T, lambda_C,w_k_grupo, ident){

  num_param <- length(beta_T)+2+length(beta_C)+2+length(alpha)+length(Sigma2)
  deriv2 <- matrix(0,num_param,num_param)  #vetor que com as derivadas de segunda ordem,
  wk = w_k_grupo[ident,]
  p <- ncol(X_T)
  q <- ncol(X_C)
  pred_T <- as.vector(exp(X_T%*%beta_T))
  pred_C <- as.vector(exp(X_C%*%beta_C))
  Delta_t <- delta_T%*%t(rep(1,ncol(wk)))
  Delta_c <- delta_C%*%t(rep(1,ncol(wk)))

  deriv2[1,1] <- mean(colSums(Delta_t*(-alpha_T^(-2)) - (t^(alpha_T))*(log(t)^(2))*lambda_T*pred_T*exp(wk)))

  deriv2[1,2] <- mean(colSums(- (t^(alpha_T))*(log(t))*pred_T*exp(wk)))
  deriv2[2,1] <- deriv2[1,2]

  for (i in 1:p){
    deriv2[1,2+i] <- mean(colSums(- (t^(alpha_T))*(log(t))*lambda_T*pred_T*X_T[,i]*exp(wk)))
    deriv2[2+i,1] <- deriv2[1,2+i]
  }

  deriv2[2,2] <- mean(colSums( Delta_t*(-lambda_T^(-2))))

  for (i in 1:p){
    deriv2[2,2+i] <- mean(colSums( - (t^(alpha_T))*pred_T*X_T[,i]*exp(wk)))
    deriv2[2+i,2] <- deriv2[2,2+i]
  }

  for (i in 1:p){
    deriv2[2+i,2+i] <- mean(colSums( - ((t^(alpha_T))*lambda_T*pred_T*(X_T[,i]^(2))*exp(wk))))
  }

  for (j in 1:(p-1)){
    for (i in 1:(p-1)){
      if (ncol(t(t(X_T))) == 1){
        next
      }
      deriv2[2+j,3+i] <- mean(colSums(- (t^(alpha_T))*lambda_T*pred_T*X_T[,1+i]*X_T[,j]*exp(wk)))
      deriv2[3+i,2+j] <- deriv2[2+j,3+i]
    }
  }

  deriv2[3+p,3+p] <- mean(colSums(Delta_c*(-alpha_C^(-2)) - (t^(alpha_C))*(log(t)^(2))*lambda_C*pred_C*exp(alpha*wk)))

  deriv2[3+p,4+p] <- mean(colSums( - (t^(alpha_C))*log(t)*pred_C*exp(alpha*wk)))
  deriv2[4+p,3+p] <- deriv2[3+p,4+p]

  for (i in 1:q){
    deriv2[3+p,(4+p+i)] <- mean(colSums( - (t^(alpha_C))*log(t)*lambda_C*pred_C*X_C[,i]*exp(alpha*wk)))
    deriv2[(4+p+i),3+p] <- deriv2[3+p,(4+p+i)]
  }

  deriv2[3+p,(4+p+q+1)] <- mean(colSums( - (t^(alpha_C))*log(t)*lambda_C*pred_C*exp(alpha*wk)*wk))
  deriv2[(4+p+q+1),3+p] <- deriv2[3+p,(4+p+q+1)]

  deriv2[4+p,4+p] <- mean(colSums(Delta_c*(-lambda_C^(-2)) ))

  for (i in 1:q){
    deriv2[4+p,(4+p+i)] <- mean(colSums( - (t^(alpha_C))*pred_C*X_C[,i]*exp(alpha*wk)))
    deriv2[(4+p+i),4+p] <- deriv2[4+p,(4+p+i)]
  }

  deriv2[4+p,(4+p+q+1)] <- mean(colSums( - (t^(alpha_C))*pred_C*exp(alpha*wk)*wk))
  deriv2[(4+p+q+1),4+p] <- deriv2[4+p,(4+p+q+1)]

  for (i in 1:q){
    deriv2[4+p+i,4+p+i] <- mean(colSums(-((t^(alpha_C))*lambda_C*pred_C*(X_C[,i]^(2))*exp(alpha*wk))))
  }

  for (j in 1:(q-1)){
    for (i in 1:(q-1)){
      if (ncol(t(t(X_T))) == 1){
        next
      }
      deriv2[4+p+j,4+p+1+i] <- mean(colSums( - (t^(alpha_C))*lambda_C*pred_C*X_C[,1+i]*X_C[,j]*exp(alpha*wk)))
      deriv2[4+p+1+i,4+p+j] <- deriv2[4+p+j,4+p+1+i]
    }
  }

  for (i in 1:q){
    deriv2[4+p+i,(4+p+q+1)] <- mean(colSums(- t^(alpha_C)*lambda_C*pred_C*X_C[,i]*exp(alpha*wk)*wk))
    deriv2[(4+p+q+1),4+p+i] <- deriv2[4+p+i,(4+p+q+1)]
  }

  deriv2[(4+p+q+1),(4+p+q+1)] <- mean(colSums( - (t^(alpha_C))*lambda_C*pred_C*exp(alpha*wk)*wk^(2)))

  deriv2[(4+p+q+2),(4+p+q+2)] <- mean(colSums( 0.5*Sigma2^(-2) - (Sigma2^(-3))*w_k_grupo^(2)))

  return((as.matrix(deriv2)))
}

###---------------------------------------------------------------------------------------------------
# Funcao para ajustar o modelo completo
###---------------------------------------------------------------------------------------------------

model_Weibull_dep <-  function(time, delta_t, delta_c, X_T, X_C, ident){

  ## Erros mais comuns
  if(any(time < 0 | time == 0)) stop('time must be greater than 0')

  ###----------------------------------------------------------------------------------------------------
  # chute inicial para os betas_T, alpha_T, lambda_Te sigma2

  q <- ncol(X_C)
  p <- ncol(X_T)
  n <- nrow(dados)
  m <- max(ident) #numero de grupos
  ###----------------------------------------------------------------------------------------------------
  # chute inicial para os betas_T,betas_C, betas_R, alpha2,alpha3  e sigma2
  # control = coxph.control()
  # ajuste_coxph_T <- coxph(Surv(time, delta_t) ~  X_T[,1] + X_T[,2] + X_T[,3] + X_T[,4], method="breslow")
  # risco_a_T <- basehaz(ajuste_coxph_T, centered=FALSE)$hazard  #cumulative hazard
  # beta_T <- ajuste_coxph_T$coef
  beta_T <- rep(0.1,p)

  # ajuste_coxph_C <- coxph(Surv(time, delta_c) ~  X_C[,1] + X_C[,2] + X_C[,3] + X_C[,4], method="breslow")
  # risco_a_C <- basehaz(ajuste_coxph_C, centered=FALSE)$hazard  #cumulative hazard
  # beta_C <- ajuste_coxph_C$coef
  beta_C <- rep(0.1,q)

  alpha <- 0
  sigma2 <- 1

  alpha_T <- 1
  lambda_T <- 1
  alpha_C <- 1
  lambda_C <- 1

  param <- c(beta_T,alpha_T,lambda_T,beta_C,alpha,alpha_C,lambda_C,sigma2)

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
    out[s,] =  c(beta_T,alpha_T,lambda_T,beta_C,alpha,alpha_C,lambda_C,sigma2)
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
    # estimando os par?metros da dist. dos tempos de falha

    S_T <- multiroot(f = modelo_T_Weibull, start = rep(0.1, 2+p),t=time, X_T=X_T, delta.t=delta_t,bi=bi)
    param_T <- S_T$root
    alpha_T <- exp(param_T[1])
    lambda_T <- exp(param_T[2])
    beta_T <- param_T[3:(2+p)]
    risco_a_T <- (time^alpha_T)*(lambda_T)

    ###-----------------------------------------------------------------------------------
    #para censura

    S_C <- multiroot(f = modelo_C_Weibull, start = rep(0.1,3+q),t=time, X_C=X_C, delta.c=delta_c,bi=bi)
    param_C <- S_C$root
    alpha_C <- exp(param_C[1])
    lambda_C <- exp(param_C[2])
    alpha <- param_C[3]
    beta_C <- param_C[4:(3+q)]
    risco_a_C <- (time^alpha_C)*(lambda_C)


    ###-----------------------------------------------------------------------------------
    # criterio de parada
    out[s+1,]=  c(beta_T,alpha_T,lambda_T,beta_C,alpha,alpha_C,lambda_C,sigma2)
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

  param_est <- c((out[s,]+out[s-1,]+out[s-2,])/3,s)


  Esp_deriv_ordem1 <- Esp_DerivPrimOrdem(t=time,delta_T=delta_t, delta_C=delta_c,X_T=X_T, X_C=X_C,beta_T=beta_T, beta_C=beta_C,
                                         alpha=alpha, Sigma2=sigma2, alpha_T=alpha_T, alpha_C=alpha_C,
                                         lambda_T=lambda_T,lambda_C=lambda_C, w_k_grupo=w_chapeu_grupo, ident=ident)


  Esp_deriv_ordem2 <- Esp_DerivParciais(t=time,delta_T=delta_t, delta_C=delta_c, X_T=X_T, X_C=X_C,beta_T=beta_T, beta_C=beta_C,
                                        alpha=alpha, Sigma2=sigma2, alpha_T=alpha_T, alpha_C=alpha_C,lambda_T=lambda_T, lambda_C=lambda_C, w_k_grupo=w_chapeu_grupo, ident=ident)


  InfFisher <- (-Esp_deriv_ordem2) - Esp_deriv_ordem1
  Var <- solve(InfFisher)

  ErroPadrao <- sqrt(diag(Var))

  return(list(param_est, ErroPadrao))

}



