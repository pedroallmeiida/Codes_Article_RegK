library(robustbase)
library(maxLik)
#library(smco)

  
#  ++++++ <Generation and density for unconditional K model> ++++++


geradorK = function(N,theta,L) {
  # if(any(theta <= 0)) stop('all parameters must be > 0')
  alpha = theta[1]
  lambda = theta[2]
  x=c()
  for(i in 1:N) x[i] = rgamma(1, L, L) * rgamma(1, alpha, lambda)
  return(x)
}
# > RNG for the K distribution <
generatorK = function(N,theta,L) {
  # if(any(theta <= 0)) stop('all parameters must be > 0')
  alpha = theta[1]
  lambda = theta[2]
  
  rgamma(N, L, L) * rgamma(N, alpha, lambda)
}
# > pdf <
Kpdf <- function(y, theta,L) {
  # if(any(c(y,theta) < 0)) stop('y and theta must be > 0')
  
  alpha <- theta[1]
  lambda <- theta[2]
  phi <- L*lambda*y
  
  2/y * (phi)^((L+alpha)/2)/(gamma(L) * gamma(alpha)) *
    besselK(2 * sqrt(phi), alpha - L)
}
#  ++++++ <generation and density for conditional K model> ++++++

# > For simple regression <

# <function for BFGS>
loglink_KBR = function(Xobs, matrix_model, par,L){
  alpha=par[1]; alpha=abs(alpha);
  bet=par[2]
  # L: variavel global (numero de looks)
  # Xobs: variavel global (vetor de dados)
  # matrix_model: variavel global (matriz modelo)
  # Xobs: variavel global (vetor de respostas)
  auxi = matrix_model*bet;
  
  res <- log(2) +
    0.5 * (L + alpha) * log(L) -
    lgamma(alpha) - lgamma(L) +
    0.5 * (L + alpha) * log(alpha) -
    0.5 * (L + alpha) * bet * mean(matrix_model) +
    0.5 * (L + alpha-2) * mean(log(Xobs)) +
    mean(
      log(
        besselK(
          2 * sqrt( L * alpha * Xobs / exp(auxi) ),
          alpha - L)
      ),na.rm=T
    )
  
  return(res)  
}

# <function for PSO/...>
loglink_KBR_0 = function(Xobs, matrix_model, par,L) - loglink_KBR(Xobs, matrix_model, par,L)

Chute.inicial <- function(Y,X,L){
  alpha = seq(1,10,0.5)
  b1 = solve( t(X)%*%X   )%*%t(X)%*%log(Y)
  if(b1 > 10){
    b1 = 10
  }else if( b1 < -10 ){
    b1 = -10
  }  
  par = cbind( alpha, rep(b1,length(alpha)) )
  
  lk=c()
  for(i in 1:length(alpha)) lk[i] = loglink_KBR( Y, X, par[i,], L  )
  
  a = cbind( alpha, lk)
  chute.a = subset( a, a[,2] == max(a[,2])  )[1]
  return( chute = c(alpha = chute.a, b1=b1  ) )
  
}

cml.K <- function(Xobs, matrix_model, par, L, method){
  # L: variavel global (numero de looks)
  # matrix_model: variavel global (matriz modelo)
  # Xobs: variavel global (vetor de respostas)
  # par: chute inicial
  #method: metodo de otimizacao: BFGS ou SMCO
  A = matrix(c(1,0,-1,0,0,1,0,-1),4,2, byrow = T)
  B = c(-0.1, 100,20,20)
  
  lk.new = function(par)  loglink_KBR(Xobs,matrix_model, par,L)
  lk.new_0 = function(par)  loglink_KBR_0(Xobs,matrix_model, par,L)
  
  if( method == "BFGS" ){
  EST = maxLik(lk.new,start=par,method="BFGS", constraints=list(ineqA = A, ineqB = B))
  Estimate = EST$estimate; Code = EST$code

  }else if( method == "smco"){
    EST=smco(par = par,
             LB = c(0.1,-10), UB = c(50,10),
             fn = lk.new_0,
             maxiter = 10000, trc = FALSE)
    Estimate = EST$par; Code = 0
  }
  
  beta=Estimate[2]
  mu.hat=exp( beta%*%matrix_model)
  
  return(list( estimate = Estimate, code = Code, mu.hat=mu.hat  ))
  
}


R2.k <- function(Y,X,par,L){
  N=length(Y)
  lk.est = loglink_KBR(Y, X, par,L)
  lk0 = loglink_KBR(Y, X, c(par[1],0),L)
  R2 = 1 - exp(  -(2/N)*(lk.est - lk0)  )
  return( R2 = R2 )
}

fR2 <- function(Y,Y.hat){
  k=1
  N=length(Y)
  SQT = sum( (Y - mean(Y))^2 )
  SQR = sum( (Y - Y.hat)^2 )

  R2 = 1 - SQR/SQT
  R2.ajust = 1 - ((N-1)/(N - (k)))*(1-R2)
  return( list(R2 = R2, R2.ajust = R2.ajust) )
}

DistK = function( alpha, L, mu, y ){
  X =  ( 1/gamma(L)*gamma(alpha) ) * meijerg( alpha, L, mu, y  )[,1]
  return(X)
}

