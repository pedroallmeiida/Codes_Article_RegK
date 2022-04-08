source("func_MEIJER.R")
source("Functions_RegK.R")


library(nortest)

geradorK = function(N,x,par,L) {
  # if(any(theta <= 0)) stop('all parameters must be > 0')
  alpha = par[1]
  beta = par[2:3]
  a = (  gamma(alpha+1)*gamma( L+1 )   )/(  gamma(alpha)*gamma(L)   )
  mu = exp( x%*%beta )
  lambda = (1/(mu*L))*a
  z=c()
  for(i in 1:N) z[i] = rgamma(1, shape=L, scale=1/L) * rgamma(1, shape=alpha, scale=1/lambda[i]  )
  return(z)
}
loglink_KBR = function(Xobs,X,par,L){
  alpha=par[1]; alpha=abs(alpha);
  bet0=par[2];
  bet1=par[3];
  N=length(Xobs)
  matrix_model = cbind(rep(1,length(X)),X)
  # L: variavel global (numero de looks)
  # Xobs: variavel global (vetor de dados)
  # matrix_model: variavel global (matriz modelo)
  # Xobs: variavel global (vetor de respostas)
  auxi = (matrix_model)%*%c(bet0,bet1);
  
  res <- N * log(2) +
    N * 0.5 * (L + alpha) * log(L) -
    N * lgamma(alpha) - N * lgamma(L) +
    N * 0.5 * (L + alpha) * log(alpha) -
    N * 0.5 * (L + alpha) * bet0 -
    0.5 * (L + alpha) * bet1 * sum(matrix_model[,2]) +
    0.5 * (L + alpha-2) * sum(log(Xobs)) +
    sum(
      log(
        besselK(
          2 * sqrt( L * alpha * Xobs / exp(auxi) ),
          alpha - L)
      ),na.rm=T
    )
  
  return(res)  
}
loglink_KBRb0 = function(Xobs,X,par,L){
  alpha=par[1]; alpha=abs(alpha);
  bet0=par[2];
  N=length(Xobs)
  matrix_model = cbind(rep(1,length(X)))
  # L: variavel global (numero de looks)
  # Xobs: variavel global (vetor de dados)
  # matrix_model: variavel global (matriz modelo)
  # Xobs: variavel global (vetor de respostas)
  auxi = (matrix_model)%*%c(bet0);
  
  res <- N * log(2) +
    N * 0.5 * (L + alpha) * log(L) -
    N * lgamma(alpha) - N * lgamma(L) +
    N * 0.5 * (L + alpha) * log(alpha) -
    N * 0.5 * (L + alpha) * bet0 +
    0.5 * (L + alpha-2) * sum(log(Xobs)) +
    sum(
      log(
        besselK(
          2 * sqrt( L * alpha * Xobs / exp(auxi) ),
          alpha - L)
      ),na.rm=T
    )
  
  return(res)  
}
loglink_KBRb1 = function(Xobs,X,par,L){
  # L: variavel global (numero de looks)
  # Xobs: variavel global (vetor de dados)
  # matrix_model: variavel global (matriz modelo)
  # Xobs: variavel global (vetor de respostas)
  
  alpha = par[1]; alpha = abs(alpha);
  bet1 = par[2];
  N = length(Xobs)
  matrix_model = X
  auxi = bet1%*%matrix_model;
  
  res <- N * log(2) +
    N * 0.5 * (L + alpha) * log(L) -
    N * lgamma(alpha) - N * lgamma(L) +
    N * 0.5 * (L + alpha) * log(alpha) -
    0.5 * (L + alpha) * bet1 * sum(matrix_model) +
    0.5 * (L + alpha-2) * sum(log(Xobs)) +
    sum(
      log(
        besselK(
          2 * sqrt( L * alpha * Xobs / exp(auxi) ),
          alpha - L)
      ),na.rm=T
    )
  return(res)  
}
cml.Kb0 <- function(z, x, par, L, method){
  # L: variavel global (numero de looks)
  # x: variavel global (matriz modelo)
  # z: variavel global (vetor de respostas)
  # par: chute inicial
  #method: metodo de otimizacao: BFGS ou SMCO
  #A = matrix(c(1,0),1,2, byrow = T)
  #B = c(-0.001)
  
  A = matrix(c(1,0,-1,0,0,1,0,-1),4,2, byrow = T)
  B = c(-0.01, 100,100,100)
  
  
  lk.new = function(par)  loglink_KBRb0(z,x, par,L)
  
  
    EST = maxLik(lk.new,start=par[1:2],method=method,constraints=list(ineqA = A, ineqB = B))
    Estimate = EST$estimate; Code = EST$code
  
  
  return(list( estimate = Estimate, code = Code  ))
  
}
cml.Kb1 <- function(z, x, par, L, method){
  # L: variavel global (numero de looks)
  # x: variavel global (matriz modelo)
  # z: variavel global (vetor de respostas)
  # par: chute inicial
  #method: metodo de otimizacao: BFGS ou SMCO
  #A = matrix(c(1,0),1,2, byrow = T)
  #B = c(-0.001)
  
  A = matrix(c(1,0,-1,0,0,1,0,-1),4,2, byrow = T)
  B = c(-0.01, 100,600,600)
  
  
  lk.new = function(par)  loglink_KBRb1(z,x, par,L)
  
  
    EST = maxLik(lk.new,start=par,method=method,constraints=list(ineqA = A, ineqB = B))
    Estimate = EST$estimate; Code = EST$code
  
  return(list( estimate = Estimate, code = Code  ))
  
}
loglink_KBR_0 = function(Xobs,X,par,L){
  alpha=abs(par[1]);
  N=length(Xobs)
  # L: variavel global (numero de looks)
  # Xobs: variavel global (vetor de dados)
  # matrix_model: variavel global (matriz modelo)
  # Xobs: variavel global (vetor de respostas)
  matrix_model = cbind(rep(1,length(X)),X)
  
  res <- N * log(2) +
    N * 0.5 * (L + alpha) * log(L) -
    N * lgamma(alpha) - N * lgamma(L) +
    N * 0.5 * (L + alpha) * log(alpha) +
    #0.5 * (L + alpha) * bet * sum(matrix_model) +
    0.5 * (L + alpha-2) * sum(log(Xobs)) +
    sum(
      log(
        besselK(
          2 * sqrt( L * alpha * Xobs ),
          alpha - L)
      ),na.rm=T
    )
  
  return(res)  
}
Chute.inicial <- function(Y,x,L){
  X = cbind(rep(1,length(x)),x)
  
  
  alpha = seq(1,15,0.5)
  b1 = solve( t(X)%*%X   )%*%t(X)%*%log(Y)
  if(b1[2] > 20){
    b1[2] = 20
  }else if( b1[2] < -20 ){
    b1[2] = -20
  }  
  
  if(b1[1] > 30){
    b1[1] = 30
  }else if( b1[1] < -30 ){
    b1[1] = -30
  }  
  
  par = cbind( alpha, rep(  b1[1],length(alpha)  ),rep(  b1[2],length(alpha)  ) )
  
  lk=c()
  for(i in 1:length(alpha)) lk[i] = loglink_KBR( Y, x, par[i,], L=4  )
  
  a = cbind( alpha, lk)
  chute.a = subset( a, a[,2] == max(a[,2])  )[1]
  return( chute = c(alpha = chute.a, b1=b1  ) )
  
}
Chute.inicial_0 <- function(Y,x,L){
  X = cbind(rep(1,length(x)),x)
  
  
  alpha = seq(1,15,0.5)
  
  lk=c()
  for(i in 1:length(alpha)) lk[i] = loglink_KBR( Y, x, alpha[i], L=4  )
  
  a = cbind( alpha, lk)
  chute.a = subset( a, a[,2] == max(a[,2])  )[1]
  return( chute = c(alpha = chute.a  ) )
  
}
cml.K <- function(z, x, par, L, method){
  # L: variavel global (numero de looks)
  # x: variavel global (matriz modelo)
  # z: variavel global (vetor de respostas)
  # par: chute inicial
  #method: metodo de otimizacao: BFGS ou SMCO
  #A = matrix(c(1,0),1,2, byrow = T)
  #B = c(-0.001)
  n = length(z)
  
  A = matrix(c(1,0,0,-1,0,0,0,1,0,0,-1,0,0,0,1,0,0,-1),6,3, byrow = T)
  B = c(-0.01, 100,200,200,900,900)
  matrix_model = cbind(rep(1,length(x)),x)
  
  lk.new = function(par)  loglink_KBR(z,x, par,L)
  
  EST = maxLik(lk.new,start=par,method=method,constraints=list(ineqA = A, ineqB = B))
  Estimate = EST$estimate; Code = EST$code;eta = (matrix_model)%*%c(Estimate[2],Estimate[3]) ; mu.hat = exp( eta )

  AIC = -2*loglink_KBR(z, x, Estimate, L) + 2*length(Estimate) 
  BIC = -2*loglink_KBR(z, x, Estimate, L) + 2*log(n)
  
  return( list( estimate = Estimate, code = Code, mu.hat = as.numeric(mu.hat), AIC=AIC, BIC=BIC  ) )
  
}
cml.K_0 <- function(z, x, par, L, method){
  # L: variavel global (numero de looks)
  # x: variavel global (matriz modelo)
  # z: variavel global (vetor de respostas)
  # par: chute inicial
  #method: metodo de otimizacao: BFGS ou SMCO
  #A = matrix(c(1,0),1,2, byrow = T)
  #B = c(-0.001)
  
  A = matrix(1,1,1, byrow = T)
  B = c(0)
  
  lk.new = function(par)  loglink_KBR_0(z,x, par,L)
  
  
    EST = maxLik(lk.new,start=par[1],method=method,constraints=list(ineqA = A, ineqB = B))
    Estimate = EST$estimate; Code = EST$code
    
  
  return(list( estimate = Estimate, code = Code  ))
  
}
R2.k <- function(Y,X,par,par0,L){
  
  N=length(Y)
  lk.est = loglink_KBR(Y, X, par,L)
  lk0 = loglink_KBRb0(Y, X,par0,L)
  R2 = 1 - exp(  -(2/N)*(lk.est - lk0)  )
  return( R2 = R2 )
}
DistK = function( alpha, L, mu, y ){
  X =  ( 1/( gamma(L)*gamma(alpha) ) ) * meijerg( alpha, L, c(mu), c(y)  )[,1]
  return(X)
}
geradorK = function(N, x, par, L) {
  X = cbind( rep(1, length(x)), x  )
  alpha = par[1]
  beta = par[2:3]
  a = (  gamma(alpha+1)*gamma( L+1 )   )/(  gamma(alpha)*gamma(L)   )
  mu = exp( X%*%beta )
  lambda = (1/(mu*L))*a
  z=c()
  for(i in 1:N) z[i] = rgamma(1, shape=L, scale=1/L) * rgamma(1, shape=alpha, scale=1/lambda[i]  )
  return(z)
}
envelope = function( ajuste, Y, X, n=100, L, method  ){
  # input ->
  # est: ajuste do modelo
  # Y: vetor resposta
  # X: vetor de covariavel
  # n: quantidade de reamostras
  # L: numero de looks
  # output -> 
  # gráfico do envelope simulado
  
  mu.hat = c(ajuste$mu.hat)
  alpha = ajuste$estimate[1]
  
  F.dist = DistK( alpha,L,  mu.hat, Y  )
  res = qnorm( F.dist, lower.tail =  TRUE ) 
  
  res.sort = matrix(NA, n, length(Y))
  
  for(i in 1:n){
    y.sim = geradorK( length(Y), X, as.numeric(est$estimate) , L)
    chute = Chute.inicial( y.sim , X, L=4 )
    ajuste = cml.K( y.sim, X, chute, L, method  ) 
    mu.hat = c(ajuste$mu.hat); alpha = ajuste$estimate[1]
    F.dist = DistK( alpha, L, mu.hat, y.sim  )
    res.quantile = qnorm( F.dist, lower.tail =  TRUE ) 
    res.sort[i,] = sort(res.quantile) 
  }
  
  UI = c();LI=c()
  for(j in 1:length(Y)){
    LI[j] = quantile(res.sort[,j], probs = 0.025)
    UI[j] = quantile(res.sort[,j], probs = 0.975)
  }
  MED = colMedians(res.sort)
  
  faixay <- range(res.sort, res)
  qq0 <- qqnorm(sort(res), main = "", xlab = "Theorical Quantiles", pch = 20,
                col = "blue", ylab = "Sample Quantiles" ,ylim = faixay, xlim = c(-3,3) ,cex.lab=2, cex.axis=2)
  eixox <- qq0$x
  lines(eixox, MED)
  lines(eixox, LI)
  lines(eixox, UI) 
  #lines( sort(res ), MED  )
}
trv.test <- function( lk, lk0){
  # ---- Inputs -> 
  # Y: Vetor da variavel resposta
  # X: Vetor da variavel explicativa
  # est: vetor parâmetros estimados pelo modelo completo
  # est0: vetor de parâmetros estimado pelo modelo com intercepto
  TRV = 2*(  lk -  lk0    )
  p.value = 1 - pchisq( L, df=1 ); 
  return( list( TRV = TRV, p.value = p.value  )  )
  
}
se.par <- function(B, x, par, L){
  N = length(x)
  est = matrix(NA, B, length(par))
  code = c()
  
  for( i in 1:B ){
    error = 1;
    while( error == 1 ){  
      y <- geradorK( N, x, par, L ) 
      ajuste = cml.K(y, x, par, L, method = "BFGS" )
      error = ifelse( ajuste$code == 0, 0, 1  )
    }
    est[i,] = ajuste$est
  }
  Sx = apply(est , 2, sd)/sqrt(N)
  MEAN = apply(est , 2, mean)
  return( list(erro.padrao = Sx, mean = MEAN) )
  
}

envel.gamma <- function(modelo=fit.model,iden=0,nome=seq(along = model.matrix(modelo)[,1]),sim=100,conf=.90,res="Q",mv=F,quad=T,maxit=20) {
  
  #
  # Descrição e detalhes:
  # A saída será o gráfico de probabilidade normal com envelopes simulados para um ajuste da distribuição gama.
  #
  # A opção res="C" faz o gráfico de probabilidade meio-normal com envelopes simulados utilizando a distância de Cook,
  # possibilitando a detecção de pontos simultaneamente aberrantes e/ou influentes.
  #
  # Atenção: a função não funcionará corretamente se o ajuste possuir offsets! Neste caso é preciso adaptá-la como foi
  # feito na função envel.pois
  #
  # Os dados devem estar disponíveis pelo comando attach( ).
  #
  # Argumentos obrigatórios:
  # modelo: deve-se informar o objeto onde está o ajuste do modelo, caso não seja informado, a função procurará
  # 	  o ajuste no objeto fit.model;
  # 
  # Argumentos opcionais:
  # iden: caso deseje, informe o número de observações que irá querer destacar. O padrão é não destacar ninguém (iden=0).
  #	Qualquer valor que não seja um inteiro positivo (por ex., negativo ou decimal) fará com que a função pergunte
  #	o número de pontos após a execução;
  # nome: esse argumento só é utilizado caso seja destacado algum ponto no gráfico. Caso não seja informado nada, os pontos
  #	identificados serão os números da ordem em que estão no banco de dados (os índices). Caso se queira, pode-se
  #	informar um vetor de nomes ou de identificações alternativas. Obrigatoriamente esse vetor deve ter o mesmo
  #	comprimento do banco de dados;
  # sim: número de simulações para gerar a banda de confiança. Atkinson sugere um mínimo de 20 simulações.
  #      O padrão é de 100;
  # conf: nível de confiança do envelope. O padrão é de 90%;
  # res: permite-se a escolha dos resíduos. As opções dos resíduos são: "Q" quantil (ver Dunn e Smyth, 1996), "D" componente
  #      do desvio, "P" Pearson padronizado, "A" Anscombe, "W" Williams e "C" distância de Cook. A opção padrão é a "D";
  # mv: o valor T (True) fará com se obtenha a estimativa de máxima verossimilhança (EMV) para o parâmetro de
  #     dispersão. O valor F (False) indicará a escolha pela estimativa consistente pelo método dos momentos. O
  #     padrão é F. Note que como a EMV é razoavelmente mais demorada para ser obtida, a função demorará mais
  #     para rodar. Para obter a EMV a biblioteca MASS deve estar presente, no entanto não requer-se que seja
  #     carregada previamente;
  # quad: o padrão (quad=T, True) faz um gráfico quadrado, enquanto quad=F (False) faz um gráfico utilizando a área máxima
  #       disponível;
  # maxit: essa opção é utilizada nos ajustes de cada simulação e indica o máximo de iterações permitidas nos ajustes.
  #	 O padrão é maxit=20.
  #
  # Autor: Frederico Zanqueta Poleto <fred@poleto.com>, arquivo disponível em http://www.poleto.com
  #
  # Referências:
  # DUNN, K. P., and SMYTH, G. K. (1996). Randomized quantile residuals. J. Comput. Graph. Statist. 5, 1-10
  #    [http://www.statsci.org/smyth/pubs/residual.html e http://www.statsci.org/smyth/pubs/residual.ps]
  # MCCULLAGH, P. e NELDER, J. A. (1989). Generalized Linear Models. 2ª ed. Chapman and Hall, London.
  # PAULA, G. A. (2003). Modelos de Regressão com apoio computacional. IME-USP, São Paulo. [Não publicado,
  #    disponível em http://www.ime.usp.br/~giapaula/Book.pdf]
  #
  # Exemplos:
  # envel.gama(ajuste,sim=1000,conf=.95,mv=T,maxit=50)
  # envel.gama(ajuste,res="C")
  #
  
  if(class(modelo)[1] != "glm") {
    stop(paste("\nA classe do objeto deveria ser glm e nao ",class(modelo),"!!!\n"))
  }
  if(modelo$family[[1]] != "Gamma") {
    stop(paste("\nA familia do objeto deveria ser Gamma e nao ",modelo$family[[1]],"!!!\n"))
  }
  
  alfa<-(1-conf)/2
  X <- model.matrix(modelo)
  n <- nrow(X)
  p <- ncol(X)
  w <- modelo$weights
  W <- diag(w)
  H <- solve(t(X)%*%W%*%X)
  H <- sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
  h <- diag(H)
  
  #para evitar divisão por 0 ao studentizar os residuos, mas tentando manter o valor exagerado da alavanca
  h[round(h,15)==1]<-0.999999999999999
  
  m<-predict(modelo,type="response")
  y<-modelo$y
  if(mv==F) {
    fi <- (n-p)/sum((resid(modelo,type="response")/m)^2)
  } else {
    library("MASS")
    fi <- 1/gamma.dispersion(modelo) #Função gamma.shape retorna phi do texto, gamma.shape$alpha=1/gamma.dispersion
  }
  
  if(res=="Q") {
    tipo<-"Resíduo Quantil"
    r<-qnorm( pgamma(y,fi,fi/m) )
  } else {
    if(res=="D") {
      tipo<-"Resíduo Componente do Desvio"
      r<-resid(modelo,type="deviance")*sqrt(fi/(1-h))
    } else {
      if(res=="P") {
        tipo<-"Resíduo de Pearson Padronizado"
        r<-resid(modelo,type="pearson")*sqrt(fi/(1-h))
      } else {
        if(res=="A") {
          tipo<-"Resíduo de Anscombe"
          r<-3*sqrt(fi)*( y^(1/3) - m^(1/3) )/(m^(1/3))
        } else {
          if(res=="W") {
            tipo<-"Resíduo de Williams"
            r<-sign(y-m)*sqrt((1-h)*(( resid(modelo,type="deviance")*sqrt(fi/(1-h)) )^2)+(h*( resid(modelo,type="pearson")*sqrt(fi/(1-h)) )^2))
          } else {
            if(res=="C") {
              tipo<-"Distância de Cook"
              r<-(h/((1-h)*p))*((resid(modelo,type="pearson")/sqrt(1-h))^2)
            } else {
              stop(paste("\nVoce nao escolheu corretamente um dos residuos disponiveis!!!\n"))
            }
          }
        }
      }
    }
  }
  
  link<-modelo$family[[2]]
  
  e <- matrix(0,n,sim)
  e1 <- numeric(n)
  e2 <- numeric(n)
  
  if (is.null(version$language) == T) {
    #No S-Plus, a opção start é para entrar com o preditor linear
    pm<-predict(modelo)
  } else {
    #No R, a opção start é para entrar com os coeficientes
    pm<-coef(modelo)
  }
  mu<-m
  for(i in 1:sim) {
    resp <- rgamma(n,fi,fi/mu)
    if ( (is.null(version$language) == T && link == "Log: log(mu)") | (is.null(version$language) == F && link == "log") ) {
      fit <- glm(resp ~ X-1,family=Gamma(link=log),maxit=maxit,start=pm)
    } else {
      if ( (is.null(version$language) == T && link == "Inverse: 1/mu") | (is.null(version$language) == F && link == "inverse") ) {
        fit <- glm(resp ~ X-1,family=Gamma,maxit=maxit,start=pm)
      } else {
        if ( (is.null(version$language) == T && link == "sqrt: mu") | (is.null(version$language) == F && link == "sqrt") ) {
          fit <- glm(resp ~ X-1,family=Gamma(link="sqrt"),maxit=maxit,start=pm)
        } else {
          stop(paste("\nEsta funcao so aceita as ligacoes: canonica, log e sqrt!!!\nLigacao ",link," desconhecida!!!\n"))
        }
      }
    }
    w <- fit$weights
    W <- diag(w)
    H <- solve(t(X)%*%W%*%X)
    H <- sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
    h <- diag(H)
    h[round(h,15)==1]<-0.999999999999999
    m <- predict(fit,type="response")
    y <- fit$y
    if(mv==F) {
      phi <- (n-p)/sum((resid(fit,type="response")/m)^2)
    } else {
      phi <- 1/gamma.dispersion(fit)
    }
    e[,i] <- 
      sort( if(res=="Q") {
        qnorm( pgamma(y/(m/phi),phi) )
      } else {
        if(res=="D") {
          resid(fit,type="deviance")*sqrt(phi/(1-h))
        } else {
          if(res=="P") {
            resid(fit,type="pearson")*sqrt(phi/(1-h))
          } else {
            if(res=="A") {
              3*sqrt(phi)*( y^(1/3) - m^(1/3) )/(m^(1/3))
            } else {
              if(res=="W") {
                sign(y-m)*sqrt((1-h)*(( resid(fit,type="deviance")*sqrt(phi/(1-h)) )^2)+(h*( resid(fit,type="pearson")*sqrt(phi/(1-h)) )^2))
              } else {
                if(res=="C") {
                  (h/((1-h)*p))*((resid(fit,type="pearson")/sqrt(1-h))^2)
                } else {
                  stop(paste("\nVoce nao escolheu corretamente um dos residuos disponiveis!!!\n"))
                }
              }
            }
          }
        }
      })
  }
  
  for(i in 1:n) {
    eo <- sort(e[i,])
    e1[i] <- quantile(eo,alfa)
    e2[i] <- quantile(eo,1-alfa)
  }
  
  med <- apply(e,1,median)
  
  if(quad==T) {
    par(pty="s")
  }
  if(res=="C") {
    #Segundo McCullagh e Nelder (1989, pág.407) e Paula (2003, pág.57) deve-se usar qnorm((n+1:n+.5)/(2*n+1.125))
    #Segundo Neter et alli (1996, pág.597) deve-se usar qnorm((n+1:n-.125)/(2*n+0.5))
    qq<-qnorm((n+1:n+.5)/(2*n+1.125))
    plot(qq,sort(r),xlab="Quantil Meio-Normal",ylab=tipo, ylim=range(r,e1,e2), pch=20, col = "blue",cex.lab=2, cex.axis=2)
  } else {
    qq<-qnorm((1:n-.375)/(n+.25))
    plot(qq,sort(r),xlab="Theorical quantiles",ylab="Samples Quantiles", ylim=range(r,e1,e2), pch=20, col = "blue",cex.lab=2, cex.axis=2)
  }
  lines(qq,e1,lty=1)
  lines(qq,e2,lty=1)
  lines(qq,med,lty=2)
  nome<-nome[order(r)]
  r<-sort(r)
  while ( (!is.numeric(iden)) || (round(iden,0) != iden) || (iden < 0) ) {
    cat("Digite o num.de pontos a ser identificado (0=nenhum) e <enter> para continuar\n")
    out <- readline()
    iden<-as.numeric(out)
  }
  if(iden>0) {identify(qq,r,n=iden,labels=nome)}
  if(quad==T) {
    par(pty="m")
  }
  cat("Banda de ",conf*100,"% de confianca, obtida por ",sim," simulacoes.\n")
}
envel.norm <- function(modelo=fit.model,iden=0,nome=seq(along = model.matrix(modelo)[,1]),sim=100,conf=.90,res=T,quad=T) {
  
  #
  # Descrição e detalhes:
  # A saída será o gráfico de probabilidade normal com envelopes simulados para um ajuste da distribuição normal.
  #
  # A opção res=F faz o gráfico de probabilidade meio-normal com envelopes simulados utilizando a distância de Cook,
  # possibilitando a detecção de pontos simultaneamente aberrantes e/ou influentes.
  #
  # Atenção: a função não funcionará corretamente se o ajuste possuir offsets! Neste caso é preciso adaptá-la como foi
  # feito na função envel.pois
  #
  # Os dados devem estar disponíveis pelo comando attach( ).
  #
  # Argumentos obrigatórios:
  # modelo: deve-se informar o objeto onde está o ajuste do modelo normal linear, caso não seja informado, a função
  # 	  procurará o ajuste no objeto fit.model;
  # 
  # Argumentos opcionais:
  # iden: caso deseje, informe o número de observações que irá querer destacar. O padrão é não destacar ninguém (iden=0).
  #	Qualquer valor que não seja um inteiro positivo (por ex., negativo ou decimal) fará com que a função pergunte
  #	o número de pontos após a execução;
  # nome: esse argumento só é utilizado caso seja destacado algum ponto no gráfico. Caso não seja informado nada, os pontos
  #	identificados serão os números da ordem em que estão no banco de dados (os índices). Caso se queira, pode-se
  #	informar um vetor de nomes ou de identificações alternativas. Obrigatoriamente esse vetor deve ter o mesmo
  #	comprimento do banco de dados;
  # sim: número de simulações para gerar a banda de confiança. Atkinson sugere um mínimo de 20 simulações.
  #      O padrão é de 100;
  # conf: nível de confiança do envelope. O padrão é de 90%;
  # res: permite-se a escolha se o gráfico será feito com os resíduos (res=T, True, padrão) ou com a distância de Cook
  #      (res=F, False);
  # quad: o padrão (quad=T, True) faz um gráfico quadrado, enquanto quad=F (False) faz um gráfico utilizando a área máxima
  #       disponível.
  #
  # Autor: Frederico Zanqueta Poleto <fred@poleto.com>, arquivo disponível em http://www.poleto.com
  #
  # Referências:
  # MCCULLAGH, P. e NELDER, J. A. (1989). Generalized Linear Models. 2ª ed. Chapman and Hall, London.
  # NETER, J., KUTNER, M. H., NACHTSHEIM, C. J. and WASSERMAN, W. (1996). Applied Linear Statistical Models. 4ª ed.
  #    Mc Graw Hill, Boston.
  # PAULA, G. A. (2003). Modelos de Regressão com apoio computacional. IME-USP, São Paulo. [Não publicado,
  #    disponível em http://www.ime.usp.br/~giapaula/Book.pdf]
  #
  # Exemplos:
  # envel.norm(ajuste,sim=10000,conf=.95)
  # envel.norm(ajuste,res=F)
  #
  
  if( class(modelo)[1]=="lm" || (class(modelo)[1]=="glm" && (modelo$family[[1]]=="Gaussian" | modelo$family[[1]]=="gaussian")) ) {
    
  } else {
    stop(paste("\nA classe do objeto deveria ser lm ou glm (com distribuicao gaussian) !!!"))
  }
  
  alfa<-(1-conf)/2
  X <- model.matrix(modelo)
  y<-predict(modelo)+resid(modelo)
  n <- nrow(X)
  p <- ncol(X)
  H <- X%*%solve(t(X)%*%X)%*%t(X)
  h <- diag(H)
  m <- fitted(modelo)
  
  #para evitar divisão por 0 ao studentizar os residuos, mas tentando manter o valor exagerado da alavanca
  h[round(h,15)==1]<-0.999999999999999
  
  si <- lm.influence(modelo)$sigma
  r <- resid(modelo)
  tsi <- r/(si*sqrt(1-h))
  sigma<- summary(modelo)$sigma
  ti <- r/(sigma*sqrt(1-h))
  di <- (1/p)*(h/(1-h))*(ti^2)
  
  e <- matrix(0,n,sim)
  e1 <- numeric(n)
  e2 <- numeric(n)
  
  for(i in 1:sim) {
    resp <- rnorm(n,m,sigma)
    fit <- lm(resp~X-1)
    ti<-resid(fit)/(summary(fit)$sigma*sqrt(1-h))
    if(res==F) {
      e[,i] <- (1/p)*(h/(1-h))*(ti^2)
    } else {
      e[,i] <- ti*sqrt( (n-p-1)/(n-p-(ti^2)) )
    }	
    e[,i] <- sort(e[,i])
  }
  
  for(i in 1:n) {
    eo <- sort(e[i,])
    e1[i] <- quantile(eo,alfa)
    e2[i] <- quantile(eo,1-alfa)
  }
  
  med <- apply(e,1,median)
  
  if(quad==T) {
    par(pty="s")
  }
  if(res==F) {
    #Segundo McCullagh e Nelder (1989, pág.407) e Paula (2003, pág.57) deve-se usar qnorm((n+1:n+.5)/(2*n+1.125))
    #Segundo Neter et alli (1996, pág.597) deve-se usar qnorm((n+1:n-.125)/(2*n+0.5))
    qq<-qnorm((n+1:n+.5)/(2*n+1.125))
    plot(qq,sort(di),xlab="Quantil Meio-Normal",ylab="Distância de Cook", ylim=range(di,e1,e2), pch=16, col = "blue",cex.lab=2, cex.axis=2)
    nome<-nome[order(di)]
    r<-sort(di)
  } else {
    qq<-qnorm((1:n-.375)/(n+.25))
    plot(qq,sort(tsi),xlab="Theorical quantiles",ylab="Samples Quantiles", ylim=range(tsi,e1,e2), pch=20, col = "blue",cex.lab=2, cex.axis=2)
    nome<-nome[order(tsi)]
    r<-sort(tsi)
  }
  lines(qq,e1,lty=1)
  lines(qq,e2,lty=1)
  lines(qq,med,lty=2)
  while ( (!is.numeric(iden)) || (round(iden,0) != iden) || (iden < 0) ) {
    cat("Digite o num.de pontos a ser identificado (0=nenhum) e <enter> para continuar\n")
    out <- readline()
    iden<-as.numeric(out)
  }
  if(iden>0) {identify(qq,r,n=iden,labels=nome)}
  if(quad==T) {
    par(pty="m")
  }
  cat("\nBanda de ",conf*100,"% de confianca, obtida por ",sim," simulacoes.\n")
}

