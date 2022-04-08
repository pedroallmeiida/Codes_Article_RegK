source("func_MEIJER.R")
source("Functions_regressaoK.R")

library(dplyr)
library(parallel)
library(pbmcapply)
library(magrittr)
library(ggplot2)
library(reshape2)

### functions
glm.K <- function( Y,X,L, method ){
  
  chute = Chute.inicial( Y, X, L=L )
  est = cml.K( Y, X, chute, L=L, method =  method  )
  est_b0 = cml.Kb0( Y, X, chute[1:2], L=L, method =  method )
  est_b1 = cml.Kb1( Y, X, chute[c(1,3)], L=L, method =  method  )
  
  alpha <- est$estimate[1]; beta0 <- est$estimate[2]; beta1 <- est$estimate[3]
  mu.hat <- est$mu.hat
  lk <- loglink_KBR( c(Y), c(X), as.numeric(est$estimate),L=L  )
  lkb0 <- loglink_KBRb0( c(Y), c(X), c(as.numeric(est_b0$estimate)[1:2]),L=L  )
  lkb1 <- loglink_KBRb1( c(Y), c(X), c(as.numeric(est_b1$estimate)[c(1,2)]),L=L  )
  
  TRV = c( 2*(lk -  lkb0), 2*(lk -  lkb1) )  
  p.value = c(1 - pchisq( TRV[1], df=1 ), 1 - pchisq( TRV[2], df=1 )   ); p.value
  
  par = c(alpha, beta0, beta1)
  AIC_Y = est$AIC;  BIC_Y = est$BIC;
  
  F.dist = DistK( alpha, L=L, mu.hat, Y  )
  res.quantile = qnorm( F.dist) 
  r2.k <- R2.k(c(Y), c(X), as.vector(est$estimate),est_b0$estimate ,  L=L)
  mae = mean( abs( Y - mu.hat ) )
  
  
  return( list(ajuste = par, estimate = est$estimate, residual = res.quantile , AIC = AIC_Y, BIC = BIC_Y, fitted.values = mu.hat, r.squared = r2.k, mae = mae) )
  
  
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
  
  mu.hat = c(ajuste$fitted.values)
  alpha = ajuste$estimate[1]
  
  F.dist = DistK( alpha,L,  mu.hat, Y  )
  res = qnorm( F.dist, lower.tail =  TRUE ) 
  
  res.sort = matrix(NA, n, length(Y))
  
  for(i in 1:n){
    y.sim = geradorK( length(Y), X, as.numeric(ajuste$estimate) , L)
    chute = Chute.inicial( y.sim , X, L=4 )
    ajuste_sim = cml.K( y.sim, X, chute, L, method  ) 
    mu.hat = c(ajuste_sim$mu.hat); alpha = ajuste_sim$estimate[1]
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
calculate.r2 <-function( Y, X, method = method, L  ){
  require( rcompanion  )
  chute = Chute.inicial( Y, X, L=L )
  est = cml.K( Y, X, chute, L=L, method = method  );est_b0 = cml.Kb0( Y, X, chute[1:2], L=L, method =  method )
  r2.k <- R2.k(c(Y), c(X), as.vector(est$estimate),est_b0$estimate ,  L=L)
  
  ajuste_norm <- lm( Y ~ X )
  ajuste_gamma <- glm( Y ~ X, family = Gamma( link = "log" ) )
  r2 = nagelkerke(ajuste_norm); r2_ajuste_norm = r2$Pseudo.R.squared.for.model.vs.null[2]
  r2 = nagelkerke(ajuste_gamma); r2_ajuste_gamma = r2$Pseudo.R.squared.for.model.vs.null[2]
  
  return( R2 = c( r2_K = r2.k,  r2_norm = r2_ajuste_norm , r2_gamma = r2_ajuste_gamma  ) )
}

### load SAR image
load("SanFranciscoImage_150_150.Rdata")


M.HH = t(apply(matrix( Re(ImgMP_Sim[,,1]),ncol=dim(ImgMP_Sim)[2] ), 2, rev))
M.HV = t(apply(matrix( Re(ImgMP_Sim[,,2]),ncol=dim(ImgMP_Sim)[2] ), 2, rev))
M.VV = t(apply(matrix( Re(ImgMP_Sim[,,3]),ncol=dim(ImgMP_Sim)[2] ), 2, rev))

## channels
HH = M.HH[-c(1:3, (dim(M.HH)[1]-2):dim(M.HH)[1]),-c(1:3, (dim(M.HH)[2]-2):dim(M.HH)[2])]
HV = M.HV[-c(1:3, (dim(M.HV)[1]-2):dim(M.HV)[1]),-c(1:3, (dim(M.HV)[2]-2):dim(M.HV)[2])]
VV = M.VV[-c(1:3, (dim(M.VV)[1]-2):dim(M.VV)[1]),-c(1:3, (dim(M.VV)[2]-2):dim(M.VV)[2])]

##### urban region


a=136:142 ; b = 42:48
urban_VV = c( t(M.VV[a,b]) )
urban_HH = c( t(M.HH[a,b]) )
urban_HV = c( t(M.HV[a,b]) )

### Regression K ----


ajuste_urban_VV_K = glm.K(urban_VV, urban_HV, method = "NM", L=4 )

mu.hat_urban_VV_K = ajuste_urban_VV_K$fitted.values
dados = data.frame( urban_VV, urban_HV ,mu.hat= mu.hat_urban_VV_K)
ggplot(data = dados, aes(x =  urban_HV, y = urban_VV)) + 
  labs(x="HV", y="VV") +
  geom_point(color='black') +
  geom_line(color='red', aes(x= urban_HV, y=mu.hat))+
  theme_minimal()+
  theme(legend.position="top",
        legend.key.size = unit(1, "lines"))+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size = 15, face = "bold"),legend.title=element_text(size=15), legend.text=element_text(size=15))


