
rm(list = ls())
library(tidyverse)
library(DataCombine)
library(xtable)
library(Benchmarking)
library(gamlss)
library(betareg)

### Directory
## setwd(choose.dir())
setwd("C:\\Users\\Bruno\\Dropbox\\Bruno Artigo\\Bruno\\codes")
getwd()

### data 

load("dados.RData")
head(dados)
dim(dados)

####  efficiency

input <- dados[,c('lxterra', 'lxtrab', 'lxresto')]
output <- dados[,'ly']

inputR  <- apply(input, 2, rank) / length(output)
outputR <- rank(output) / length(output)

e_vrs_out0 <- dea(X = inputR, 
                  Y = outputR, 
                  RTS="vrs", 
                  ORIENTATION= "out", 
                  DUAL=TRUE)

Efi_out0 <- 1/e_vrs_out0$eff

Efi_out0[which(Efi_out0 == 1)]

dados <- cbind(Efi_out0, dados)


#### Without Endogeneity

dados$regiao <- relevel(dados$regiao, ref = '5')

mu.formula0    = Efi_out0 ~ financi + social + demo + ambi +  ica460s + ginitotal
alpha.formula0 =          ~ 1
phi.formula0   =          ~ regiao 

fit00 = gamlss(formula = mu.formula0, 
               nu.formula= alpha.formula0, 
               sigma.formula= phi.formula0, 
               family= BEOI(mu.link = "logit", 
                            sigma.link = "log", 
                            nu.link = "probit"), 
               data = dados, trace = FALSE,
               method=RS())

summary(fit00)


#### Financi

End = financi ~ 
  lxterra + lxtrab + lxresto +
  social + demo + ambi + ica460s + ginitotal

## it will be used in Murphy and Topel
mmm <- model.frame(End, dados)
X_end <- model.matrix(End, mmm)
kd <- dim(X_end)[2]
##

dados0 <- dados
dados0[ which(dados$financi == 1), 'financi' ] <- 0.999999999999

dados0$regiao <- relevel(dados0$regiao, ref = '5')

gy_logit <- betareg(formula = End, 
                  data = dados0,
                  link.phi = "log",
                  link = "logit")

summary(gy_logit) 

fin_pred <- predict(gy_logit, type="response")
phi_pred <- predict(gy_logit, type="precision")


## Beta regression inflated one with the endogenous credit variable.
dadosA <- cbind(dados, fin_pred)

dadosA$regiao <- relevel(dadosA$regiao, ref = '5')

mu.formula    = Efi_out0 ~ fin_pred  + social + demo + ambi + ica460s + ginitotal
alpha.formula =          ~ 1
phi.formula   =          ~ regiao 

fit=gamlss(formula = mu.formula, 
           nu.formula= alpha.formula, 
           sigma.formula= phi.formula, 
           family=BEOI(mu.link = "logit", 
                       sigma.link = "log", 
                       nu.link = "probit"), 
           data = dadosA, trace = FALSE)

summary(fit)

coef_mu    <- coef(fit, what = 'mu')
coef_alpha <- coef(fit, what = 'nu')
coef_phi   <- coef(fit, what = 'sigma')

coef0_mu    <- coef(fit00, what = 'mu')
coef0_alpha <- coef(fit00, what = 'nu')
coef0_phi   <- coef(fit00, what = 'sigma')

coef_fit <- c(coef_mu, coef_phi, coef_alpha)
coef0_fit <- c(coef0_mu, coef0_phi, coef0_alpha)

ses <- vcov(fit)
ses0 <- vcov(fit00)

mu    <- fitted(fit, what = 'mu')
alpha <- fitted(fit, what = 'nu')
phi   <- fitted(fit, what = 'sigma')

#### Take care from now on

mu.formula    = Efi_out0 ~ fin_pred + social + demo + ambi + ica460s + ginitotal
alpha.formula =          ~ 1
phi.formula   =          ~ regiao 

m <- model.frame(mu.formula, dadosA)
X_mu <- model.matrix(mu.formula, m)
y <- model.response(m, "numeric")

position <- which(y %in% 1)
yy <- replace(y, y != 1, 0)

mm <- model.frame(phi.formula, dadosA)
X_phi <- model.matrix(phi.formula, mm)

mmm     <- model.frame(alpha.formula, dadosA)
X_alpha <- model.matrix(alpha.formula, mmm)

mu_X <- X_mu[-position,]
phi_X <- X_phi[-position,]
y_mu <- y[-position]

y_fin <-  dados$financi

b_mu    <- coef(fit, what = 'mu')
b_alpha <- coef(fit, what = 'nu')
b_phi   <- coef(fit, what = 'sigma')

k0 <- dim(X_mu)[2]
k1 <- k0 + dim(X_alpha)[2]
kk <- k1 + dim(X_phi)[2]

parms <- z <-  c(b_mu, b_alpha, b_phi)

{
  lG22 <- function(z) {
    betas <- z[1:k0]
    b_alfa_eta  <- z[(k0 + 1):k1]
    b_phi_eta <- z[(k1 + 1):kk]
    
    ###### Parte discreta
    alfa_est  <- as.vector(X_alpha   %*%   b_alfa_eta)
    alpha <- pnorm(alfa_est)
    
    lalpha <- (yy/alpha -  (1 - yy)/ (1 - alpha))*dnorm(alfa_est)
    Balpha <- X_alpha*as.vector(lalpha)
    
    ###### Parte contínua
    
    y_eta     <- as.vector(mu_X    %*%   betas)
    phi_eta   <- as.vector(phi_X   %*%   b_phi_eta)
    
    phi <- pmax(exp(phi_eta), 1e-150)
    mu <- pmax(( 1 / (1 + exp(-y_eta))), 1e-150)
    q_mu <- pmax((1 - mu), 1e-150)
    
    lmu <- phi*(digamma(q_mu*phi) - digamma(mu*phi) + log(y_mu/(1 - y_mu)))*(mu*q_mu)
    lphi <- (digamma(phi) - q_mu*digamma(q_mu*phi) - mu*digamma(mu*phi) + mu*log(y_mu/(1 - y_mu)) + log(1 - y_mu))*phi
    
    Bmu <- mu_X*as.vector(lmu)
    Bphi <- phi_X*as.vector(lphi)
    
    new_Bmu<- matrix(0, nrow= dim(X_alpha)[1], ncol= dim(mu_X)[2]) 
    new_Bmu[-position,] <- Bmu
    
    new_Bphi <- matrix(0, nrow= dim(X_alpha)[1], ncol= dim(phi_X)[2]) 
    new_Bphi[-position,] <- Bphi
    
    TE <- cbind(new_Bmu, new_Bphi, Balpha)
    colnames(TE) <- c(colnames(Bmu), colnames(Bphi), colnames(Balpha))
    
    return(TE)
  }
  
  lG21 <- function(z, par_fin, Xreg) {
    b_mu    <- z[1:k0]
    b_alfa  <- z[(k0 + 1):k1]
    b_phi   <- z[(k1 + 1):kk]
    
    prob_fin <- plogis(Xreg %*% par_fin)

    ###### Parte discreta
    alfa_eta  <- as.vector(X_alpha   %*%   b_alfa)
    mu_eta     <- as.vector(mu_X    %*%   b_mu)
    phi_eta   <- as.vector(phi_X   %*%   b_phi)
    
    phi   <- pmax(exp(phi_eta), 1e-150)
    mu    <- pmax(( 1 / (1 + exp(-mu_eta))), 1e-150)
    alpha <- pnorm(alfa_eta)
    q_mu  <- pmax((1 - mu), 1e-150)
    
    ############ Verossimilhança
    mu_a <- digamma(mu*phi) - digamma(q_mu*phi)
    y_a  <- log(y_mu/(1 - y_mu))
    
    lalpha <- (yy/alpha -  (1 - yy)/ (1 - alpha))*dnorm(alfa_eta)
    lmu <- phi*(y_a - mu_a)*(mu*q_mu)
    
    Lmu <- matrix(0, nrow= dim(X_alpha)[1], ncol= 1) 
    Lmu[-position,] <- lmu
    
    # g_fin <- Xreg * as.vector(prob_fin*(1 - prob_fin)*(lalpha*b_alfa[2] +  Lmu*b_mu[2])) ### Regressão logistica
    g_fin <- Xreg * as.vector(prob_fin*(1 - prob_fin)*(Lmu*b_mu[2])) ### Regressão logistica
    
    gra <- cbind(g_fin)
    colnames(gra) <- c(colnames(Xreg))
    
    return(gra)
  }
  
  G1 <- function(par_fin, phi, y, Xreg) {

    betas <- par_fin
    eta_mu <- as.vector(Xreg %*% betas)
    mu <- exp(eta_mu) / ( 1 + exp(eta_mu) )
    mu.u <- digamma(mu * phi) - digamma((1 - mu) * phi) 
    ystar <- log(y/ (1 - y))
    
    dev_link <- 1/(mu*(1 - mu))
    TO <- 1/dev_link
    
    L_mu <- cbind( phi * Xreg * TO * (ystar - mu.u) )
    L_phi <- cbind( 'phi' = (mu * (ystar - mu.u) + log(1-y) - digamma((1-mu)*phi) + digamma(phi)) )
    
    rval <- cbind(L_mu, L_phi )
    return(rval)
  }
}


G22 <- lG22(parms)
colSums(G22)

par_fin <- coef(gy_logit)[1:9]
par_phi <- coef(gy_logit)[10]

G21 <- lG21(z = parms, par_fin = par_fin, Xreg = X_end )
colSums(G21)

fin <- dados0$financi

G_fin <- G1(par_fin = par_fin, phi = phi_pred, y = fin, Xreg = X_end)
colSums(G_fin)

V0 <- vcov(gy_logit)

V1 <- V0[-10,-10]
V2 <- vcov(fit)
G11 <- G_fin[, -10]

#### Corrigido

C <- t(G22) %*% G21
R <- t(G22) %*% G11

V33 <- ((C %*% V1 %*% t(C)) - (R %*% V1 %*% t(C)) - (C %*% V1 %*% t(R))) 

V2_corrigida <- V2 + V2 %*% V33 %*% V2

## Teste Hausman

d <-  coef_fit[1:7] - coef0_fit[1:7]

H <- solve(V2_corrigida[1:7, 1:7] - ses0[1:7, 1:7])

Teste <- t(d) %*% H %*% d

1 - pchisq(as.numeric(Teste), 7)

#############
names0 <- colnames(V2_corrigida)
EP2 <- sqrt(diag(V2))

EP_corrigido <- sqrt(diag(V2_corrigida))
names(EP_corrigido) <- names0

DS <- data.frame(names0, EP2, EP_corrigido)
DS

b_mu    <- coef(fit, what = 'mu')
b_alpha <- coef(fit, what = 'nu')
b_phi   <- coef(fit, what = 'sigma')

mu_desvio_padrao <- EP_corrigido[1:7]
alpha_desvio_padrao <- EP_corrigido[13]
phi_desvio_padrao <- EP_corrigido[8:12]

muL_inf <- b_mu - qnorm(1-0.05/2)*mu_desvio_padrao
muL_sup <- b_mu + qnorm(1-0.05/2)*mu_desvio_padrao
z_mu <-  b_mu/mu_desvio_padrao
p_valor_mu <- 2*pnorm(-abs(z_mu))
result_mu <- round(data.frame("Estimado" = b_mu, "DesvioPadrao" = mu_desvio_padrao, "L_inferior" = muL_inf, "L_superior" = muL_sup, "P-valor" = p_valor_mu), 4)

phiL_inf <- b_phi - qnorm(1-0.05/2)*phi_desvio_padrao
phiL_sup <- b_phi + qnorm(1-0.05/2)*phi_desvio_padrao
z_phi <-  b_phi/phi_desvio_padrao
p_valor_phi <- 2*pnorm(-abs(z_phi))
result_phi <- round(data.frame("Estimado" = b_phi, "DesvioPadrao" = phi_desvio_padrao, "L_inferior" = phiL_inf, "L_superior" = phiL_sup, "P-valor" = p_valor_phi), 4)

alphaL_inf <- b_alpha - qnorm(1-0.05/2)*alpha_desvio_padrao
alphaL_sup <- b_alpha + qnorm(1-0.05/2)*alpha_desvio_padrao
z_alpha <-  b_alpha/alpha_desvio_padrao
p_valor_alpha <- 2*pnorm(-abs(z_alpha))
result_alpha <- round(data.frame("Estimado" = b_alpha, "DesvioPadrao" = alpha_desvio_padrao, "L_inferior" = alphaL_inf, "L_superior" = alphaL_sup, "P-valor" = p_valor_alpha), 4)

saida0 <- list("Regressão mu" = result_mu,"Regressão phi" = result_phi,  "Regressão alpha" = result_alpha)

saida0




























