rm(list = ls())
library(DataCombine)
library(xtable)
library(Benchmarking)
library(gamlss)
library(betareg)
library(tidyverse)

### Semente 
set.seed(2005, kind = "Marsaglia-Multicarry")

##### Directory
## setwd(choose.dir())
setwd("C:\\Users\\Bruno\\Dropbox\\Bruno Artigo\\Bruno\\codes")
getwd()

### data 
load("dados.RData")
head(dados)
dim(dados)

input <- dados[, c('lxterra', 'lxtrab', 'lxresto') ]
output <- dados[,'ly']

inputR  <- apply(input, 2, rank) / length(output)
outputR <- rank(output) / length(output)

e_vrs_out0 <- dea(X = inputR, 
                  Y = outputR, 
                  RTS="vrs", 
                  ORIENTATION= "out", 
                  DUAL=TRUE)

Efi_out0 <- ( 1 / e_vrs_out0$eff)

Efi_out0 <- Efi_out0 %>% as.numeric()

#### Financit

End = financi ~ 
  lxterra + lxtrab + lxresto +
  social + demo + ambi + ica460s + ginitotal

a <- which(dados$financi == 1)

dados[a, 'financi'] <- 0.999999999999

resul9 <- betareg(formula = End, 
                    data = dados,
                    link.phi = "log",
                    link = "logit")

resumo <- summary(resul9)[[1]] 

lapply(resumo,round,4)

confint(resul9) %>% round(4)

coef_beta <- coef(resul9)

fin_pred <- predict(resul9, type="response")
phi_pred <- predict(resul9, type="precision")

dadosO <- cbind(dados, Efi_out0, fin_pred)

dadosO$regiao <- relevel(dadosO$regiao, ref = '5')

formula    = Efi_out0 ~ fin_pred + social + demo + ambi + ica460s + ginitotal + regiao 

dadosA <- model.frame(formula, dadosO)

New_data <- dadosA

i <- 0

{
  boot_betaflacionada <- function(data, indices){
    
    data_new <- data[indices, ]
    
    data_new$regiao <- relevel(data_new$regiao, ref = '5')
    
    i <<- i + 1
    print(paste0('count ', i))
    
    mu.formula    =  Efi_out0 ~ fin_pred  + social + demo + ambi + ica460s + ginitotal
    alpha.formula =          ~ 1
    phi.formula   =          ~ regiao 
    
    fit = gamlss(formula = mu.formula, 
                   nu.formula= alpha.formula, 
                   sigma.formula= phi.formula, 
                   family= BEOI(mu.link = "logit", 
                                sigma.link = "log", 
                                nu.link = "probit"), 
                   data = data_new, trace = FALSE,
                   method=RS())
    
    #summary(fit)
    coef_mu    <- coef(fit, what = 'mu')
    coef_alpha <- coef(fit, what = 'nu')
    coef_phi   <- coef(fit, what = 'sigma')
    
    coef_fit <- c(coef_mu, coef_phi, coef_alpha)
    
    return(coef_fit)  
    
  }
}

library(boot)
duncan.boot <- boot(New_data, boot_betaflacionada, R = 2000)
























