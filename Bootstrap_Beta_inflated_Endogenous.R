rm(list = ls())

library(Benchmarking)
library(gamlss)

# Directory
setwd(choose.dir())
getwd()

# data 
load("dados.RData")
head(dados)
dim(dados)

# DEA model
input <- dados[, c('lxterra', 'lxtrab', 'lxresto') ]
output <- dados[,'ly']

inputR  <- apply(input, 2, rank) / length(output)
outputR <- rank(output) / length(output)

e_vrs_out0 <- dea(X = inputR, 
                  Y = outputR, 
                  RTS="vrs", 
                  ORIENTATION= "out", 
                  DUAL=TRUE)

Efi_out <- ( 1 / e_vrs_out0$eff) %>% as.numeric()

dados <- cbind(Efi_out, dados)

# Endogenous Variable
End = financi ~ 
  lxterra + lxtrab + lxresto +
  social + demo + ambi + ica460s + ginitotal

reg <- glm(End, data = dados, family = quasibinomial(link = "logit"))

summary(reg)$coef %>% round(4)
######################

formula    = Efi_out ~ lxterra + lxtrab + lxresto + financi + social + demo + ambi + ica460s + ginitotal + regiao 

New_data <- model.frame(formula, dados)

New_data$regiao <- relevel(New_data$regiao, ref = '5')

{
  boot_betaflacionada <- function(data, indices){
    
    data_new <- data[indices, ]
    
    data_new$regiao <- relevel(data_new$regiao, ref = '5')
    
    i <<- i + 1
    print(paste0('count ', i))
    
    End = financi ~ 
      lxterra + lxtrab + lxresto +
      social + demo + ambi + ica460s + ginitotal
    
    reg <- glm(End, data = data_new, family = quasibinomial(link = "logit"))
    fin_pred <- predict(reg, type="response")
    
    dadosO <- cbind(data_new, fin_pred)
    
    mu.formula    =  Efi_out ~ fin_pred  + social + demo + ambi + ica460s + ginitotal
    alpha.formula =          ~ 1
    phi.formula   =          ~ regiao 
    
    fit = gamlss(formula = mu.formula, 
                 nu.formula= alpha.formula, 
                 sigma.formula= phi.formula, 
                 family= BEOI(mu.link = "logit", 
                              sigma.link = "log", 
                              nu.link = "probit"), 
                 data = dadosO, trace = FALSE,
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
### Semente 
set.seed(2005)
i <- 0
duncan.boot <- boot(New_data, boot_betaflacionada, R = 2000)
duncan.boot
duncan.boot[[1]] %>% cbind() %>% round(5)



















