## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message=FALSE, warning=FALSE)
library(knitr)

## ----echo=T, results='hide'---------------------------------------------------
library(StatComp21003)
library(MASS)
library(glmnet)
library(ncpen)
library(SIS)
library(VariableScreening)
library(bestglm)
library(hdi)
ss <- seq(15, 25, by=5)
finalresult_ldp <- matrix(nrow=length(ss), ncol=6)
finalresult_ap <- matrix(nrow=length(ss), ncol=6)
n <- 200;  p <- 300;  rho <- 0.8
counter <- 1
for(k in 1:length(ss)){
  
  set.seed(99)
  s <- ss[k]
  S_0 <- c(1:s)#sample(seq(1, 500, 1), 15)
  c <- 2
  truebeta <- rep(0, p)
  truebeta[S_0] <- runif(s, 0, c)
  
  m <- 3
  cpmax_ldp <- rep(-1,m); cpall_ldp <- rep(-1, m); len_ldp <- rep(-1, m)
  cpmax_ap <- rep(-1,m); cpall_ap <- rep(-1, m); len_ap <- rep(-1, m)
  for(i in 1:m){
    
    dat <- simdata(n, p, rho, truebeta)
    X <- dat$X
    y <- dat$y
    objsis <- suppressMessages(SIS(X, y, family = 'gaussian', tune = 'bic'))
    S <- objsis$sis.ix0
    
    objLDP <- LDP(X, y)
    CI <- objLDP$CI;  len_CI <- objLDP$len_CI
    
    cpmax_ldp[i] <- length(intersect(which(truebeta[S_0]>=CI[1,S_0]), 
                                     which(truebeta[S_0]<=CI[2,S_0])))/s
    cpall_ldp[i] <- length(intersect(which(truebeta>=CI[1,]), 
                                     which(truebeta<=CI[2,])))/p
    len_ldp[i] <- mean(len_CI)
    
    objAP <- AP(X, y, S)
    CI <- objAP$CI; len_CI <- objAP$len_CI
    
    cpmax_ap[i] <- length(intersect(which(truebeta[S_0]>=CI[1,S_0]), 
                                    which(truebeta[S_0]<=CI[2,S_0])))/s
    cpall_ap[i] <- length(intersect(which(truebeta>=CI[1,]), 
                                    which(truebeta<=CI[2,])))/p
    len_ap[i] <- mean(len_CI)
  }
  
  
  finalresult_ldp[counter, ] <- c(p, s, c, mean(cpall_ldp), mean(cpmax_ldp), mean(len_ldp))
  finalresult_ap[counter, ] <- c(p, s, c, mean(cpall_ap), mean(cpmax_ap), mean(len_ap))
  counter <- counter+1
}

## -----------------------------------------------------------------------------
cpmax <- cbind(finalresult_ap[,c(2,5)], finalresult_ldp[,5])
dimnames(cpmax)[[2]] <- c("s", "CP_APE", "CP_LDPE")
kable(cpmax, caption = "Averaged coverage probability of the nonzero coefficients")

## -----------------------------------------------------------------------------
len <- cbind(finalresult_ap[,c(2,6)], finalresult_ldp[,6])
dimnames(len)[[2]] <- c("s", "length_APE", "length_LDPE")
kable(len, caption = "Averaged length of confidence intervals of all coefficients")

