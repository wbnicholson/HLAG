#### R-script Simulation Study : Measuring Forecast Accuracy for Factor Models, AR and VAR(1) benchmark ####
rm(list=ls())

#### Check packages installed ####
checkpackage<-function(U){
  if((U %in% rownames(installed.packages()))==F){
    install.packages(U)
    library(U, character.only = TRUE)
  }else{
    library(U, character.only = TRUE)
  }
}
packagelist<-list("lattice", "Rcpp", "MASS","methods", "zoo", "stats","utils","grDevices",
                  "graphics","RcppArmadillo", "RcppEigen", "R.matlab", "vars", "bigtime")
lapply(packagelist,checkpackage)

#### Source Functions ####
source("factorfunctions.R") # factor functions

  oldw <- getOption("warn")
  options(warn = -1)

##############################################################
#### Simulation Scenario 1 : COmponentwise HLag Structure ####
##############################################################

#### Setting ####
p=5 # p the maximal lag order
k=45 # k is the number of series
n=100 # n is the time series length. We ran the simulations for both n=100 and n=200
Nsim=500
#### Forecast Performance ####
Sim1MSFE <- matrix(NA, ncol=3, nrow=Nsim) # Note : VAR1 is already contained in the other file 
colnames(Sim1MSFE) <- c("DFM", "FAVAR", "AR")
datasim1 <- read.table("sim1dat.txt")
library(bigtime) # for AR model 

for(r in 1:Nsim){
  # Data
  Y <- datasim1[((n)*(r-1)+1):(r*n), ]
  
  # DFM
  SFMfit <- SFM(Y = as.matrix(Y[-nrow(Y), ]), horizon = 1)
  DFMfit <- DFM(Y = as.matrix(Y[-nrow(Y), ]), f = SFMfit$f, rank = SFMfit$rank, horizon = 1, 
                lag.max = p, Yhat_static = SFMfit$Yhat_static,  decomp = SFMfit$decomp) 
  MSFEs_DFM <- (Y[nrow(Y), ] - DFMfit$Yhat_dynamic_AIC)^2
  
  # FAVAR
  FAVARfit <- FAVAR(Y = as.matrix(Y[-nrow(Y), ]), horizon = 1, lag.max = p)
  MSFEs_FAVAR <- (Y[nrow(Y), ] - FAVARfit$YhatsAIC)^2
  
  # AR
  Y <- matrix(as.matrix(Y), ncol=ncol(Y), nrow=nrow(Y))
  MSFEs_ar <- matrix(NA, ncol = ncol(Y), nrow = 1)
  for(i in 1:ncol(Y)){
    # AR
    ourar <- sparseVAR(Y[-nrow(Y), i], p = p); 
    MSFEs_ar[1, i] <- (Y[nrow(Y), i] - directforecast(fit = ourar, model ="VAR", h=1))^2
  }
  
  Sim1MSFE[r, 1] <- mean(MSFEs_DFM)
  Sim1MSFE[r, 2] <- mean(MSFEs_FAVAR)
  Sim1MSFE[r, 3] <- mean(MSFEs_ar)

}
apply(Sim1MSFE, 2, mean)

##########################################################
#### Simulation Scenario 2 : Own-Other HLag Structure ####
##########################################################

#### Setting ####
k=45;p=2;n=100
Nsim <- 500 # Number of simulations

#### Forecast Performance ####
Sim2MSFE <- matrix(NA, ncol=3, nrow=Nsim) # Note : VAR1 is already contained in the other file 
colnames(Sim2MSFE) <- c("DFM", "FAVAR", "AR")
datasim2 <- read.table("sim2dat.txt")
library(bigtime)
for(r in 1:Nsim){
  # Data
  Y <- datasim2[((n)*(r-1)+1):(r*n), ]
  
  # DFM
  SFMfit <- SFM(Y = as.matrix(Y[-nrow(Y), ]), horizon = 1)
  DFMfit <- DFM(Y = as.matrix(Y[-nrow(Y), ]), f = SFMfit$f, rank = SFMfit$rank, horizon = 1, 
                lag.max = p, Yhat_static = SFMfit$Yhat_static,  decomp = SFMfit$decomp) 
  MSFEs_DFM <- (Y[nrow(Y), ] - DFMfit$Yhat_dynamic_AIC)^2
  
  # FAVAR
  FAVARfit <- FAVAR(Y = as.matrix(Y[-nrow(Y), ]), horizon = 1, lag.max = p)
  MSFEs_FAVAR <- (Y[nrow(Y), ] - FAVARfit$YhatsAIC)^2
  
  # AR
  Y <- matrix(as.matrix(Y), ncol=ncol(Y), nrow=nrow(Y))
  MSFEs_ar <- matrix(NA, ncol = ncol(Y), nrow = 1)
  for(i in 1:ncol(Y)){
    # AR
    ourar <- sparseVAR(Y[-nrow(Y), i], p = p); 
    MSFEs_ar[1, i] <- (Y[nrow(Y), i] - directforecast(fit = ourar, model ="VAR", h=1))^2
  }
  
  Sim2MSFE[r, 1] <- mean(MSFEs_DFM)
  Sim2MSFE[r, 2] <- mean(MSFEs_FAVAR)
  Sim2MSFE[r, 3] <- mean(MSFEs_ar)
}
apply(Sim2MSFE, 2, mean)


############################################################
#### Simulation Scenario 3 : Elementwise HLag Structure ####
############################################################

#### Setting ####
k=45;p=4;n=100
Nsim <- 500 # Number of simulations
#### Forecast Performance ####
Sim3MSFE <- matrix(NA, ncol=3, nrow=Nsim) # Note : VAR1 is already contained in the other file 
colnames(Sim3MSFE) <- c("DFM", "FAVAR", "AR")
datasim3 <- read.table("sim3dat.txt")
library(bigtime)
for(r in 1:Nsim){
  # Data
  Y <- datasim3[((n)*(r-1)+1):(r*n), ]
  
  # DFM
  SFMfit <- SFM(Y = as.matrix(Y[-nrow(Y), ]), horizon = 1)
  DFMfit <- DFM(Y = as.matrix(Y[-nrow(Y), ]), f = SFMfit$f, rank = SFMfit$rank, horizon = 1, 
                lag.max = p, Yhat_static = SFMfit$Yhat_static,  decomp = SFMfit$decomp) 
  MSFEs_DFM <- (Y[nrow(Y), ] - DFMfit$Yhat_dynamic_AIC)^2
  
  # FAVAR
  FAVARfit <- FAVAR(Y = as.matrix(Y[-nrow(Y), ]), horizon = 1, lag.max = p)
  MSFEs_FAVAR <- (Y[nrow(Y), ] - FAVARfit$YhatsAIC)^2
  
  # AR
  Y <- matrix(as.matrix(Y), ncol=ncol(Y), nrow=nrow(Y))
  MSFEs_ar <- matrix(NA, ncol = ncol(Y), nrow = 1)
  for(i in 1:ncol(Y)){
    # AR
    ourar <- sparseVAR(Y[-nrow(Y), i], p = p); 
    MSFEs_ar[1, i] <- (Y[nrow(Y), i] - directforecast(fit = ourar, model ="VAR", h=1))^2
  }
  
  Sim3MSFE[r, 1] <- mean(MSFEs_DFM)
  Sim3MSFE[r, 2] <- mean(MSFEs_FAVAR)
  Sim3MSFE[r, 3] <- mean(MSFEs_ar)
}
apply(Sim3MSFE, 2, mean)

######################################################
#### Simulation Scenario 4 : Data-based Structure ####
######################################################

#### Setting ####
k <- 40; p <- 4; n=195
Nsim <- 500 # Number of simulations

#### Forecast Performance ####
Sim4MSFE <- matrix(NA, ncol=4, nrow=Nsim) # Note : VAR1 is already contained in the other file 
colnames(Sim4MSFE) <- c("DFM", "FAVAR", "AR", "VAR")
datasim4 <- read.table("sim4dat.txt")
library(bigtime)

for(r in 1:Nsim){
  # Data
  Y <- datasim4[((n)*(r-1)+1):(r*n), ]
  
  # DFM
  SFMfit <- SFM(Y = as.matrix(Y[-nrow(Y), ]), horizon = 1)
  DFMfit <- DFM(Y = as.matrix(Y[-nrow(Y), ]), f = SFMfit$f, rank = SFMfit$rank, horizon = 1, 
                lag.max = p, Yhat_static = SFMfit$Yhat_static,  decomp = SFMfit$decomp) 
  MSFEs_DFM <- (Y[nrow(Y), ] - DFMfit$Yhat_dynamic_AIC)^2
  
  # FAVAR
  FAVARfit <- FAVAR(Y = as.matrix(Y[-nrow(Y), ]), horizon = 1, lag.max = p)
  MSFEs_FAVAR <- (Y[nrow(Y), ] - FAVARfit$YhatsAIC)^2
  
  # AR
  Y <- matrix(as.matrix(Y), ncol=ncol(Y), nrow=nrow(Y))
  MSFEs_ar <- matrix(NA, ncol = ncol(Y), nrow = 1)
  for(i in 1:ncol(Y)){
    # AR
    ourar <- sparseVAR(Y[-nrow(Y), i], p = p); 
    MSFEs_ar[1, i] <- (Y[nrow(Y), i] - directforecast(fit = ourar, model ="VAR", h=1))^2
  }

  # VAR estimation
  VARfit <- VAR(y = Y[-nrow(Y),], type = "none")
  VARpredict <- predict(VARfit, n.ahead = 1)
  collectresults <- matrix(unlist(VARpredict$fcst), nrow = k, ncol = 4, byrow = T)
  VARforecast <- collectresults[,1]
  MSFEs_var <- (Y[nrow(Y), ] - VARforecast)^2
  
  Sim4MSFE[r, 1] <- mean(MSFEs_DFM)
  Sim4MSFE[r, 2] <- mean(MSFEs_FAVAR)
  Sim4MSFE[r, 3] <- mean(MSFEs_ar)
  Sim4MSFE[r, 4] <- mean(MSFEs_var)
}
apply(Sim4MSFE, 2, mean)
