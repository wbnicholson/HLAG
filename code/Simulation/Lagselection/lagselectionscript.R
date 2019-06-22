#### R-script Simulation Study : Measuring Lag order selection performance ####
# Set your working directory to Code\Simulation\Lagselection
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
                  "graphics","RcppArmadillo", "RcppEigen", "R.matlab")
lapply(packagelist,checkpackage)

#### Source Functions ####
source("mainfunctionsinR.R") # Main Functions in R 
sourceCpp('mainfunctionsinC.cpp') # Main Functions in C
sourceCpp('auxftc.cpp') # Auxiliary Functions

#####################################################
#### Simulation Scenario 5 : Robustness Scenario ####
#####################################################

#### Load objects to reproduce sim results ####
load(file="TrueZerosS5.RData")
load(file="YYS5.RData")
#### Details on how objects YYS5 and TrueZerosS5 were created ####
# p=5;k=10;n=200
# Fp = matrix(1,k*p,k*p)
# while(max(Mod(eigen(Fp)$values)>1))
# {
#   Aex <- matrix(0,k,k*p)
#   
#   for(i in 1:k)
#   {
#     if(i<=4)
#     {
#       bb <- runif(k,-.4,.4)
#       Aex[i,1:(k)] <- bb
#       Aex[i,(k+1):(2*k)] <- -.9*bb
#       
#     }
#     
#     if(i>=5 & i<=7)
#     {
#       
#       bb <- runif(k,-.6,.6)
#       Aex[i,1:(k)] <- bb
#       Aex[i,(k+1):(2*k)] <- -.8*bb
#       Aex[i,(2*k+1):(3*k)] <- .6*bb
#       Aex[i,(3*k+1):(4*k)] <- -.5*bb
#       Aex[i,(4*k+1):(5*k)] <- .45*bb
#     }
#     
#   }
#   
#   Fp = matrix(0,k*p,k*p)
#   Fp[1:k,] = Aex
#   Fp[-(1:k),1:(k*(p-1))] = diag(k*(p-1))
#   print(max(Mod(eigen(Fp)$values)))
#   
# }
# TrueZerosS5 <- LagMatrix(Fp[1:k,],k,p,10e-4)
# save(TrueZerosS5, file="TrueZerosS5.RData")

# Nsim <- 500 # Number of simulations
# YYS5 <- vector("list", Nsim)
# for(i in 1:Nsim){
#   YYS5[[i]] <- MultVarSim(k, A1=Fp, p, .01*diag(k), n)
# }
# save(YYS5, file="YYS5.RData")

#### Lag order selection Performance ####
Nsim <- 500
k=10; p=12; n=200
Sim5Lag <- matrix(NA, ncol=11, nrow=Nsim)
colnames(Sim5Lag) <- c("HLagC", "HLagOO", "HLagElem", "Lasso", "LagLasso", "LS AIC", "LS BIC", "Mean", "RW", "BGR", "GLP")

# Methods that do not perform lag selection performance
TrueZeros <- TrueZerosS5
Sim5Lag[,8] <- (sum(abs(TrueZeros-matrix(0, ncol=k, nrow=k))))/sum(abs(TrueZeros))
Sim5Lag[,9] <- (sum(abs(TrueZeros-diag(1,k))))/sum(abs(TrueZeros))
Sim5Lag[,10:11] <- (sum(abs(TrueZeros-matrix(12, ncol=k, nrow=k))))/sum(abs(TrueZeros))

for(r in 1:Nsim){
  
  # Data
  Y <- YYS5[[r]]
  
  A <- constructModel(Y, p=12, "Basic", gran=c(5000,10),verbose=F,ONESE=TRUE,RVAR=TRUE,
                      lagselect = TRUE,tol=1e-6, T1=floor(0.67*n), T2=n-2)
  
  # Lasso
  resH <- cv.BigVAR(A)
  Sim5Lag[r, 4] <- resH@SR[length(resH@SR)]

  # AIC and BIC 
  AICp <- resH@AICpvec[1]
  Sim5Lag[r, 6] <- (sum(abs(TrueZeros-matrix(AICp, ncol=k, nrow=k))))/sum(abs(TrueZeros))
  BICp <- resH@BICpvec[1]
  Sim5Lag[r, 7] <- (sum(abs(TrueZeros-matrix(BICp, ncol=k, nrow=k))))/sum(abs(TrueZeros))
  
  # Componentwise
  A@Structure <- "HVARC"
  resH <- cv.BigVAR(A)
  Sim5Lag[r, 1] <- resH@SR[length(resH@SR)]
  
  # Own-other
  A@Structure <- "HVAROO"
  resH <- cv.BigVAR(A)
  Sim5Lag[r, 2] <- resH@SR[length(resH@SR)]
  
  # Elementwise
  A@Structure <- "HVARELEM"
  resH <- cv.BigVAR(A)
  Sim5Lag[r, 3] <- resH@SR[length(resH@SR)]
  
  # Lag-Weighted Lasso
  A@Structure <- "Tapered"
  resH <- cv.BigVAR(A)
  Sim5Lag[r, 5] <- resH@SR[length(resH@SR)]


}
apply(Sim5Lag, 2, mean)
