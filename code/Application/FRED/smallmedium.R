#### Macro-Economic Application smallmedium VAR model ####
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
source("factorfunctions.R") # Functions for factor-based models
sourceCpp('mainfunctionsinC.cpp') # Main Functions in C
sourceCpp('auxftc.cpp') # Auxiliary Functions

#### Import the data ####
load("data_FRED_QD_current_all.RData")
data <- data_FRED_QD_current_all
smallmediumset <- c(which(colnames(data)=="GDPC1"),
                    which(colnames(data)=="CPIAUCSL"),
                    which(colnames(data)=="FEDFUNDS"),
                    which(colnames(data)=="PPIACO"),
                    which(colnames(data)=="CMRMTSPLx"),
                    which(colnames(data)=="TABSHNOx"),
                    which(colnames(data)=="M2REAL"),
                    which(colnames(data)=="PCECC96"),
                    which(colnames(data)=="INDPRO"),
                    which(colnames(data)=="CUMFNS")) 
FREDsmallmedium <- data[, smallmediumset]

#### Forecast Accuracy ####
h <- 1# Specify Forecast Horizon, manuscript considers h = 1, h = 4, h = 8 

#### Set dimensions ####
p <- 4 # Number of lags
recursive <- T
k <- ncol(FREDsmallmedium)
Y <- FREDsmallmedium
Ysmallmedium <- Y

# Set dimensions forecast exercise 
T1 <- floor(nrow(Y)/3)+p
T2 <- floor(2*nrow(Y)/3)+p

#### Forecast Exercise ####
names_methods <- c("Componentwise", "Own-other", "Elementwise",
                   "Lasso", "Lag-weighted Lasso", "AIC", "BIC",
                   "BGR","Sample mean", "Random walk", 
                   "DFM", "FAVAR", 
                   "AR", "VAR1")
Yhats <- array(NA, c(nrow(Y) - T2, k, 14),
               dimnames=list(paste0("t=",1:(nrow(Y) - T2)), colnames(Y), 
                             names_methods))

standardize <- T
train_roll <- T
test_roll <- T
A <- constructModel(Y,p,"Basic", c(50,10), RVAR=FALSE,h, "Rolling", MN=FALSE, verbose=TRUE, T1=T1, T2=T2, 
                    recursive=recursive, standardize = standardize)

# Lasso
resLasso <- cv.BigVAR(A)
Yhats[, , 4] <- resLasso@preds
Ytest <- resLasso@Ytest_stand

# AIC
Yhats[, , 6] <- resLasso@AICPreds

# BIC
Yhats[, , 7] <- resLasso@BICPreds

# Mean
Yhats[, , 9] <- resLasso@MeanPreds

# Random Walk
Yhats[, , 10] <- resLasso@RWPreds

# BGR
A <- constructModel(Y,p,"BGR",c(50,10),RVAR=FALSE,h,"Rolling",MN=F,verbose=TRUE,T1=T1,T2=T2, 
                    recursive=recursive, standardize = standardize)
resBGR <- cv.BigVAR(A)
Yhats[, , 8] <- resBGR@preds

# Own Other
A <- constructModel(Y,p,"HVAROO",c(50,10),RVAR=FALSE,h,"Rolling",MN=FALSE,verbose=TRUE,T1=T1,T2=T2, 
                    recursive=recursive, standardize = standardize)
resHOO <- cv.BigVAR(A)
Yhats[, , 2] <- resHOO@preds

# Elementwise
A <- constructModel(Y,p,"HVARELEM",c(50,10),RVAR=FALSE,h,"Rolling",MN=FALSE,verbose=TRUE,T1=T1,T2=T2, 
                    recursive=recursive, standardize = standardize)
resHELEM <- cv.BigVAR(A)
Yhats[, , 3] <- resHELEM@preds


# Componentwise
A <- constructModel(Y,p,"HVARC",c(50,10),RVAR=FALSE,h,"Rolling",MN=FALSE,verbose=TRUE,T1=T1,T2=T2, recursive=recursive)
resHC <- cv.BigVAR(A)
Yhats[, , 1] <- resHC@preds

# Tapered
A <- constructModel(Y,p,"Tapered",c(50,10),RVAR=FALSE,h,"Rolling",MN=FALSE,verbose=TRUE,T1=T1,T2=T2, 
                    recursive=recursive, standardize = standardize)
resLW <- cv.BigVAR(A)
Yhats[, , 5] <- resLW@preds

VAR_LS_AICs <- resHC@AICpvec
VAR_LS_BICs <- resHC@BICpvec

# Factor Model and Simple Benchmarks
library(bigtime)
library(vars)
MSFEs_others <- matrix(NA, 76, 4)
MSFEs_others_array <- array(NA, c(76, ncol(Ysmallmedium), 4))
for(it in 1:76){
  Yestim <- resHELEM@Ytrain_rol[[it]]
  
  # Factor-based models
  SFMfit <- SFM(Y = as.matrix(Yestim), r.max = min(ncol(Yestim), 5), horizon = h)
  DFMfit <- DFM(Y = as.matrix(Yestim), f = SFMfit$f, rank = SFMfit$rank, horizon = h, 
                lag.max = min(p, floor(nrow(Yestim)/SFMfit$rank)), Yhat_static = SFMfit$Yhat_static, decomp = SFMfit$decomp)
  FAVARfit <- FAVAR(Y = as.matrix(Yestim), horizon = h, lag.max = p, r.max = min(ncol(Yestim) - 1, 5))
  
  MSFEs_DFM <- MSFEs_FAVAR <- matrix(NA, ncol = ncol(Y), nrow = 1)
  MSFEs_DFM <- (resHELEM@Ytest_stand[it,] - DFMfit$Yhat_dynamic_AIC[h,])^2
  MSFEs_FAVAR <- (resHELEM@Ytest_stand[it,] - FAVARfit$YhatsAIC[h,])^2
  
  MSFEs_others_array[it, , 1] <- MSFEs_DFM
  MSFEs_others_array[it, , 2] <- MSFEs_FAVAR
  MSFEs_others[it, 1] <- mean(MSFEs_DFM)
  MSFEs_others[it, 2] <- mean(MSFEs_FAVAR)
  
  # VAR and AR
  MSFEs_ar <- MSFEs_var <- matrix(NA, ncol = ncol(Y), nrow = 1)
  VARfit <- VAR(y = as.data.frame(Yestim), type = "none", p = 1)
  VARpredict <- predict(VARfit, n.ahead = h)
  for(i in 1:ncol(Yestim)){
    arp4 <- sparseVAR(Yestim[, i],p = 4, h=h)
    pAR <- arp4$p
    for(i.h in 1:h){
      if(i.h ==1){
        Ylastc <- Ylast <- Yestim[nrow(Yestim) : ( nrow(Yestim) - pAR + 1), i]
      }
      
      Yhatc <- arp4$phi0hat + arp4$Phihat%*%Ylast[1:pAR]
      Yhat <- arp4$Phihat%*%Ylast[1:pAR]
      Ylastc <- c(Yhatc, Ylastc)
      Ylast <- c(Yhat, Ylast)
    }
    MSFEs_ar[1, i] <- (resHELEM@Ytest_stand[it, i] - Yhatc)^2
    MSFEs_var[1, i] <- (resHELEM@Ytest_stand[it, i] - VARpredict$fcst[[i]][h, 1])^2
  }
  
  MSFEs_others_array[it, , 3] <- MSFEs_ar
  MSFEs_others_array[it, , 4] <- MSFEs_var
  MSFEs_others[it, 3] <- mean(MSFEs_ar)
  MSFEs_others[it, 4] <- mean(MSFEs_var)
}

###############
#### MSFEs ####
###############
MSFEsmallmedium <- matrix(NA, ncol=1, nrow= 14)
colnames(MSFEsmallmedium) <- c("MSFE")
rownames(MSFEsmallmedium) <- names_methods
MSFEsmallmedium[, 1] <- c(mean(resHC@OOSMSFE), mean(resHOO@OOSMSFE), mean(resHELEM@OOSMSFE),
                          mean(resLasso@OOSMSFE), mean(resLW@OOSMSFE), mean(resHELEM@AICMSFE), mean(resHELEM@BICMSFE),
                          mean(resBGR@OOSMSFE), mean(resHELEM@MeanMSFE), mean(resHELEM@RWMSFE),
                          apply(MSFEs_others, 2, mean)) 
MSFEsmallmedium[1:10, 1] <- MSFEsmallmedium[1:10, 1]/ncol(Y)

##########################
#### Individual MSFEs ####
##########################
MSFEindivsmallmedium <- array(NA, dim(Yhats))
dimnames(MSFEindivsmallmedium) <- dimnames(Yhats)
for(i.ts in 1:10){
  MSFEindivsmallmedium[,,i.ts] <- (Yhats[,,i.ts] - Ytest)^2
}
MSFEindivsmallmedium[,,11:14] <- MSFEs_others_array

GDP <- apply(MSFEindivsmallmedium[,1, c(1:5,8,11:12)], 2, mean); round(t(GDP), 3)
CPI <- apply(MSFEindivsmallmedium[,2, c(1:5,8,11:12)], 2, mean); round(t(CPI), 3)
FYFF <- apply(MSFEindivsmallmedium[,3, c(1:5,8,11:12)], 2, mean); round(t(FYFF), 3)

#################
#### wMSFE  #####
#################
varTS <- apply(Ytest, 2, var)
wMSFEindivsmallmedium <- MSFEindivsmallmedium
for(i in 1:76){ # Time Point
  for(j in 1:ncol(Ysmallmedium)){ # Variable 
    for(k in 1:14){# Methods
      wMSFEindivsmallmedium[i, j, k] <- MSFEindivsmallmedium[i, j, k]/varTS[j]
    }
  }
}

round(apply(wMSFEindivsmallmedium, 3, mean), 3)
