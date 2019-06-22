#### Macro-Economic Application large VAR model ####
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
FREDlarge <- data

#### Forecast Accuracy ####
h <- 1# Specify Forecast Horizon, manuscript considers h = 1, h = 4, h = 8 

#### Set dimensions ####
p <- 4 # Number of lags
recursive <- T
k <- ncol(FREDlarge)
Y <- FREDlarge
Ylarge <- Y

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
MSFEs_others_array <- array(NA, c(76, ncol(Ylarge), 4))
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
    MSFEs_var[1, i] <- NA
  }
  
  MSFEs_others_array[it, , 3] <- MSFEs_ar
  MSFEs_others_array[it, , 4] <- MSFEs_var
  MSFEs_others[it, 3] <- mean(MSFEs_ar)
  MSFEs_others[it, 4] <- NA
}

###############
#### MSFEs ####
###############
MSFElarge <- matrix(NA, ncol=1, nrow= 14)
colnames(MSFElarge) <- c("MSFE")
rownames(MSFElarge) <- names_methods
MSFElarge[, 1] <- c(mean(resHC@OOSMSFE), mean(resHOO@OOSMSFE), mean(resHELEM@OOSMSFE),
                          mean(resLasso@OOSMSFE), mean(resLW@OOSMSFE), mean(resHELEM@AICMSFE), mean(resHELEM@BICMSFE),
                          mean(resBGR@OOSMSFE), mean(resHELEM@MeanMSFE), mean(resHELEM@RWMSFE),
                          apply(MSFEs_others, 2, mean)) 
MSFElarge[1:10, 1] <- MSFElarge[1:10, 1]/ncol(Y)

##########################
#### Individual MSFEs ####
##########################
MSFEindivlarge <- array(NA, dim(Yhats))
dimnames(MSFEindivlarge) <- dimnames(Yhats)
for(i.ts in 1:10){
  MSFEindivlarge[,,i.ts] <- (Yhats[,,i.ts] - Ytest)^2
}
MSFEindivlarge[,,11:14] <- MSFEs_others_array

GDPC1 <- apply(MSFEindivlarge[,which(colnames(data)=="GDPC1"), c(1:5,8,11:12)], 2, mean); round(t(GDP), 3)
CPIAUCSL <- apply(MSFEindivlarge[,which(colnames(data)=="CPIAUCSL"), c(1:5,8,11:12)], 2, mean); round(t(CPIAUCSL), 3)
FEDFUNDS <- apply(MSFEindivlarge[,which(colnames(data)=="FEDFUNDS"), c(1:5,8,11:12)], 2, mean); round(t(FEDFUNDS), 3)

#################
#### wMSFE  #####
#################
varTS <- apply(Ytest, 2, var)
wMSFEindivlarge <- MSFEindivlarge
for(i in 1:76){ # Time Point
  for(j in 1:ncol(Ylarge)){ # Variable 
    for(k in 1:14){# Methods
      wMSFEindivlarge[i, j, k] <- MSFEindivlarge[i, j, k]/varTS[j]
    }
  }
}

round(apply(wMSFEindivlarge, 3, mean), 3)

###################################
#### MACRO-ECONOMIC CATEGORIES ####
###################################
FREDgroups <- read.table("FRED_groups.txt", header = T)
attach(FREDgroups)

TABLE_GROUPS_1 <- matrix(NA, 14, 7)
rownames(TABLE_GROUPS_1) <- c("Componentwise", "Own-other", "Elementwise",
                              "Lasso", "Lag-weighted Lasso", "AIC", "BIC",
                              "BGR","Sample mean", "Random walk", 
                              "DFM", "FAVAR", 
                              "AR", "VAR1")
colnames(TABLE_GROUPS_1) <- names(FREDgroups)[-1][1:7]

FREDred <- FREDgroups[1:dim(wMSFEindivlarge)[2], ]
MSFEs <-array(NA, c(76, dim(wMSFEindivlarge)[2], 14))
MSFEs[,,1:14] <- wMSFEindivlarge
names_methods <- c("Componentwise", "Own-other", "Elementwise",
                   "Lasso", "Lag-weighted Lasso", "AIC", "BIC",
                   "BGR","Sample mean", "Random walk", 
                   "DFM", "FAVAR", 
                   "AR", "VAR1")


##### NIPA  #####
series <- which(FREDred$NIPA==1)
MSFE_series <- MSFEs[, series,]
MSFE_mseries <- apply(MSFE_series, c(1,3), mean)
colnames(MSFE_mseries) <- names_methods
TABLE_GROUPS_1[, 1] <- apply(MSFE_mseries, 2, mean)

##### IP #####
series <- which(FREDred$IP==1)
MSFE_series <- MSFEs[, series,]
MSFE_mseries <- apply(MSFE_series, c(1,3), mean)
colnames(MSFE_mseries) <- names_methods
TABLE_GROUPS_1[, 2] <- apply(MSFE_mseries, 2, mean)

#### Employment####
series <- which(FREDred$Emp==1)
MSFE_series <- MSFEs[, series,]
MSFE_mseries <- apply(MSFE_series, c(1,3), mean)
colnames(MSFE_mseries) <- names_methods
TABLE_GROUPS_1[, 3] <- apply(MSFE_mseries, 2, mean)

#### Housing  ####
series <- which(FREDred$Housing==1)
MSFE_series <- MSFEs[, series,]
MSFE_mseries <- apply(MSFE_series, c(1,3), mean)
colnames(MSFE_mseries) <- names_methods
TABLE_GROUPS_1[, 4] <- apply(MSFE_mseries, 2, mean)

#### IOS####
series <- which(FREDred$IOS==1)
MSFE_series <- MSFEs[, series,]
MSFE_mseries <- apply(MSFE_series, c(1,3), mean)
colnames(MSFE_mseries) <- names_methods
TABLE_GROUPS_1[, 5] <- apply(MSFE_mseries, 2, mean)

#### Prices  ####
series <- which(FREDred$Prices==1)
MSFE_series <- MSFEs[, series,]
MSFE_mseries <- apply(MSFE_series, c(1,3), mean)
colnames(MSFE_mseries) <- names_methods
TABLE_GROUPS_1[, 6] <- apply(MSFE_mseries, 2, mean)

#### Earnings  ####
series <- which(FREDred$Earnings==1)
MSFE_series <- MSFEs[, series,]
MSFE_mseries <- apply(MSFE_series, c(1,3), mean)
colnames(MSFE_mseries) <- names_methods
TABLE_GROUPS_1[, 7] <- apply(MSFE_mseries, 2, mean)

round(TABLE_GROUPS_1, 3)

TABLE_GROUPS_2 <- matrix(NA, dim(TABLE_GROUPS_1)[1], dim(TABLE_GROUPS_1)[2]-1)
rownames(TABLE_GROUPS_2) <- rownames(TABLE_GROUPS_1)
colnames(TABLE_GROUPS_2) <- names(FREDred)[-c(1:8)][-5]

#### Interest Rates ####
series <- which(FREDred$Interest==1)
MSFE_series <- MSFEs[, series,]
MSFE_mseries <- apply(MSFE_series, c(1,3), mean)
colnames(MSFE_mseries) <- names_methods
TABLE_GROUPS_2[, 1] <- apply(MSFE_mseries, 2, mean)

#### Money ####
series <- which(FREDred$Money==1)
MSFE_series <- MSFEs[, series,]
MSFE_mseries <- apply(MSFE_series, c(1,3), mean)
colnames(MSFE_mseries) <- names_methods
TABLE_GROUPS_2[, 2] <- apply(MSFE_mseries, 2, mean)

#### HBS ####
series <- which(FREDred$HBS==1)
MSFE_series <- MSFEs[, series,]
MSFE_mseries <- apply(MSFE_series, c(1,3), mean)
colnames(MSFE_mseries) <- names_methods
TABLE_GROUPS_2[, 3] <- apply(MSFE_mseries, 2, mean)

#### Exchange Rates  ####
series <- which(FREDred$ER==1)
MSFE_series <- MSFEs[, series,]
MSFE_mseries <- apply(MSFE_series, c(1,3), mean)
colnames(MSFE_mseries) <- names_methods
TABLE_GROUPS_2[, 4] <- apply(MSFE_mseries, 2, mean)

#### Stock Market  ####
series <- which(FREDred$SM==1)
MSFE_series <- MSFEs[, series,]
MSFE_mseries <- apply(MSFE_series, c(1,3), mean)
colnames(MSFE_mseries) <- names_methods
TABLE_GROUPS_2[, 5] <- apply(MSFE_mseries, 2, mean)

### NHBS ####
series <- which(FREDred$NHBS==1)
MSFE_series <- MSFEs[, series,]
MSFE_mseries <- apply(MSFE_series, c(1,3), mean)
colnames(MSFE_mseries) <- names_methods
TABLE_GROUPS_2[, 6] <- apply(MSFE_mseries, 2, mean)

round(TABLE_GROUPS_2, 3)