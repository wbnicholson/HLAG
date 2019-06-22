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
koopact <- read.csv('SW.csv',head=TRUE)
koopact <- na.omit(koopact)
attach(koopact)
kooplarge <- as.matrix(cbind(GDP251,CPIAUCSL,FYFF,PSCCOMR,FMRNBA,FMRRA,FM2,GDP252,IPS10,UTL11,LHUR,HSFR,PWFSA,GDP273,CES275R,
                             FM1,FSPIN,FYGT10,EXRUS,CES002,Sfygt10,HHSNTN,PMI,PMDEL,PMCP,GDP256,LBOUT,PMNV,GDP263,GDP264,GDP265,
                             LBMNU,PMNO,CCINRV,BUSLOANS,PMP,GDP276_1,GDP270,GDP253,LHEL,GDP254,GDP255,GDP257,GDP258,GDP259,GDP260,
                             GDP261,GDP266,GDP267,GDP268,GDP269,GDP271,GDP272,GDP274,GDP275,GDP276,GDP277,GDP278,GDP279,GDP280,
                             GDP281,GDP282,GDP284,GDP285,GDP286,GDP287,GDP288,GDP289,GDP290,GDP291,GDP292,LBPUR7,LBLCPU,GDP274_1,
                             GDP274_2,GDP274_3,GDP275_1,GDP275_2,GDP275_3,GDP275_4,GDP276_2,GDP276_3,GDP276_4,GDP276_5,GDP276_6,
                             GDP276_7,GDP276_8,GDP284_1,GDP284_2,GDP285_1,GDP285_2,IPS11,IPS299,IPS12,IPS13,IPS18,IPS25,IPS32,IPS34,
                             IPS38,IPS43,IPS307,IPS306,CES275,CES277,CES278,CES277R,CES278.R,CES003,CES006,CES011,CES015,CES017,
                             CES033,CES046,CES048,CES049,CES053,CES088,CES140,LHELX,LHEM,LHNAG,LHU680,LHU5,LHU14,LHU15,LHU26,LHU27,
                             CES151,CES155,HSBR,HSNE,HSMW,HSSOU,HSWST,FYGM3,FYGM6,FYGT1,FYGT5,FYGT10,FYAAAC,FYBAAC,Sfygm6,Sfygt1,
                             sFYAAAC,sFYBAAC,MZMSL,FMFBA,CPILFESL,PCEPILFE,PWFCSA,PWIMSA,PWCMSA,PWCMSAR,PSCCOM,PW561,PW561R,EXRSW,
                             EXRJAN,EXRUK,EXRCAN,FSPCOM,FSDXP,FSPXE,FSDJ,MOCMQ,MSONDQ))

#### Forecast Accuracy ####
h <- 1# Specify Forecast Horizon, manuscript considers h = 1, h = 4, h = 8 

#### Set dimensions ####
p <- 4 # Number of lags
recursive <- T
k <- ncol(kooplarge)
Y <- kooplarge
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
MSFEs_others <- matrix(NA, 61, 4)
MSFEs_others_array <- array(NA, c(61, ncol(Ylarge), 4))
for(it in 1:61){
  cat("start", it, "\n")
  Yestim <- resHELEM@Ytrain_rol[[it]]
  
  # Factor-based models
  SFMfit <- SFM(Y = as.matrix(Yestim), r.max = min(ncol(Yestim), 10), horizon = h)
  DFMfit <- DFM(Y = as.matrix(Yestim), f = SFMfit$f, rank = SFMfit$rank, horizon = h, 
                lag.max = min(p, floor(nrow(Yestim)/SFMfit$rank)), Yhat_static = SFMfit$Yhat_static, decomp = SFMfit$decomp)
  FAVARfit <- FAVAR(Y = as.matrix(Yestim), horizon = h, lag.max = p)
  
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

GDP <- apply(MSFEindivlarge[,1, c(1:5,8,11:12)], 2, mean); round(t(GDP), 3)
CPI <- apply(MSFEindivlarge[,2, c(1:5,8,11:12)], 2, mean); round(t(CPI), 3)
FYFF <- apply(MSFEindivlarge[,3, c(1:5,8,11:12)], 2, mean); round(t(FYFF), 3)

#################
#### wMSFE  #####
#################
varTS <- apply(Ytest, 2, var)
wMSFEindivlarge <- MSFEindivlarge
for(i in 1:61){ # Time Point
  for(j in 1:ncol(Ylarge)){ # Variable 
    for(k in 1:14){# Methods
      wMSFEindivlarge[i, j, k] <- MSFEindivlarge[i, j, k]/varTS[j]
    }
  }
}

round(apply(wMSFEindivlarge, 3, mean), 3)

############################################
#### GROUPS OF MACRO_ECONOMIC VARIABLES ####
############################################
SWgroups <- read.table("SWgrouping.txt", header = T)
attach(SWgroups)

TABLE_GROUPS_1 <- matrix(NA, 14, 7)
rownames(TABLE_GROUPS_1) <- c("Componentwise", "Own-other", "Elementwise",
                              "Lasso", "Lag-weighted Lasso", "AIC", "BIC",
                              "BGR","Sample mean", "Random walk", 
                              "DFM", "FAVAR", 
                              "AR", "VAR1")
colnames(TABLE_GROUPS_1) <- names(SWgroups)[-1][1:7]

SWred <- SWgroups[1:dim(wMSFEindivlarge)[2], ]
MSFEs <-array(NA, c(61, dim(wMSFEindivlarge)[2], 14))
MSFEs[,,1:14] <- wMSFEindivlarge
names_methods <- c("Componentwise", "Own-other", "Elementwise",
                   "Lasso", "Lag-weighted Lasso", "AIC", "BIC",
                   "BGR","Sample mean", "Random walk", 
                   "DFM", "FAVAR", 
                   "AR", "VAR1")

##### GDP #####
GDPseries <- which(SWred$GDP==1)
MSFE_GDP_series <- MSFEs[, GDPseries,]
MSFE_GDP <- apply(MSFE_GDP_series, c(1,3), mean)
colnames(MSFE_GDP) <- names_methods
TABLE_GROUPS_1[, 1] <- apply(MSFE_GDP, 2, mean)

##### IP #####
IPseries <- which(SWred$IP==1)
MSFE_IPseries <- MSFEs[, IPseries,]
MSFE_IP <- apply(MSFE_IPseries, c(1,3), mean)
colnames(MSFE_IP) <- names_methods
TABLE_GROUPS_1[, 2] <- apply(MSFE_IP, 2, mean)

#### Employment ####
Emplseries <- which(SWred$Employment==1)
MSFE_Emplseries <- MSFEs[, Emplseries, ]
MSFE_Employment <- apply(MSFE_Emplseries, c(1,3), mean)
colnames(MSFE_Employment) <- names_methods
TABLE_GROUPS_1[, 3] <- apply(MSFE_Employment, 2, mean)

#### Unemployment ####
Unemplseries <- which(SWred$Unemployment==1)
MSFE_Unemplseries <- MSFEs[, Unemplseries, ]
MSFE_Unemployment <- apply(MSFE_Unemplseries, c(1,3), mean)
colnames(MSFE_Unemployment) <- names_methods
TABLE_GROUPS_1[, 4] <- apply(MSFE_Unemployment, 2, mean)

#### Housing  ####
Housingseries <- which(SWred$Housing==1)
MSFE_Housingseries <- MSFEs[, Housingseries, ]
MSFE_Housing <- apply(MSFE_Housingseries, c(1,3), mean)
colnames(MSFE_Housing) <- names_methods
TABLE_GROUPS_1[, 5] <- apply(MSFE_Housing, 2, mean)

#### Inventories ####
Invseries <- which(SWred$Inventories==1)
MSFEInvseries <- MSFEs[, Invseries, ]
MSFE_Inventory <- apply(MSFEInvseries, c(1,3), mean)
colnames(MSFE_Inventory) <- names_methods
TABLE_GROUPS_1[, 6] <- apply(MSFE_Inventory, 2, mean)

#### Prices  ####
Priceseries <- which(SWred$Prices==1)
MSFE_Priceseries <- MSFEs[, Priceseries, ]
MSFE_Prices <- apply(MSFE_Priceseries, c(1,3), mean)
colnames(MSFE_Prices) <- names_methods
TABLE_GROUPS_1[, 7] <- apply(MSFE_Prices, 2, mean)

round(TABLE_GROUPS_1, 3)

TABLE_GROUPS_2 <- matrix(NA, dim(TABLE_GROUPS_1)[1], dim(TABLE_GROUPS_1)[2]-1)
rownames(TABLE_GROUPS_2) <- rownames(TABLE_GROUPS_1)
colnames(TABLE_GROUPS_2) <- names(SWred)[-c(1:8)]

#### Wages  ####
Wagesseries <- which(SWred$Wages==1)
MSFE_Wagesseries <- MSFEs[, Wagesseries, ]
MSFE_Wages <- apply(MSFE_Wagesseries, c(1,3), mean)
colnames(MSFE_Wages) <- names_methods
TABLE_GROUPS_2[, 1] <- apply(MSFE_Wages, 2, mean)

#### Interest Rates ####
IRseries <- which(SWred$InterestRates==1)
MSFE_IRseries <- MSFEs[, IRseries, ]
MSFE_IR <- apply(MSFE_IRseries, c(1,3), mean)
TABLE_GROUPS_2[, 2] <- apply(MSFE_IR, 2, mean)

#### Money  ####
Money_series <- which(SWred$Money==1)
MSFE_Moneyseries <- MSFEs[, Money_series, ]
MSFE_Money <- apply(MSFE_Moneyseries, c(1,3), mean)
colnames(MSFE_Money) <- names_methods
TABLE_GROUPS_2[, 3] <- apply(MSFE_Money, 2, mean)

#### Exchange Rates  ####
Exchangeseries <- which(SWred$ExchangeRate==1)
MSFE_Exchangeseries <- MSFEs[, Exchangeseries, ]
MSFE_Exchange <- apply(MSFE_Exchangeseries, c(1,3), mean)
colnames(MSFE_Exchange) <- names_methods
TABLE_GROUPS_2[, 4] <- apply(MSFE_Exchange, 2, mean)

#### Stock Price ####
StockPriceseries <- which(SWred$StockPrice==1)
MSFE_StockPriceseries <- MSFEs[, StockPriceseries, ]
MSFE_StockPrice <- apply(MSFE_StockPriceseries, c(1,3), mean)
colnames(MSFE_StockPrice) <- names_methods
TABLE_GROUPS_2[, 5] <- apply(MSFE_StockPrice, 2, mean)

### Consumer Expencations  ####
ConsumerExpseries <- which(SWred$ConsumerExpectations==1)
MSFE_ConsumerExpseries <- MSFEs[, ConsumerExpseries, ]
MSFE_ConsumerExp <- MSFE_ConsumerExpseries
colnames(MSFE_ConsumerExp) <- names_methods
TABLE_GROUPS_2[, 6] <- apply(MSFE_ConsumerExp, 2, mean)

round(TABLE_GROUPS_2, 3)


