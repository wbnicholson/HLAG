#### Macro-Economic Application mediumlarge VAR model ####
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
koopmediumlarge <- as.matrix(cbind(GDP251,CPIAUCSL,FYFF,PSCCOMR,FMRNBA,FMRRA,FM2,GDP252,IPS10,UTL11,LHUR,HSFR,PWFSA,GDP273,CES275R,FM1,FSPIN,FYGT10,EXRUS,CES002,Sfygt10,HHSNTN,PMI,PMDEL,PMCP,GDP256,LBOUT,PMNV,GDP263,GDP264,GDP265,LBMNU,PMNO,CCINRV,BUSLOANS,PMP,GDP276_1,GDP270,GDP253,LHEL))


#### Forecast Accuracy ####
h <- 1# Specify Forecast Horizon, manuscript considers h = 1, h = 4, h = 8 

#### Set dimensions ####
p <- 4 # Number of lags
recursive <- T
k <- ncol(koopmediumlarge)
Y <- koopmediumlarge
Ymediumlarge <- Y

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
A <- constructModel(Y,p,"Basic", c(25,10), RVAR=FALSE,h, "Rolling", MN=FALSE, verbose=TRUE, T1=T1, T2=T2, 
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
A <- constructModel(Y,p,"BGR",c(25,10),RVAR=FALSE,h,"Rolling",MN=F,verbose=TRUE,T1=T1,T2=T2, 
                    recursive=recursive, standardize = standardize)
resBGR <- cv.BigVAR(A)
Yhats[, , 8] <- resBGR@preds

# Own Other
A <- constructModel(Y,p,"HVAROO",c(25,10),RVAR=FALSE,h,"Rolling",MN=FALSE,verbose=TRUE,T1=T1,T2=T2, 
                    recursive=recursive, standardize = standardize)
resHOO <- cv.BigVAR(A)
Yhats[, , 2] <- resHOO@preds

# Elementwise
A <- constructModel(Y,p,"HVARELEM",c(25,10),RVAR=FALSE,h,"Rolling",MN=FALSE,verbose=TRUE,T1=T1,T2=T2, 
                    recursive=recursive, standardize = standardize)
resHELEM <- cv.BigVAR(A)
Yhats[, , 3] <- resHELEM@preds


# Componentwise
A <- constructModel(Y,p,"HVARC",c(25,10),RVAR=FALSE,h,"Rolling",MN=FALSE,verbose=TRUE,T1=T1,T2=T2, recursive=recursive)
resHC <- cv.BigVAR(A)
Yhats[, , 1] <- resHC@preds

# Tapered
A <- constructModel(Y,p,"Tapered",c(25,10),RVAR=FALSE,h,"Rolling",MN=FALSE,verbose=TRUE,T1=T1,T2=T2, 
                    recursive=recursive, standardize = standardize)
resLW <- cv.BigVAR(A)
Yhats[, , 5] <- resLW@preds

VAR_LS_AICs <- resHC@AICpvec
VAR_LS_BICs <- resHC@BICpvec

# Factor Model and Simple Benchmarks
library(bigtime)
library(vars)
MSFEs_others <- matrix(NA, 61, 4)
MSFEs_others_array <- array(NA, c(61, ncol(Ymediumlarge), 4))
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
MSFEmediumlarge <- matrix(NA, ncol=1, nrow= 14)
colnames(MSFEmediumlarge) <- c("MSFE")
rownames(MSFEmediumlarge) <- names_methods
MSFEmediumlarge[, 1] <- c(mean(resHC@OOSMSFE), mean(resHOO@OOSMSFE), mean(resHELEM@OOSMSFE),
                     mean(resLasso@OOSMSFE), mean(resLW@OOSMSFE), mean(resHELEM@AICMSFE), mean(resHELEM@BICMSFE),
                     mean(resBGR@OOSMSFE), mean(resHELEM@MeanMSFE), mean(resHELEM@RWMSFE),
                     apply(MSFEs_others, 2, mean)) 
MSFEmediumlarge[1:10, 1] <- MSFEmediumlarge[1:10, 1]/ncol(Y)

##########################
#### Individual MSFEs ####
##########################
MSFEindivmediumlarge <- array(NA, dim(Yhats))
dimnames(MSFEindivmediumlarge) <- dimnames(Yhats)
for(i.ts in 1:10){
  MSFEindivmediumlarge[,,i.ts] <- (Yhats[,,i.ts] - Ytest)^2
}
MSFEindivmediumlarge[,,11:14] <- MSFEs_others_array

GDP <- apply(MSFEindivmediumlarge[,1, c(1:5,8,11:12)], 2, mean); round(t(GDP), 3)
CPI <- apply(MSFEindivmediumlarge[,2, c(1:5,8,11:12)], 2, mean); round(t(CPI), 3)
FYFF <- apply(MSFEindivmediumlarge[,3, c(1:5,8,11:12)], 2, mean); round(t(FYFF), 3)

#################
#### wMSFE  #####
#################
varTS <- apply(Ytest, 2, var)
wMSFEindivmediumlarge <- MSFEindivmediumlarge
for(i in 1:61){ # Time Point
  for(j in 1:ncol(Ymediumlarge)){ # Variable 
    for(k in 1:14){# Methods
      wMSFEindivmediumlarge[i, j, k] <- MSFEindivmediumlarge[i, j, k]/varTS[j]
    }
  }
}

round(apply(wMSFEindivmediumlarge, 3, mean), 3)

##########################################
#### LAG ORDER STABILITY ACROSS TIME  ####
##########################################
betaHELEM_mediumlarge_h1_new <- resHELEM@betalist
names_medium_large <- colnames(koopmediumlarge)

# Lhat matrix
p <- 4
LhatHELEM <- array(NA, c(dim(betaHELEM_mediumlarge_h1_new[[1]])[1], dim(betaHELEM_mediumlarge_h1_new[[1]])[1], 
                         length(betaHELEM_mediumlarge_h1_new)),
                   dimnames = list(names_medium_large, names_medium_large, paste0("t=", 1:61)))
for(i.t in 1:length(betaHELEM_mediumlarge_h1_new)){
  for(i.eq in 1:nrow(betaHELEM_mediumlarge_h1_new[[1]])){
    for(i.var in 1:nrow(betaHELEM_mediumlarge_h1_new[[1]])){
      betas <- betaHELEM_mediumlarge_h1_new[[i.t]][i.eq,seq(from = i.var+1, length = p, by = dim(betaHELEM_mediumlarge_h1_new[[1]])[1])]
      LhatHELEM[i.eq, i.var, i.t] <- length(which(betas!=0))
    }
  }
}

# Series according to macro-economic categories 
SWgroups <- read.table("SWgrouping.txt", header = T)
attach(SWgroups)
SWred <- SWgroups[1:40, ]

GDP <- which(SWred$GDP==1)
IP <- which(SWred$IP==1)
Employment <- which(SWred$Employment==1)
Unemployment <- which(SWred$Unemployment==1)
Housing <- which(SWred$Housing==1)
Inventories <- which(SWred$Inventories==1)
Prices <- which(SWred$Prices==1)
Wages <- which(SWred$Wages==1)
InterestRates <- which(SWred$InterestRates==1)
Money <- which(SWred$Money==1)
ExchangeRate <- which(SWred$ExchangeRate==1)
StockPrice <- which(SWred$StockPrice==1)
ConsumerExpectations <- which(SWred$ConsumerExpectations==1)

mediumlarge_groups <- c(1:13)
HELEMres <- array(NA, c(3, ncol(SWred), 61), 
                  dimnames = list(c("GDP", "CPI", "FYFF"), c("total", names(SWred)[-c(1)]), paste0("t=", 1:61)))

for(it in 1:61){
  for(iv in 1:3){
    datagroup <- LhatHELEM[iv,,it]!=0
    HELEMres[which(iv==c(1:3)), 1, it] <- sum(datagroup)
    for(ig2 in 1:13){
      HELEMres[which(iv==c(1:3)), 1 + which(ig2==mediumlarge_groups), it] <- sum(datagroup[which(SWred[, 1 + ig2]==1)])
    }
  }
}

library(RColorBrewer)
coul = colorRampPalette(c("blue", "white"))( 4 ) ## (n)

names(SWred)[2] <- "GDP components"
names(SWred)[5] <- "Unemployment rate"
names(SWred)[10] <- "Interest rates"
names(SWred)[12] <- "Exchange rates"
names(SWred)[13] <- "Stock prices"
names(SWred)[14] <- "Consumer expectations"
tpoints <- c("1992Q3",	"1992Q4",	"1993Q1",	"1993Q2",	"1993Q3",	"1993Q4",	"1994Q1",	"1994Q2",	"1994Q3",	"1994Q4",	
             "1995Q1",	"1995Q2",	"1995Q3",	"1995Q4",	"1996Q1",	"1996Q2",	"1996Q3",	"1996Q4",	"1997Q1",	"1997Q2",	
             "1997Q3",	"1997Q4",	"1998Q1",	"1998Q2",	"1998Q3",	"1998Q4",	"1999Q1",	"1999Q2",	"1999Q3",	"1999Q4",	
             "2000Q1",	"2000Q2",	"2000Q3",	"2000Q4",	"2001Q1",	"2001Q2",	"2001Q3",	"2001Q4",	"2002Q1",	"2002Q2",	
             "2002Q3",	"2002Q4",	"2003Q1",	"2003Q2",	"2003Q3",	"2003Q4",	"2004Q1",	"2004Q2",	"2004Q3",	"2004Q4",	
             "2005Q1",	"2005Q2",	"2005Q3",	"2005Q4",	"2006Q1",	"2006Q2",	"2006Q3",	"2006Q4",	"2007Q1",	"2007Q2",	"2007Q3")
ARRAY_HELEM <- array(NA, c(3, ncol(SWred)-1, 61), dimnames = list(c("GDP", "CPI", "FYFF"),
                                                                  names(SWred)[-c(1)], tpoints))

for(it in 1:61){
  for(ig1 in 1:3){
    ARRAY_HELEM[ig1, , it] <- HELEMres[ig1, -1, it]/sum(HELEMres[ig1, -1, it])
  }
}

library(gplots)
coul = colorRampPalette(c("blue", "white"))(10) ## (n)

# LAG STABILITY GDP
heatmap(ARRAY_HELEM[1,13:1,], Colv = NA, Rowv = NA,
        scale = "none",
        col = coul[length(coul):1],
        main = "GDP251",
        xlab = "End Point Rolling Window",
        ylab = " ",
        margins = c(5, 5))

# LAG STABILITY CPI
heatmap(ARRAY_HELEM[2,13:1,], Colv = NA, Rowv = NA,
        scale = "none",
        col = coul[length(coul):1],
        main = "CPIAUSL",
        xlab = "End Point Rolling Window",
        ylab = " ",
        margins = c(5, 5))

# LAG STABILITY FYFF
heatmap(ARRAY_HELEM[3,13:1,], Colv = NA, Rowv = NA,
        scale = "none",
        col = coul[length(coul):1],
        main = "FYFF",
        xlab = "End Point Rolling Window",
        ylab = " ",
        margins = c(5, 5))

##############################################
#### LAG ORDER SELECTION PLOT FULL SAMPLE ####
##############################################
source("fullsample.R")
h <- 1
p <- 4 # Number of lags
recursive <- T
k <- ncol(koopmediumlarge)
Y <- koopmediumlarge
# transform to zero mean, unit variance
for (i in 1:k) {
  Y[, i] <- Y[, i]/apply(Y, 2, sd)[i]
}
Y <- t(Y)
Y <- Y - c(apply(Y, 1, mean)) %*% t(c(rep(1, ncol(Y))))
Y <- t(Y)
Ymedium <- Y

# Set dimensions forecast exercise 
T1 <- floor(nrow(Y)/3)+p
T2 <- floor(2*nrow(Y)/3)+p



# Elementwise
A <- constructModel(Y,p,"HVARELEM",c(25,10), RVAR=F,h,"Rolling",
                    MN=FALSE,verbose=TRUE,
                    ONESE=F, recursive=T) #  full sample
resHELEM <- cv.BigVAR(A)

library(grid)
B1 <- resHELEM@betaPred[,2:ncol(resHELEM@betaPred)]

Names <- names(as.data.frame(koopmediumlarge))
names(koopmediumlarge)
library(reshape2)
library(ggplot2)
k=40;p=4

LagMaxElem <-   LagMatrix(B1,k,p,10e-8)
dat <- LagMaxElem[1:3,]

dat <- melt(dat)
dat$Var1 <- ifelse(dat$Var1==1,.25,dat$Var1)
dat$Var1 <- ifelse(dat$Var1==2,.5,dat$Var1)
dat$Var1 <- ifelse(dat$Var1==3,.75,dat$Var1)
dat$Var2 <- ifelse(dat$Var2==0,"",dat$Var2)

N2 <- rep("",40)
N2[1:3] <- Names[1:3]
N2[29] <- Names[29]
N2[18] <- Names[18]
N2[14] <- Names[14]
N2[8] <- Names[8]
N2[30] <- Names[30]
N2[37] <- Names[37]
N2[38] <- Names[38]
N2[39] <- Names[39]
N2[26] <- Names[26]

breaks=c(.25,.5,.75)

P <- ggplot(data =  dat, aes(y = rev(Var1), x = Var2)) +
  geom_tile(aes(fill =ifelse(value>0,as.numeric(round(value)),0)), colour = "white") +scale_y_continuous("",breaks=c(.25,.5,.75)
                                                                                                         , labels=rev(N2[1:3]),expand=c(0.001,0.001))+ scale_x_continuous("",breaks=c(1:40),labels=N2,expand=c(0.001,0.001))+
  geom_text(aes(label = ifelse(value>0,as.numeric(round(value)),"")), vjust = 1)+scale_fill_gradient(low = "white", high = "steelblue")+coord_fixed(ratio=10)+theme(axis.text.x=element_text(angle=60, size=10, vjust=0.5),legend.position="none")

print(P)