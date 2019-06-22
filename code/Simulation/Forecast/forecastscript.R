#### R-script Simulation Study : Measuring Forecast Accuracy ####
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


##############################################################
#### Simulation Scenario 1 : COmponentwise HLag Structure ####
##############################################################

#### Create autoregressive matrices ####
p=5 # p the maximal lag order
k=45 # k is the number of series
n=100 # n is the time series length. We ran the simulations for both n=100 and n=200
Fp = matrix(1,k*p,k*p)
set.seed(1500)
while(max(Mod(eigen(Fp)$values)>1))
{

  Aex <- matrix(0,k,k*p)

  for(i in 1:k)
  {
    if(i<10)
    {
      Aex[i,1:(k)] <- runif(k,-.2,.2)

    }
    if(i>=10 & i<19)
    {
      Aex[i,1:(2*k)] <- runif(2*k,-.15,.15)

    }
    if(i>=19 & i<28)
    {

      Aex[i,1:(3*k)] <- runif(3*k,-.15,.15)


    }
    if(i>=28 & i<37)
    {


      Aex[i,1:(4*k)] <- runif(4*k,-.15,.15)


    }

    if(i>=37)
    {

      Aex[i,1:(5*k)] <- runif(5*k,-.125,.125)


    }



  }

  Fp = matrix(0,k*p,k*p)
  Fp[1:k,] = Aex
  Fp[-(1:k),1:(k*(p-1))] = diag(k*(p-1))
  print(max(Mod(eigen(Fp)$values)))

}

#### Figure autoregressive matrices ####
print(SparsityPlotLarge(Fp[1:k,1:(5*k)],5,k)) # Figure 4, panel (1)

#### Simulated data ####
Nsim <- 500 # Number of simulations
YY <- vector("list", Nsim)
set.seed(1989)
for(i in 1:Nsim){
  YY[[i]] <-  MultVarSim(k, A1=Fp, p, .01*diag(k), n)
  write(t(YY[[i]]),file=paste0("sim1dat.txt"), ncolumns=ncol(YY[[i]]), sep="\t", append=T)
}

#### Forecast Performance ####
Sim1MSFE <- matrix(NA, ncol=9, nrow=Nsim)
colnames(Sim1MSFE) <- c("HLagC", "HLagOO", "HLagElem", 
                        "Lasso", "LagLasso", "LS", 
                        "Mean", "RW", "BGR")

datasim1 <- read.table("sim1dat.txt")
for(r in 1:Nsim){

  # Data
  Y <- datasim1[((n)*(r-1)+1):(r*n), ]
  Y <- matrix(as.matrix(Y), ncol=ncol(Y), nrow=nrow(Y))
  
  A <- constructModel(Y, p=5, "Basic", gran=c(500,10),verbose=FALSE, 
                      T1=floor(0.67*nrow(Y)), T2=nrow(Y)-1) 
  
  # Lasso
  resH <- cv.BigVAR(A)
  Sim1MSFE[r, 4] <- mean(resH@OOSMSFE)

  # BGR
  A@Structure <- "BGR"
  resBGR <- cv.BigVAR(A)
  Sim1MSFE[r, 9] <- mean(resBGR@OOSMSFE)
  
  # Componentwise
  A@Structure <- "HVARC"
  resH <- cv.BigVAR(A)
  Sim1MSFE[r, 1] <- mean(resH@OOSMSFE)
  
  # Own-other
  A@Structure <- "HVAROO"
  resH <- cv.BigVAR(A)
  Sim1MSFE[r, 2] <- mean(resH@OOSMSFE)
  
  # Elementwise
  A@Structure <- "HVARELEM"
  resH <- cv.BigVAR(A)
  Sim1MSFE[r, 3] <- mean(resH@OOSMSFE)
  
  # Lag-Weighted Lasso
  A@Structure <- "Tapered"
  resH <- cv.BigVAR(A)
  Sim1MSFE[r, 5] <- mean(resH@OOSMSFE)
  
  # VAR(1), Mean and Random Walk
  Sim1MSFE[r, 6] <- mean(resH@AICMSFE)
  Sim1MSFE[r, 7] <- mean(resH@MeanMSFE)
  Sim1MSFE[r, 8] <- mean(resH@RWMSFE)

}
apply(Sim1MSFE, 2, mean)/k



##########################################################
#### Simulation Scenario 2 : Own-Other HLag Structure ####
##########################################################

#### Create autoregressive matrices ####
.vecoovars<-function(p,k,k1)
{
  vv=list()
  vv[[1]]=1:(p*k)
  vv[[2]]=vv[[1]][-k1]
  q1=1
  for(i in 3:(2*p))
  {
    if(i%%2!=0)
    {
      vv[[i]]=(q1*k+1):(k*p)
      q1=q1+1
    }
    else{
      vv[[i]]=vv[[i-1]][-k1]
    }
  }
  return(vv)
}

set.seed(1901)
k=45;p=2;n=100
Fp = matrix(1,k*p,k*p)
while(max(Mod(eigen(Fp)$values)>.99))
{
  Aex <- matrix(0,k,k*p)
  for(i in 1:k)
  {
    kk <- .vecoovars(1,k,i)
    if(i>15 & i<=30)
    {
      aa <- length(intersect(kk[[1]],kk[[2]]))
      aaa <- setdiff(kk[[1]],kk[[2]])
      cc <- length(setdiff(kk[[2]],kk[[3]]))
      bb <- runif(aa,-.25,.25)
      dd <- 1.25*max(bb)
      Aex[i,setdiff(kk[[1]],kk[[2]])] <- dd
      Aex[i,intersect(kk[[1]],kk[[2]])] <- bb
    }
    Aex[16:30,61:75] <- diag(runif(15,-.3,.3))
    kk <- .vecoovars(2,k,i)
    if(i>30)
    {
      bb <- runif(aa,-.2,.2)
      dd <- 1.25*max(bb)
      i1 <- setdiff(kk[[1]],kk[[2]])
      ij <- 1:k
      ij <- ij[-i1]
      i2 <- setdiff(kk[[3]],kk[[4]])
      ij2 <- (k +1):(2*k)
      ij2 <- ij2[(ij2!=i2)]
      bb <- runif(length(ij),-.25,.25)
      dd <- 1.25*max(bb)
      Aex[i,i1] <- dd
      Aex[i,ij] <- bb
      Aex[i,i2] <- .2*dd
      Aex[i,ij2] <- .2*bb



    }
  }
  Fp = matrix(0,k*p,k*p)
  Fp[1:k,] = Aex
  Fp[1:15,1:15] <- diag(runif(15,-.15,.15))
  Fp[-(1:k),1:(k*(p-1))] = diag(k*(p-1))
  print(max(Mod(eigen(Fp)$values)))
}

#### Figure autoregressive matrices ####
print(SparsityPlotLarge(Fp[1:k,1:(2*k)],2,k)) # Figure 4, panel (2)

#### Simulated data ####
Nsim <- 500 # Number of simulations
YY <- vector("list", Nsim)
set.seed(1990)
for(i in 1:Nsim){
  YY[[i]] <-  MultVarSim(k, A1=Fp, p, .01*diag(k), n)
  write(t(YY[[i]]),file=paste0("sim2dat.txt"), ncolumns=ncol(YY[[i]]), sep="\t", append=T)
}

#### Forecast Performance ####
Sim2MSFE <- matrix(NA, ncol=9, nrow=Nsim)
colnames(Sim2MSFE) <- c("HLagC", "HLagOO", "HLagElem", 
                        "Lasso", "LagLasso", "LS", 
                        "Mean", "RW", "BGR")

datasim2 <- read.table("sim2dat.txt")
for(r in 1:Nsim){
  
  # Data
  Y <- datasim2[((n)*(r-1)+1):(r*n), ]
  Y <- matrix(as.matrix(Y), ncol=ncol(Y), nrow=nrow(Y))
  
  A <- constructModel(Y, p=2, "Basic", gran=c(1000,10),verbose=FALSE, 
                      T1=floor(0.67*nrow(Y)), T2=nrow(Y)-1) 
  
  # Lasso
  resH <- cv.BigVAR(A)
  Sim2MSFE[r, 4] <- mean(resH@OOSMSFE)
  
  # BGR
  A@Structure <- "BGR"
  resBGR <- cv.BigVAR(A)
  Sim2MSFE[r, 9] <- mean(resBGR@OOSMSFE)
  
  # Componentwise
  A@Structure <- "HVARC"
  resH <- cv.BigVAR(A)
  Sim2MSFE[r, 1] <- mean(resH@OOSMSFE)
  
  # Own-other
  A@Structure <- "HVAROO"
  resH <- cv.BigVAR(A)
  Sim2MSFE[r, 2] <- mean(resH@OOSMSFE)
  
  # Elementwise
  A@Structure <- "HVARELEM"
  resH <- cv.BigVAR(A)
  Sim2MSFE[r, 3] <- mean(resH@OOSMSFE)
  
  # Lag-Weighted Lasso
  A@Structure <- "Tapered"
  resH <- cv.BigVAR(A)
  Sim2MSFE[r, 5] <- mean(resH@OOSMSFE)
  
  # VAR(1), Mean and Random Walk
  Sim2MSFE[r, 6] <- mean(resH@AICMSFE)
  Sim2MSFE[r, 7] <- mean(resH@MeanMSFE)
  Sim2MSFE[r, 8] <- mean(resH@RWMSFE)
  
}
apply(Sim2MSFE, 2, mean)/k



############################################################
#### Simulation Scenario 3 : Elementwise HLag Structure ####
############################################################

#### Create autoregressive matrices ####
k=45;p=4
set.seed(1000)
Fp = matrix(1,k*p,k*p)
while(max(Mod(eigen(Fp)$values)>1))
{
  Aex <- matrix(0,k,k*p)
  for(j in 1:10)
  {
    for(i in 1:30)
    {
      test <- runif(1)
      if(test<.1){
        bb <- seq(i,4*k,k)
        mag=.4
      }else{


        bb<- seq(i,k,k)
        mag=.1
      }


      Aex[j,bb] <- runif(1,-mag,mag)
    }
    for(i in 31:45)
    {

      test <- runif(1)
      if(test<.1){
        bb <- seq(i,4*k,k)
        mag <- .4
      }else{

        bb <- seq(i,k,k)
        mag=.1
      }

      Aex[j,bb] <- runif(length(bb),-mag,mag)
    }

  }

  for(j in 11:30)
  {
    for(i in 1:30)
    {
      test <- runif(1)
      if(test<.1){
        bb <- seq(i,4*k,k)
        mag <- .4

      }else{

        bb <- seq(i,k,k)
        mag=.1
      }
      Aex[j,bb] <- runif(length(bb),-mag,mag)
    }

    for(i in 31:45)
    {
      test <- runif(1)
      if(test<.1){
        bb <- seq(i,4*k,k)
        mag=.4
      }else{

        bb <- seq(i,k,k)
        mag=.1
      }
      Aex[j,bb] <- runif(length(bb),-mag,mag)
    }
  }

  for(j in 31:45)
  {
    for(i in 1:30)
    {

      test <- runif(1)
      if(test<.1){
        bb <- seq(i,4*k,k)
        mag=.3
      }else{

        bb <- seq(i,k,k)
        mag=.1
      }

      Aex[j,bb] <- runif(length(bb),-mag,mag)
    }
    for(i in 31:45)
    {
      test <- runif(1)
      if(test<.1){
        bb <- seq(i,4*k,k)
        mag=.3
      }else{

        bb <- seq(i,k,k)
        mag=.1
      }


      Aex[j,bb] <- runif(length(bb),-mag,mag)
    }
  }

  Fp = matrix(0,k*p,k*p)
  Fp[1:k,] = Aex
  Fp[-(1:k),1:(k*(p-1))] = diag(k*(p-1))
  print(max(Mod(eigen(Fp)$values)))
}

#### Figure autoregressive matrices ####
print(SparsityPlotLarge(Fp[1:k,1:(4*k)],4,k)) # Figure 4, panel (3)

#### Simulated data ####
Nsim <- 500 # Number of simulations
YY <- vector("list", Nsim)
set.seed(1900)
for(i in 1:Nsim){
  YY[[i]] <-  MultVarSim(k, A1=Fp, p, .01*diag(k), n)
  write(t(YY[[i]]),file=paste0("sim3dat.txt"), ncolumns=ncol(YY[[i]]), sep="\t", append=T)
}

#### Forecast Performance ####
Sim3MSFE <- matrix(NA, ncol=9, nrow=Nsim)
colnames(Sim3MSFE) <- c("HLagC", "HLagOO", "HLagElem", 
                        "Lasso", "LagLasso", "LS", 
                        "Mean", "RW", "BGR")

datasim3 <- read.table("sim3dat.txt")
for(r in 1:Nsim){
  
  # Data
  Y <- datasim3[((n)*(r-1)+1):(r*n), ]
  Y <- matrix(as.matrix(Y), ncol=ncol(Y), nrow=nrow(Y))
  
  A <- constructModel(Y, p=4, "Basic", gran=c(2500,10),verbose=FALSE, 
                      T1=floor(0.67*nrow(Y)), T2=nrow(Y)-1) 
  
  # Lasso
  resH <- cv.BigVAR(A)
  Sim3MSFE[r, 4] <- mean(resH@OOSMSFE)
  
  # BGR
  A@Structure <- "BGR"
  resBGR <- cv.BigVAR(A)
  Sim3MSFE[r, 9] <- mean(resBGR@OOSMSFE)
  
  # Componentwise
  A@Structure <- "HVARC"
  resH <- cv.BigVAR(A)
  Sim3MSFE[r, 1] <- mean(resH@OOSMSFE)
  
  # Own-other
  A@Structure <- "HVAROO"
  resH <- cv.BigVAR(A)
  Sim3MSFE[r, 2] <- mean(resH@OOSMSFE)
  
  # Elementwise
  A@Structure <- "HVARELEM"
  resH <- cv.BigVAR(A)
  Sim3MSFE[r, 3] <- mean(resH@OOSMSFE)
  
  # Lag-Weighted Lasso
  A@Structure <- "Tapered"
  resH <- cv.BigVAR(A)
  Sim3MSFE[r, 5] <- mean(resH@OOSMSFE)
  
  # VAR(1), Mean and Random Walk
  Sim3MSFE[r, 6] <- mean(resH@AICMSFE)
  Sim3MSFE[r, 7] <- mean(resH@MeanMSFE)
  Sim3MSFE[r, 8] <- mean(resH@RWMSFE)
}
apply(Sim3MSFE, 2, mean)




######################################################
#### Simulation Scenario 4 : Data-based Structure ####
######################################################

#### Import residuals and estimated coefficients from the macro data set ####
load("GLPresid.RData")
load("GLPbeta.RData")

#### Figure autoregressive matrices ####
k <- 40; p <- 4; n=195
B <- GLPbeta[,2:ncol(GLPbeta)]
BB <- VarptoVar1MC(B,4,k)
print(SparsityPlotLarge(BB[1:k,1:(4*k)],4,k)) # Figure 4, panel (4)

#### Simulated data ####
Nsim <- 500 # Number of simulations
YY <- vector("list", Nsim)
for(i in 1:Nsim){
  YY[[i]] <- SimBoot(k=k, A1=BB, p=p, Resid=GLPresid, T=n)
  write(t(YY[[i]]),file=paste0("sim4dat.txt"), ncolumns=ncol(YY[[i]]), sep="\t", append=T)
}

#### Forecast Performance ####
Sim4MSFE <- matrix(NA, ncol=10, nrow=Nsim)
colnames(Sim4MSFE) <- c("HLagC", "HLagOO", "HLagElem", 
                        "Lasso", "LagLasso", "LS AIC", "LS BIC", 
                        "Mean", "RW", "BGR")

datasim4 <- read.table("sim4dat.txt")
for(r in 1:Nsim){
  
  # Data
  Y <- datasim4[((n)*(r-1)+1):(r*n), ]
  Y <- matrix(as.matrix(Y), ncol=ncol(Y), nrow=nrow(Y))
  
  A <- constructModel(Y, p=4, "Basic", gran=c(2500,10),verbose=FALSE, 
                      T1=floor(0.67*nrow(Y)), T2=nrow(Y)-1) 
  
  # Lasso
  resH <- cv.BigVAR(A)
  Sim4MSFE[r, 4] <- mean(resH@OOSMSFE)
  
  # BGR
  A@Structure <- "BGR"
  resBGR <- cv.BigVAR(A)
  Sim4MSFE[r, 10] <- mean(resBGR@OOSMSFE)
  
  # Componentwise
  A@Structure <- "HVARC"
  resH <- cv.BigVAR(A)
  Sim4MSFE[r, 1] <- mean(resH@OOSMSFE)
  
  # Own-other
  A@Structure <- "HVAROO"
  resH <- cv.BigVAR(A)
  Sim4MSFE[r, 2] <- mean(resH@OOSMSFE)
  
  # Elementwise
  A@Structure <- "HVARELEM"
  resH <- cv.BigVAR(A)
  Sim4MSFE[r, 3] <- mean(resH@OOSMSFE)
  
  # Lag-Weighted Lasso
  A@Structure <- "Tapered"
  resH <- cv.BigVAR(A)
  Sim4MSFE[r, 5] <- mean(resH@OOSMSFE)
  
  # VAR(1), Mean and Random Walk
  Sim4MSFE[r, 6] <- mean(resH@AICMSFE)
  Sim4MSFE[r, 7] <- mean(resH@BICMSFE)
  Sim4MSFE[r, 8] <- mean(resH@MeanMSFE)
  Sim4MSFE[r, 9] <- mean(resH@RWMSFE)
}
apply(Sim4MSFE, 2, mean)/k



#####################################################
#### Simulation Scenario 5 : Robustness Scenario ####
#####################################################

#### Create autoregressive matrices ####
set.seed(1234)
p=5 # We ran the simulations for p=1,5,12,25,50
k=10
n=200
Fp = matrix(1,k*p,k*p)
while(max(Mod(eigen(Fp)$values)>1))
{
  Aex <- matrix(0,k,k*p)

  for(i in 1:k)
  {
    if(i<=4)
    {
      bb <- runif(k,-.4,.4)
      Aex[i,1:(k)] <- bb
      Aex[i,(k+1):(2*k)] <- -.9*bb

    }

    if(i>=5 & i<=7)
    {

      bb <- runif(k,-.6,.6)
      Aex[i,1:(k)] <- bb
      Aex[i,(k+1):(2*k)] <- -.8*bb
      Aex[i,(2*k+1):(3*k)] <- .6*bb
      Aex[i,(3*k+1):(4*k)] <- -.5*bb
      Aex[i,(4*k+1):(5*k)] <- .45*bb
    }

  }

  Fp = matrix(0,k*p,k*p)
  Fp[1:k,] = Aex
  Fp[-(1:k),1:(k*(p-1))] = diag(k*(p-1))
  print(max(Mod(eigen(Fp)$values)))

}

#### Figure autoregressive matrices ####
print(SparsityPlotLarge(Fp[1:k,1:(4*k)],4,k)) # Figure 4, panel 5

#### Simulated data ####
Nsim <- 500 # Number of simulations
YY <- vector("list", Nsim)
for(i in 1:Nsim){
  YY[[i]] <-  MultVarSim(k, A1=Fp, p, .01*diag(k), n)
  write(t(YY[[i]]),file=paste0("sim5dat.txt"), ncolumns=ncol(YY[[i]]), sep="\t", append=T)
}

#### Forecast Performance ####
Sim5MSFE <- matrix(NA, ncol=10, nrow=Nsim)
colnames(Sim5MSFE) <- c("HLagC", "HLagOO", "HLagElem", "Lasso", "LagLasso", "LS AIC", "LS BIC", "Mean", "RW", "BGR")

datasim5 <- read.table("sim5dat.txt")
for(r in 1:Nsim){
  
  # Data
  Y <- datasim5[((n)*(r-1)+1):(r*n), ]
  Y <- matrix(as.matrix(Y), ncol=ncol(Y), nrow=nrow(Y))
  
  A <- constructModel(Y, p=5, "Basic", gran=c(1000,10),verbose=FALSE, 
                      T1=floor(0.67*nrow(Y)), T2=nrow(Y)-1,tol=1e-5) 
  
  # Lasso
  resH <- cv.BigVAR(A)
  Sim5MSFE[r, 4] <- mean(resH@OOSMSFE)
  
  # BGR
  A@Structure <- "BGR"
  resBGR <- cv.BigVAR(A)
  Sim5MSFE[r, 10] <- mean(resBGR@OOSMSFE)
  
  # Componentwise
  A@Structure <- "HVARC"
  resH <- cv.BigVAR(A)
  Sim5MSFE[r, 1] <- mean(resH@OOSMSFE)
  
  # Own-other
  A@Structure <- "HVAROO"
  resH <- cv.BigVAR(A)
  Sim5MSFE[r, 2] <- mean(resH@OOSMSFE)
  
  # Elementwise
  A@Structure <- "HVARELEM"
  resH <- cv.BigVAR(A)
  Sim5MSFE[r, 3] <- mean(resH@OOSMSFE)
  
  # Lag-Weighted Lasso
  A@Structure <- "Tapered"
  resH <- cv.BigVAR(A)
  Sim5MSFE[r, 5] <- mean(resH@OOSMSFE)
  
  # VAR(1), Mean and Random Walk
  Sim5MSFE[r, 6] <- mean(resH@AICMSFE)
  Sim5MSFE[r, 7] <- mean(resH@BICMSFE)
  Sim5MSFE[r, 8] <- mean(resH@MeanMSFE)
  Sim5MSFE[r, 9] <- mean(resH@RWMSFE)
}
apply(Sim5MSFE, 2, mean)/k
