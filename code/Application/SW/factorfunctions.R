SFM <- function(Y, r.max = min(ncol(Y), 10), horizon = 1){
  
  ptm1 <- proc.time()
  # PCA
  decomp <- eigen(cov(Y))
  factor <- Y%*%decomp$vector[,1:r.max]
  
  # Number of factors
  PCp <- rep(NA, r.max)
  residmax <- Y - factor[, r.max]
  s2hat <- mean(diag((t(residmax)%*%residmax))/nrow(Y)) 
  for(i.r in 1:r.max){
    resid <- (Y - factor[, i.r])
    PCp[i.r] <- mean(diag((t(resid)%*%resid))/nrow(Y)) + i.r*s2hat*((ncol(Y) + nrow(Y))/(ncol(Y)*nrow(Y)))*log(min(ncol(Y), nrow(Y)))
  }
  rank = which.min(PCp)
  
  if(rank==1){
    f <- matrix(factor[,1], nrow(Y), 1)
  }else{
    f <- factor[, 1:rank]
  }
  
  
  ##################
  # Static version #
  ##################
  
  # Estimates
  Ydata <- Y[-1, ]
  Xdata <- f[-dim(Y)[1], ]
  if(rank==1){
    Xdata <- matrix(Xdata, dim(Y)[1] - 1, 1)
  }
  beta <- solve(t(Xdata)%*%Xdata)%*%t(Xdata)%*%Ydata
  
  # Forecast
  Yhat_static <- matrix(NA, horizon, ncol(Y))
  for(i.h in 1:horizon){
    if(i.h==1){
      Yhat_static[i.h, ] <- Y[nrow(Y),]%*%decomp$vector[,1:rank]%*%beta
    }else{
      Yhat_static[i.h, ] <- Yhat_static[i.h - 1, ]%*%decomp$vector[,1:rank]%*%beta
    }
  }
  
  timeSFM <- proc.time() - ptm1
  out <- list( "rank" = rank, "eigenvector" = decomp$vector, "factor" = factor, "f" = f,
               "Yhat_static" = Yhat_static, "runtime" = timeSFM, "decomp" = decomp)
  
}

DFM <- function(Y, f, rank, horizon = 1, lag.max = 12, Yhat_static, decomp){
  
  ptm1 <- proc.time()
  
  ###################
  # Dynamic Version #
  ###################
  flags <- embed(f, lag.max) 
  
  Ydata <- Y[-(1:lag.max), ]
  Xdata <- flags[-nrow(flags), ]
  
  ### BIC ###
  BIClag <- AIClag <- rep(NA, lag.max)
  for(i.l in 1:lag.max){
    YBIC <- c(Ydata)
    XBIC <- kronecker(diag(1, ncol(Ydata)), Xdata[, (1:(rank*i.l))])
    Xinv <- kronecker(solve(t(diag(1, ncol(Ydata)))%*%diag(1, ncol(Ydata))), 
                      solve(t(Xdata[, (1:(rank*i.l))])%*%Xdata[, (1:(rank*i.l))]))
    betaBIC <- kronecker(solve(t(diag(1, ncol(Ydata)))%*%diag(1, ncol(Ydata)))%*%t(diag(1, ncol(Ydata))), 
                         solve(t(Xdata[, (1:(rank*i.l))])%*%Xdata[, (1:(rank*i.l))])%*%t(Xdata[, (1:(rank*i.l))]))%*%YBIC
    s2BIC <- mean((YBIC - XBIC%*%betaBIC)^2)
    BIClag[i.l] <- nrow(Ydata)%*%log(s2BIC) + (length(betaBIC)/ncol(Ydata))*log(nrow(Ydata))
    AIClag[i.l] <- nrow(Ydata)%*%log(s2BIC) + (length(betaBIC)/ncol(Ydata))*2
  }
  lag.opt.BIC <- which.min(BIClag)
  
  if(lag.opt.BIC==1){
    Yhat_dynamic_BIC = Yhat_static
  }else{
    Yopt <- Ydata
    Xopt <- Xdata[, (1:(rank*lag.opt.BIC))]  
    betaDFM <- solve(t(Xopt)%*%Xopt)%*%t(Xopt)%*%Yopt
    
    Ylast <- Y[nrow(Y) : (nrow(Y) - lag.opt.BIC + 1), ] 
    Yhat_dynamic_BIC <- matrix(NA, horizon, ncol(Y))
    for(i.h in 1:horizon){
      if(i.h==1){
        Yhat_dynamic_BIC[i.h, ] <- c(t(Ylast[1:lag.opt.BIC, ]%*%decomp$vector[,1:rank]))%*%betaDFM
        Ylast <- rbind(Yhat_dynamic_BIC[i.h, ], Ylast)
      }else{
        Yhat_dynamic_BIC[i.h, ] <- c(t(Ylast[1:lag.opt.BIC, ]%*%decomp$vector[,1:rank]))%*%betaDFM
        Ylast <- rbind(Yhat_dynamic_BIC[i.h, ], Ylast)
      }
    }
  }
  
  lag.opt.AIC <- which.min(AIClag)
  
  if(lag.opt.AIC==1){
    Yhat_dynamic_AIC = Yhat_static
  }else{
    
    if(lag.opt.AIC==lag.opt.BIC){
      Yhat_dynamic_AIC = Yhat_dynamic_BIC
    }else{
      Yopt <- Ydata
      Xopt <- Xdata[, (1:(rank*lag.opt.AIC))]  
      betaDFM <- solve(t(Xopt)%*%Xopt)%*%t(Xopt)%*%Yopt
      
      Ylast <- Y[nrow(Y) : (nrow(Y) - lag.opt.AIC + 1), ] # ROW1: Y1t-1 Y2t-1, ROW2: Y1t-2 Y2t-2, ....
      Yhat_dynamic_AIC <- matrix(NA, horizon, ncol(Y))
      for(i.h in 1:horizon){
        if(i.h==1){
          Yhat_dynamic_AIC[i.h, ] <- c(t(Ylast[1:lag.opt.AIC, ]%*%decomp$vector[,1:rank]))%*%betaDFM
          Ylast <- rbind(Yhat_dynamic_AIC[i.h, ], Ylast)
        }else{
          Yhat_dynamic_AIC[i.h, ] <- c(t(Ylast[1:lag.opt.AIC, ]%*%decomp$vector[,1:rank]))%*%betaDFM
          Ylast <- rbind(Yhat_dynamic_AIC[i.h, ], Ylast)
        }
      }
    }
    
  }
  
  runtime <- proc.time() - ptm1
  
  out <- list("lagsBIC" = lag.opt.BIC, "Yhat_dynamic_BIC" = Yhat_dynamic_BIC,
              "lag.opt.AIC" = lag.opt.AIC, "Yhat_dynamic_AIC" = Yhat_dynamic_AIC, "runtime" = runtime)
  
}

FAVAR <- function(Y, r.max = min(ncol(Y) - 1, 10), horizon = 1, lag.max = 12){
  
  ptm1 <- proc.time()
  ranks <- rep(NA, ncol(Y))
  lagsAIC <- rep(NA, ncol(Y))
  lagsBIC <- rep(NA, ncol(Y))
  Yhat_dynamic_BIC <- Yhat_dynamic_AIC <- matrix(NA, horizon, ncol(Y))
  
  
  for(i.h in 1:horizon){
    
    if(i.h==1){
      YlastBIC <- YlastAIC <- Y[nrow(Y) : (nrow(Y) - lag.max + 1), ]
      if(lag.max==1){
        YlastBIC <- matrix(YlastBIC, nrow = 1)
        YlastAIC <- matrix(YlastAIC, nrow = 1)
      }
    }
    
    for(i.k in 1:ncol(Y)){
      
      # PCA
      decomp <- eigen(cov(Y[,-i.k]))
      factor <- Y[, -i.k]%*%decomp$vector[,1:r.max]
      
      # Determine the number of factors
      PCp <- rep(NA, r.max)
      residmax <- Y[,-i.k] - factor[, r.max]
      s2hat <- mean(diag((t(residmax)%*%residmax))/nrow(Y)) 
      for(i.r in 1:r.max){
        resid <- (Y[,-i.k] - factor[, i.r])
        PCp[i.r] <- mean(diag((t(resid)%*%resid))/nrow(Y)) + i.r*s2hat*((ncol(Y[,-i.k]) + nrow(Y[,-i.k]))/(ncol(Y[,-i.k])*nrow(Y[,-i.k])))*log(min(ncol(Y[,-i.k]), nrow(Y[,-i.k])))
      }
      
      # Number of Factors
      ranks[i.k] = which.min(PCp)
      
      if(ranks[i.k]==1){
        f <- matrix(factor[,1], nrow(Y[,-i.k]), 1)
      }else{
        f <- factor[, 1:ranks[i.k]]
      }
      
      # Estimate FAVAR
      data <- cbind(Y[, i.k], f)
      lag.max <- min(lag.max, floor(nrow(Y)/(ranks[i.k])))
      datalags <- embed(data, lag.max + 1) 
      
      BIClag <- AIClag <- rep(NA, lag.max)
      for(i.l in 1:lag.max){
        YBIC <- c(datalags[, 1:ncol(data)]) 
        XBIC <- kronecker(diag(1, ncol(data)), datalags[, -c(1:ncol(data))][, 1:(ncol(data)*i.l)])
        Xinv <- kronecker(solve(t(diag(1, ncol(data)))%*%diag(1, ncol(data))), 
                          solve(t(datalags[, -c(1:ncol(data))][, 1:(ncol(data)*i.l)])%*%datalags[, -c(1:ncol(data))][, 1:(ncol(data)*i.l)]))
        betaBIC <- kronecker(solve(t(diag(1, ncol(data)))%*%diag(1, ncol(data)))%*%t(diag(1, ncol(data))), 
                             solve(t(datalags[, -c(1:ncol(data))][, 1:(ncol(data)*i.l)])%*%datalags[, -c(1:ncol(data))][, 1:(ncol(data)*i.l)])%*%t(datalags[, -c(1:ncol(data))][, 1:(ncol(data)*i.l)]))%*%YBIC
        
        Xeq1 <- datalags[, -c(1:ncol(data))][, 1:(ncol(data)*i.l)]
        betaeq1 <- betaBIC[1:dim(Xeq1)[2]]
        s2BIC <- mean((datalags[, 1] - Xeq1%*%betaeq1)^2) 
        
        BIClag[i.l] <- nrow(datalags)%*%log(s2BIC) + length(betaeq1)*log(nrow(datalags))
        AIClag[i.l] <- nrow(datalags)%*%log(s2BIC) + length(betaeq1)*2
      }

      lag.opt.BIC <- which.min(BIClag)
      lag.opt.AIC <- which.min(AIClag)
      lagsBIC[i.k] <- lag.opt.BIC
      lagsAIC[i.k] <- lag.opt.AIC

      # BIC
      Yopt <- datalags[, 1:ncol(data)]
      Xopt <- datalags[, -c(1:ncol(data))][, (1:(ncol(data)*lag.opt.BIC))]  
      betaFAVARBIC <- solve(t(Xopt)%*%Xopt)%*%t(Xopt)%*%Yopt
      
      coefs <- rep(NA, nrow(betaFAVARBIC))
      coefs[seq(from = 1, by = ncol(data), length = lag.opt.BIC)] <- YlastBIC[1 : (lag.opt.BIC), i.k]
      coefsf <- YlastBIC[1 : (lag.opt.BIC), -i.k]%*%decomp$vector[,1:ranks[i.k]]
      for(i.r in 1:ranks[i.k]){
        coefs[seq(from = 1+ i.r, by = ncol(data), length = lag.opt.BIC)] <- coefsf[, i.r]
      }
      Yhats <- coefs%*%betaFAVARBIC
      Yhat_dynamic_BIC[i.h, i.k] <- Yhats[1]
      
      # AIC
      Yopt <- datalags[, 1:ncol(data)]
      Xopt <- datalags[, -c(1:ncol(data))][, (1:(ncol(data)*lag.opt.AIC))]  
      betaFAVARAIC<- solve(t(Xopt)%*%Xopt)%*%t(Xopt)%*%Yopt
      
      
      coefs <- rep(NA, nrow(betaFAVARAIC))
      coefs[seq(from = 1, by = ncol(data), length = lag.opt.AIC)] <- YlastAIC[1 : (lag.opt.AIC), i.k]
      coefsf <- YlastAIC[1 : (lag.opt.AIC), -i.k]%*%decomp$vector[,1:ranks[i.k]]
      for(i.r in 1:ranks[i.k]){
        coefs[seq(from = 1+ i.r, by = ncol(data), length = lag.opt.AIC)] <- coefsf[, i.r]
      }
      Yhats <- coefs%*%betaFAVARAIC
      Yhat_dynamic_AIC[i.h, i.k] <- Yhats[1]
    
    }  
    
    if(i.h!=horizon){
      YlastBIC <- rbind(Yhat_dynamic_BIC[i.h, ], YlastBIC)
      YlastAIC <- rbind(Yhat_dynamic_AIC[i.h, ], YlastAIC)
    }
    
  }
  
  runtime <- proc.time() - ptm1
  out <- list("ranks" = ranks, "YhatsBIC" = Yhat_dynamic_BIC, "YhatsAIC" = Yhat_dynamic_AIC,
              "lagBIC" = lagsBIC, "lagAIC" = lagsAIC, "runtime" = runtime)
}