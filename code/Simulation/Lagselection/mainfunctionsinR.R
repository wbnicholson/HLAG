# Additional functions outside of BigVAR that are required for lag selection simulations
findLagLength <- function(B,k=60,p=12,thresh=1e-8)
    {
        if(max(abs(B))<thresh){return(0)}
        else{   
        B2 <- max(which(abs(B)>thresh))
            return(floor(B2/k))
        }
    }

 LagMatrix <- function(B,k,p,thresh)
    {
        L <- matrix(0,nrow=k,ncol=k)
        for(i in 1:k)
            {
         for(j in 1:k)       
                {
                 kp <- seq(j,k*p,k)
                 ML1 <- which(abs(B[i,kp])>thresh)                 
                 ML <- ifelse(length(ML1==0),max(ML1),0)
                 L[i,j] <- ML
                    


                }
         }

return(L)


     }   



.EvalLag <- function(ZFull,gamopt,k,p,group,palpha,T1,T2,tol,intercept)
{
    gran2=1	
    gamm <- gamopt
    Y <- ZFull$Y
    Z <- ZFull$Z
    ## T1 <- floor(2 * nrow(Y)/3)
    ## T2 <- nrow(Y)
    MSFE <- c()
    alpha=1/(ncol(Y)+1)	
    b2=array(0,dim=c(k,k*p+1,1))
    MN=FALSE

    ## v <- T1+1
    ## v <- T2
       ##  trainY <- ZFull$Y[1:(v-p-1), ]
       ## trainZ <- ZFull$Z[,1:(v-p-1)]
    PercentageSupport <- c()
    for(i in (T1+1):T2)
    {
        trainY <- ZFull$Y[1:(i-p-1), ]
       trainZ <- ZFull$Z[,1:(i-p-1)]

        
    if (group == "Basic") {
        b2 <- .lassoVARFist(b2, trainZ, trainY,gamm, tol,p,MN,C1=NULL,intercept)
        b2temp <- b2[,2:(k*p+1),1]
        }
                 if (group=="HVARC")
           {
		b2 <- .HVARCAlg(b2,trainY,trainZ,gamm,tol,p,MN,C1=NULL,intercept)
        b2temp <- b2[,2:(k*p+1),1]
             }
        if(group=="HVAROO")
          {
              b2 <- .HVAROOAlg(b2,trainY,trainZ,gamm,tol,p,MN,C1=NULL,intercept)
              b2temp <- b2[,2:(k*p+1),1]
          }
	if(group=="HVARELEM")
	{
	    b2<-.HVARElemAlg(b2,trainY,trainZ,gamm,tol,p,MN,C1=NULL,intercept)
        b2temp <- b2[,2:(k*p+1),1]
    }
    if(group=="Tapered")
        {
        b2 <- .lassoVARTL(b2, trainZ, trainY,gamm, tol,p,MN,palpha,C1=NULL,intercept)
        b2temp <- b2[,2:(k*p+1),1]
            }
 ## browser()
EstZeros <- LagMatrix(b2temp,k,p,1e-3)
    PercentageSupport[i] <- (sum(abs(TrueZeros-EstZeros)))/sum(abs(TrueZeros))
}

return(PercentageSupport)
    }    


LagMatrix<-function(B,k,p,thresh)
    {
        L <- matrix(0,nrow=k,ncol=k)
        for(i in 1:k)
            {
         for(j in 1:k)       
         {
                 kp <- seq(j,k*p,k)
                 ML1 <- which(abs(B[i,kp])>thresh)                 
                 ML <- ifelse(length(ML1==0),max(ML1),0)
                 L[i,j] <- ML
                    


                }
         }

return(L)


}


# Ensures that the created BigVAR object is valid
check.BigVAR <- function(object){
  
  errors <- character()
  
  VARX <- object@VARX
  Y <- object@Data
  
  if(any(is.na(Y))){msg <- c("Remove NA values before running ConstructModel")
  
  errors <- c(errors,msg)
  }      
  if(dim(Y)[2]>dim(Y)[1] & length(VARX)==0){msg <- paste("k is",ncol(Y),"which is greater than T, is Y formatted correctly (k x T)?")}      
  if(object@lagmax<0){msg <- c("Maximal lag order must be at least 0")
  errors <- c(errors,msg)
  }
  if(object@lagmax==0& object@Structure!="Basic"){
    msg <- c("Only Basic VARX-L supports a transfer function")
    errors <- c(errors,msg)
  }
  structures=c("Basic","Lag","SparseLag","OwnOther","SparseOO","HVARC","HVAROO","HVARELEM","Tapered","EFX","BGR")
  cond1=object@Structure%in% structures
  if(cond1==FALSE){
    msg <- paste("struct must be one of",structures)
    errors <- c(errors,msg)      
  }
  if(object@horizon<1){msg <- paste("Forecast Horizon is ",object@horizon, " must be at least 1")
  
  }
  if(object@crossval!="Rolling" & object@crossval!="LOO"){msg <- c("Cross-Validation type must be one of Rolling or LOO")
  errors <- c(errors,msg)
  }
  if(length(object@Granularity)!=2&object@ownlambdas==FALSE){msg("Granularity must have two parameters")
    errors <- c(errors,msg)
    
  }
  if(any(object@Granularity<=0)){
    msg <- c("Granularity parameters must be positive")
    errors <- c(errors,msg)
  }
  structure2 <- c("Basic","Lag","HVARC")
  cond2=object@Structure %in% structure2
  k1=0
  if(length(VARX)!=0){
    
    k1 <- VARX$k
    if(k1>ncol(Y)){msg <- c("k is greater than the number of columns in Y")
    
    errors <- c(errors,msg)
    }        
  }else{k=0}
  m <- ncol(Y)-k1
  
  if(object@tf & object@lagmax!=0){
    msg <- c("p must be 0 if fitting a transfer function")
    errors <- c(errors,msg)
  }
  nseries <- ncol(Y)-ifelse(m<ncol(Y),m,0)    
  if(nseries==1 & cond2==FALSE ){
    msg <- c("Univariate support is only available for Basic VARX-L, Lag Group VARX-L, and Componentwise HVAR")
    
    errors <- c(errors,msg)
    
  }
  if(length(VARX)==0 & object@Structure=="EFX"){
    
    msg <- c("EFX is only supported in the VARX framework")
    
    errors <- c(errors,msg)          
    
  }
  
  if(is.list(VARX) & length(VARX)>0 & !(exists('k',where=VARX) & exists('s',where=VARX)))
  {
    
    msg <- c("VARX Specifications entered incorrectly")
    
    errors <- c(errors,msg)
  }
  
  
  if(object@Structure=="EFX" & !is.null(VARX$contemp)){
    if(VARX$contemp){
      msg <- c("EFX does not support contemporaneous dependence")
      errors <- c(errors,msg)
    }
    
  }
  structs=c("HVARC","HVAROO","HVARELEM")
  if(length(VARX)!=0& object@Structure %in% structs){
    msg <- c("EFX is the only nested model supported in the VARX framework")
    
    errors <- c(errors,msg)
    
  }
  if(object@T1>nrow(Y) | object@T2>nrow(Y) |object@T2<object@T1){
    msg <- c("Training dates exceed series length")
    errors <- c(errors,msg)
    
  }
  
  if(any(object@alpha<0) || any(object@alpha>1)){
    msg <- c("alpha must be between zero and 1")
    errors <- c(errors,msg)
  }
  
  if(length(errors)==0) TRUE else errors
  
}

#' BigVAR Object Class
#'
#' An object class to be used with cv.BigVAR
#' 
#' @slot Data a \eqn{T \times k} multivariate time Series
#' @slot lagmax Maximal lag order for modeled series
#' @slot intercept Indicator as to whether an intercept should be included 
#' @slot Structure Penalty Structure
#' @slot Relaxed Indicator for relaxed VAR
#' @slot Granularity Granularity of Penalty Grid
#' @slot horizon Desired Forecast Horizon
#' @slot crossval Cross-Validation Procedure
#' @slot Minnesota Minnesota Prior Indicator
#' @slot verbose Indicator for Verbose output
#' @slot dates dates extracted from an xts object 
#' @slot ic Indicator for including AIC and BIC benchmarks
#' @slot VARX VARX Model Specifications
#' @slot T1 Index of time series in which to start cross validation
#' @slot T2  Index of times series in which to start forecast evaluation
#' @slot ONESE Indicator for "One Standard Error Heuristic"
#' @slot ownlambdas Indicator for user-supplied lambdas
#' @slot tf Indicator for transfer function
#' @slot alpha Grid of candidate alpha values (applies only to Sparse VARX-L models)
#' @slot recursive Indicator as to whether recursive multi-step forecasts are used (applies only to multiple horizon VAR models)
#' @slot constvec vector indicating variables to shrink toward a random walk instead of toward zero (valid only if Minnesota is \code{TRUE})
#' @slot tol optimization tolerance
#' @slot lagselect lag selection indicator
#' @details To construct an object of class BigVAR, use the function \code{\link{constructModel}}
#' @seealso \code{\link{constructModel}}
#' @export
setClass(
  Class="BigVAR",
  representation(
    Data="matrix",
    lagmax="numeric",
    Structure="character",
    Relaxed="logical",
    Granularity="numeric",
    intercept="logical",
    Minnesota="logical",  
    horizon="numeric",
    verbose="logical",  
    crossval="character",
    ic="logical",
    VARX="list",
    T1="numeric",
    T2="numeric",
    ONESE="logical",
    ownlambdas="logical",
    tf="logical",
    alpha="numeric",
    recursive="logical",
    dates="character",
    constvec="numeric",
    tol="numeric",
    lagselect="logical"
  ),validity=check.BigVAR
)


#' Construct an object of class BigVAR
#' @param Y \eqn{T \times k} multivariate time series or Y \eqn{T \times (k+m)} endogenous and exogenous series, respectively 
#' @param p Predetermined maximal lag order (for modeled series)
#' @param struct The choice of penalty structure (see details).
#' @param gran vector of penalty parameter specifications.
#' @param intercept True or False: option to fit an intercept
#' @param RVAR True or False: option to refit based upon the support selected using the Relaxed-VAR procedure
#' @param h Desired forecast horizon
#' @param cv Cross-validation approach, either "Rolling" for rolling cross-validation or "LOO" for leave-one-out cross-validation.
#' @param MN Minnesota Prior Indicator
#' @param verbose Verbose output while estimating
#' @param IC True or False: whether to include AIC and BIC benchmarks
#' @param VARX List containing VARX model specifications. 
#' @param T1 Index of time series in which to start cross validation
#' @param T2  Index of times series in which to start forecast evaluation
#' @param ONESE True or False: whether to use the "One Standard Error Heuristic"
#' @param ownlambdas True or False: Indicator for user-supplied penalty parameters
#' @param alpha grid of candidate parameters for the alpha in the Sparse Lag and Sparse Own/Other VARX-L 
#' @param recursive True or False: Indicator as to whether iterative multi-step predictions are desired in the VAR context if the forecast horizon is greater than 1
#' @param C vector of coefficients to shrink toward a random walk (if \code{MN} is \code{TRUE})
#' @param tol optimization tolerance (default 1e-4)
#' @param dates optional vector of dates corresponding to \eqn{Y}
#'
#' 
#'  @details The choices for "struct" are as follows
#' \itemize{
#' \item{  "Basic" (Basic VARX-L)}
#' \item{  "Lag" (Lag Group VARX-L)} 
#' \item{  "SparseLag" (Lag Sparse Group VARX-L)} 
#' \item{  "OwnOther" (Own/Other Group VARX-L) }
#' \item{  "SparseOO" (Own/Other Sparse Group VARX-L) }
#' \item{  "EFX" (Endogenous First VARX-L)}
#' \item{  "HVARC" (Componentwise HVAR) }
#' \item{  "HVAROO" (Own/Other HVAR) }
#' \item{  "HVARELEM" (Elementwise HVAR)}
#' \item{  "Tapered" (Lag weighted Lasso VAR)}
#' \item{  "BGR" (Bayesian Ridge Regression (cf. Banbura et al))}
#' }
#'
#' The first number in the vector "gran" specifies how deep to construct the penalty grid and the second specifies how many penalty parameters to use  If ownlambas is set to TRUE, gran should contain the user-supplied penalty parameters.
#' 
#' VARX specifications consist of a list with entry k denoting the series that are to be modeled and entry s to denote the maximal lag order for exogenous series.
#'
#' The argument alpha is ignored unless the structure choice is "SparseLag" or "Lag."  By default "alpha" is set to \code{NULL} and will be initialized as 1/(k+1) in \code{cv.BigVAR} and \code{BigVAR.est}.  Any user supplied values must be between 0 and 1.  

#' @note The specifications "Basic", "Lag," "SparseLag," "SparseOO," and "OwnOther" can accommodate both VAR and VARX models.  EFX only applies to VARX models.  "HVARC," "HVAROO," "HVARELEM," and "Tapered" can only be used with VAR models.
#'
#' @seealso \code{\link{cv.BigVAR}},\code{\link{BigVAR.est}}
#' 
#' @references  William B Nicholson, Jacob Bien, and David S Matteson. "High Dimensional Forecasting via Interpretable Vector Autoregression." arXiv preprint arXiv:1412.5250, 2016.
#' William B Nicholson, David S. Matteson, and Jacob Bien (2015), "VARX-L Structured regularization for large vector autoregressions with exogenous variables," arXiv preprint arXiv:1508.07497, 2016.
#' William B Nicholson, David S. Matteson, and Jacob Bien (2016), "BigVAR: Dimension Reduction Reduction Methods for Multivariate Time Series," \url{http://www.wbnicholson.com/BigVAR.pdf}.
#'
#' Bańbura, Marta, Domenico Giannone, and Lucrezia Reichlin. "Large Bayesian vector auto regressions." Journal of Applied Econometrics 25.1 (2010): 71-92.
#' @examples
#' # VARX Example
#' # Create a Basic VARX-L with k=2, m=1, s=2, p=4
#' VARX=list()
#' VARX$k=2 # indicates that the first two series are modeled
#' VARX$s=2 # sets 2 as the maximal lag order for exogenous series
#' data(Y)
#' T1=floor(nrow(Y)/3)
#' T2=floor(2*nrow(Y)/3)
#' Model1=constructModel(Y,p=4,struct="Basic",gran=c(50,10),verbose=FALSE,VARX=VARX,T1=T1,T2=T2)
#' @export
constructModel <- function(Y,p,struct,gran,RVAR=FALSE,h=1,cv="Rolling",MN=FALSE,verbose=TRUE,IC=TRUE,VARX=list(),T1=floor(nrow(Y)/3),T2=floor(2*nrow(Y)/3),ONESE=FALSE,ownlambdas=FALSE,alpha=as.double(NULL),recursive=FALSE,C=as.double(NULL),dates=as.character(NULL),intercept=TRUE,tol=1e-4,lagselect=FALSE)
{
  if(any(is.na(Y))){stop("Remove NA values before running constructModel")}      
  if(dim(Y)[2]>dim(Y)[1] & length(VARX)==0){warning("k is greater than T, is Y formatted correctly (k x T)?")}      
  if(p<0){stop("Maximal lag order must be at least 0")}
  if(p==0& struct!="Basic"){stop("Only Basic VARX-L supports a transfer function")}
  oldnames <- c("None","Diag","SparseDiag")
  if(struct%in%oldnames) stop("Naming Convention for these structures has changed. Use Basic, OwnOther, and SparseOO.")
  structures=c("Basic","Lag","SparseLag","OwnOther","SparseOO","HVARC","HVAROO","HVARELEM","Tapered","EFX","BGR")
  cond1=struct %in% structures
  if(!cond1){stop(cat("struct must be one of",structures))}
  if(h<1){stop("Forecast Horizon must be at least 1")}
  if(cv!="Rolling" & cv!="LOO"){stop("Cross-Validation type must be one of Rolling or LOO")}
  if(length(gran)!=2&ownlambdas==FALSE){stop("Granularity must have two parameters")}
  if(any(gran<=0)){stop("Granularity parameters must be positive")}
  if(tol<0 | tol>1e-1){stop("Tolerance must be positive")}
  structure2 <- c("Basic","Lag","HVARC")
  cond2=struct %in% structure2
  ## k <- ncol(Y)
  
  if(length(VARX)!=0){
    
    k <- VARX$k
    if(k>ncol(Y)){stop("k is greater than the number of columns in Y")}
  }else{k=ncol(Y)}
  m <- ncol(Y)-k
  nseries <- ncol(Y)-ifelse(m<ncol(Y),m,0)
  if(p==0){tf=TRUE
  }else{
    tf=FALSE
  }
  if(nseries==1 & cond2==FALSE ){stop("Univariate support is only available for Lasso, Lag Group, and Componentwise HVAR")}
  if(length(VARX)==0 & struct=="EFX"){stop("EFX is only supported in the VARX framework")}
  if(struct=="EFX" & !is.null(VARX$contemp)){
    if(VARX$contemp){
      stop("EFX does not support contemporaneous dependence")}
  }
  structs=c("HVARC","HVAROO","HVARELEM")
  if(length(VARX)!=0& struct %in% structs){stop("EFX is the only nested model supported in the VARX framework")}
  if(T1>nrow(Y) | T2>nrow(Y) |T2<T1){stop("Training dates exceed series length")}
  
  if(is.list(VARX) & length(VARX)>0 & !(exists('k',where=VARX) & exists('s',where=VARX)))
  {
    
    stop("VARX Specifications entered incorrectly")
    
  }
  
  if(!is.null(alpha)){
    if(any(alpha<0) || any(alpha>1)){stop("alpha must be between 0 and 1")}
  }
  if(length(C)!=0){
    
    if(length(C)!=k){stop("C must have length k")}
    if(!all(C%in%c(0,1))){stop("Values of C must be either 0 or 1")}
    
  }else{
    
    C <- rep(1,k)
    
  }
  ## if("xts"%in%class(Y)){
  ##     ind <- as.character(index(Y))
  ##     Y <- as.matrix(Y)
  
  if(length(dates)!=0){
    
    ind <- dates
    
  }else{
    ind <- as.character(NULL)
  }
  # Can't have a class named C
  (BV1 <- new(
    "BigVAR",
    Data=Y,
    lagmax=p,
    Structure=struct,
    Relaxed=RVAR,
    Granularity=gran,
    Minnesota=MN,
    verbose=verbose,
    horizon=h,
    crossval=cv,
    ic=IC,
    VARX=VARX,
    T1=T1,
    T2=T2,
    ONESE=ONESE,
    ownlambdas=ownlambdas,
    tf=tf,
    alpha=alpha,
    recursive=recursive,
    dates=ind,
    constvec=C,
    intercept=intercept,
    tol=tol,
    lagselect=lagselect
  ))
  
  return(BV1)
  
}




# show-default method to show an object when its name is printed in the console.
#' Default show method for an object of class BigVAR
#'
#' @param object \code{BigVAR} object created from \code{ConstructModel}
#' @return Displays the following information about the BigVAR object:
#' \itemize{
#' \item{Prints the first 5 rows of \code{Y}}
#' \item{ Penalty Structure}
#' \item{ Relaxed Least Squares Indicator}
#' \item{Maximum lag order} 
#' \item{ VARX Specifications (if applicable)}
#' \item{Start, end of cross validation period}
#' }
#' @seealso \code{\link{constructModel}} 
#' @name show.BigVAR
#' @aliases show,BigVAR-method
#' @docType methods
#' @rdname show-methods
#' @export
setMethod("show","BigVAR",
          function(object)
          {
            
            
            T1P <- ifelse(length(object@dates)!=0,object@dates[object@T1],object@T1)
            
            T2P <- ifelse(length(object@dates)!=0,object@dates[object@T2],object@T2)
            
            nrowShow <- min(5,nrow(object@Data))
            cat("*** BIGVAR MODEL *** \n")
            cat("Data (First 5 Observations):\n")
            print(formatC(object@Data[1:nrowShow,]),quote=FALSE)
            cat("Structure\n") ;print(object@Structure)
            cat("Forecast Horizon \n") ;print(object@horizon)
            cat("Relaxed VAR \n") ;print(object@Relaxed)
            cat("Minnesota Prior \n") ;print(object@Minnesota)
            cat("Maximum Lag Order \n") ;print(object@lagmax)
            if(length(object@VARX)!=0){
              cat("VARX Specs \n") ;print(object@VARX)}
            cat("Start of Cross Validation Period \n") ;print(T1P)
            cat("End of Cross Validation Period \n") ;print(T2P)
            
          }
)

#' Plot a BigVAR object
#' 
#' @param x BigVAR object created from \code{ConstructModel}
#' @param y NULL
#' @param ... additional plot arguments
#' @return NA, side effect is graph
#' @details Uses plot.zoo to plot each individual series of \code{Y} on a single plot
#' @name plot.BigVAR
#' @import methods
#' @seealso \code{\link{constructModel}}
#' @aliases plot,BigVAR-method
#' @aliases plot-methods
#' @docType methods
#' @method plot
#' @rdname plot.BigVAR-methods
#' @export
#' @importFrom zoo plot.zoo
#' @importFrom zoo as.zoo
#' @importFrom zoo zoo
#' @importFrom zoo as.yearqtr
#' @importFrom zoo index
#' @importFrom graphics abline
#' @importFrom graphics legend
setMethod(f="plot",signature="BigVAR",
          definition= function(x,y=NULL,...)
          {
            
            T1P <- ifelse(length(x@dates)!=0,x@dates[x@T1],x@T1)
            
            T2P <- ifelse(length(x@dates)!=0,x@dates[x@T2],x@T2)
            
            g=ncol(x@Data)
            names <- ifelse(rep(!is.null(colnames(x@Data)),ncol(x@Data)),colnames(x@Data),as.character(1:g))
            if(length(x@dates)!=0){
              
              dates <- as.yearqtr(x@dates)
            }else{
              
              dates <- 1:nrow(x@Data)
            }
            
            Yzoo <- zoo(as.matrix(x@Data),order.by=dates)
            plot.zoo(Yzoo,plot.type="single",col=1:g)
            legend('topright',names,lty=1,col=1:g)
            
            abline(v=index(Yzoo[as.yearqtr(T1P)]))
            abline(v=index(Yzoo[as.yearqtr(T2P)]))
            
          }
)

#' Cross Validation for BigVAR
#' 
#' @usage cv.BigVAR(object)
#' @param object BigVAR object created from \code{ConstructModel}
#' @details The main function of the BigVAR package. Performs cross validation to select penalty parameters over a training sample (as the minimizer of in-sample MSFE), then evaluates them over a test set.  Compares against sample mean, random walk, AIC, and BIC benchmarks.  Creates an object of class \code{BigVAR.results}
#' @return An object of class \code{BigVAR.results}.
#' @seealso \code{\link{constructModel}}, \code{\link{BigVAR.results}},\code{\link{BigVAR.est}} 
#' @name cv.BigVAR
#' @aliases cv.BigVAR,BigVAR-method
#' @docType methods
#' @rdname cv.BigVAR-methods
#' @examples
#' data(Y)
#' # Fit a Basic VARX-L with rolling cross validation 
#' Model1=constructModel(Y,p=4,struct="Basic",gran=c(50,10))
#' results=cv.BigVAR(Model1)
#' @export
setGeneric(
  
  name="cv.BigVAR",
  def=function(object)
  {
    standardGeneric("cv.BigVAR")
    
  }
  
  
)
# Cross-validation and evaluation function
setMethod(
  
  f="cv.BigVAR",
  signature="BigVAR",
  definition=function(object){
    p <- object@lagmax
    s1 <- 0
    Y <- object@Data
    k <- ncol(Y)
    alpha <- object@alpha
    RVAR <- object@Relaxed
    group <- object@Structure
    cvtype <- object@crossval
    intercept=object@intercept
    recursive <- object@recursive
    VARX <- object@VARX
    tol=object@tol
    lagselect=object@lagselect
    if(length(alpha)==0){
      
      if(length(VARX)>0){    
        alpha <- 1/(VARX$k+1)
        
      }else{
        
        alpha <- 1/(k+1)
        
      }
    }
    
    C <- object@constvec
    
    if(length(alpha)>1 & group%in%c("SparseLag","SparseOO"))
    {
      dual <- TRUE
      
    }else{
      
      dual <- FALSE
    }
    
    MN <- object@Minnesota
    h <- object@horizon
    jj <- 0
    
    if(!"matrix"%in%class(Y)){Y=matrix(Y,ncol=1)}
    
    if(object@crossval=="Rolling"){
      T1 <- object@T1
      
    }else{
      
      T1 <- p+2    
      
    }
    T2 <- object@T2
    s <- ifelse(length(object@VARX)!=0,object@VARX$s,0)
    
    if(object@ownlambdas){
      gamm <- object@Granularity
      gran2 <- length(gamm)
    }     
    
    ONESE <- object@ONESE
    
    
    
    # Adjust T1, T2 by maximum lag order to account for initialization
    if(object@crossval=="Rolling"){
      T1 <- T1-max(p,s)
      T2 <- T2-max(p,s)
    }
    if(!object@ownlambdas){
      
      gran2 <- object@Granularity[2]
      gran1 <- object@Granularity[1]
      
    }
    
    
    
    # Constructing lag matrix in VARX setting
    if(length(VARX)!=0){
      
      VARX <- TRUE
      k1 <- object@VARX$k
      s <- object@VARX$s
      
      if(!is.null(object@VARX$contemp)){
        
        contemp <- TRUE
        s1 <- 1
        
      }else{
        
        contemp <- FALSE
        s1 <- 0
      }
      
      m <- k-k1
      Y1 <-matrix(Y[,1:k1],ncol=k1)
      X <- matrix(Y[,(ncol(Y)-m+1):ncol(Y)],ncol=m)
      
      if(!object@tf){
        trainZ <- VARXCons(Y1,X,k1,p,m,s,contemp=contemp)
        
      }else{
        
        trainZ <- VARXCons(matrix(0,ncol=1,nrow=nrow(X)),matrix(X,ncol=m),k=0,p=0,m=m,s=s,contemp=contemp,oos=FALSE)
        
      }
      
      trainZ <- trainZ[2:nrow(trainZ),]
      
      trainY <- matrix(Y[(max(c(p,s))+1):nrow(Y),1:k1],ncol=k1)
      
      if(group=="Lag"|group=="SparseLag"){
        
        # Lag based groupings
        jj <- groupfunVARXLG(p,k,k1,s+s1)
        
        # initialize activeset as empty
        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                         gran2*length(alpha))
        
      }
      if(group=="OwnOther"|group=="SparseOO"){
        # Own other based groupings
        jj <- diaggroupfunVARXLG(p,k,k1,s+s1)
        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                         gran2)
        
      }
      if(object@ownlambdas==FALSE){
        if(dual){
          # Constructs penalty grid if both alpha and lambda are selected
          gamm <- .LambdaGridXDual(gran1, gran2, jj, trainY, trainZ,group,p,k1,s,m,k,MN,alpha,C,intercept,tol)
          
        }else{
          
          # Penalty parameter grid for just lambda
          gamm <- .LambdaGridX(gran1, gran2, jj, as.matrix(trainY[1:T2,]), trainZ[,1:T2],group,p,k1,s+s1,m,k,MN,alpha,C,intercept,tol)
        }
      }
      
      # Coefficient matrices
      ## if(!dual){
      
      ##     beta <- array(0,dim=c(k1,k1*p+(k-k1)*(s+s1)+1,gran2))
      
      ## }else{
      
      beta <- array(0,dim=c(k1,k1*p+(k-k1)*(s+s1)+1,gran2*length(alpha)))
      
      ## }
      # Groupings in accordance with C++ indexing standards
      if (group == "Lag") {
        jj <- groupfunVARX(p,k,k1,s+s1)
        jjcomp <- groupfunVARXcomp(p,k,k1,s+s1)
        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                         gran2)
      }
      if (group == "SparseLag") {
        
        jj <- groupfunVARX(p, k,k1,s+s1)
        q1a <- list()
        # Initializing warm start vectors for the power method
        for (i in 1:(p+s+s1)) {
          q1a[[i]] <- matrix(runif(length(jj[[i]]), -1, 1), ncol = 1)
        }
        
        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                         gran2*length(alpha))
        
      }
      if (group == "OwnOther") {
        kk <- diaggroupfunVARX(p,k,k1,s+s1)
        activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                         gran2)
      }
      if (group == "SparseOO") {
        kk <- diaggroupfunVARX(p,k,k1,s+s1)
        activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                         gran2*length(alpha))
        q1a <- list()
        
        for (i in 1:(2*p+s+s1)) {
          q1a[[i]] <- matrix(runif(length(jj[[i]]), -1, 1), ncol = 1)
        }
        
        
      }
    }else{
      # VAR estimation
      contemp <- FALSE
      if(group=="Lag"|group=="SparseLag")
      {
        jj <- .groupfun(p,k)
      }else{
        jj <- .lfunction3(p,k)
      }
      Z1 <- VARXCons(Y,matrix(0,nrow=nrow(Y)),k,p,0,0) 
      
      trainZ <- Z1[2:nrow(Z1),]   
      
      trainY <- matrix(Y[(p+1):nrow(Y),],ncol=k)          
      
      GY <- matrix(trainY[1:T2,],ncol=k)
      
      # We only use training period data to construct the penalty grid
      GZ <- trainZ[,1:T2]
      
      if(object@ownlambdas==FALSE){
        
        
        # Find starting values for penalty grid
        
        if(dual){
          
          gamm <- .LambdaGridEDual(gran1, gran2, jj, GY, GZ,group,p,k,MN,alpha,C,intercept,tol)
          
        }else{
          
          ## gamm <- .LambdaGridE(gran1, gran2, jj, GY, GZ,group,p,k,MN,alpha,C,intercept,tol)
          if(group!="BGR"){
            gamm <- .LambdaGridE(gran1, gran2, jj, GY, GZ,group,p,k,MN,alpha,C,intercept,tol)
            
          }else{
            gamm <- seq(1,5,by=.025)
            gamm <- gamm*sqrt(k*p)
            
          }
          
        }
        
      }
      VARX <- FALSE
      k1 <- k
      s <- 0   
      if (group == "Lag") {
        jj <- .groupfuncpp(p, k)
        jjcomp <- .groupfuncomp(p,k)
        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                         gran2)
        k1=k
      }
      if (group == "SparseLag") {
        jj <- .groupfuncpp(p, k)
        q1a <- list()
        for (i in 1:p) {
          q1a[[i]] <- matrix(runif(k, -1, 1), ncol = 1)
        }
        
        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                         gran2*length(alpha))
        
      }
      if (group == "OwnOther") {
        kk <- .lfunction3cpp(p, k)
        activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                         gran2)
      }
      if (group == "SparseOO") {
        kk <- .lfunction3cpp(p, k)
        jjcomp <- .lfunctioncomp(p,k)
        jj <- .lfunction3(p,k)
        activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                         gran2*length(alpha))
        q1a <- list()
        for (i in 1:(2*p)) {
          q1a[[i]] <- matrix(runif(length(kk[[i]]), -1, 1), ncol = 1)
        }
        
        
      }
      
      beta <- array(0,dim=c(k,k*p+1,gran2*length(alpha)))
    }           
    h <- object@horizon
    verbose <- object@verbose         
    ZFull <- list()
    
    if(!is.matrix(trainZ)){trainZ <- matrix(trainZ,ncol=1)}
    
    if(!is.matrix(trainY)){trainY <- matrix(trainY,ncol=1)}
    
    ZFull$Z <- trainZ
    ZFull$Y <- trainY
    T <- nrow(trainY)
    
    if(object@ownlambdas){gamm <- object@Granularity}    
    
    # Constructs grid for lag weighted lasso
    if(group=="Tapered")
    {
      
      palpha <- seq(0,1,length=10)
      palpha <- rev(palpha)
      gran2 <- length(gamm)*length(palpha)
      beta <- array(0,dim=c(k,k*p+1,gran2))
      
    }
    
    if(class(ZFull$Y)!="matrix" ){
      ZFull$Y <- matrix(ZFull$Y,ncol=1)
    }
    
    if(!dual){
      MSFE <- matrix(0, nrow = T2 - T1+1, ncol = gran2)
    }else{
      
      nalpha <- length(alpha)
      MSFE <- matrix(0, nrow = T2 - T1+1, ncol = gran2*nalpha)
      
    }
    if(verbose){
      pb <- txtProgressBar(min = T1, max = T2, style = 3)
      cat("Cross-Validation Stage:",group)}
    YT <- Y[1:T2,]
    
    # Start of penalty parameter selection     
    
    for (v in (T1-h+1):T2) {
      
      if(cvtype=="Rolling")
        
      {
        
        if(v+h-1>T2){
          break
        }
        
        if(h>1 & !recursive){
          
          trainY <- ZFull$Y[(h):(v-1), ]
          
          trainZ <- ZFull$Z[, 1:(v-h)]
          
        }else{
          
          trainY <- ZFull$Y[1:(v-1), ]
          
          
          trainZ <- ZFull$Z[, 1:(v-1)]         
        }
      }else{
        
        if(VARX)
          
        {
          
          YT2 <- YT[-v,]
          
          Y1 <- matrix(YT2[,1:k1],ncol=k1)
          
          X <- matrix(YT2[,(ncol(YT2)-m+1):ncol(YT2)],ncol=m)
          
          trainZ <- VARXCons(Y1,X,k1,p,m,s,contemp=contemp)
          
          trainZ <- trainZ[2:nrow(trainZ),]
          
          trainY <- matrix(YT2[(max(c(p,s))+1):nrow(YT2),1:k1],ncol=k1)
          
        }else{
          
          
          YT2 <- YT[-v,]
          
          Z1 <- VARXCons(YT2,matrix(0,nrow=nrow(YT2)),k,p,0,0) 
          
          trainZ <- Z1[2:nrow(Z1),]        
          
          trainY <- matrix(YT2[(p+1):nrow(YT2),],ncol=k)                                  
          
        }
        
      }
      
      
      if (group == "Basic") {
        
        if(VARX){
          
          beta <- .lassoVARFistX(beta, trainZ, trainY,gamm, tol,p,MN,k,k1,s+s1,m,C,intercept)
          
        }else{
          
          beta <- .lassoVARFist(beta, trainZ, trainY,gamm, tol,p,MN,C,intercept)
        }
        
      }
      
      if (group == "Lag") {
        
        GG <- .GroupLassoVAR1(beta,jj,jjcomp,trainY,trainZ,gamm,activeset,tol,p,MN,k,k1,s+s1,C,intercept)
        
        beta <- GG$beta
        
        activeset <- GG$active
        
      }
      
      if (group == "SparseLag") {
        
        if(VARX){
          
          if(!dual){
            
            GG <- .SparseGroupLassoVARX(beta, jj, trainY, trainZ, 
                                        gamm, alpha, INIactive = activeset, tol, q1a,p,MN,k,s+s1,k1,C,intercept)
            
          }else{
            
            GG <- .SparseGroupLassoVARXDual(beta, jj, trainY, trainZ, 
                                            gamm, alpha, INIactive = activeset, tol, q1a,p,MN,k,s+s1,k1,C,intercept)
            
          }
        }else{
          
          if(!dual){
            GG <- .SparseGroupLassoVAR(beta, jj, trainY, trainZ, 
                                       gamm, alpha, INIactive = activeset, tol, q1a,p,MN,C,intercept)
            
          }else{
            GG <- .SparseGroupLassoVARDual(beta, jj, trainY, trainZ, 
                                           gamm, alpha, INIactive = activeset, tol, q1a,p,MN,C,intercept)
            
            
          }
        }
        
        beta <- GG$beta
        
        activeset <- GG$active
        
        q1a <- GG$q1
        
      }
      
      if (group == "OwnOther") {
        
        if(VARX){
          
          GG <- .GroupLassoOOX(beta, kk, trainY, trainZ, gamm, 
                               activeset, tol,p,MN,k,k1,s+s1,C,intercept)
          
        }else{
          
          
          GG <- .GroupLassoOO(beta, kk, trainY, trainZ, gamm, 
                              activeset, tol,p,MN,C,intercept)
        }
        
        beta <- GG$beta
        
        activeset <- GG$active
        
      }
      
      if (group == "SparseOO") {
        if(VARX){
          
          
          GG <- .SparseGroupLassoVAROOX(beta, kk, trainY, trainZ, 
                                        gamm, alpha, INIactive = activeset, tol,p,MN,k1,s+s1,k,dual,C,intercept)
          
        }else{
          
          GG <- .SparseGroupLassoVAROO(beta, kk, trainY, trainZ, 
                                       gamm, alpha, INIactive = activeset, tol,q1a,p,MN,dual,C,intercept)
          
          q1a <- GG$q1
          
        }
        
        beta <- GG$beta
        
        activeset <- GG$active
        
      }
      
      if(group=="Tapered")
      {
        
        beta <- .lassoVARTL(beta,trainZ,trainY,gamm,tol,p,MN,palpha,C,intercept)    
      }
      
      if(group=="EFX")
      {
        
        beta <- .EFVARX(beta,trainY,trainZ,gamm,tol,MN,k1,s,m,p,C,intercept)
        
      }
      
      if(group=="HVARC")
      {
        beta <- .HVARCAlg(beta,trainY,trainZ,gamm,tol,p,MN,C,intercept)
        
      }
      
      if(group=="HVAROO")
      {
        beta <- .HVAROOAlg(beta,trainY,trainZ,gamm,tol,p,MN,C,intercept)
      }
      
      if(group=="HVARELEM")
      {
        
        beta <- .HVARElemAlg(beta,trainY,trainZ,gamm,tol,p,MN,C,intercept)
        
      }
      
      eZ <- c(1,ZFull$Z[,v])
      
      
      if(group!="BGR"){
        if(!dual)
        {
          
          
          # Calculate h-step MSFE for each penalty parameter
          for (ii in 1:gran2) {
            
            if (RVAR)
            {
              # Relaxed Least Squares (intercept ignored)
              beta[,,ii] <- RelaxedLS(cbind(t(trainZ),trainY),beta[,,ii],k,p,k1,s+s1)
            }
            if (MN){
              # shrink toward vector random walk
              
              pred <- beta[,2:dim(beta)[2],ii] %*% eZ[2:length(eZ)]
              
              if(h>1 & recursive){
                pred <- matrix(pred,nrow=1)
                pred <- predictMS(pred,trainY,h-1,beta[,2:dim(beta)[2],ii],p,MN)
              }
              
              MSFE[v - (T1 - h), ii] <- norm2(ZFull$Y[v+h-1,1:k1] - pred)^2
              # Subtract one from diagonal for warm start purposes
              
              diag(beta[,2:(k1+1),ii]) <- diag(beta[,2:(k1+1),ii])-C
              
            }else{
              
              if(object@crossval=="Rolling"){
                
                pred <- beta[,,ii] %*% eZ
                
                if(h>1 & recursive){
                  pred <- matrix(pred,nrow=1)
                  pred <- predictMS(pred,trainY,h-1,beta[,,ii],p)
                }
                MSFE[v - (T1 -h), ii] <- norm2(ZFull$Y[v+h-1,1:k1] - pred)^2
                
                
              }else{
                if(VARX){          
                  
                  eZ <- VARXCons(matrix(Y[(v-p):(v),1:k1],ncol=k1),Y[(v-p):(v),(ncol(Y)-m+1):(ncol(Y))],k1,p
                                 ,m,s,contemp=contemp)
                  
                  pred <- beta[,,ii] %*% eZ
                  
                  
                }else{
                  
                  eZ<- VARXCons(Y[(v-p):(v),1:k1],matrix(0,nrow=length((v-p):v)),k1,p,0,0)
                  
                  pred <- beta[,,ii] %*% eZ
                  
                }
                MSFE[v - (T1 - h), ii] <- norm2(Y[v+h-1,1:k1] - pred)^2     
                
                
              }                    
              
            }
            
          }
        }else{
          
          # If alpha and lambda are jointly fit, calculate MSFE for each alpha, lambda combination
          for (ii in 1:gran2) {
            for(jj in 1:length(alpha)){
              if (RVAR) {
                
                # Relaxed Least Squares (intercept ignored)
                beta[,,(ii-1)*nalpha+jj] <- RelaxedLS(cbind(t(trainZ),trainY),beta[,,(ii-1)*nalpha+jj],k,p,k1,s+s1)
              }
              
              if(MN){
                
                pred <- beta[,2:dim(beta)[2],(ii-1)*nalpha+jj] %*% eZ[2:length(eZ)]                            
                if(h>1 & recursive){
                  pred <- matrix(pred,nrow=1)
                  pred <- predictMS(pred,trainY,h-1,beta[,2:dim(beta)[2],(ii-1)*nalpha+jj],p,MN)
                }
                
                
                MSFE[v - (T1 - h), (ii-1)*nalpha+jj] <- norm2(ZFull$Y[v+h-1,1:k1] - beta[,2:dim(beta)[2],(ii-1)*nalpha+jj] %*% eZ[2:length(eZ)])^2
                
              }else{
                pred <- beta[,,(ii-1)*nalpha+jj] %*% eZ
                
                if(h>1 & recursive){
                  pred <- matrix(pred,nrow=1)
                  pred <- predictMS(pred,trainY,h-1,beta[,,(ii-1)*nalpha+jj],p)
                }
                
                if(object@crossval=="Rolling"){
                  MSFE[v - (T1 - h), (ii-1)*nalpha+jj] <- norm2(ZFull$Y[v+h-1,1:k1] - beta[,,(ii-1)*nalpha+jj] %*% eZ)^2
                  err <- norm2(ZFull$Y[v,1:k1] - beta[,,ii] %*% eZ)^2
                  
                }else{
                  if(VARX){          
                    eZ<- VARXCons(matrix(Y[(v-p):(v),1:k1],ncol=k1),Y[(v-p):(v),(ncol(Y)-m+1):(ncol(Y))],k1,p,m,s,contemp=contemp)
                  }else{
                    eZ<- VARXCons(Y[(v-p):(v),1:k1],matrix(0,nrow=length((v-p):v)),k1,p,0,0)
                  }
                  MSFE[v - (T1 - h), (ii-1)*nalpha+jj] <- norm2(Y[v,1:k1] - beta[,,(ii-1)*nalpha+jj] %*% eZ)^2
                  
                  
                  
                  
                }
              }
              
              
            }
          }
        }
        
      }else{
        for (ii in 1:ncol(MSFE)) {
          MSFE[v - (T1 - h), ii] <- norm2(Y[v+h-1,1:k1] - beta[,,ii])^2
        }
      }
      
      
      if(verbose){
        setTxtProgressBar(pb, v)}
    }
    
    # Sort out indexing for 2-d gridsearch
    
    if(group=="Tapered")
      
    {
      
      indopt <- which.min(colMeans(MSFE))
      
      
      if(indopt<length(gamm))
      {
        
        alphaind <- 1
        
        alphaopt <- palpha[1]
        
      }
      
      else{
        
        alphaopt <- palpha[floor(indopt/length(gamm))]
        
        alphaind <- floor(indopt/length(gamm))
        
      }
      
      if(alphaind==1)
      {
        
        gamopt <- gamm[indopt]
        
      }else if(indopt %% length(gamm)==0)
      {
        
        gamopt <- gamm[length(gamm)]
        
      }else{
        gamind <- indopt-length(gamm)*alphaind
        gamopt <- gamm[gamind]
        
      }
      
      palpha<-alphaopt
      optind <- indopt
    }
    
    
    # one standard error correction     
    if(ONESE & !dual){
      MSFE2 <- MSFE 
      G2 <- colMeans(na.omit(MSFE2))
      G3 <- sd(na.omit(MSFE2))/sqrt(nrow(na.omit(MSFE2)))
      optind <- min(which(G2<(min(G2)+G3)))
      gamopt <- gamm[optind]
    }else{
      
      if(group!="Tapered" & !dual){
        # in rare cases in which MSFE is equal, the smaller penalty parameter is chosen.
        # This prevents extremely sparse solutions
        optind <- max(which(colMeans(na.omit(MSFE))==min(colMeans(na.omit(MSFE)))))
        gamopt <- gamm[optind]
      }else if(dual){
        
        if(!ONESE){
          
          optind <- max(which(colMeans(na.omit(MSFE))==min(colMeans(na.omit(MSFE)))))
          
          inds <- findind(optind,gamm[,1],alpha)                                   
        }else{
          
          G2 <- colMeans(na.omit(MSFE))
          
          G3 <- sd(na.omit(MSFE))/sqrt(nrow(na.omit(MSFE)))
          
          optind <- min(which(G2<(min(G2)+G3)))
          inds <- findind(optind,gamm[,1],alpha)
          
        }
        gamopt <- gamm[inds[1],inds[2]]
        gamm <- gamm[,inds[2]]
        alphaopt <- alpha[inds[2]]
        optind <- inds
        ## print(paste("alphaopt",alphaopt))
      }
    }
    if(!dual){
      alphaopt <- alpha
    }
    
    if(VARX){
      
      # Out of sample forecast evaluation: VARX
      OOSEval <- .BigVAREVALX(ZFull,gamopt,k,p,group,h,MN,verbose,RVAR,palpha,T2,T,k1,s,m,contemp,alphaopt,C,intercept,tol)
      
    }else{
      # Out of sample evaluation for VAR
      if(!lagselect){
        OOSEval <- .BigVAREVAL(ZFull,gamopt,k,p,group,h,MN,verbose,RVAR,palpha,T2,T,alphaopt,recursive,C,intercept,tol)
        PercentageSupport <- 0
      }else{
        PercentageSupport <- .EvalLag(ZFull,gamopt,k,p,group,palpha,T2,T,tol,intercept)
        OOSEval <- .BigVAREVAL(ZFull,gamopt,k,p,group,h,MN,verbose,RVAR,palpha,T2,T2+2,alphaopt,recursive,C,intercept,tol)
      }
    }
    MSFEOOSAgg <- na.omit(OOSEval$MSFE)
    betaPred <- OOSEval$betaPred
    preds <- OOSEval$predictions
    Y <- object@Data # preserve Y for BigVAR.results object
    if(VARX){
      
      # Construct lag matrix for OOS predictions
      # need to be careful with OOS predictions if contemporaneous exogenous covariates are included
      if(contemp){OOS <- FALSE}else{OOS <- TRUE}
      # special case of transfer function                                      
      if(!object@tf){
        Zvals <- VARXCons(matrix(Y[,1:k1],ncol=k1),matrix(Y[,(ncol(Y)-m+1):ncol(Y)],ncol=m),k1,p,m,s,oos=OOS,contemp=contemp)
      }else{
        
        Zvals <- VARXCons(matrix(0,ncol=1,nrow=nrow(Y)),matrix(Y[,(ncol(Y)-m+1):ncol(Y)],ncol=m),0,0,m,s,oos=FALSE,contemp=contemp)
      }
      
    }else{
      m <- 0;s <- 0
      Zvals <- VARXCons(matrix(Y[,1:k1],ncol=k1),matrix(0,nrow=nrow(Y)),k1,p,m,s,oos=TRUE)
    }
    
    Zvals <- matrix(Zvals[,ncol(Zvals)],ncol=1)
    if(ncol(Y)==1| k1==1){betaPred <- matrix(betaPred,nrow=1)}
    lagmatrix <- rbind(rep(1,ncol(ZFull$Z)),ZFull$Z)
    
    fitted <- t(betaPred%*%lagmatrix)
    #Residuals
    resids <- ((ZFull$Y)-fitted)
    
    ## lagmatrix <- ZFull$Z
    
    MSFEOOS<-mean(na.omit(MSFEOOSAgg))
    
    seoos <- sd(na.omit(MSFEOOSAgg))/sqrt(length(na.omit(MSFEOOSAgg)))
    
    if(!VARX){k1 <- k}
    
    
    # naive benchmarks     
    meanbench <- .evalMean(ZFull$Y[,1:k1],T2,T,h=h)
    RWbench <- .evalRW(ZFull$Y[,1:k1],T2,T,h=h)
    
    
    if(object@ic==FALSE|object@tf){
      
      AICbench <- list()
      AICbench$Mean <- as.double(NA)
      AICbench$SD <- as.double(NA)
      AICbench$preds <- as.matrix(NA)
      AICbench$pvec <- as.double(NA)
      AICbench$svec <- as.double(NA)                          
      
      BICbench <- list()
      BICbench$Mean <- as.double(NA)
      BICbench$SD <- as.double(NA)                          
      BICbench$preds <- as.matrix(NA)                          
      
      BICbench$pvec <- as.double(NA)
      BICbench$svec <- as.double(NA)                          
      
      # Information Criterion Benchmarks    
      
    }else{
      
      if(!VARX){
        
        X <- matrix(0,nrow=nrow(Y),ncol=k)
        
        AICbench1 <- VARXForecastEval(matrix(ZFull$Y,ncol=k),X,p,0,T2,T,"AIC",h)
        
        AICbench <- list()
        
        AICbench$Mean <- mean(AICbench1$MSFE)
        
        AICbench$SD <- sd(AICbench1$MSFE)/sqrt(length(AICbench1$MSFE))
        AICbench$preds <- AICbench1$pred
        AICbench$pvec <- AICbench1$p
        AICbench$svec <- AICbench1$s
        
        BICbench1 <- VARXForecastEval(matrix(ZFull$Y,ncol=k),X,p,0,T2,T,"BIC",h)
        
        BICbench <- list()
        
        BICbench$Mean <- mean(BICbench1$MSFE)
        
        BICbench$SD <- sd(BICbench1$MSFE)/sqrt(length(BICbench1$MSFE))
        BICbench$preds <- BICbench1$pred
        BICbench$pvec <- BICbench1$p
        BICbench$svec <- BICbench1$s
        
      }else{
        
        offset <- max(c(p,s))
        
        X <- matrix(Y[(offset+1):nrow(Y),(k1+1):ncol(Y)],ncol=m)
        
        AICbench1 <- VARXForecastEval(matrix(ZFull$Y[,1:k1],ncol=k1),as.matrix(X),p,s,T2,T,"AIC",h=h)
        
        AICbench <- list()
        
        AICbench$Mean <- mean(AICbench1$MSFE)
        
        AICbench$SD <- sd(AICbench1$MSFE)/sqrt(length(AICbench1$MSFE))
        AICbench$preds <- AICbench1$pred
        AICbench$pvec <- AICbench1$p
        AICbench$svec <- AICbench1$s
        
        BICbench1 <- VARXForecastEval(matrix(ZFull$Y[,1:k1],ncol=k1),X,p,s,T2,T,"BIC",h=h)
        
        BICbench <- list()
        
        BICbench$Mean <- mean(BICbench1$MSFE)
        
        BICbench$SD <- sd(BICbench1$MSFE)/sqrt(length(BICbench1$MSFE))  
        BICbench$preds <- BICbench1$pred
        BICbench$pvec <- BICbench1$p
        BICbench$svec <- BICbench1$s
        
      }
      
    }
    
    if(!VARX){contemp=FALSE}
    
    if(VARX & contemp){
      VARXL <- list(k=k1,s=s,contemp=contemp)
    }else if(VARX & !contemp){
      
      VARXL <- list(k=k1,s=s,contemp=FALSE)
      
    }else{
      VARXL <- list()
    }
    # Create a new BigVAR.Results Object
    results <- new("BigVAR.results",InSampMSFE=colMeans(MSFE),InSampSD=apply(MSFE,2,sd)/sqrt(nrow(MSFE)),LambdaGrid=gamm,index=optind,OptimalLambda=gamopt,OOSMSFE=MSFEOOSAgg,seoosmsfe=seoos,MeanMSFE=meanbench$Mean,AICMSFE=AICbench$Mean,AICpvec=AICbench$pvec,AICsvec=AICbench$svec,AICPreds=AICbench$preds,BICpvec=BICbench$pvec,BICsvec=BICbench$svec,BICPreds=BICbench$preds,RWMSFE=RWbench$Mean,RWPreds=RWbench$preds,MeanSD=meanbench$SD,MeanPreds=meanbench$preds,AICSD=AICbench$SD,BICMSFE=BICbench$Mean,BICSD=BICbench$SD,RWSD=RWbench$SD,Data=object@Data,lagmax=object@lagmax,Structure=object@Structure,Minnesota=object@Minnesota,Relaxed=object@Relaxed,Granularity=object@Granularity,horizon=object@horizon,betaPred=betaPred,Zvals=Zvals,resids=resids,VARXI=VARX,VARX=VARXL,preds=preds,T1=T1,T2=T2,dual=dual,alpha=alphaopt,crossval=object@crossval,ownlambdas=object@ownlambdas,tf=object@tf,recursive=recursive,constvec=C,intercept=intercept,tol=tol,fitted=fitted,lagmatrix=lagmatrix,SR=PercentageSupport)
    
    return(results)
  }
)


#' BigVAR Estimation
#' @description
#' Fit a BigVAR object with a structured penalty (VARX-L or HVAR).
#' @usage BigVAR.est(object)
#' @param object BigVAR object created from \code{ConstructModel}
#' @details Fits HVAR or VARX-L model on a BigVAR object.  Does not perform cross-validation.  This method allows the user to construct their own penalty parameter selection procedure.
#' @return An array of \eqn{k \times kp \times n} or \eqn{k\times kp+ms \times n} coefficient matrices; one for each of the n values of lambda.
#' @seealso \code{\link{constructModel}}, \code{\link{BigVAR.results}},\code{\link{cv.BigVAR}} 
#' @name BigVAR.est
#' @aliases BigVAR.est,BigVAR-method
#' @docType methods
#' @rdname BigVAR.est-methods
#' @examples
#' data(Y)
#' Y=Y[1:100,]
#' #construct a Basic VAR-L
#' Model1=constructModel(Y,p=4,struct="Basic",gran=c(50,10))
#' BigVAR.est(Model1)
#' @export
setGeneric(
  name="BigVAR.est",
  def=function(object)
  {
    standardGeneric("BigVAR.est")
    
  }
)
setMethod(
  f="BigVAR.est",
  signature="BigVAR",
  definition=function(object){
    p=object@lagmax
    tol=object@tol
    if(object@ownlambdas==TRUE){
      gamm=object@Granularity
      gran2 <- length(gamm)
      
    }      
    
    C <- object@constvec
    group <- object@Structure
    Y <- object@Data
    k <- ncol(Y)
    VARX <- object@VARX
    intercept <- object@intercept
    if(length(object@alpha)==0){
      
      if(length(VARX)>0){    
        alpha <- 1/(VARX$k+1)
        
      }else{
        
        alpha <- 1/(k+1)
        
      }
    }else{
      
      alpha <- object@alpha
      
    }
    
    
    if(length(alpha)>1 & group%in%c("SparseLag","SparseOO"))
    {
      dual <- TRUE
      
    }else{
      
      dual <- FALSE
    }
    
    RVAR <- object@Relaxed
    group <- object@Structure
    MN <- object@Minnesota
    jj <- 0
    if(!"matrix"%in%class(Y)){Y=matrix(Y,ncol=1)}
    s <- ifelse(length(object@VARX)!=0,object@VARX$s,0)
    T <- nrow(Y)-max(p,s)
    VARX <- object@VARX
    
    if(!object@ownlambdas){
      gran2 <- object@Granularity[2]
      gran1 <- object@Granularity[1]}
    
    if(length(VARX)!=0){
      
      VARX <- TRUE
      k1 <- object@VARX$k
      s <- object@VARX$s
      
      if(!is.null(object@VARX$contemp)){
        
        contemp <- TRUE
        s1 <- 1
        
      }else{
        
        contemp <- FALSE
        s1 <- 0
      }
      
      m <- k-k1
      Y1 <- matrix(Y[,1:k1],ncol=k1)
      X <- matrix(Y[,(ncol(Y)-m+1):ncol(Y)],ncol=m)
      
      if(!object@tf){
        
        trainZ <- VARXCons(Y1,X,k1,p,m,s,contemp=contemp)
        
      }else{
        
        trainZ <- VARXCons(matrix(0,ncol=1,nrow=nrow(X)),matrix(X,ncol=m),k=0,p=0,m=m,s=s,contemp=contemp,oos=FALSE)
        
      }
      
      
      trainZ <- trainZ[2:nrow(trainZ),]
      
      trainY <- matrix(Y[(max(c(p,s))+1):nrow(Y),1:k1],ncol=k1)
      
      
      
      if(group=="Lag"|group=="SparseLag"){
        
        
        jj=groupfunVARXLG(p,k,k1,s)
        
        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                         gran2*length(alpha))
        
        
      }
      
      if(group=="OwnOther"|group=="SparseOO"){
        
        jj=diaggroupfunVARXLG(p,k,k1,s)
        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                         gran2)
        
        
      }
      
      
      if(group=="BGR"){
        Grid <- seq(1,5,by=.025)
        grid <- Grid*sqrt(k*p)
        MSFE <- matrix(0, nrow = T2 - T1+1, ncol = length(grid))
      }
      
      
      if(!object@ownlambdas){
        
        
        if(dual){
          # Constructs penalty grid if both alpha and lambda are selected
          gamm <- .LambdaGridXDual(gran1, gran2, jj, trainY, trainZ,group,p,k1,s,m,k,MN,alpha,C,intercept,tol)
          
        }else{
          
          # Penalty parameter grid for just lambda
          ## gamm <- .LambdaGridX(gran1, gran2, jj, as.matrix(trainY), trainZ,group,p,k1,s+s1,m,k,MN,alpha,C,intercept,tol)
          if(group!="BGR"){
            gamm <- .LambdaGridX(gran1, gran2, jj, as.matrix(trainY[1:T2,]), trainZ[,1:T2],group,p,k1,s+s1,m,k,MN,alpha,C,intercept,tol)
          }else{
            gamm <- seq(1,5,by=.025)
            gamm <- gamm*sqrt(k*p)
            
          }
        }
      }else{
        
        gamm <- object@Granularity
        
      }
      
      # Initial Coefficient Matrix         
      beta=array(0,dim=c(k1,k1*p+(k-k1)*s+1,gran2*length(alpha)))
      
      # Initialize groups, active sets, power method calculations, etc
      if (group == "Lag") {
        
        jj=groupfunVARX(p,k,k1,s+s1)
        
        jjcomp <- groupfunVARXcomp(p,k,k1,s+s1)
        
        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                         gran2)
        
      }
      
      if (group == "SparseLag") {
        
        jj <- groupfunVARX(p, k,k1,s+s1)
        
        q1a <- list()
        
        
        for (i in 1:(p+s+s1)) {
          
          q1a[[i]] <- matrix(runif(length(jj[[i]]), -1, 1), ncol = 1)
          
        }
        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                         gran2*length(alpha))
        
      }
      if (group == "OwnOther") {
        
        kk <- diaggroupfunVARX(p,k,k1,s+s1)
        
        activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                         gran2)
        
      }
      if (group == "SparseOO") {
        
        
        kk <- diaggroupfunVARX(p,k,k1,s+s1)
        
        activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                         gran2*length(alpha))
        
        q1a <- list()
        
        for (i in 1:(2*p+s+s1)) {
          
          q1a[[i]] <- matrix(runif(length(jj[[i]]), -1, 1), ncol = 1)
          
        }
        
      }
    }else{
      
      s=p
      
      if (group=="Lag"|group=="SparseLag")
        
      {
        
        jj=.groupfun(p,k)
        
      }else{
        
        jj <- .lfunction3(p,k)
        
      }
      
      
      Z1 <- VARXCons(Y,matrix(0,nrow=nrow(Y)),k,p,0,0) 
      
      
      trainZ <- Z1[2:nrow(Z1),]   
      
      
      trainY <- matrix(Y[(p+1):nrow(Y),],ncol=k)          
      
      
      if(!object@ownlambdas){
        
        if(dual){
          
          gamm <- .LambdaGridEDual(gran1, gran2, jj, trainY, trainZ,group,p,k,MN,alpha,C,intercept,tol)
          
          
        }else{
          gamm <- .LambdaGridE(gran1, gran2, jj, trainY, trainZ,group,p,k,MN,alpha,C,intercept,tol)
        }
      }else{
        
        gamm <- object@Granularity
        
      }
      
      VARX <- FALSE
      k1 <- k
      s <- 0   
      
      if (group == "Lag") {
        
        jj <- .groupfuncpp(p, k)
        
        jjcomp <- .groupfuncomp(p,k)
        
        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                         gran2)
        
        k1 <- k
        
      }
      
      if (group == "SparseLag") {
        
        jj <- .groupfuncpp(p, k)
        
        q1a <- list()
        
        for (i in 1:p) {
          
          q1a[[i]] <- matrix(runif(k, -1, 1), ncol = 1)
        }
        
        
        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                         gran2)
        
        
      }
      
      if (group == "OwnOther") {
        
        kk <- .lfunction3cpp(p, k)
        
        activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                         gran2)
        
      }
      
      if (group == "SparseOO") {
        
        kk <- .lfunction3cpp(p, k)
        jjcomp <- .lfunctioncomp(p,k)
        jj <- .lfunction3(p,k)
        activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                         gran2)
        q1a <- list()
        
        for (i in 1:(2*p)) {
          
          q1a[[i]] <- matrix(runif(length(kk[[i]]), -1, 1), ncol = 1)
        }
        
      }
      beta <- array(0,dim=c(k,k*p+1,gran2*length(alpha)))
      
    }           
    
    
    
    if(group=="Tapered")
      
    {
      palpha <- seq(0,1,length=10)
      palpha <- rev(palpha)
      gran2 <- length(gamm)*length(palpha)
      beta <- array(0,dim=c(k,k*p+1,gran2))
    }
    
    
    if (group == "BGR") {
      
      trainZ <- rbind(1,trainZ)
      beta <- BGRGridSearch(trainY,trainZ,p,gamm,as.numeric(MN))
    }
    
    if (group == "Basic") {
      
      if(VARX){
        beta <- .lassoVARFistX(beta, trainZ, trainY,gamm, tol,p,MN,k,k1,s,m,C,intercept)}
      else{
        beta <- .lassoVARFist(beta, trainZ, trainY,gamm, tol,p,MN,C,intercept)
      }
    }
    
    if (group == "Lag") {
      
      GG <- .GroupLassoVAR1(beta,jj,jjcomp,trainY,trainZ,gamm,activeset,tol,p,MN,k,k1,s,C,intercept)
      beta <- GG$beta
      activeset <- GG$active
    }
    
    if (group == "SparseLag") {
      
      if(VARX){
        
        if(!dual){
          
          GG <- .SparseGroupLassoVARX(beta, jj, trainY, trainZ, 
                                      gamm, alpha, INIactive = activeset, tol, q1a,p,MN,k,s+s1,k1,C,intercept)
          
        }else{
          
          GG <- .SparseGroupLassoVARXDual(beta, jj, trainY, trainZ, 
                                          gamm, alpha, INIactive = activeset, tol, q1a,p,MN,k,s+s1,k1,C,intercept)
          
        }
      }else{
        
        if(!dual){
          GG <- .SparseGroupLassoVAR(beta, jj, trainY, trainZ, 
                                     gamm, alpha, INIactive = activeset, tol, q1a,p,MN,C,intercept)
          
        }else{
          GG <- .SparseGroupLassoVARDual(beta, jj, trainY, trainZ, 
                                         gamm, alpha, INIactive = activeset, tol, q1a,p,MN,C,intercept)
          
          
        }
      }
      
      beta <- GG$beta
      
      activeset <- GG$active
      
      q1a <- GG$q1
      
    }
    
    if (group == "OwnOther") {
      if(VARX){
        GG <- .GroupLassoOOX(beta, kk, trainY, trainZ, gamm, 
                             activeset, tol,p,MN,k,k1,s,C,intercept)
      }else{
        
        GG <- .GroupLassoOO(beta, kk, trainY, trainZ, gamm, 
                            activeset, tol,p,MN,C,intercept)
      }
      beta <- GG$beta
      activeset <- GG$active
    }
    
    if (group == "SparseOO") {
      if(VARX){
        
        
        GG <- .SparseGroupLassoVAROOX(beta, kk, trainY, trainZ, 
                                      gamm, alpha, INIactive = activeset, tol,p,MN,k1,s+s1,k,dual,C,intercept)
        
      }else{
        
        GG <- .SparseGroupLassoVAROO(beta, kk, trainY, trainZ, 
                                     gamm, alpha, INIactive = activeset, tol,q1a,p,MN,dual,C,intercept)
        
        q1a <- GG$q1
        
      }
      
      beta <- GG$beta
      
      activeset <- GG$active
      
    }
    
    if(group=="Tapered")
    {
      
      beta <- .lassoVARTL(beta,trainZ,trainY,gamm,tol,p,MN,palpha,C,intercept)
      
    }
    
    if(group=="EFX")
    {
      
      beta <- .EFVARX(beta,trainY,trainZ,gamm,tol,MN,k1,s,m,p,C,intercept)
      
    }
    
    if(group=="HVARC")
    {
      
      beta <- .HVARCAlg(beta,trainY,trainZ,gamm,tol,p,MN,C,intercept)
      
    }
    if(group=="HVAROO")
    {
      
      beta <- .HVAROOAlg(beta,trainY,trainZ,gamm,tol,p,MN,C,intercept)
      
    }
    
    if(group=="HVARELEM")
    {
      
      beta <- .HVARElemAlg(beta,trainY,trainZ,gamm,tol,p,MN,C,intercept)
    }      
    
    return(list(B=beta,lambdas=gamm))
    
  }
  
)


## New object class: bigVAR results, inherits class bigVAR, prints results from cv.bigVAR

#' BigVAR.results
#' This class contains the results from cv.BigVAR.
#'
#' It inherits the class BigVAR, but contains substantially more information. 
#' 
#' @field InSampMSFE In-sample MSFE from optimal value of lambda
#' @field LambdaGrid Grid of candidate lambda values
#' @field index Rank of optimal lambda value 
#' @field OptimalLambda Value of lambda that minimizes MSFE
#' @field OOSMSFE Average Out of sample MSFE of BigVAR model with optimal lambda
#' @field seoosfmsfe Standard error of out of sample MSFE of BigVAR model with optimal lambda
#' @field MeanMSFE Average out of sample MSFE of unconditional mean forecast
#' @field MeanSD Standard error of out of sample MSFE of unconditional mean forecast
#' @field MeanPreds predictions from conditional mean model
#' @field RWMSFE Average out of sample MSFE of random walk forecast
#' @field RWPreds Predictions from random walk model
#' @field RWSD Standard error of out of sample MSFE of random walk forecast
#' @field AICMSFE Average out of sample MSFE of AIC forecast
#' @field AICSD Standard error of out of sample MSFE of AIC forecast
#' @ield AICPreds Predictions from AIC VAR/VARX model
#' @field AICpvec Lag orders selected from AIC VAR model
#' @field AICpvec Lag orders selected from AIC VARX model
#' @field BICMSFE Average out of sample MSFE of BIC forecast
#' @field BICSD Standard error of out of sample MSFE of BIC forecast
#' @field BICPreds Predictions from BIC VAR/VARX model
#' @field BICpvec Lag orders selected from BIC VAR model
#' @field BICpvec Lag orders selected from BIC VARX model
#' @field betaPred The final estimated \eqn{k\times kp+ms+1} coefficient matrix, to be used for prediction
#' @field Zvals The final lagged values of \code{Y}, to be used for prediction
#' @field fitted fitted values obtained from betaPred
#' @field resids residuals obtained from betaPred
#' @field Data a \eqn{T \times k} or \eqn{T\times k + m} multivariate time Series
#' @field lagmax Maximal lag order
#' @field Structure Penalty structure
#' @field Relaxed Indicator for relaxed VAR
#' @field Granularity Granularity of penalty grid
#' @field horizon Desired forecast horizon
#' @field crossval Cross-Validation procedure
#' @field alpha additional penalty parameter for Sparse Lag Group or Sparse Own/Other methods. Will contain either the heuristic choice of \eqn{1/(k+1)} or the value selected by cross validation if the argument \code{dual} is set to \code{TRUE}
#' @field VARXI VARX Indicator 
#' @field Minnesota Minnesota Prior Indicator
#' @field verbose  verbose indicator
#' @field dual indicator as to whether dual cross validation was conducted
#' @field contemp indicator if contemporaneous exogenous predictors are used
#' @field lagmatrix matrix of lagged values used to compute residuals (of which Zvals is the final column)

#'
#' @note One can also access any object of class BigVAR from BigVAR.results
#' @name BigVAR.results 
#' @rdname BigVAR.results
#' @aliases BigVAR.results-class
#' @exportClass BigVAR.results
#' @author Will Nicholson
#' @export
setClass("BigVAR.results",
         representation(InSampMSFE="numeric",InSampSD="numeric",LambdaGrid="numeric",index="numeric",OptimalLambda="numeric",OOSMSFE="numeric",seoosmsfe="numeric",MeanMSFE="numeric",AICMSFE="numeric",AICPreds="matrix",BICMSFE="numeric",BICpvec="numeric",BICsvec="numeric",AICpvec="numeric",AICsvec="numeric",BICSD="numeric",BICPreds="matrix",RWMSFE="numeric",RWPreds="matrix",MeanSD="numeric",MeanPreds="matrix",AICSD="numeric",RWSD="numeric",betaPred="matrix",Zvals="matrix",VARXI="logical",resids="matrix",preds="matrix",dual="logical",contemp="logical",fitted="matrix",lagmatrix="matrix",SR='numeric'),
         contains="BigVAR"
)


#' Plot an object of class BigVAR.results
#' 
#' @param x BigVAR.results object created from \code{cv.BigVAR}
#' @param y NULL
#' @param ... additional arguments
#' @details Plots the in sample MSFE of all values of lambda with the optimal value highlighted.
#' @name plot
#' @import methods
#' @aliases plot,BigVAR.results-method
#' @aliases plot-methods
#' @docType methods
#' @method plot
#' @rdname BigVAR.results-plot-methods
#' @importFrom graphics abline
#' @export
setMethod(f="plot",signature="BigVAR.results",
          definition= function(x,y=NULL,...)
          {
            
            plot(x@LambdaGrid,x@InSampMSFE,type="o",xlab="Value of Lambda",ylab="MSFE",log="x")
            
            abline(v=x@OptimalLambda,col="green")
            
          }
)

#' Default show method for an object of class BigVAR.results
#' 
#' @param object BigVAR.results object created from \code{cv.BigVAR}
#' @details prints forecast results and additional diagnostic information as well as comparisons with mean, random walk, and AIC, and BIC benchmarks
#' @seealso \code{\link{cv.BigVAR}},\code{\link{BigVAR.results}} 
#' @name show
#' @aliases show,BigVAR.results-method
#' @docType methods
#' @method show
#' @rdname show-methods-BigVAR.results
#' @export
setMethod("show","BigVAR.results",
          function(object)
          {
            cat("*** BIGVAR MODEL Results *** \n")
            cat("Structure\n") ;print(object@Structure)
            if(object@Relaxed==TRUE){
              cat("Relaxed VAR \n") ;print(object@Relaxed)}
            
            cat("Forecast Horizon \n") ;print(object@horizon)
            cat("Minnesota VAR\n") ;print(object@Minnesota)
            
            if(object@VARXI){
              cat("VARX Specs \n") ;print(object@VARX)}
            cat("Maximum Lag Order \n") ;print(object@lagmax)
            cat("Optimal Lambda \n"); print(round(object@OptimalLambda,digits=4))
            if(object@dual){
              
              cat("Optimal Alpha \n"); print(round(object@alpha,digits=2))
              
            }        
            cat("Grid Depth \n") ;print(object@Granularity[1])
            cat("Index of Optimal Lambda \n");print(object@index)
            cat("In-Sample MSFE\n");print(round(min(object@InSampMSFE),digits=3))
            cat("BigVAR Out of Sample MSFE\n");print(round(mean(object@OOSMSFE),digits=3))
            cat("*** Benchmark Results *** \n")
            cat("Conditional Mean Out of Sample MSFE\n");print(round(object@MeanMSFE,digits=3))
            cat("AIC Out of Sample MSFE\n");print(round(object@AICMSFE,digits=3))
            cat("BIC Out of Sample MSFE\n");print(round(object@BICMSFE,digits=3))
            cat("RW Out of Sample MSFE\n");print(round(object@RWMSFE,digits=3))
          }
)



#' Forecast using a BigVAR.results object
#' 
#' @usage predict(object,...)
#' @param object BigVAR.results object from \code{cv.BigVAR}
#' @param ... additional arguments affecting the predictions produced (e.g. \code{n.ahead})
#' @details Provides \code{n.ahead} step forecasts using the model produced by cv.BigVAR. 
#' @seealso \code{\link{cv.BigVAR}} 
#' @name predict
#' @aliases predict,BigVAR.results-method
#' @docType methods
#' @method predict
#' @rdname predict-methods-BigVAR.results
#' @examples
#' data(Y)
#' Y=Y[1:100,]
#' Model1=constructModel(Y,p=4,struct="Basic",gran=c(50,10),verbose=FALSE)
#' results=cv.BigVAR(Model1)
#' predict(results,n.ahead=1)
#' @export
setMethod("predict","BigVAR.results",
          function(object,n.ahead,newxreg=NULL,...)
          {
            # MN option removes intercept
            MN <- object@Minnesota
            eZ <- object@Zvals
            betaPred <- object@betaPred
            Y <- object@Data
            k <-object@VARX$k
            m <- ncol(object@Data)-k
            p <- object@lagmax
            s <- object@VARX$s
            VARX <- object@VARXI
            contemp <- object@VARX$contemp
            s1 <- 0
            if(VARX){
              if(contemp){
                s1 <- 1
              }
            }else{
              s1=0
            }
            fcst <- betaPred%*%eZ
            
            if(n.ahead==1)
            {
              return(fcst)
            }else{
              if(!VARX){
                # iterative multistep forecasts
                fcst <- predictMS(matrix(fcst,nrow=1),Y[(nrow(Y)-p+1):nrow(Y),],n.ahead-1,betaPred,p,MN)
                
              }else{
                
                # experimental VARX multistep forecasts
                if(is.null(newxreg))
                {
                  stop("Need new data for multi-step VARX forecasts.  Re-run with new data in newxreg")
                  
                }else{
                  if(nrow(newxreg)<n.ahead-1){
                    stop(paste("Need at least ",n.ahead-1,"rows of new data"))
                  }
                  C <- max(p,s)
                  
                  fcst <- predictMSX(matrix(fcst,nrow=1),as.matrix(Y[(nrow(Y)-C+1):nrow(Y),1:(k)]),n.ahead-1,betaPred,p,newxreg,matrix(Y[(nrow(Y)-C+1):nrow(Y),(ncol(Y)-m+1):ncol(Y)],ncol=m),m,s,1,MN)
                }
              }
              
            }
            
            return(fcst)
          }
)

#' Sparsity Plot of a BigVAR.results object 
#'
#' @param object BigVAR.results object
#' @return NA, side effect is graph
#' @details Uses \code{levelplot} from the \code{lattice} package to plot the magnitude of each coefficient in the last coefficient estimated by \code{cv.BigVAR}.
#' @name SparsityPlot.BigVAR.results
#' @aliases SparsityPlot.BigVAR.results,BigVAR.results-method
#' @seealso \code{\link{cv.BigVAR}}, \code{\link{BigVAR.results}}
#' @docType methods
#' @rdname SparsityPlot.BigVAR.results-methods
#' @examples
#' data(Y)
#' Y <- Y[1:100,]
#' Model1 <- constructModel(Y,p=4,struct="Basic",gran=c(50,10),verbose=FALSE)
#' SparsityPlot.BigVAR.results(cv.BigVAR(Model1))
#' @export
#' @importFrom lattice levelplot
#' @importFrom lattice panel.abline
#' @importFrom lattice panel.levelplot
#' @importFrom grDevices colorRampPalette

setGeneric(
  name="SparsityPlot.BigVAR.results",
  def=function(object)
  {
    standardGeneric("SparsityPlot.BigVAR.results")
  }
)
setMethod(
  f="SparsityPlot.BigVAR.results",
  signature="BigVAR.results",
  definition=function(object){
    B <- object@betaPred
    if(nrow(B)==1){
      B <- matrix(B[,2:ncol(B)],nrow=1)
    }else{
      B <- B[,2:ncol(B)]
    }
    k <- nrow(B)
    p <- object@lagmax
    s1 <- 0
    if(length(object@VARX!=0)){
      s <- object@VARX$s
      m <- ncol(object@Data)-object@VARX$k
      contemp <- object@VARX$contemp
      if(!is.null(contemp)){
        if(contemp){
          
          s1 <- 1
        }            
      }else{
        s1 <- 0
        
      }
    }
    else{
      m <- 0;s <- 0
    }
    
    s <- s+s1
    
    text <- c()
    
    for (i in 1:p) {
      
      text1 <- as.expression(bquote(bold(Phi)^(.(i))))
      
      text <- append(text, text1)
      
    }
    
    if(m>0){
      
      for (i in (p+1):(p+s+1)) {
        
        text1 <- as.expression(bquote(bold(beta)^(.(i-p-s1))))
        
        text <- append(text, text1)
        
      }
      
    }
    
    f <- function(m) t(m)[, nrow(m):1]
    
    rgb.palette <- colorRampPalette(c("white", "blue" ),space = "Lab")
    
    at <- seq(k/2 + 0.5, p * (k)+ 0.5, by = k)
    
    if(m>0){
      
      at2 <- seq(p*k+m/2+.5,p*k+s*m+.5,by=m)
      
    }else{
      at2=c()
      
    }
    
    at <- c(at,at2)
    
    se2 = seq(1.75, by = k, length = k)
    
    L2 <- levelplot(as.matrix(f(abs(B))), col.regions = rgb.palette, colorkey = NULL, 
                    xlab = NULL, ylab = NULL, main = list(label = "Sparsity Pattern Generated by BigVAR", 
                                                          cex = 1), panel = function(...) {
                                                            panel.levelplot(...)
                                                            panel.abline(a = NULL, b = 1, h = seq(1.5, m*s+p* k + 
                                                                                                    0.5, by = 1), v = seq(1.5, by = 1, length = p*k+m*s),lwd=.5)
                                                            bl1 <- seq(k + 0.5,p*k + 0.5, by = k)
                                                            b23 <- seq(p*k + 0.5, p * 
                                                                         k + 0.5+s*m, by = m)
                                                            b1 <- c(bl1,b23)
                                                            panel.abline(a = NULL, b = 1, v = p*k+.5, lwd = 3)
                                                            
                                                            panel.abline(a = NULL, b = 1, v = b1, lwd = 2)
                                                            
                                                            
                                                          }, scales = list(x = list(alternating = 1, labels = text, 
                                                                                    cex = 1, at = at, tck = c(0, 0)), y = list(alternating = 0,tck
                                                                                                                               = c(0, 0))))
    
    return(L2)
    
  }
  
)


### Support functions for BigVAR Package
### These are mostly utility functions that will not be seen by the user


# mean benchmark
.evalMean <- function(Y,T1,T2,h=1)
{
  
  ypredF <- NULL
  if(!"matrix"%in%class(Y)){
    Y <- matrix(Y,ncol=1)
  }
  
  MSFE <- c()
  
  k <- ncol(Y)
  
  for (u in (T1-h+2):T2) {
    
    if(h+u-1>T2){break}
    
    trainY1 <- Y[1:(u-1), ]
    
    if(k>1){
      
      ypred <- colMeans(trainY1)
      ypredF <- rbind(ypredF,ypred)
    }else{
      ypred <- mean(trainY1)
      ypredF <- c(ypredF,ypred)
    }
    
    uhat <- matrix(Y[u+h-1, ] - ypred, 
                   ncol = k)
    
    MSFE <- c(MSFE,norm2(uhat)^2)
    
  }
  return(list(Mean=mean(na.omit(MSFE)),SD=sd(na.omit(MSFE))/sqrt(length(na.omit(MSFE))),preds=as.matrix(ypredF)))
}

# random walk benchmark
.evalRW <- function(Y,T1,T2,h=1)
{
  
  if(!"matrix"%in%class(Y)){
    Y <- matrix(Y,ncol=1)
  }
  ypredF <- NULL
  MSFE <- c()
  
  k <- ncol(Y)
  
  for (u in (T1-h+2):T2) {
    
    if(h+u-1>T2){break}
    
    trainY1 <- Y[u-1, ]
    ypredF <- rbind(ypredF,trainY1)
    uhat <- matrix(Y[u+h-1, ] - trainY1, 
                   ncol = k)
    
    MSFE <- c(MSFE,norm2(uhat)^2)
  }
  return(list(Mean=mean(na.omit(MSFE)),SD=sd(na.omit(MSFE))/sqrt(length(na.omit(MSFE))),preds=as.matrix(ypredF)))
}

# Constructs the Grid of Lambda Values: VAR
.LambdaGridE<- function (gran1, gran2, jj = jj, Y, Z, group,p,k,MN,alpha,C,intercept,tol) 
{
  
  if (group == "Lag") {
    mat <- list()
    for (i in 1:length(jj)) {
      if(k>1){
        mat[[i]] <- norm2(Z[jj[[i]], ]%*%Y)
      }
      
      else{
        mat[[i]] <- norm2(t(Y)%*%Z[jj[[i]],])
      }
    }
    gamstart <- max(unlist(mat))
  }
  if (group == "Basic"|group=="Tapered") {
    
    gamstart <- max(t(Y) %*% t(Z))
    
  }
  
  if (group == "SparseLag") {
    
    mat <- list()
    
    if(alpha>0){
      for (i in 1:length(jj)) {
        
        if(k>1){
          
          mat[[i]] <- norm2(Z[jj[[i]], ]%*%Y)*(1/(alpha))
          
        }
        
        else{mat[[i]] <- norm2(t(Y)%*%Z[jj[[i]],])
        
        }
        
      }
      
      gamstart <- max(unlist(mat))
      
    }else{
      gamstart <- max(t(Y) %*% t(Z))
      
    }
  }
  
  if (group == "OwnOther") {
    
    mat <- list()
    
    ZZ <- kronecker(t(Z), diag(k))
    
    for (i in 1:length(jj)) {
      
      mat[[i]] <- norm(as.vector(t(Y)) %*% ZZ[, jj[[i]]]/sqrt(length(jj[[i]])), 
                       "F")
      
    }
    
    gamstart <- max(unlist(mat))
    
  }
  if (group == "SparseOO") {
    
    mat <- list()
    
    ZZ <- kronecker(t(Z), diag(k))
    
    if(alpha>0){
      for (i in 1:length(jj)) {
        
        mat[[i]] <- norm(1/(k + 1) * as.vector(t(Y)) %*% 
                           ZZ[, jj[[i]]]/sqrt(length(jj[[i]])), "F")
        
      }
      
      gamstart <- max(unlist(mat))
      
    }else{
      gamstart <- max(t(Y) %*% t(Z))
      
    }
  }
  
  if (group == "HVARC"|group=="HVAROO"|group=="HVARELEM") {
    
    gmax <- c()
    
    for (i in 1:k) {
      
      gmax[i]  <- norm2(Z %*% Y[, i])
      
    }
    
    gamstart <- max(gmax)
    
  }
  
  if(group=="Tapered"){
    
    beta <- array(0,dim=c(k,k*p+1,10))
    
  }else{
    beta <- array(0,dim=c(k,k*p+1,1))
  }
  
  gamstart <- LGSearch(gamstart,Y,Z,beta,group,k,p,jj,MN,alpha,C,intercept,tol)
  
  gamm <- exp(seq(from = log(gamstart), to = log(gamstart/gran1), 
                  length = gran2))
  
  return(gamm)
  
}
# Constructs penalty grid for each value of alpha in case of dual cv
.LambdaGridEDual <- function(gran1,gran2,jj,Y,Z,group,p,k,MN,alpha,C,intercept,tol){
  
  Lambda <- matrix(0,nrow=gran2,ncol=length(alpha))
  
  for(i in 1:length(alpha)){
    
    Lambda[,i] <- .LambdaGridE(gran1,gran2,jj,Y,Z,group,p,k,MN,alpha[i],C,intercept,tol)
    
  }
  return(Lambda)
  
  
  
  
}

# Construct Lambda Grid: VARX
.LambdaGridX<- function (gran1, gran2, jj = jj, Y, Z, group,p,k1,s,m,k,MN,alpha,C,intercept,tol) 
{
  
  if (group == "Lag") {
    
    mat  <- list()
    
    for (i in 1:length(jj)) {
      
      if(k>1){
        
        mat[[i]] <- norm2(Z[jj[[i]], ]%*%Y)
        
      }
      
      else{
        mat[[i]] <- norm2(t(Y)%*%Z[jj[[i]],])
      }
      
    }
    
    gamstart <- max(unlist(mat))
    
  }
  
  if (group == "Basic") {
    
    gamstart <- max(t(Y)%*%t(Z))
    
  }
  
  if (group == "SparseLag") {
    
    mat  <- list()
    
    if(alpha>0){
      
      for (i in 1:length(jj)) {
        
        if(k>1){
          
          mat[[i]] <- norm2(Z[jj[[i]], ]%*%Y*(1/(alpha)))
          
        }
        
        else{
          
          mat[[i]] <- norm2(t(Y)%*%Z[jj[[i]],])
          
        }
        
      }
      
      gamstart <- max(unlist(mat))
      
    }else{
      
      gamstart <- max(t(Y)%*%t(Z))
      
      
      
      
    }
    
    
  }
  
  if (group == "OwnOther") {
    
    mat <- list()
    
    ZZ <- kronecker(t(Z), diag(k1))
    
    for (i in 1:length(jj)) {
      
      mat[[i]] <- norm(as.vector(t(Y)) %*% ZZ[, jj[[i]]]/sqrt(length(jj[[i]])), 
                       "F")
      
    }
    
    gamstart <- max(unlist(mat))
    
  }
  
  if (group == "SparseOO") {
    
    mat <- list()
    
    ZZ <- kronecker(t(Z), diag(k1))
    
    if(alpha>0){
      for (i in 1:length(jj)) {
        
        mat[[i]] <- norm(1/(alpha) * as.vector(t(Y)) %*% 
                           ZZ[, jj[[i]]]/sqrt(length(jj[[i]])), "F")
      }
      
      gamstart <- max(unlist(mat))
      
    }else{
      
      gamstart <- max(t(Y)%*%t(Z))
      
    }
    
  }
  
  if (group=="EFX") {
    
    gmax <- c()
    
    for (i in 1:k1) {
      
      gmax[i] <- norm2(Z %*% Y[, i])/sqrt(k*p)
      
    }
    
    gamstart <- max(gmax)
    
  }
  
  
  beta <- array(0,dim=c(k1,k1*p+s*m+1,1))
  gamstart <- LGSearchX(gamstart,Y,Z,beta,group,k1,p,s,m,jj,k,MN,alpha,C,intercept,tol)
  gamm <- exp(seq(from = log(gamstart), to = log(gamstart/gran1), 
                  length = gran2))
  
  return(gamm)
  
}

.LambdaGridXDual <- function(gran1,gran2,jj,Y,Z,group,p,k1,s,m,k,MN,alpha,C,intercept,tol){
  
  Lambda <- matrix(0,nrow=gran2,ncol=length(alpha))
  for(i in 1:length(alpha)){
    
    Lambda[,i] <- .LambdaGridX(gran1,gran2,jj,Y,Z,group,p,k1,s,m,k,MN,alpha[i],C,intercept,tol)
    
  }
  return(Lambda)
  
  
  
}

# Forecast evaluation: VARX (called in cv.bigvar)
.BigVAREVALX <- function(ZFull,gamopt,k,p,group,h,MN,verbose,RVAR,palpha,T2,T,k1,s,m,contemp,alpha,C,intercept,tol)
{
  
  if(contemp){
    s1 <- 1
    
  }else{
    
    s1 <- 0
    
  }
  
  preds <- matrix(NA,nrow=length((T2+1):T),ncol=k1)
  
  gran2 <- 1	
  
  gamm <- gamopt
  
  Y <- ZFull$Y
  
  MSFE <- rep(NA,length((T2+1):T))
  
  beta <- array(0,dim=c(k1,k1*p+(k-k1)*(s+s1)+1,1))  
  
  if (group == "Lag") {
    
    jj <- groupfunVARX(p,k,k1,s+s1)
    
    jjcomp <- groupfunVARXcomp(p,k,k1,s+s1)
    
    activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                     gran2)
    
  }
  
  if (group == "SparseLag") {
    
    jj <- groupfunVARX(p, k,k1,s+s1)
    
    q1a <- list()
    
    for (i in 1:(p+s+s1)) {
      
      q1a[[i]] <- matrix(runif(k1, -1, 1), ncol = 1)
      
    }
    
    activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                     gran2)
    
    
  }
  
  if (group == "OwnOther") {
    
    kk <- diaggroupfunVARX(p, k,k1,s+s1)
    activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                     gran2)
  }
  
  if (group == "SparseOO") {
    
    kk <- diaggroupfunVARX(p, k,k1,s+s1)
    
    activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                     gran2)
  }
  
  
  if(verbose){
    
    print("Evaluation Stage")
    
    pb <- txtProgressBar(min = T2-h+2, max = T, style = 3)
  }
  
  for (v in (T2-h+2):T) {
    
    if(v+h-1>T){
      break
    }
    
    trainY <- ZFull$Y[h:(v-1), ]
    
    trainZ <- ZFull$Z[,1:(v-h)]
    
    if (group == "Basic") {
      
      beta <- .lassoVARFistX(beta, trainZ, trainY,gamm, tol,p,MN,k,k1,s+s1,m,C,intercept)
    }
    
    if(group=="BGR"){
      trainZ <- rbind(1,trainZ)
      beta <- BGRGridSearch(trainY,trainZ,p,gamm,as.numeric(MN))
      
    }
    
    
    if (group == "Lag") {
      
      GG <- .GroupLassoVAR1(beta,jj,jjcomp,trainY,trainZ,gamm,activeset,tol,p,MN,k,k1,s+s1,C,intercept)
      
      
      beta <- GG$beta
      
      activeset <- GG$active
      
    }
    
    if (group == "SparseLag") {
      
      GG <- .SparseGroupLassoVARX(beta, jj, trainY, trainZ, 
                                  gamm, alpha, INIactive = activeset, tol, q1a,p,MN,k,s+s1,k1,C,intercept)
      
      beta <- GG$beta
      
      activeset <- GG$active
      
      q1a <- GG$q1
      
    }
    
    if (group == "OwnOther") {
      
      GG <- .GroupLassoOOX(beta, kk, trainY, trainZ, gamm, 
                           activeset, tol,p,MN,k,k1,s+s1,C,intercept)
      
      beta <- GG$beta
      
      activeset <- GG$active
      
    }
    if (group == "SparseOO") {
      
      
      GG <- .SparseGroupLassoVAROOX(beta, kk, trainY, trainZ, 
                                    gamm, alpha, INIactive = activeset, tol,p,MN,k1,s+s1,k,FALSE,C,intercept)
      
      beta <- GG$beta
      
      activeset <- GG$active
      
    }
    
    
    if(group=="EFX"){
      
      beta <- .EFVARX(beta,trainY,trainZ,gamm,tol,MN,k1,s,m,p,C,intercept)
    }
    
    betaEVAL <- matrix(beta[,,1],nrow=k1,ncol=(k1*p+(k-k1)*(s+s1)+1))
    
    if (RVAR) {
      
      betaEVAL <- RelaxedLS(cbind(t(trainZ),trainY),betaEVAL,k,p,k1,s+s1)
      
    }
    
    
    if(MN){
      
      eZ <- matrix(ZFull$Z[,v],ncol=1)
      
    }else{
      
      eZ <- matrix(c(1,ZFull$Z[,v]),ncol=1)
      
    }
    
    if(MN){
      
      MSFE[v-T2+h-1] <- norm2(ZFull$Y[v+h-1,1:k1] - betaEVAL[,2:ncol(betaEVAL)] %*% eZ)^2
      
      preds[v-T2+h-1,] <-  betaEVAL[,2:ncol(betaEVAL)] %*% eZ
      
      diag(beta[,2:(k1+1),1]) <- diag(beta[,2:(k1+1),1])-C # subtract one for warm start purposes 
      
    }else{
      MSFE[v-T2+h-1] <- norm2(ZFull$Y[v+h-1,1:k1] - betaEVAL %*% eZ)^2
      
      preds[v-T2+h-1,] <- betaEVAL %*% eZ
      
    }
    
    if(verbose){
      
      setTxtProgressBar(pb, v)
      
    }
  }
  # Parameter estimates using all available data
  if (group == "Basic") {
    
    betaPred <- .lassoVARFistX(beta, ZFull$Z, ZFull$Y,gamm, tol,p,MN,k,k1,s+s1,m,C,intercept)            
  }
  if (group == "Lag") {
    
    GG <- .GroupLassoVAR1(beta,jj,jjcomp,ZFull$Y,ZFull$Z,gamm,activeset,tol,p,MN,k,k1,s+s1,C,intercept)
    
    betaPred <- GG$beta
    
  }
  
  if (group == "SparseLag") {
    
    GG <- .SparseGroupLassoVARX(beta, jj, ZFull$Y, ZFull$Z, 
                                gamm, alpha, INIactive = activeset, tol, q1a,p,MN,k,s+s1,k1,C,intercept)
    
    betaPred <- GG$beta
    
  }
  
  if (group == "OwnOther") {
    
    GG <- .GroupLassoOOX(beta, kk, ZFull$Y, ZFull$Z, gamm, 
                         activeset, tol,p,MN,k,k1,s+s1,C,intercept)
    
    betaPred <- GG$beta
    
  }
  if (group == "SparseOO") {
    
    GG <- .SparseGroupLassoVAROOX(beta, kk, ZFull$Y, ZFull$Z, 
                                  gamm, alpha, INIactive = activeset, tol,p,MN,k1,s+s1,k,FALSE,C,intercept)
    
    betaPred <- GG$beta
    
    
  }
  
  if(group=="EFX")
    
  {
    
    
    betaPred <- .EFVARX(beta,ZFull$Y,ZFull$Z,gamm,tol,MN,k1,s+s1,m,p,C,intercept)
    
    
  }
  
  
  betaPred <- as.matrix(betaPred[,,1])
  
  
  return(list(MSFE=MSFE,betaPred=betaPred,predictions=preds))
  
  
}


# Forecast evaluation: VAR (called in cv.bigvar)
.BigVAREVAL <- function(ZFull,gamopt,k,p,group,h,MN,verbose,RVAR,palpha,T1,T2,alpha,recursive,C,intercept,tol)
{
  
  
  gran2 <- 1
  
  preds <- matrix(NA,nrow=length((T1+1):T2),ncol=k)
  
  
  gamm <- gamopt
  
  Y <- ZFull$Y
  
  s <- p
  k1 <- k
  
  MSFE <- rep(NA,length((T1+1):T2))
  
  
  beta <- array(0,dim=c(k,p*k+1,1))
  
  if (group == "Lag")
  {
    
    jj <- .groupfuncpp(p,k)
    
    jjcomp <- .groupfuncomp(p,k)
    
    activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                     gran2)
    
  }
  
  if (group == "SparseLag")
  {
    
    jj <- .groupfun(p, k)
    
    q1a <- list()
    
    for (i in 1:(p)) {
      
      q1a[[i]] <- matrix(runif(k, -1, 1), ncol = 1)
      
    }
    
    activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                     gran2)
    
    
    
  }
  
  
  if (group == "OwnOther")
  {
    
    
    kk <- .lfunction3cpp(p, k)
    
    activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                     gran2)
    
  }
  
  if (group == "SparseOO") {
    
    kk <- .lfunction3cpp(p, k)
    
    jjcomp <- .lfunctioncomp(p,k)
    
    jj=.lfunction3(p,k)
    
    activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                     gran2)
    
    q1a <- list()
    
    for (i in 1:(2*p))
    {
      
      q1a[[i]] <- matrix(runif(length(kk[[i]]), -1, 1), ncol = 1)
      
    }
    
  }
  
  
  if(verbose)
  {
    
    print("Evaluation Stage")
    
    pb <- txtProgressBar(min = T1-h+2, max = T2, style = 3)
    
  }
  
  for (v in (T1-h+2):T2)
  {
    
    if(h>1 & ! recursive){
      
      
      trainY <- ZFull$Y[h:(v-1), ]
      
      
      trainZ <- ZFull$Z[,1:(v-h)]
      
    }else{
      
      
      trainY <- ZFull$Y[1:(v-1), ]
      
      trainZ <- ZFull$Z[,1:(v-1)]
    }
    if(v+h-1>T2){
      break
    }
    if (group == "Basic") {
      
      beta <- .lassoVARFist(beta, trainZ, trainY,gamm, tol,p,MN,C,intercept)
      
    }
    
    if (group == "Lag") {
      
      GG <- .GroupLassoVAR1(beta,jj,jjcomp,trainY,trainZ,gamm,activeset,tol,p,MN,k,k1,s,C,intercept)
      
      beta <- GG$beta
      
      activeset <- GG$active
    }
    
    
    if (group == "SparseLag") {
      
      GG <- .SparseGroupLassoVAR(beta, jj, trainY, trainZ, 
                                 gamm, alpha, INIactive = activeset, tol, q1a,p,MN,C,intercept)
      
      beta <- GG$beta
      
      activeset <- GG$active
      
      q1a <- GG$q1
      
    }
    
    
    if (group == "OwnOther") {
      
      GG <- .GroupLassoOO(beta, kk, trainY, trainZ, gamm, 
                          activeset, tol,p,MN,C,intercept)
      
      beta <- GG$beta
      
      activeset <- GG$active
      
    }
    
    if (group == "SparseOO") {
      
      GG <- .SparseGroupLassoVAROO(beta, kk, trainY, trainZ, 
                                   gamm, alpha, INIactive = activeset, tol,q1a,p,MN,FALSE,C,intercept)
      
      beta <- GG$beta
      
      activeset <- GG$active
      
      q1a <- GG$q1
      
    }
    
    if (group=="HVARC"){
      
      beta <- .HVARCAlg(beta,trainY,trainZ,gamm,tol,p,MN,C,intercept)
      
    }
    
    
    if(group=="HVAROO"){
      
      beta <- .HVAROOAlg(beta,trainY,trainZ,gamm,tol,p,MN,C,intercept)
      
    }
    
    if(group=="HVARELEM"){
      
      beta <- .HVARElemAlg(beta,trainY,trainZ,gamm,tol,p,MN,C,intercept) 	
      
      
    }
    
    if(group=="Tapered"){
      
      beta <- .lassoVARTL(beta, trainZ, trainY,gamm, tol,p,MN,palpha,C,intercept)
    }
    
    if(group=="BGR"){
      trainZ <- rbind(1,trainZ)
      beta <- BGRGridSearch(trainY,trainZ,p,gamm,as.numeric(MN))
      
    }
    
    
    ## betaEVAL <- matrix(beta[,,1],nrow=k,ncol=(k*p+1))
    ## if(group!="BGR"){
    betaEVAL <- matrix(beta[,,1],nrow=k,ncol=(k*p+1))
    
    if (RVAR) {
      
      betaEVAL <- RelaxedLS(cbind(t(trainZ),trainY),betaEVAL,k,p,k1,s)      
    }
    
    if(MN){
      eZ <- matrix(ZFull$Z[,v],ncol=1)
      
    }else{
      
      
      eZ <- matrix(c(1,ZFull$Z[,v]),ncol=1)
      
      
    }
    
    # We don't consider an intercept for the MN lasso
    if(MN){
      
      preds[v-T1+h-1,] <- betaEVAL[,2:ncol(betaEVAL)] %*% eZ
      
      if(h>1 & recursive){
        
        pred <- matrix(preds[v-T1,],nrow=1)
        
        preds[v-T1+h-1,] <- predictMS(pred,trainY,h-1,betaEVAL[,2:ncol(betaEVAL)],p,TRUE)
        
      }else{
        
        MSFE[v-T1+h-1] <- norm2(ZFull$Y[v+h-1, ] - preds[v-T1-h+1,])^2
        
        diag(beta[,2:(k+1),1]) <- diag(beta[,2:(k+1),1])-C # subtract one for warm start purposes 
        
      }
      
    }else{
      preds[v-T1+h-1,] <- betaEVAL %*% eZ                       
      if(h>1 & recursive){
        
        pred <- matrix(preds[v-T1+h-1,],nrow=1)
        
        preds[v-T1+h-1,] <- predictMS(pred,trainY,h-1,betaEVAL,p,FALSE)
        
      }
      
      MSFE[v-T1+h-1] <- norm2(ZFull$Y[v+h-1,] - preds[v-T1+h-1,])^2                    
      
    }
    ## }else{
    
    ##     ## browser()
    ##     MSFE[v - (T1 - h)] <- norm2(Y[v+h-1,1:k1] - beta[,,1])^2
    
    ##     }
    if(verbose){
      
      
      setTxtProgressBar(pb, v)
      
    }
  }
  
  if (group == "Basic") {
    betaPred <- .lassoVARFist(beta, ZFull$Z, ZFull$Y,gamm, tol,p,MN,C,intercept)
    
  }
  
  if (group == "Lag") {
    
    GG <- .GroupLassoVAR1(beta,jj,jjcomp,ZFull$Y,ZFull$Z,gamm,activeset,tol,p,MN,k,k,s,C,intercept)
    
    
    betaPred <- GG$beta
    
  }
  
  if (group == "SparseLag") {
    
    GG <- .SparseGroupLassoVAR(beta, jj, ZFull$Y, ZFull$Z, 
                               gamm, alpha, INIactive = activeset, tol, q1a,p,MN,C,intercept)
    
    betaPred <- GG$beta
    
  }
  if(group=="BGR"){
    ZF <- rbind(1,ZFull$Z)
    betaPred <- BGRGridSearch(ZFull$Y,ZF,p,gamm,as.numeric(MN))
    
  }
  
  
  
  if (group == "OwnOther") {
    
    ## beta <- array(0,dim=c(k,p*k+1,1))
    ##             activeset <- rep(list(rep(rep(list(0), length(kk)))), 
    ##                  gran2)
    
    
    GG <- .GroupLassoOO(beta, kk, ZFull$Y, ZFull$Z, gamm, 
                        activeset, tol,p,MN,C,intercept)
    
    betaPred <- GG$beta
    
  }
  
  if (group == "SparseOO") {
    
    GG <- .SparseGroupLassoVAROO(beta, kk, ZFull$Y, ZFull$Z, 
                                 gamm, alpha, INIactive = activeset, tol,q1a,p,MN,FALSE,C,intercept)
    
    betaPred <- GG$beta
    
    
  }
  
  if (group=="HVARC")
  {
    
    betaPred <- .HVARCAlg(beta,ZFull$Y,ZFull$Z,gamm,tol,p,MN,C,intercept)
    
  }
  
  if(group=="HVAROO")
  {
    
    betaPred <- .HVAROOAlg(beta,ZFull$Y,ZFull$Z,gamm,tol,p,MN,C,intercept)
    
    
  }
  
  if(group=="HVARELEM")
  {
    
    betaPred<-.HVARElemAlg(beta,ZFull$Y,ZFull$Z,gamm,tol,p,MN,C,intercept) 	
    
  }
  
  if(group=="Tapered")
  {
    betaPred <- .lassoVARTL(beta,ZFull$Z,ZFull$Y,gamm,tol,p,MN,palpha,C,intercept)            
  }    
  
  ## if(group!="BGR"){
  betaPred <- as.matrix(betaPred[,,1])
  ## browser()
  
  betaPred <- matrix(betaPred,nrow=k)
  ## if(group=="BGR"){
  ##     betaPred <- t(betaPred)
  ##     }
  ## }else{
  ##     betaPred <- matrix(0,nrow=1,ncol=1)
  ##     }
  
  ## betaPred <- as.matrix(betaPred[,,1])
  
  ## betaPred <- matrix(betaPred,nrow=k)
  
  
  
  return(list(MSFE=MSFE,betaPred=betaPred,predictions=preds))
  
  
}



#' Converts a VAR coefficient matrix of order p to multiple companion form
#' 
#' @param B a \eqn{k \times kp} coefficient matrix
#' @param p Lag order
#' @param k Number of Series
#' @return Returns a \eqn{kp \times kp} coefficient matrix representing all coefficient matrices contained in Ai as a VAR(1).
#' @references See page 15 of Lutkepohl, "A New Introduction to Multiple Time Series Analysis"
#' @seealso \code{\link{MultVarSim}}
#' @examples
#' k=3;p=6
#' B=matrix(0,nrow=k,ncol=p*k)
#' A1<- matrix(c(.4,-.02,.01,-.02,.3,.02,.01,.04,.3),ncol=3,nrow=3)
#' A2 <- matrix(c(.2,0,0,0,.3,0,0,0,.13),ncol=3,nrow=3)
#' B[,1:k]=A1
#' B[,(4*k+1):(5*k)]=A2
#' A <- VarptoVar1MC(B,p,k)
#' @export
VarptoVar1MC <- function(B,p,k)
{
  
  Fp=matrix(0,nrow=k*p,ncol=k*p)      
  
  Fp[1:k,] = B
  
  Fp[-(1:k),1:(k*(p-1))] = diag(k*(p-1))
  # We require that the coefficient matrix generates a stationary VAR
  if(max(Mod(eigen(Fp)$values))>1){warning("Coefficient Matrix is not stationary")}
  
  return(Fp)
  
}


#' Simulate a VAR
#' 
#' @param k Number of Series
#' @param A1 Either a \eqn{k \times k} coefficient matrix or a \eqn{kp \times kp} matrix created using \code{\link{VarptoVar1MC}}. 
#' @param p Maximum Lag Order
#' @param Sigma Residual Covariance Matrix of dimension \eqn{k\times k}
#' @param T Number of simulations
#' @return Returns a \eqn{T \times k} of realizations from a VAR.
#' @references Lutkepohl, "A New Introduction to Multiple Time Series Analysis"
#' @seealso \code{\link{VarptoVar1MC}}
#' @examples
#' k=3;p=6
#' B=matrix(0,nrow=k,ncol=p*k)
#' A1<- matrix(c(.4,-.02,.01,-.02,.3,.02,.01,.04,.3),ncol=3,nrow=3)
#' A2 <- matrix(c(.2,0,0,0,.3,0,0,0,.13),ncol=3,nrow=3)
#' B[,1:k]=A1
#' B[,(4*k+1):(5*k)]=A2
#' A <- VarptoVar1MC(B,p,k)
#' Y <-MultVarSim(k,A,p,.1*diag(k),100)
#' @export
#' @importFrom MASS mvrnorm
MultVarSim <- function (k, A1, p, Sigma, T) 
{
  
  if(max(Mod(eigen(A1)$values))>1){stop("Error: Generator Matrix is not stationary")}
  
  # add 500 observations for initialization purposes
  
  Y <- matrix(0, nrow = T+500+p , ncol = k)
  
  YY <- as.vector(Y)
  
  for (i in seq(from = (k * p + 1), to = (nrow(Y) * k - 1), 
                by = k)) {
    
    u <- as.vector(c(mvrnorm(1, rep(0, k), Sigma), rep(0, 
                                                       k * p - k)))
    
    YY[(i + k):(i - k * p + 1 + k)] <- A1 %*% YY[(i):(i - 
                                                        k * p + 1)] + as.matrix(u, ncol = 1)
    
  }
  
  YY <- YY[YY!=0]
  
  Y <- matrix(YY, ncol = k, byrow = TRUE)
  
  Y <- Y[,c(ncol(Y):1)]
  
  Y <- Y[501:nrow(Y), ]
  
  return(Y)
}

# function to create subsets for lag group VARX-L
.groupfun <- function(p,k)
{
  
  jjj <- list()
  jjj[[1]] <- 1:k
  
  if(p>1){
    
    for(i in 2:p){
      jjj[[i]] <- jjj[[i-1]]+k
      
    }
    
  }
  
  return(jjj)
  
}

#C++ groupings to account for 0 indexing
.groupfuncpp <- function(p,k)
{
  jjj <- list()
  
  jjj[[1]] <- 0:(k-1)
  
  if(p>1)
  {
    for(i in 2:p){
      
      jjj[[i]] <- jjj[[i-1]]+k
    }
  }
  return(jjj)
}

# subsetting complement of groups in rcpp
.groupfuncomp <- function(p,k)
{
  
  
  ownoth <- .groupfuncpp(p,k)
  
  kk2 <- list()
  
  pmax <- max(unlist(ownoth))
  
  to <- 0:(pmax)
  
  for(i in 1:length(ownoth))
  {
    
    kk2[[i]] <- to[is.na(pmatch(to,ownoth[[i]]))]
  }
  
  return(kk2)
  
}    

# Group indexing for own/other VARX-L
.lfunction2 <- function(p,k)
{
  
  kk <- list()
  
  kk[[1]] <- 1:(k^2)
  
  if(p>1)
  {
    for(i in 2:p)
    {
      kk[[i]] <- 1:(k^2)+tail(kk[[i-1]],1)
      
    }
  }    
  return(kk)
}

.lfunction2cpp <- function(p,k)
{
  
  kk <- list()
  
  kk[[1]] <- 0:(k^2-1)
  
  if(p>1)
    
  {
    
    for(i in 2:p)
      
    {
      
      kk[[i]] <- 0:(k^2-1)+tail(kk[[i-1]],1)+1
      
      
    }
    
  }    
  
  return(kk)
  
  
}

.lfunction3 <- function(p,k)
{
  
  kk <- .lfunction2(p,k)
  
  oo <- list()
  
  pp <- list()
  
  for(i in 1:length(kk))
  {
    j <- 0
    oo[[i]] <- kk[[i]][(seq(1,length(kk[[1]]),k+1)+(j*k^2))]
    pp[[i]] <- kk[[i]][-(seq(1,length(kk[[1]]),k+1)+(j*k^2))]
    j <- j+1
    
  }
  
  ownoth <- c(oo,pp)
  return(ownoth)
}    

.lfunction3cpp <- function(p,k)
{
  
  kk <- .lfunction2cpp(p,k)
  oo <- list()
  pp <- list()
  
  for(i in 1:length(kk))
  {
    j <- 0
    oo[[i]] <- kk[[i]][(seq(1,length(kk[[1]]),k+1)+(j*k^2))]
    pp[[i]] <- kk[[i]][-(seq(1,length(kk[[1]]),k+1)+(j*k^2))]
    j <- j+1
  }
  
  ownoth <- c(oo,pp)
  
  return(ownoth)
  
}    


.lfunctioncomp <- function(p,k)
{
  
  ownoth <- .lfunction3cpp(p,k)
  kk2 <- list()
  pmax <- max(unlist(ownoth))
  to <- 0:(pmax)
  for(i in 1:length(ownoth))
  {
    kk2[[i]] <- to[is.na(pmatch(to,ownoth[[i]]))]
  }
  
  return(kk2)
}    


# This function should work for arbitrary groups
.lfunction <- function(groups,p)
{
  
  H <- as.vector(do.call('cbind',groups))
  kk <- list()
  kk[[1]] <- H
  if(p>1)
  {
    for(i in 2:p){
      kk[[i]] <-as.vector(do.call('cbind',groups))+tail(kk[[i-1]],1) }
  }
  return(kk)
  
}

# Indexing for HVAR
.vsubs <- function(p,k)
{
  vi <- list()
  
  for(i in p:1)
  {
    g <- max(k*i-k,0)
    vi[[i]] <- g:(k*(p)-1)
  }
  return(vi)   
  
}

#indexing for HVAR OO
.oofun <- function(p,k)
{
  
  kk <- .lfunction2(p, k)
  oo <- list()
  pp <- list()
  
  for (i in 1:length(kk)) {
    j <- 0
    oo[[i]] <- kk[[i]][(seq(1, length(kk[[1]]), k + 1) + 
                          (j * k^2))]
    pp[[i]] <- kk[[i]][-(seq(1, length(kk[[1]]), k + 1) + 
                           (j * k^2))]
    j <- j + 1
  }
  ownoth <- list()
  jj <- .lfunction3(p,k)
  
  for(i in 1:length(jj))
    
  {
    
    if(i==1)
      
    {
      
      ownoth[[i]] <- jj[[i]]
      oo[[1]] <- NULL
      
    }
    
    if(i==length(jj))
    {
      ownoth[[i]] <- tail(jj,1)
      pp[[1]] <- NULL
    }
    
    if(i!=1& i%%2!=0)
    {
      
      ownoth[[i]] <- head(oo,1)
      oo[[1]] <- NULL
    }
    
    if(i!=length(jj)&i%%2==0)
    {
      ownoth[[i]] <- head(pp,1)
      pp[[1]] <- NULL
    }
    
    
  }
  
  return(rev(ownoth))
  
}

.oocumfun <-function(p,k)
{
  
  kk <- rev(.oofun(p,k))
  oogroups <- list()
  oogroups[[1]] <- unlist(kk)
  for(i in 2:length(kk))
  {
    oogroups[[i]] <- unlist(kk[-(1:(i-1))])
    
  }
  
  return(oogroups)
}


# indexing function for Own/Other HVAR
.vecoovars<-function(p,k,k1)
{
  
  vv <- list()
  
  vv[[1]] <- 1:(p*k)
  
  vv[[2]] <- vv[[1]][-k1]
  
  q1 <- 1
  
  if(p>1){
    for(i in 3:(2*p))
    {
      if(i%%2!=0)
      {
        vv[[i]] <- (q1*k+1):(k*p)
        q1 <- q1+1
      }else{
        vv[[i]] <- vv[[i-1]][-k1]
      }
    }
  }
  return(vv)
  
}

# indexing to start at zero for use within rcpp
.vecoovarscpp<-function(p,k,k1)
{
  
  vv <- list()
  vv[[1]] <- 0:(p*k-1)
  vv[[2]] <- vv[[1]][-(k1)]
  q1 <- 1
  
  if(p>1){
    for(i in 3:(2*p))
    {
      if(i%%2!=0)
      {
        vv[[i]] <- (q1*k):(k*p-1)
        q1 <- q1+1
      }else{
        
        vv[[i]] <- vv[[i-1]][-(k1)]
        
      }
      
    }
    
  }
  
  return(vv)
  
}


#VARX Lag Group function
groupfunVARX <- function(p,k,k1,s)
{
  
  jj <- list()
  m <- k-k1
  jj <- .groupfuncpp(p, k1)
  kp <- k1*p+m*s-1
  jj2 <- list()
  startjj <- max(unlist(jj))+1
  for(i in seq(startjj,kp,by=1))
  {
    jj[[i]] <- i                
  }
  jj[sapply(jj, is.null)] <- NULL
  
  return(jj)
}

groupfunVARXcomp <- function(p,k,k1,s)
{
  ownoth <- groupfunVARX(p,k,k1,s)
  kk2 <- list()
  pmax <- max(unlist(ownoth))
  to <- 0:(pmax)
  for(i in 1:length(ownoth))
  {
    kk2[[i]] <- to[is.na(pmatch(to,ownoth[[i]]))]
  }
  return(kk2)
  
}


# indexing starts at 1 for use in R        
groupfunVARXLG <- function(p,k,k1,s)
{
  
  jj <- list()
  jj <- .groupfun(p, k1)
  m <- k-k1
  kp <- k1*p+m*s
  jj2 <- list()
  startjj <- max(unlist(jj))+1
  for(i in seq(startjj,kp,by=1))
  {
    jj[[i]] <- i
    
  }
  jj[sapply(jj, is.null)] <- NULL
  
  return(jj)
}



diaggroupfunVARX <- function(p,k,k1,s)
{
  m <- k-k1
  jj <- list()
  jj <- .lfunction3cpp(p, k1)
  kp <-k1*(p*k1+s*m)-1
  jj2 <- list()
  startjj <- max(unlist(jj))+1
  
  for(i in seq(startjj,kp,by=k1))
  {
    jj[[i]] <- i:(i+k1-1)
    
  }
  
  jj[sapply(jj, is.null)] <- NULL
  
  
  return(jj)
  
}

diaggroupfunVARXcomp <- function(p,k,k1,s)
{
  
  ownoth <- diaggroupfunVARX(p,k,k1,s)
  kk2 <- list()
  pmax <- max(unlist(ownoth))
  to <- 0:(pmax)
  
  for(i in 1:length(ownoth))
  {
    kk2[[i]] <- to[is.na(pmatch(to,ownoth[[i]]))]
  }
  
  return(kk2)
}

diaggroupfunVARXLG <- function(p,k,k1,s)
{
  
  m <- k-k1
  jj <- list()
  jj <- .lfunction3(p, k1)
  kp <- k1*(p*k1+s*m)
  jj2 <- list()
  startjj <- max(unlist(jj))+1
  for(i in seq(startjj,kp,by=k1))
  {
    jj[[i]] <- i:(i+(k1-1))                
  }
  
  jj[sapply(jj, is.null)] <- NULL
  
  return(jj)
  
}

diaggroupfunVARXLGL <- function(p,k,k1)
{
  
  jj <- list()
  jj <- .lfunction3(p, k1)
  kp <- k1*p*k
  jj2 <- list()
  startjj <- max(unlist(jj))+1
  
  for(i in seq(startjj,kp,by=1))
  {
    jj[[i]] <- i                
  }
  
  jj[sapply(jj, is.null)] <- NULL
  
  return(jj)
}


diaggroupfunVARXL <- function(p,k,k1)
{
  jj <- list()
  jj <- .lfunction3cpp(p, k1)
  kp <- k1*p*k-1
  jj2 <- list()
  startjj <- max(unlist(jj))+1
  
  for(i in seq(startjj,kp,by=1))
  {
    jj[[i]] <- i
    
  }
  
  jj[sapply(jj, is.null)] <- NULL
  return(jj)
  
}

diaggroupfunVARXcompL <- function(p,k,k1)
{
  
  ownoth <- diaggroupfunVARXL(p,k,k1)
  kk2 <- list()
  pmax <- max(unlist(ownoth))
  
  to <- 0:(pmax)
  
  for(i in 1:length(ownoth))
    
  {
    
    kk2[[i]] <- to[is.na(pmatch(to,ownoth[[i]]))]
    
  }
  
  return(kk2)
  
}


# iterative procedure to find a less coarse bound for lambda starting value via binary search
LGSearchX <- function(gstart,Y,Z,BOLD,group,k1,p,s,m,gs,k,MN,alpha,C,intercept,tol)
{
  
  tk <- 1/max(Mod(eigen(Z%*%t(Z))$values))
  lambdah <- gstart
  lambdal <- 0
  activeset <- list(rep(rep(list(0), length(gs))))
  
  while(max(lambdah-lambdal)>.00001)
  {
    
    lambda <- (lambdah+lambdal)/2
    if(group=="EFX"){
      BOLD <- .EFVARX(BOLD,Y,Z,lambda,tol,MN,k1,s,m,p,C,intercept)
      param <- BOLD[,2:(k1*p+m*s+1),1]
    }
    
    if(group=="Basic"){
      param <- .lassoVARFistX(BOLD,Z,Y[,1:k1],lambda,tol,p,MN,k1+m,k1,s,m,C,intercept)[,2:(k1*p+m*s+1),]
    }
    
    if(group=="Lag"){
      
      jj <- groupfunVARX(p,k,k1,s)
      jjcomp <- groupfunVARXcomp(p,k,k1,s)
      BB <- .GroupLassoVAR1(BOLD,jj,jjcomp,Y[,1:k1],Z,lambda,activeset,tol,p,MN,k,k1,s,C,intercept)
      BOLD <- BB$beta
      param <- BB$beta[,2:(k1*p+m*s+1),]
      activeset <- BB$active
    }
    
    
    if(group=="OwnOther")
    {
      
      kk <- diaggroupfunVARX(p, k,k1,s)
      BB <- .GroupLassoOOX(BOLD, kk, Y, Z, lambda,activeset, tol,p,MN,k,k1,s,C,intercept)
      param <- BB$beta[,2:(k1*p+m*s+1),]
      BOLD <- BB$beta
      activeset <- BB$active
    }
    
    if(group=="SparseOO")
    {
      
      kk <- diaggroupfunVARX(p, k,k1,s)
      BB <- .SparseGroupLassoVAROOX(BOLD, kk, Y[,1:k1], Z, lambda,alpha,activeset, tol,p,MN,k1,s,k,FALSE,C,intercept)
      param <- BB$beta[,2:(k1*p+m*s+1),]
      BOLD <- BB$beta
      activeset <- BB$active
      
    }
    
    
    if(group=="SparseLag"){
      
      jj <- groupfunVARX(p, k,k1,s)
      jjcomp <- groupfunVARXcomp(p,k,k1,s)
      q1a=list()
      
      for (i in 1:(p+s)) {
        
        q1a[[i]] <- matrix(runif(length(jj[[i]]), -1, 1), ncol = 1)
        
      }
      
      BB <- .SparseGroupLassoVARX(BOLD,jj,Y[,1:k1],Z,lambda,alpha,activeset,tol,q1a,p,MN,k,s,k1,C,intercept)
      param <- BB$beta[,2:(k1*p+m*s+1),]
      BOLD <- BB$beta
      activeset <- BB$active
    }
    
    if(MN){
      ## diag(param[1:k1,1:k1]) <- 0
      ## diag(BOLD[1:(k1),2:(k1+1),1]) <- 0
      
      
      diag(param[1:k1,1:k1]) <- ifelse(C==0,diag(param[1:k1,1:k1]),diag(param[1:k1,1:k1])-C)
      diag(BOLD[,2:(k1*p+1+m*s),]) <- ifelse(C==0,diag(BOLD[,2:(k1*p+1+m*s),]),diag(BOLD[,2:(k1*p+1+m*s),])-C)
      
    }
    
    if(max(abs(param))<sqrt(.Machine$double.eps))
    {
      lambdah <- lambda
    }else{
      
      lambdal <- lambda
      
    }
    
  }
  
  lambdah
}

# Same as above, but for the VAR
LGSearch <- function(gstart,Y,Z,BOLD,group,k,p,gs,MN,alpha,C,intercept,tol)
{
  
  tk <- 1/max(Mod(eigen(Z%*%t(Z))$values))
  lambdah <- gstart
  lambdal <- 0
  activeset <- list(rep(rep(list(0), length(gs))))
  
  if(group=="SparseOO")
  {
    kk <- .lfunction3cpp(p, k)
    activeset <- rep(list(rep(rep(list(0), length(kk)))),1)
    q1a=list()
    for (i in 1:(2*p)) {
      q1a[[i]] <- matrix(runif(length(gs[[i]]), -1, 1), ncol = 1)
    }
    
  }
  
  while(max(abs(lambdah-lambdal))>.00001)
  {
    
    lambda <- (lambdah+lambdal)/2
    if(group=="Basic"){
      param <- .lassoVARFist(BOLD,Z,Y,lambda,tol,p,MN,C,intercept)[,2:(k*p+1),]
    }
    
    if(group=="Tapered"){
      
      param <- .lassoVARTL(BOLD,Z,Y,lambda,tol,p,MN,rev(seq(0,1,length=10)),C,intercept)[,2:(k*p+1),]
      
      param <- param[,,1]                    
    }
    
    if(group=="Lag"){
      
      jj <- .groupfuncpp(p, k)
      jjcomp <- .groupfuncomp(p,k)
      BB <- .GroupLassoVAR1(BOLD,jj,jjcomp,Y,Z,lambda,activeset,tol,p,MN,k,k,p,C,intercept)
      BOLD <- BB$beta
      param <- BB$beta[,2:(k*p+1),]
      activeset <- BB$active                  
    }
    
    if(group=="SparseLag"){
      
      jj <- .groupfuncpp(p, k)
      jjcomp <- .groupfuncomp(p,k)
      q1a=list()
      for (i in 1:p) {
        q1a[[i]] <- matrix(runif(k, -1, 1), ncol = 1)
      }
      BB <- .SparseGroupLassoVAR(BOLD,jj,Y,Z,lambda,alpha,activeset,tol,q1a,p,MN,C,intercept)
      param <- BB$beta[,2:(k*p+1),]
      BOLD <- BB$beta
      activeset <- BB$active
    }
    
    
    if(group=="OwnOther")
    {
      kk <- .lfunction3cpp(p, k)
      BB <- .GroupLassoOO(BOLD, kk, Y, Z, lambda,activeset, tol,p,MN,C,intercept)
      param <- BB$beta[,2:(k*p+1),]
      BOLD <- BB$beta
      activeset <- BB$active
      
    }
    
    
    if(group=="SparseOO")
    {
      BB <- .SparseGroupLassoVAROO(BOLD, kk, Y, Z, lambda,alpha,activeset, tol,q1a,p,MN,FALSE,C,intercept)
      param <- BB$beta[,2:(k*p+1),]
      BOLD <- BB$beta
      activeset <- BB$active                  
    }
    
    if(group=="HVARC")
    {
      BOLD <- .HVARCAlg(BOLD,Y,Z,lambda,tol,p,MN,C,intercept)                  
      param <- BOLD[,2:(k*p+1),]
    }
    if(group=="HVAROO")
    {
      BOLD <- .HVAROOAlg(BOLD,Y,Z,lambda,tol,p,MN,C,intercept)
      param <- BOLD[,2:(k*p+1),]
    }
    
    if(group=="HVARELEM")
    {
      BOLD <- .HVARElemAlg(BOLD,Y,Z,lambda,tol,p,MN,C,intercept)
      param <- BOLD[,2:(k*p+1),]
    }
    
    
    if(MN){
      
      diag(param[1:k,1:k]) <- ifelse(C==0,diag(param[1:k,1:k]),0)
      if(group!="Tapered"){
        diag(BOLD[,2:(k*p+1),]) <- ifelse(C==0,diag(BOLD[,2:(k*p+1),]),0)
      }
    }
    ## browser()
    if(max(abs(param))< tol)
    {
      lambdah <- lambda
      
    }else{
      lambdal <- lambda
      
      
    }
    
    
  }
  
  lambdah
  
}


#' Evaluate forecasts from a VAR or VARX with lag orders selected by AIC/BIC
#' 
#' @param Y a \eqn{T \times k} multivariate time series 
#' @param X a \eqn{T \times m} multivariate time series of unmodeled exogenous variables
#' @param p maximum lag order for endogenous series
#' @param s maximum lag order for exogenous series
#' @param T1 start of forecast evaluation period.
#' @param T2 end of forecast evaluation period
#' @param IC specifies whether to select lag order according to "AIC" or "BIC"
#' @param h desired forecast horizon
#' @param iterated indicator as to whether to use iterated or direct multistep forecasts (if applicable, VAR context only)
#' @return Returns the one-step ahead MSFE as well as the forecasts over the evaluation period and lag order selected.
#' @details This function evaluates the one-step ahead forecasts of a VAR or VARX fit by least squares over an evaluation period.  At every point in time, lag orders for the endogenous and exogenous series are selected according to AIC or BIC.  This function is run automatically when \code{\link{cv.BigVAR}} is called unless \code{ic} is set to \code{FALSE} in \code{\link{constructModel}}.      
#' @references Neumaier, Arnold, and Tapio Schneider. "Estimation of parameters and eigenmodes of multivariate autoregressive models." ACM Transactions on Mathematical Software (TOMS) 27.1 (2001): 27-57.
#' @seealso \code{\link{VARXFit}},\code{\link{constructModel}}, \code{\link{cv.BigVAR}}
#' @examples
#' data(Y)
#'
#' # Evaluate the performance of a VAR with lags selected by BIC.
#' p <- 4
#' T1 <- floor(nrow(Y))/3
#' T2 <- floor(2*nrow(Y))/3
#' # Matrix of zeros for X
#' X <- matrix(0,nrow=nrow(Y),ncol=ncol(Y))
#' BICMSFE <- VARXForecastEval(Y,X,p,0,T1,T2,"BIC",1)
#' 
#' @export
VARXForecastEval <- function(Y,X,p,s,T1,T2,IC,h,iterated=FALSE)
{
  
  
  if(T1>nrow(Y) | T2>nrow(Y) |T2<T1){stop("Training dates exceed series length")}
  
  if(!IC%in%c("AIC","BIC") )
  {
    
    stop("IC must either be AIC or BIC")
    
  }
  
  MSFE <- c()
  predF <- NULL
  pvec <- NULL
  svec <- NULL
  k <- ncol(Y)
  m <- ifelse(s!=0,ncol(X),0)
  for(i in (T1-h+2):T2){
    
    if(h+i-1>T2){break}
    
    testY <- as.matrix(Y[1:(i-1),])
    testX <- as.matrix(X[1:(i-1),])
    
    if(!iterated){
      hd=h
      
    }else{
      
      hd=1
      
    }
    if(IC=="BIC"){
      popt <- ICX(testY,testX,k,p,s,m,"BIC",h=hd)
    }
    if(IC=="AIC"){
      
      popt <- ICX(testY,testX,k,p,s,m,"AIC",h=hd) 
    }
    B1 <- popt$B
    
    if(popt$p==0 &popt$s==0){
      
      eZ <- matrix(rep(1,1),ncol=1)
      
      pred <- B1%*%eZ
      
    }else{
      
      
      C <- max(popt$p,popt$s)
      
      ## print(C)
      ## if(C==1){
      
      ## # possibly memory leak in VARX lag matrix construction in Eigen if maxlag is 1.
      ## # to be on the safe side, we will perform it in R
      ## if(popt$s==0){
      ## eZ <- c(1,Y[i-1,])
      ## eZ <- matrix(eZ,ncol=1)
      ## }else{
      ##     eZ <- c(1,Y[i-1,],X[i-1,])               
      ##     }
      ## }else{
      eZ <- VARXCons(as.matrix(Y[(i-C):(i),]),as.matrix(X[(i-C):(i),]),k,popt$p,m,popt$s)
      
      ## }
      
      pred <- B1%*%eZ
      
      # iterated multistep forecasts (if VAR and horizon greater than 1)
      if(h>1 & s==0 & iterated ){
        
        pred <- predictMS(matrix(pred,nrow=1),Y,h-1,B1,C,FALSE)
        
      }
      
    }
    ## browser()
    predF <- rbind(predF,t(pred))
    MSFEi <- norm2(Y[i+h-1,]-pred)^2
    MSFE <- c(MSFE,MSFEi)
    svec <- c(popt$s,s)
    pvec <- c(popt$p,p)
  }
  
  return(list(MSFE=MSFE,pred=as.matrix(predF),p=pvec,s=svec))
  
  
}

#' Fit a VAR or VARX model by least squares
#' 
#' @param Y a \eqn{t \times k} multivariate time series
#' @param p maximum lag order
#' @param IC Information criterion indicator, if set to \code{NULL}, it will fit a least squares VAR(X) of orders p and s.  Otherwise, if set to "AIC" or "BIC" it return the model with lag orders that minimize the given IC. 
#' @param VARX a list of VARX specifications (as in \code{\link{constructModel}} (or NULL )
#' @return Returns a list with four entries:
#' \itemize{
#' \item{"Bhat"}{Estimated \eqn{k\times kp+ms} coefficient matrix}
#' \item{"SigmaU}{Estimated \eqn{k\times k} residual covariance matrix}
#' \item{"phat"}{Selected lag order for VAR component}
#' \item{"shat"}{Selected lag order for VARX component}
#' }
#' @details This function uses a modified form of the least squares technique proposed by Neumaier and Schneider (2001).  It fits a least squares VAR or VARX via a QR decomposition that does not require explicit matrix inversion.  This results in improved computational performance as well as numerical stability over the conventional least squares approach. 
#' @references Neumaier, Arnold, and Tapio Schneider. "Estimation of parameters and eigenmodes of multivariate autoregressive models." ACM Transactions on Mathematical Software (TOMS) 27.1 (2001): 27-57.
#' @seealso \code{\link{constructModel}}, \code{\link{cv.BigVAR}}
#' @examples
#' data(Y)
#' # fit a VAR_3(3)
#' mod <- VARXFit(Y,3,NULL,NULL)
#' # fit a VAR_3 with p= 6 and lag selected according to AIC
#' modAIC <- VARXFit(Y,6,"AIC",NULL)
#' # Fit a VARX_{2,1} with p=6, s=4 and lags selected by BIC
#' modXBIC <- VARXFit(Y,6,"BIC",list(k=1,s=4))
#' 
#' @export
VARXFit <- function(Y,p,IC,VARX=NULL)
{
  
  if(!is.null(VARX)){
    
    if(is.list(VARX) & !(exists('k',where=VARX) & exists('s',where=VARX)))
    {
      
      stop("VARX Specifications entered incorrectly")
      
    }
    
  }
  if(is.list(VARX) & (length(VARX)!=0)){
    
    k1 <- VARX$k
    s <- VARX$s
    Y1 <- matrix(Y[,1:k1],ncol=k1)
    m <- ncol(Y)-k1
    X <- matrix(Y[,(k1+1):ncol(Y)],ncol=m)
    Z <- VARXCons(Y1,X,k1,p,m,s)
    offset <- max(p,s)+1
    YT <- matrix(Y1[offset:nrow(Y),],ncol=k1)
    X <- matrix(X[offset:nrow(X),],ncol=m)
    
  }else{
    
    k <- ncol(Y)
    k1 <- k
    s=0;m=0
    offset <- p+1
    X <- matrix(0,nrow=nrow(Y))
    
    Z <- VARXCons(Y,X,k,p,m,s)
    YT <- matrix(Y[(offset):nrow(Y),],ncol=ncol(Y))
    
  }
  if(is.null(IC)){
    
    Res <- ARFitVARXR(cbind(t(Z),YT),k1,p,m,s)
    
    shat <- s
    phat <- p
    
  }else{
    if(!IC%in%c("AIC","BIC") )
    {
      
      stop("IC must either be AIC,BIC, or set to NULL")
      
    }
    
    Res <- ICX(YT,X,k1,p,s,m,IC)
    
    shat <- Res$s
    phat <- Res$p
    
  }
  
  list(Bhat=Res$B,SigmaU=Res$SigmaU,phat=phat,shat=shat)
  
}

# Recursive multi-step predictions


predictMS <- function(pred,Y,n.ahead,B,p,MN=FALSE){
  
  # Augment Y with predictions, create lag matrix (no intercept if MN)
  Y <- rbind(Y,pred)
  
  # Can't call this function directly from R due to assert errors
  ## Z <- ZmatF(Y,p,ncol(Y),oos=TRUE,intercept=!MN)
  
  Z <- VARXCons(Y,matrix(0,nrow=nrow(Y),ncol=1),ncol(Y),p,0,0,oos=TRUE)
  if(MN){
    Z <- Z[2:nrow(Z),]
  }
  Z <- Z[,ncol(Z)]
  
  pred <- matrix(B%*%Z,ncol=ncol(Y),nrow=1)
  
  if(n.ahead==1){return(pred)}
  
  predictMS(pred,Y,n.ahead-1,B,p,MN)
  
}

# Multi-step VARX with new data.
predictMSX <- function(pred,Y,n.ahead,B,p,newxreg,X,m,s,cumulative,MN){
  
  Y <- rbind(Y,pred)
  X <- rbind(X,matrix(newxreg[cumulative,],ncol=m))
  
  if(nrow(Y)!=nrow(X)){stop("error, dimension issue")}
  Z <- VARXCons(as.matrix(Y),X,ncol(Y),p,m,s,TRUE)
  Z <- Z[,ncol(Z)]
  if(MN){
    
    Z <- as.matrix(Z[2:nrow(Z),])
  }    
  pred <- matrix(B%*%Z,ncol=ncol(Y),nrow=1)
  
  if(n.ahead==1){return(pred)}
  
  predictMSX(pred,Y,n.ahead-1,B,p,newxreg,X,m,s,cumulative+1,MN)    
  
}




# Find optimal values in 2-d gridsearch
findind <- function(opt,lambda1,lambda2)
{
  if(opt<length(lambda2)){
    lambda1ind <- 1
  }else{
    lambda1ind <- ceiling(opt/length(lambda2))
  }
  if(lambda1ind==1)
  {
    jind <- opt
  }else{
    jind <- opt-(length(lambda2))*(lambda1ind-1)
  }
  return(c(lambda1ind,jind))
}




BVARLitterman <- function(Y,Z,p,tau,mu,H,iRW)
{
  T <- nrow(Y); k <- ncol(Y)
  # prior covariance based on univariate AR models
  sigmas <- c()
  for(i in 1:k){
    Z1 <- VARXCons(matrix(Y[,i],ncol=1),matrix(0,nrow=nrow(Y),ncol=1),1,p,0,0)
    ## Y1 <- matrix(Y[(p+1):nrow(Y),i],ncol=1)
    # get the prior cov
    K <- cbind(t(Z1),Y[(p+1):nrow(Y),i])
    sigmas[i] <- sqrt(ARFitV2(K,1,p)$SigmaU)    
    ## print(sigmas[i])
  }
  ## browser()
  MMO <- colMeans(Y)
  
  
  
  ## Z <- VARX(Y,matrix(0,nrow=Y,ncol=1),k,p,0,0)
  
  ## Y1 <- matrix(Y[(p+1):nrow(Y),],ncol=k)
  
  # create prior random walk dummy
  Yrw1 <- diag(sigmas*iRW)
  Yrw2 <- matrix(0,nrow=k*(p-1),ncol=k)
  Yrw<- tau*(rbind(Yrw1,Yrw2))
  Zrw <- tau*cbind(kronecker(diag(1:p),diag(sigmas)),matrix(0,nrow=k*p,ncol=1))
  
  
  # create dummy for intercept
  epsilon=1e-5
  Ycs <- 1e-5*matrix(0,nrow=1,ncol=k)
  Zcs <- epsilon*cbind(matrix(0,ncol=k*p,nrow=1),1)
  
  # dummy on the sums of coefficients
  Ylr <- mu*diag(MMO*iRW)
  Zlr1 <- kronecker(matrix(1,nrow=1,ncol=p),diag(MMO)*iRW)
  Zlr <- mu*(cbind(Zlr1,matrix(0,nrow=k,ncol=1)))
  
  
  # Dummy for residual covariance matrix
  Ycv <- diag(sigmas)
  Zcv <- matrix(0,nrow=k,ncol=k*p+1)
  
  Yprior <- rbind(Yrw,Ylr,Ycv,Ycs)
  Zprior <- rbind(Zrw,Zlr,Zcv,Zcs)
  
  ## Zp2 <- rev(Zprior[(nrow(Zprior)),])
  ## Zp3 <- Zprior[1:(nrow(Zprior)-1),]
  ## Zp3 <- rbind(Zp2,Zp3)
  ## ZZI2 <- solve(t(Zp3)%*%Zp3+t(Z)%*%Z)
  ## dim(ZZI2)
  Tstar <- nrow(Yprior)
  ## dim(Zprior)
  ## # posterior
  ## dim(Z)
  ## browser()
  Z <- t(Z)
  Z <- cbind(Z[,2:ncol(Z)],1)
  ## dim(Yprior)
  ## browser()        
  
  ## kappa(crossprod(Z))
  ZZinv <- solve(t(Zprior)%*%Zprior+t(Z)%*%Z)
  ZY <- t(Zprior)%*%Yprior+t(Z)%*%Y
  beta <- ZZinv%*%ZY
  ## browser()
  ## preds <- Z[nrow(Z),]%*%beta
  
  return(t(beta))
  
}


BGRGridSearch <- function(Y,Z,p,grid,RWIND)
{
  preds <- list()
  for(i in 1:length(grid))
  {
    pi <- grid[i]
    mu=pi*.1 # used in BGR paper
    preds[[i]] <- BVARLitterman(Y,Z,p,pi,mu,-1,RWIND)
    
  }
  preds <- array(unlist(preds), dim = c(nrow(preds[[1]]), ncol(preds[[1]]), length(preds)))
  
  return(preds)
  
}

# BigVAR Algorithms
# Most of the computationally expensive portions of the code have been exported to C++

# Sparse Own/Other (VAR)
.SparseGroupLassoVAROO<-function(beta,groups,Y,Z,gamm,alpha,INIactive,eps,q1a,p,MN,dual=FALSE,C1,intercept) 
{
  ## browser()
  
  
  m <- 0
  k <- ncol(Y)
  if(MN)
  {
    
    C <- matrix(0,nrow=k,ncol=k*p)        
    diag(C) <- C1
    Y <- t(Y)
    Y <- Y-C%*%Z
    Y <- t(Y)
  }
  
  Y <- t(Y)
  YOLD <- Y
  ZOLD <- Z
  if(intercept){
    YMean <- c(apply(Y, 1, mean))
    ZMean <- c(apply(Z, 1, mean))    
    Y <- Y - c(apply(Y, 1, mean)) %*% t(c(rep(1, ncol(Y))))
    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
  }else{
    YMean <- rep(0,nrow(Y))
    ZMean <- rep(0,nrow(Z))
  }
  ZZ <- kronecker(t(Z),diag(k))  
  M1f <- list()
  M2f <- list()
  jj <- .lfunction3(p,k)
  eigs <- c()
  q1 <- list()
  # get step size from inverse of max eigenvalue via power method
  for (j in 1:length(jj)) {
    M1f[[j]] <- ZZ[,jj[[j]] ]
    M2f[[j]] <- crossprod(ZZ[,jj[[j]]])
    gg1 <- powermethod(M2f[[j]], q1a[[j]])
    eigs[j] <- gg1$lambda
    q1[[j]] <- gg1$q1
  }    
  jj <- .lfunction3cpp(p,k)
  jjfull <- jj
  jjcomp <- .lfunctioncomp(p, k)
  dims <- dim(beta)
  beta <- array(beta[,2:ncol(beta[,,1]),],dim=c(dims[1],dims[2]-1,dims[3]))
  if(!dual){
    
    BB <- GamLoopSGLOO(beta,INIactive,gamm,alpha,Y,ZZ,jj,jj,jjcomp,eps,YMean,ZMean,k,p*k,M2f,eigs,m)
    
    
  }else{
    
    ## browser()
    BB <- GamLoopSGLOODP(beta,INIactive,gamm,alpha,Y,ZZ,jj,jj,jjcomp,eps,YMean,ZMean,k,p*k,M2f,eigs,m)
    
  }
  
  BB$q1 <- q1
  
  if (MN)
  {
    for(i in 1:(dim(BB$beta)[3]))
      
    {
      
      BB$beta[,2:ncol(BB$beta[,,i]),i] <- BB$beta[,2:ncol(BB$beta[,,i]),i]+C
    }
  }
  return(BB)
}



# Sparse Lag (VAR)
.SparseGroupLassoVAR<-function(beta,groups,Y,Z,gamm,alpha,INIactive,eps,q1a,p,MN,C1,intercept) 
{
  k <- ncol(Y)
  
  if (MN)
  {
    
    C <- matrix(0,nrow=k,ncol=k*p)        
    diag(C) <- C1
    Y <- t(Y)
    Y <- Y-C%*%Z
    Y <- t(Y)
  }
  
  Y <- t(Y)
  YOLD <- Y
  ZOLD <- Z
  if(intercept){
    YMean <- c(apply(Y, 1, mean))
    ZMean <- c(apply(Z, 1, mean))    
    Y <- Y - c(apply(Y, 1, mean)) %*% t(c(rep(1, ncol(Y))))
    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
  }else{
    YMean <- rep(0,nrow(Y))
    ZMean <- rep(0,nrow(Y))
  }
  M1f <- list()
  M2f <- list()
  eigs <- c()
  q1 <- list()
  jj <- .groupfun(p,k)
  
  
  # get step size from inverse of max eigenvalue via power method
  for (j in 1:length(jj)) {
    M1f[[j]] <- Z[jj[[j]], ]
    M2f[[j]] <- Z[jj[[j]], ] %*% t(Z[jj[[j]], ])
    gg1 <- powermethod(M2f[[j]], q1a[[j]])
    eigs[j] <- gg1$lambda
    q1[[j]] <- gg1$q1
  }
  
  jj <- .groupfuncpp(p,k)
  jjfull <- jj
  jjcomp <- .groupfuncomp(p, k)
  dims <- dim(beta)
  beta <- array(beta[,2:ncol(beta[,,1]),],dim=c(dims[1],dims[2]-1,dims[3]))
  ## browser()    
  ## beta <- beta[,2:ncol(beta[,,1]),]
  
  BB <- GamLoopSGL(beta,INIactive,gamm,alpha,Y,Z,jj,jjfull,jjcomp,eps,YMean,ZMean,k,p*k,M1f,M2f,eigs)
  BB$q1 <- q1
  
  if (MN)
  {
    for(i in 1:(dim(BB$beta)[3]))
    {
      
      BB$beta[,2:ncol(BB$beta[,,i]),i] <- BB$beta[,2:ncol(BB$beta[,,i]),i]+C
    }
    
  }
  
  return(BB)
}

# Sparse Lag (VAR) Dual Search
.SparseGroupLassoVARDual<-function(beta,groups,Y,Z,gamm,alpha,INIactive,eps,q1a,p,MN,C1,intercept) 
{
  k <- ncol(Y)
  if (MN)
  {
    
    C <- matrix(0,nrow=k,ncol=k*p)
    diag(C) <- C1
    ## diag(C) <- rep(1,k)
    Y <- t(Y)
    Y <- Y-C%*%Z
    Y <- t(Y)
  }
  
  Y <- t(Y)
  YOLD <- Y
  ZOLD <- Z
  if(intercept){
    YMean <- c(apply(Y, 1, mean))
    ZMean <- c(apply(Z, 1, mean))    
    Y <- Y - c(apply(Y, 1, mean)) %*% t(c(rep(1, ncol(Y))))
    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
  }else{
    YMean <- rep(0,nrow(Y))
    ZMean <- rep(0,nrow(Z))
  }
  M1f <- list()
  M2f <- list()
  eigs <- c()
  q1 <- list()
  jj <- .groupfun(p,k)
  # get step size from inverse of max eigenvalue via power method
  for (j in 1:length(jj)) {
    M1f[[j]] <- Z[jj[[j]], ]
    M2f[[j]] <- Z[jj[[j]], ] %*% t(Z[jj[[j]], ])
    gg1 <- powermethod(M2f[[j]], q1a[[j]])
    eigs[j] <- gg1$lambda
    q1[[j]] <- gg1$q1
  }
  
  jj <- .groupfuncpp(p,k)
  jjfull <- jj
  jjcomp <- .groupfuncomp(p, k)
  ngp <- length(alpha)*length(gamm)
  dims <- dim(beta)
  beta <- array(beta[,2:ncol(beta[,,1]),],dim=c(dims[1],dims[2]-1,dims[3]))
  ## browser()
  BB <- GamLoopSGLDP(beta,INIactive,gamm,alpha,Y,Z,jj,jjfull,jjcomp,eps,YMean,ZMean,k,p*k,M1f,M2f,eigs)
  BB$q1 <- q1
  
  if (MN)
  {
    for(i in 1:(dim(BB$beta)[3]))
    {
      
      BB$beta[,2:ncol(BB$beta[,,i]),i] <- BB$beta[,2:ncol(BB$beta[,,i]),i]+C
    }
    
  }
  
  return(BB)
}


# Lag Group (VAR/VARX-L)
.GroupLassoVAR1 <- 
  function (beta, groups,jjcomp, Y, Z, gamm, INIactive, eps,p,MN,k,k1,s,C1,intercept) 
  {
    
    
    
    if(!"matrix"%in%class(Y))
    {
      Y <- matrix(Y,ncol=1)
    }
    
    m <- k-k1
    
    if(MN)
    {
      C <- matrix(0,nrow=k1,ncol=k1*p+s*m)        
      ## diag(C) <- rep(1,k1)
      diag(C) <- C1
      Y <- t(Y)
      Y <- Y-C%*%Z
      Y <- t(Y)
    }
    
    Y <- t(Y)
    YOLD <- Y
    ZOLD <- Z
    if(intercept){
      YMean <- c(apply(Y, 1, mean))
      ZMean <- c(apply(Z, 1, mean))    
      Y <- Y - c(apply(Y, 1, mean)) %*% t(c(rep(1, ncol(Y))))
      Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
    }else{
      YMean <- rep(0,nrow(Y))
      ZMean <- rep(0,nrow(Z))
    }
    
    if(k>1){    
      
      Y <- Y - c(apply(Y, 1, mean)) %*% t(c(rep(1, ncol(Y))))}
    
    else{Y <- Y-mean(Y)}    
    
    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
    jj <- groups
    jjfull <- groups
    
    Eigsys <- Eigencomp(Z,jj,length(jj),k1)
    M2 <- Eigsys$M3
    eigvals <- Eigsys$eigval
    eigvecs <- Eigsys$eigvec
    
    beta <- array(beta[,2:ncol(as.matrix(beta[,,1])),],dim=c(k1,(k1)*p+s*m,length(gamm)))
    BB <- GamLoopGL2(beta,INIactive,gamm,Y,Z,jj,jjfull,jjcomp,eps,YMean,ZMean,k1,p*k1+m*s,M2,eigvals,eigvecs)
    
    if(MN)
    {
      for(i in 1:(dim(BB$beta)[3]))
      {
        
        BB$beta[,2:ncol(BB$beta[,,i]),i] <- BB$beta[,2:ncol(BB$beta[,,i]),i]+C
      }
      
    }
    return(BB)
  }



# Group Lasso Own/Other (VARXL)
.GroupLassoOOX <- function (beta, groups, Y, Z, gamm, INIactive, eps,p,MN,k,k1,s,C1,intercept) 
{
  m <- k-k1
  
  if(MN)
  {
    C <- matrix(0,nrow=k1,ncol=k1*p+s*m)        
    ## diag(C) <- rep(1,k1)
    diag(C) <- C1
    Y <- t(Y)
    Y <- Y-C%*%Z
    Y <- t(Y)
  }
  
  
  Y <- t(Y)
  YOLD <- Y
  ZOLD <- Z
  
  if(intercept){
    if(k1>1){
      
      YMEAN <- apply(Y, 1, mean)
      
    }else{
      YMEAN <- mean(Y)
    }         
    
    ZMEAN <- apply(Z,1,mean)
  }else{
    
    YMEAN <- rep(0,nrow(Y))
    ZMEAN <- rep(0,nrow(Z))
    
  }
  if(intercept){
    if(k>1){    
      Y <- Y - c(apply(Y, 1, mean)) %*% t(c(rep(1, ncol(Y))))
    }else{Y <- Y-mean(Y)}    
    
    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
  }
  kk=diaggroupfunVARX(p,k,k1,s)
  kkcomp <- diaggroupfunVARXcomp(p,k,k1,s)
  
  ZZ <- kronecker(t(Z), diag(k1))
  
  Eigsys <- EigencompOO(ZZ,kk,length(kk),k1)
  M2 <- Eigsys$M3
  eigvals <- Eigsys$eigval
  eigvecs <- Eigsys$eigvec
  
  kkfull <- kk
  
  beta <- array(beta[,2:ncol(as.matrix(beta[,,1])),],dim=c(k1,(k1)*p+m*s,length(gamm)))
  
  BB <- GamLoopGLOO(beta,INIactive,gamm,Y,ZZ,kk,kkfull,kkcomp,eps,YMEAN,ZMEAN,k1,p*(k1)+m*s,M2,eigvals,eigvecs,k1)
  
  if(MN)
  {
    for(i in 1:(dim(BB$beta)[3]))
    {
      BB$beta[,2:ncol(BB$beta[,,i]),i] <- BB$beta[,2:ncol(BB$beta[,,i]),i]+C
    }
  }
  
  return(BB)
}


# Own/Other Group VAR-L
.GroupLassoOO <- function (beta, groups, Y, Z, gamm, INIactive, eps,p,MN,C1,intercept) 
{
  
  if(!"matrix"%in%class(Y))
  {
    Y <- matrix(Y,ncol=1)
  }
  
  k <- ncol(Y)
  
  if(MN)
  {
    C <- matrix(0,nrow=k,ncol=k*p)        
    ## diag(C) <- rep(1,k)
    diag(C) <- C1
    Y <- t(Y)
    Y <- Y-C%*%Z
    Y <- t(Y)
  }
  
  
  Y <- t(Y)
  YOLD <- Y
  ZOLD <- Z
  
  if(intercept){
    if(k>1){
      
      YMEAN <- apply(Y, 1, mean)
      
    }else{
      YMEAN <- mean(Y)
    }         
    
    ZMEAN <- apply(Z,1,mean)
  }else{
    
    YMEAN <- rep(0,nrow(Y))
    ZMEAN <- rep(0,nrow(Z))
    
  }
  if(intercept){
    if(k>1){    
      Y <- Y - c(apply(Y, 1, mean)) %*% t(c(rep(1, ncol(Y))))
    }else{Y <- Y-mean(Y)}    
    
    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
  }
  kk <- groups
  kkfull <- groups
  kkcomp <- .lfunctioncomp(p, k)
  
  ZZ <- kronecker(t(Z), diag(k))
  
  Eigsys <- EigencompOO(ZZ,kk,length(kk),k)
  M2 <- Eigsys$M3
  eigvals <- Eigsys$eigval
  eigvecs <- Eigsys$eigvec
  
  beta <- array(beta[,2:ncol(as.matrix(beta[,,1])),],dim=c(k,k*p,length(gamm)))
  
  
  BB <- GamLoopGLOO(beta,INIactive,gamm,Y,ZZ,kk,kkfull,kkcomp,eps,YMEAN,ZMEAN,k,p*k,M2,eigvals,eigvecs,k)
  
  if(MN)
  {
    for(i in 1:(dim(BB$beta)[3]))
    {
      BB$beta[,2:ncol(BB$beta[,,i]),i] <- BB$beta[,2:ncol(BB$beta[,,i]),i]+C
    }
    
  }
  
  return(BB)
}



# Elementwise HVAR
.HVARElemAlg <- function (beta, Y, Z, lambda, eps,p,MN,C1,intercept) 
{
  
  k <- ncol(Y)
  if(MN)
  {
    C <- matrix(0,nrow=k,ncol=k*p)        
    ## diag(C) <- rep(1,k)
    diag(C) <- C1
    Y <- t(Y)
    Y <- Y-C%*%Z
    Y <- t(Y)
  }
  
  k <- ncol(Y)
  betafin <- beta
  YOLD <- Y
  ZOLD <- Z
  if(intercept){
    YMean <- c(apply(Y, 2, mean))
    ZMean <- c(apply(Z, 1, mean))
    
    Y <- Y - matrix(c(rep(1, nrow(Y))), ncol = 1) %*% matrix(c(apply(Y, 
                                                                     2, mean)), nrow = 1)
    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
  }else{
    YMean <- rep(0,ncol(Y))
    ZMean <- rep(0,nrow(Z))
  }
  
  tk <- 1/max(Mod(eigen(Z %*% t(Z))$values))
  
  betaini <- array(beta[, 2:ncol(beta[, , 1]), ],dim=c(k,k*p,length(lambda)))
  betafin <- gamloopElem(betaini,Y,Z,lambda,eps,YMean,ZMean,as.matrix(betaini[,,1]),k,p)
  
  if(MN)
  {
    for(i in 1:(dim(betafin)[3]))
    {
      betafin[,2:ncol(betafin[,,i]),i] <- betafin[,2:ncol(betafin[,,i]),i]+C
    }
    
  }
  
  return(betafin)
}

# Basic-VARX-L Fista Implementation
.lassoVARFistX <- function (B, Z, Y, gamm, eps,p,MN,k,k1,s,m,C1,intercept) 
{
  
  if(!"matrix"%in%class(Y))
  {
    Y <- matrix(Y,ncol=1)
  }
  
  
  if(MN)
  {
    C <- matrix(0,nrow=k1,ncol=k1*p+s*m)        
    diag(C) <- C1
    Y <- t(Y)
    Y <- Y-C%*%Z
    Y <- t(Y)
  }
  
  Y <- t(Y)
  if(intercept){
    YMean <- c(apply(Y, 1, mean))
    ZMean <- c(apply(Z, 1, mean))
    Y <- Y - YMean %*% t(c(rep(1, ncol(Y))))
    Z <- Z - ZMean %*% t(c(rep(1, ncol(Z))))
  }else{
    YMean <- rep(0,nrow(Y))
    ZMean <- rep(0,nrow(Z))
  }
  Y <- t(Y)
  
  tk <- 1/max(Mod(eigen(Z%*%t(Z),only.values=TRUE)$values))
  
  BFOO1 <- matrix(B[, 2:dim(B)[2], 1],ncol=k1)
  BFOO <- array(B[,2:ncol(as.matrix(B[,,1])),],dim=c(k1,k1*p+(k-k1)*s,length(gamm)))
  beta <- gamloopFista(BFOO, Y, Z, gamm, eps, 
                       as.matrix(YMean), as.matrix(ZMean), BFOO1,k,p,tk,k1,s)
  
  if(MN)
  {
    for(i in 1:(dim(beta)[3]))             
      beta[,2:ncol(beta[,,i]),i] <- beta[,2:ncol(beta[,,i]),i]+C
  }
  
  return(beta)
  
}



#Basic VAR Fista Implementation
.lassoVARFist <- function (B, Z, Y, gamm, eps,p,MN,C1,intercept) 
{
  
  if(!"matrix"%in%class(Y))
  {
    Y <- matrix(Y,ncol=1)
  }
  
  k <- ncol(Y)
  
  if(MN)
  {
    C <- matrix(0,nrow=k,ncol=k*p)        
    diag(C) <- C1
    Y <- t(Y)
    Y <- Y-C%*%Z
    Y <- t(Y)
  }
  
  Y <- t(Y)
  if(intercept){
    YMean <- c(apply(Y, 1, mean))
    ZMean <- c(apply(Z, 1, mean))
    Y <- Y - YMean %*% t(c(rep(1, ncol(Y))))
    Z <- Z - ZMean %*% t(c(rep(1, ncol(Z))))
  }else{
    YMean <- rep(0,nrow(Y))
    ZMean <- rep(0,nrow(Z))
  }
  Y <- t(Y)
  
  tk <- 1/max(Mod(eigen(Z%*%t(Z),only.values=TRUE)$values))
  
  BFOO1 <- matrix(B[, 2:dim(B)[2], 1],nrow=k,ncol=k*p)
  BFOO <- array(B[,2:ncol(as.matrix(B[,,1])),],dim=c(k,k*p,length(gamm)))
  
  beta <- gamloopFista(BFOO, Y, Z, gamm, eps, 
                       as.matrix(YMean), as.matrix(ZMean), BFOO1,k,p,tk,k,p)
  
  if(MN)
  {
    for(i in 1:(dim(beta)[3]))             
      beta[,2:ncol(beta[,,i]),i] <- beta[,2:ncol(beta[,,i]),i]+C
  }
  return(beta)
}

# Componentwise HVAR
.HVARCAlg <- function(beta,Y,Z,lambda,eps,p,MN,C1,intercept)
{
  if(!"matrix"%in%class(Y))
  {
    Y <- matrix(Y,ncol=1)
  }
  
  k <- ncol(Y)
  
  if(MN)
  {
    C <- matrix(0,nrow=k,ncol=k*p)        
    ## diag(C) <- rep(1,k)
    diag(C) <- C1
    Y <- t(Y)
    Y <- Y-C%*%Z
    Y <- t(Y)
  }
  
  
  betafin <- beta
  YOLD <- Y
  ZOLD <- Z
  if(intercept){
    YMean <- c(apply(Y, 2, mean))
    ZMean <- c(apply(Z, 1, mean))
  }else{
    YMean <- rep(0,ncol(Y))
    ZMean <- rep(0,nrow(Z))
  }
  if(intercept){
    if(k>1){
      Y <- Y -  matrix(c(rep(1, nrow(Y))),ncol=1)%*%matrix(c(apply(Y, 2, mean)),nrow=1)
    }else{
      Y <- Y-mean(Y)
    }
    
    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
  }
  tk <- 1/max(Mod(eigen(Z%*%t(Z))$values))
  
  betaini <- array(beta[,2:ncol(as.matrix(beta[,,1])),],dim=c(k,k*p,length(lambda)))
  betafin <- gamloopHVAR(betaini,Y,Z,lambda,eps,YMean,ZMean,as.matrix(betaini[,,1]),k,p)
  
  if(MN)
  {
    for(i in 1:(dim(betafin)[3]))
    {
      betafin[,2:ncol(betafin[,,i]),i] <- betafin[,2:ncol(betafin[,,i]),i]+C
    }
  }    
  return(betafin)
}

# Endogenous First VARX-L
.EFVARX <- function(beta,Y,Z,lambda,eps,MN,k1,s,m,p,C1,intercept)
{        
  
  
  
  if(MN)
  {
    
    
    C <- matrix(0,nrow=k1,ncol=k1*p+s*m)        
    ## diag(C) <- rep(1,k1)
    diag(C) <- C1
    Y <- t(Y)
    Y <- Y-C%*%Z
    Y <- t(Y)
  }
  
  betafin <- beta  
  Y <- t(Y)
  YOLD <- Y
  ZOLD <- Z
  if(intercept){  
    YMEAN <- apply(Y, 1, mean)
    ZMEAN <- apply(Z, 1, mean)
    Y <- Y - c(apply(Y, 1, mean)) %*% t(c(rep(1, ncol(Y))))
    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
  }else{
    YMEAN <- rep(0,nrow(Y))
    ZMEAN <- rep(0,nrow(Z))
  }
  Y <- t(Y)
  tk <- 1/max(Mod(eigen(Z%*%t(Z))$values))
  
  
  if (k1==1)
  {
    
    betaini <- array(beta[,2:(k1*p+m*s+1),],dim=c(1,k1*p+m*s,dim(beta)[3]))
    
  }else{
    
    betaini <- beta[,2:ncol(beta[,,1]),]
    
  }
  
  for(i in 1:length(lambda))
  {
    
    if(dim(beta)[3]>1){
      
      betaF <- fistaX(Y,Z,matrix(betaini[,,i],nrow=k1),p,k1,lambda[i],eps,tk,m,s)
    }
    
    else{ 
      
      betaF <- fistaX(Y,Z,matrix(betaini,nrow=k1),p,k1,lambda[i],eps,tk,m,s)
      
    }
    
    nu <- YMEAN - betaF %*% ZMEAN
    
    betafin[,,i] <- cbind(nu, betaF)
    
  }
  
  if(MN)
  {
    
    for(i in 1:(dim(betafin)[3]))
    {
      
      betafin[,2:ncol(betafin[,,i]),i] <- betafin[,2:ncol(betafin[,,i]),i]+C
      
    }
    
  }
  
  return(betafin)
}

# HVAR Own/Other
.HVAROOAlg <- function(beta,Y,Z,lambda,eps,p,MN,C1,intercept)
{
  
  k <- ncol(Y)
  betafin <- beta
  YOLD <- Y
  ZOLD <- Z
  
  if(MN)
  {
    C <- matrix(0,nrow=k,ncol=k*p)        
    ## diag(C) <- rep(1,k)
    diag(C) <- C1
    Y <- t(Y)
    Y <- Y-C%*%Z
    Y <- t(Y)
  }
  
  weights <- sqrt(c(rep(c(1,k-1),length=2*p)))
  if(intercept){
    YMEAN <- c(apply(Y, 2, mean))
    ZMEAN <- c(apply(Z, 1, mean))
    Y <- Y -  matrix(c(rep(1, nrow(Y))),ncol=1)%*%matrix(c(apply(Y, 2, mean)),nrow=1)
    Z <- Z -  c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
  }else{
    YMEAN <- rep(0,ncol(Y))
    ZMEAN <- rep(0,nrow(Z))
  }
  groups <- list()
  for(i in 1:k)
  {
    
    groups[[i]] <- .vecoovarscpp(p,k,i)
    
  }
  
  betaini <- array(beta[,2:ncol(as.matrix(beta[,,1])),],dim=c(k,k*p,length(lambda)))
  
  betafin <- gamloopOO(betaini,Y,Z,lambda,eps,YMEAN,ZMEAN,as.matrix(betaini[,,1]),k,p,weights,groups)
  
  
  if(MN)
  {
    for(i in 1:(dim(betafin)[3]))
    {
      
      betafin[,2:ncol(betafin[,,i]),i] <- betafin[,2:ncol(betafin[,,i]),i]+C
      
    }
    
  }
  
  
  return(betafin)
  
}


# indexing for efx
vxsubs <- function(i,k,m,p,s)
{
  
  
  vv <- c(((i-1)*k+1):((i-1)*k+k),((i-1)*m+k*p+1):((i-1)*m+k*p+m))
  
  vv
  
}

prox2 <- function(v,lambda,k,p,m,s)
{
  
  for(i in 1:p)
  {
    
    if(i<=s){
      
      vv <- vxsubs(i,k,m,p,s)
      F1 <- 0
      
    }
    
    if(i>s){
      vv <- ((i-1)*k+1):((i-1)*k+k)
      F1 <- 1
    }
    
    v2 <- proxvx2(v[vv],p,lambda,m,k,F1)
    v[vv] <- v2
    
    
  }
  
  v
  
}


fistaX <- function(Y,Z,beta,p,k1,lambda,eps,tk,m,s)
{
  
  for(i in 1:k1)
  {
    
    phiOLD <- beta[i,]
    phiOLDOLD <- beta[i,]
    j <- 1
    thresh <- 10*eps
    
    while(thresh>eps)
    {
      
      v <- matrix(phiOLD+((j-2)/(j+1))*(phiOLD-phiOLDOLD),nrow=1)
      phiR <- prox2(v+tk*as.vector((Y[,i]-v%*%Z)%*%t(Z)),tk*lambda,k1,p,m,s)
      thresh <- max(abs(phiR-v))
      phiOLDOLD <- phiOLD
      phiOLD <- phiR
      j <- j+1
    }
    
    beta[i,] <- phiR
    
  }
  
  beta
}

# Sparse Lag VARX-L
.SparseGroupLassoVARX<-function(beta,groups,Y,Z,gamm,alpha,INIactive,eps,q1a,p,MN,k,s,k1,C1,intercept) 
{
  
  m <- k-k1
  
  if(MN)
  {
    C <- matrix(0,nrow=k1,ncol=k1*p+s*m)        
    diag(C) <- C1
    Y <- t(Y)
    Y <- Y-C%*%Z
    Y <- t(Y)
  }
  
  Y <- t(Y)
  YOLD <- Y
  ZOLD <- Z
  if(intercept){
    YMean <- c(apply(Y, 1, mean))
    ZMean <- c(apply(Z, 1, mean))
    Y <- Y - c(apply(Y, 1, mean)) %*% t(c(rep(1, ncol(Y))))
    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
  }else{
    YMean <- rep(0,nrow(Y))
    ZMean <- rep(0,nrow(Z))
  }
  M1f <- list()
  M2f <- list()
  eigs <- c()
  q1 <- list()
  jj <- groupfunVARXLG(p,k,k1,s)
  
  for (j in 1:length(jj)) {
    M2f[[j]] <- Z[jj[[j]], ] %*% t(Z[jj[[j]], ])
    
    if(j<=p){
      
      gg1 <- powermethod(M2f[[j]], q1a[[j]])
      eigs[j] <- gg1$lambda
      q1[[j]] <- gg1$q1
    }
    
    else{
      
      M2f[[j]] <- as.vector(Z[jj[[j]], ]) %*% as.vector(t(Z[jj[[j]], ]))
      eigs[j] <- M2f[[j]]
      
    }
    
  }
  
  
  jj <- groupfunVARX(p,k,k1,s)    
  jjfull <- jj
  jjcomp <- groupfunVARXcomp(p,k,k1,s)
  
  beta <- array(beta[,2:ncol(as.matrix(beta[,,1])),],dim=c(k1,(k1)*p+(s*m),length(gamm)))
  BB <- GamLoopSGLX(beta,INIactive,gamm,alpha,Y,Z,jj,jjfull,jjcomp,eps,YMean,ZMean,k1,(k1)*p+(s*m),M2f,eigs,k1)
  
  BB$q1 <- q1
  
  if(MN)
  {
    
    for(i in 1:(dim(BB$beta)[3]))
    {
      
      BB$beta[,2:ncol(BB$beta[,,i]),i] <- BB$beta[,2:ncol(BB$beta[,,i]),i]+C
      
    }
    
  }
  
  return(BB)
}



.SparseGroupLassoVARXDual<-function(beta,groups,Y,Z,gamm,alpha,INIactive,eps,q1a,p,MN,k,s,k1,C1,intercept) 
{
  
  m <- k-k1
  
  if(MN)
  {
    C <- matrix(0,nrow=k1,ncol=k1*p+s*m)        
    ## diag(C) <- rep(1,k1)
    diag(C) <- C1
    Y <- t(Y)
    Y <- Y-C%*%Z
    Y <- t(Y)
  }
  
  Y <- t(Y)
  YOLD <- Y
  ZOLD <- Z
  if(intercept){
    YMean <- c(apply(Y, 1, mean))
    ZMean <- c(apply(Z, 1, mean))
    Y <- Y - c(apply(Y, 1, mean)) %*% t(c(rep(1, ncol(Y))))
    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
  }else{
    YMean <- rep(0,nrow(Y))
    ZMean <- rep(0,nrow(Z))
  }
  M1f <- list()
  M2f <- list()
  eigs <- c()
  q1 <- list()
  jj <- groupfunVARXLG(p,k,k1,s)
  
  for (j in 1:length(jj)) {
    
    M2f[[j]] <- Z[jj[[j]], ] %*% t(Z[jj[[j]], ])
    
    if(j<=p){
      
      gg1 <- powermethod(M2f[[j]], q1a[[j]])
      eigs[j] <- gg1$lambda
      q1[[j]] <- gg1$q1
    }
    
    else{
      
      M2f[[j]] <- as.vector(Z[jj[[j]], ]) %*% as.vector(t(Z[jj[[j]], ]))
      eigs[j] <- M2f[[j]]
      
    }
    
  }
  
  
  jj <- groupfunVARX(p,k,k1,s)    
  jjfull <- jj
  jjcomp <- groupfunVARXcomp(p,k,k1,s)
  
  beta <- array(beta[,2:ncol(as.matrix(beta[,,1])),],dim=c(k1,(k1)*p+(s*m),nrow(gamm)*length(alpha)))
  
  BB <- GamLoopSGLXDP(beta,INIactive,gamm,alpha,Y,Z,jj,jjfull,jjcomp,eps,YMean,ZMean,k1,(k1)*p+(s*m),M2f,eigs,k1)
  
  BB$q1 <- q1
  
  if(MN)
  {
    
    for(i in 1:(dim(BB$beta)[3]))
    {
      
      BB$beta[,2:ncol(BB$beta[,,i]),i] <- BB$beta[,2:ncol(BB$beta[,,i]),i]+C
      
    }
    
  }
  
  return(BB)
}


# Sparse Own/Other (VARX)
.SparseGroupLassoVAROOX<-function(beta,groups,Y,Z,gamm,alpha,INIactive,eps,p,MN,k1,s,k,dual=FALSE,C1,intercept) 
{
  
  m <- k-k1
  if(MN)
  {
    
    C <- matrix(0,nrow=k1,ncol=k1*p+s*m)        
    diag(C) <- C1
    Y <- t(Y)
    Y <- Y-C%*%Z
    Y <- t(Y)
  }
  
  Y <- t(Y)
  YOLD <- Y
  ZOLD <- Z
  if(intercept){
    YMean <- c(apply(Y, 1, mean))
    ZMean <- c(apply(Z, 1, mean))
    Y <- Y - c(apply(Y, 1, mean)) %*% t(c(rep(1, ncol(Y))))
    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
  }else{
    YMean <- rep(0,nrow(Y))
    ZMean <- rep(0,nrow(Z))
  }
  M1f <- list()
  M2f <- list()
  eigs <- c()
  q1 <- list()
  
  jj <- diaggroupfunVARXLG(p,k,k1,s)
  kk <- diaggroupfunVARX(p,k,k1,s)
  jjcomp <- diaggroupfunVARXcomp(p,k,k1,s)
  
  ZZ <- kronecker(t(Z),diag(k1))  
  
  
  for (j in 1:length(jj)) {
    
    M2f[[j]] <- crossprod(ZZ[,jj[[j]]])
    
    eigs[j] <- max(Mod(eigen(M2f[[j]],only.values=TRUE)$values))
  }
  
  jjfull <- kk
  
  if(!dual){
    
    gran2 <- length(gamm)
    
  }else{
    
    gran2 <- nrow(gamm)*ncol(gamm)
    
  }
  
  beta <- array(beta[,2:ncol(as.matrix(beta[,,1])),],dim=c(k1,k1*p+m*s,gran2))
  
  if(dual){
    BB <- GamLoopSGLOODP(beta,INIactive,gamm,alpha,Y,ZZ,kk,jjfull,jjcomp,eps,YMean,ZMean,k1,p*k1+m*s,M2f,eigs,m)
  }else{
    
    BB <- GamLoopSGLOO(beta,INIactive,gamm,alpha,Y,ZZ,kk,jjfull,jjcomp,eps,YMean,ZMean,k1,p*k1+m*s,M2f,eigs,m)
    
  }
  if(MN)
  {
    for(i in 1:(dim(BB$beta)[3]))
    {
      
      BB$beta[,2:ncol(BB$beta[,,i]),i] <- BB$beta[,2:ncol(BB$beta[,,i]),i]+C
      
    }
    
  }
  
  
  return(BB)
}



# Lag weighted lasso: VAR only

.lassoVARTL <- function (B, Z, Y, gamm, eps,p,MN,alpha,C1,intercept) 
{
  
  if(!"matrix"%in%class(Y))
  {
    Y <- matrix(Y,ncol=1)
  }
  
  k <- ncol(Y)
  
  if(MN)
  {
    C <- matrix(0,nrow=k,ncol=k*p)        
    diag(C) <- C1
    Y <- t(Y)
    Y <- Y-C%*%Z
    Y <- t(Y)
    
  }
  
  Y <- t(Y)
  if(intercept){
    YMean <- c(apply(Y, 1, mean))
    ZMean <- c(apply(Z, 1, mean))
    Y <- Y - YMean %*% t(c(rep(1, ncol(Y))))
    Z <- Z - ZMean %*% t(c(rep(1, ncol(Z))))
  }else{
    YMean <- rep(0,nrow(Y))
    ZMean <- rep(0,nrow(Z))
  }
  Y <- t(Y)
  tk <- 1/max(Mod(eigen(Z%*%t(Z),only.values=TRUE)$values))
  BFOO1 <- as.matrix(B[, 2:ncol(B[, , 1]), 1])
  BFOO <- array(B[,2:ncol(as.matrix(B[,,1])),],dim=c(k,k*p,length(gamm)*length(alpha)))
  p2 <- 1:p
  
  gran2 <- length(gamm)
  
  for(i in 1:length(alpha))
  {
    
    W <- rep(p2^(-alpha[i]),each=k)
    ZADJ <- diag(W)%*%Z
    
    B[,,(1+(i-1)*gran2):(i*length(gamm))] <- gamloopFista(array(BFOO[,,(1+(i-1)*gran2):(i*length(gamm))],dim=c(k,k*p,length(gamm))), Y, ZADJ, gamm, eps,as.matrix(YMean), as.matrix(ZMean), BFOO1,k,p,tk,k,p)
    
    
    for(j in (1+(i-1)*gran2):(i*length(gamm)))
    {
      
      B[,2:(k*p+1),j] <- B[,2:(k*p+1),j]%*%diag(W)
      
    }               
    
    if(MN)
    {
      for(i in 1:(dim(B)[3]))
        
        B[,2:ncol(B[,,i]),i] <- B[,2:ncol(B[,,i]),i]+C
      
    }
    
  }
  
  return(B)
}

SimBoot <- function (k, A1, p, Resid, T) 
{
  
  if(max(Mod(eigen(A1)$values))>1){stop("Error: Generator Matrix is not stationary")}
  
  
  # add 500 observations for initialization purposes
  Resid <- scale(Resid, center=TRUE, scale=FALSE)
  
  Y <- matrix(0, nrow = T+500+p , ncol = k)
  YY <- as.vector(Y)
  
  for (i in seq(from = (k * p + 1), to = (nrow(Y) * k - 1), 
                by = k)) {
    
    bootresid <- sample(1:nrow(Resid), 1)
    u <- as.vector(c(Resid[bootresid,], rep(0, 
                                            k * p - k)))
    
    YY[(i + k):(i - k * p + 1 + k)] <- A1 %*% YY[(i):(i - 
                                                        k * p + 1)] + as.matrix(u, ncol = 1)
    
  }
  
  YY <- YY[YY!=0]
  
  Y <- matrix(YY, ncol = k, byrow = TRUE)
  
  Y <- Y[,c(ncol(Y):1)]
  
  Y <- Y[501:nrow(Y), ]
  
  return(Y)
}
