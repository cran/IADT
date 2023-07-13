#' Partial dependence plot with specific numerical precision
#' 
#' Estimates the partial dependence plot (PDP) curve given specified numerical 
#' precision.
#' 
#' @param colInd Index of columns of covariates to specify the null hypothesis 
#' set \emph{s} (integer vector). 
#' @param object Prediction model object (class flexible). 
#' @param predictfun Prediction function to be evaluated (class function). The
#' prediction function needs to be specified with two arguments 
#' \emph{predictfun(object, X)}. The argument \emph{object} is the prediction 
#' model and \emph{X} the data on which the partial dependence functions
#' are evaluated. 
#' @param X Data on that the partial dependence function is evaluated 
#' (class matrix or data.frame). The structure of the data depends on the 
#' specified argument \emph{predictfun}.
#' @param centering Should the resulting values be mean centered? 
#' (logical scalar). Default corresponds to output original values.
#' @param nCores Number of cores used in standard parallel computation
#' setup based on R \emph{parallel} package. The default value of one uses 
#' serial processing across observations.
#' @param precBits Numerical precision that are used in computation after the 
#' calculation of the predictions from the estimated model. Default is defined 
#' to be double the amount of the 53 Bits usually used in R.
#' @return Vector of estimated the PDP curve values for each sample in \emph{X}. 
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @seealso \code{\link{testIAD_mpfr}}
#' @references 
#' \insertRef{friedmanInteract}{IADT}
#' @keywords prediction
#' @examples
#' 
#' #####################
#' # Simulation example
#' 
#' # Simulate covariates from multivariate standard normal distribution
#' set.seed(-72498)
#' library(mvnfast)
#' X <- mvnfast::rmvn(n=1e2, mu=rep(0, 2), sigma=diag(2))
#' 
#' # Response generation
#' y <- X[, 1]^2 + rnorm(n=1e2, mean=0, sd=0.5)
#' trainDat <- data.frame(X, y=y)
#' 
#' # Estimate generalized additive model
#' library(mgcv)
#' gamFit <- gam(formula=y~s(X1)+s(X2), data=trainDat, 
#' family=gaussian())
#' 
#' # Estimate PDP function
#' pdpEst1 <- pdpEst_mpfr(colInd=1, object=gamFit, 
#' predictfun=function(object, X){
#' predict(object=object, newdata=X, type="response")
#' }, X=trainDat, 
#' centering=FALSE, nCores=1, precBits=53*2)
#' 
#' # Convert to standard precision and order in sequence of observations
#' pdpEst1 <- as.numeric(pdpEst1)
#' ordInd <- order(X[, 1])
#' pdpEst1 <- pdpEst1[ordInd]
#' 
#' # Plot: PDP curve vs. true effect
#' plot(x=X[ordInd, 1], y=pdpEst1, type="l")
#' lines(x=X[ordInd, 1], y=X[ordInd, 1]^2, lty=2, col="red")
#' # -> Both curves are similiar
#' 
#' @export pdpEst_mpfr
pdpEst_mpfr <- function(colInd, object, predictfun, X, 
                        centering=FALSE, nCores=1, precBits=53*2){
  
  if( all(colInd==0) ){
    
    preds <- mpfr(predictfun(object=object, X=X), precBits=precBits)
    
    if(centering){
      return(0)
    } else{
      return(mean(preds))
    }
  }
  
  if( nCores==1 ){
    
    nObs <- dim(X)[1]

    # Predictions
    preds <- vector("list", nObs)
    Xtemp <- X
    for( i in 1:nObs ){

      for( k in 1:length(colInd) ){
        Xtemp[, colInd[k] ] <- X[i, colInd[k] ]
      }

      preds[[i]] <- mean(mpfr(predictfun(object=object, X=Xtemp), 
                              precBits = precBits))
      
    }
    
    preds <- new("mpfr", unlist(preds))
    
  } else{
    
    # Definitions of parallel environment
    packageToLoad <- names(sessionInfo()$otherPkgs)
    clust0 <- makeCluster(nCores)
    clusterExport(clust0, varlist=c("predictfun", "object", "X",
                                    "packageToLoad", "colInd"),
                  envir=environment())
    invisible(clusterEvalQ(clust0, expr={
      for( i in 1:length(packageToLoad) ){
        library(packageToLoad[i], character.only = TRUE)
      }
    }))
    
    # Evaluation function
    evalFunc <- function(i){
      Xtemp <- X
      
      for( k in 1:length(colInd) ){
        Xtemp[, colInd[k] ] <- X[i, colInd[k] ]
      }
  
      pdp1 <- mean(mpfr(predictfun(object=object, X=Xtemp), 
                        precBits = precBits))
      rm(Xtemp)
      gc()
      return(pdp1)
    }
    
    # Run parallel processing
    preds <- parLapplyLB(cl=clust0, X=1:dim(X)[1], fun=evalFunc)
    stopCluster(clust0)
    preds <- unlist(preds)
  }
  
  # Aggregation
  if(centering){
      return(preds-mean(preds))
  } else{
    return(preds)
  }
  
}
