#' Model-agnostic interaction difference test for prediction models 
#' 
#' Tests if a given set \emph{s} of covariates contributes to interaction
#' effects in the prediction model.
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
#' @param nCoresPDP Number of cores used in standard parallel computation
#' setup based on R \emph{parallel} package to compute partial 
#' dependence function. 
#' The default value of one uses serial processing across observations.
#' @param precBits Numerical precision that are used in computation after the 
#' calculation of the predictions from the estimated model. Default is defined 
#' to be double the amount of the 53 Bits usually used in R.
#' @return Vector of estimated the PDP curve values for each sample in \emph{X}. 
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @details The data set used to evaluate the interaction test coulbe be training 
#' or test data. The proof of the asymptotic distribution only works in case of
#' test data. Therefore we recommend usage of test data for evaluation.
#' @note Numerical precision does have a big influence on Type I error
#' percentages. Under the null hypothesis values of the test statistic can be
#' very small near zero. To lessen the impact of random rounding errors increased
#' numerical precision is important.
#' @seealso \code{\link{pdpEst_mpfr}}
#' @keywords prediction
#' @examples
#' 
#' #####################
#' # Simulation example
#' #####################
#' 
#' library(mvnfast)
#' library(mgcv)
#' 
#' ############################################################################
#' # H0: Covariate X1 does not contribute to interaction effects in the 
#' # prediction model.
#' 
#' # Train data
#' # Simulate covariates from multivariate standard normal distribution
#' set.seed(-72498)
#' nSize <- 100
#' X <- mvnfast::rmvn(n=nSize, mu=rep(0, 2), sigma=diag(2))
#' 
#' # Response generation
#' y <- X[, 1]^2 + rnorm(n=nSize, mean=0, sd=0.25)
#' 
#' # Complete data frame
#' trainDat <- data.frame(X, y=y)
#' 
#' # Test data
#' # Simulate covariates from multivariate standard normal distribution
#' testX <- mvnfast::rmvn(n=nSize, mu=rep(0, 2), sigma=diag(2))
#' 
#' # Response generation
#' testY <- testX[, 1]^2 + rnorm(n=nSize, mean=0, sd=0.25)
#' testDat <- data.frame(testX, y=testY)
#' 
#' # Estimate generalized additive model with training data
#' library(mgcv)
#' gamFit <- gam(formula=y~s(X1)+s(X2), data=trainDat,
#'               family=gaussian())
#' 
#' # Test interaction
#' testIAD_mpfr1 <- testIAD_mpfr(colInd=1, object=gamFit,
#'                               predictfun=function(object, X){
#'                                 predict(object=object, newdata=X, type="response")
#'                               }, X=testDat)
#' testIAD_mpfr1$pValue
#' # -> H0 is not rejected with alpha=0.05
#' 
#' ###############################################################
#' # H1: X1 does contribute to at least one interaction effect
#' # in the prediction model.
#' 
#' # Train data
#' # Simulate covariates from multivariate standard normal distribution
#' set.seed(-72498)
#' nSize <- 150
#' X <- mvnfast::rmvn(n=nSize, mu=rep(0, 2), sigma=diag(2))
#' 
#' # Response generation
#' y <- X[, 1]^2 * X[, 2] + rnorm(n=nSize, mean=0, sd=0.25)
#' 
#' # Complete data frame
#' trainDat <- data.frame(X, y=y)
#' 
#' # Test data
#' # Simulate covariates from multivariate standard normal distribution
#' testX <- mvnfast::rmvn(n=nSize, mu=rep(0, 2), sigma=diag(2))
#' 
#' # Response generation
#' testY <- testX[, 1]^2 * testX[, 2] + rnorm(n=nSize, mean=0, sd=0.25)
#' testDat <- data.frame(testX, y=testY)
#' 
#' # Estimate generalized additive model with training data
#' library(mgcv)
#' gamFit <- gam(formula=y~s(X1, X2, k=25), data=trainDat,
#'               family=gaussian())
#' 
#' # Test interaction
#' testIAD_mpfr1 <- testIAD_mpfr(colInd=1, object=gamFit,
#'                               predictfun=function(object, X){
#'                                 predict(object=object, newdata=X, type="response")
#'                               }, X=testDat)
#' testIAD_mpfr1$pValue
#' # -> H0 is rejected with alpha=0.05
#' 
#' @export testIAD_mpfr
testIAD_mpfr <- function(colInd, object, predictfun, X, centering=FALSE, 
                         nCoresPDP=1, precBits=53*2){
  
  # 1. Calculate predictions
  IAD_f <- mpfr(predictfun(object=object, X=X), precBits=precBits)

  # 2. PD function S \ s
  if( all(colInd==0) ){
    calcSet <- 1:dim(X)[2]
  } else{
    calcSet <- setdiff(1:dim(X)[2], colInd)
  }
  pdpSetDiff <- pdpEst_mpfr(colInd=calcSet, 
                            object=object, X=X, 
                            predictfun=predictfun, centering=centering,
                            nCores=nCoresPDP,
                            precBits = precBits)
  
  # 3. Sum of PDP functions over each variable
  sumPDP <- rep(0, dim(X)[1])
  for( k in 1:length(colInd) ){
    
    sumPDP <- sumPDP + pdpEst_mpfr(colInd=colInd[k], object=object, X=X,
                                   predictfun=predictfun, centering=centering,
                                   nCores=nCoresPDP, precBits = precBits)
    
  }
  
  # 4. Sum of 2. and 3.
  IAD_PD <- pdpSetDiff + sumPDP
  
  # 5. Statistical hypothesis test
  x1 <- IAD_f + IAD_PD
  x2 <- IAD_f - IAD_PD
  z3 <- (x1-mean(x1))*(x2-mean(x2))
  z3_mean <- mean(z3)
  z3_sd <- sd(z3)
  
  # 6. Check zero variance
  if( !identical(x=z3_sd, y=0) ){
    
    testStat <- z3_mean / z3_sd * sqrt(length(z3))
    pVal <- 2*(1-pnorm(abs(testStat)))
    outputPrep <- list(testStat=testStat,
                       pValue=pVal,
                       IAD_f=IAD_f,
                       IAD_PD=IAD_PD,
                       z3=z3)
    
  } else{
    outputPrep <- list(testStat=0,
                       pValue=1,
                       IAD_f=IAD_f,
                       IAD_PD=IAD_PD,
                       z3=z3)
  }
  
  # Output
  return(outputPrep) 
}
