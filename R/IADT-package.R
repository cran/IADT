#' Interaction Difference Test for Prediction Models
#' 
#' Provides functions to conduct a model-agnostic asymptotic 
#' hypothesis test for the identification of interaction effects in 
#' black-box machine learning models. The null hypothesis assumes that a given 
#' set of covariates does not contribute to interaction effects in the prediction 
#' model. The test statistic is based on the difference of variances of 
#' partial dependence functions with respect to the original 
#' black-box predictions and the predictions under the null hypothesis. 
#' The hypothesis test can be applied to any black-box prediction model, 
#' and the null hypothesis of the test can be flexibly specified according 
#' to the research question of interest. Furthermore, the test is computationally 
#' fast to apply as the null distribution does not require resampling or 
#' refitting black-box prediction models.
#' 
#' \tabular{ll}{ Package: \tab IADT \cr Type: \tab Package \cr Version:
#' \tab 1.2.1 \cr Date: \tab 2024-05-14 \cr License: \tab GPL-3 \cr }
#' 
#' @name IADT-package
#' @aliases IADT-package
#' @docType package
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @references
#' \insertRef{welchowskiIntroIML}{IADT} \cr\cr
#' \insertRef{friedmanInteract}{IADT}
#' @keywords package
NULL

# Namespace code
#' @importFrom Rmpfr mpfr mean pnorm
#' @importFrom methods new
#' @importFrom parallel makeCluster clusterExport clusterEvalQ parLapplyLB stopCluster
#' @importFrom stats predict sd
#' @importFrom mgcv gam
#' @importFrom Rdpack reprompt
#' @importFrom utils sessionInfo
#' @importFrom mvnfast rmvn
NULL
