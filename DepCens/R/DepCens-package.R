#' The 'DepCens' package.
#'
#' @description Dependent censoring regression models for survival multivariate data. These models are based on extensions of the frailty models, capable to accommodating the dependence between failure and censoring times, with Weibull and piecewise exponential marginal distributions. Theoretical details regarding the models implemented in the package can be found in Schneider et al. (2019) <doi: 10.1002/bimj.201800391>.
#' @docType package
#' @name DepCens-package
#' @aliases DepCens
#' @import Formula
#' @import survival
#' @importFrom rootSolve multiroot
#' @importFrom dlm arms
#' @importFrom matrixStats colProds
#' @importFrom stats pnorm lowess
#' @importFrom graphics lines
#'
#' @references
#' Schneider, S.; Demarqui, F. N.; Colosimo, E. A.; Mayrink, V. D. (2020). An approach to model clustered survival data with dependent censoring. Biometrical Journal, v.62, n.1, 157--174.
#'
#' Louis, T. A. (1982). Finding the observed information matrix when using the EM algorithm. Journal of the Royal Statistical Society, B44, 226–233.
#'
NULL
