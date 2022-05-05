#---------------------------------------------

#' Dependent Censoring model
#' @aliases dependent.censoring
#' @export
#' @description dependent.censoring is used to fit a survival analysis with dependent censoring.
#' @param formula an object of class "formula" (or one that can be coerced to that class): should be used as 'time ~ failure covariates | informative covariates'.
#' @param data a data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
#' @param delta_t Indicator function of the event of interest.
#' @param delta_c Indicator function of the dependent censoring.
#' @param ident Cluster indicator variable.
#' @param Num_intervals Number of intervals of the time grid (mep only).
#' @param approach approach to be used to fit the model (must be either weibull or mep).
#' @details function estimates the parameters of the Piecewise exponential model (approach = "mep") or Weibull model (approach="weibull") with dependent censoring, considering the frailty model to capture the dependence between failure and dependent censoring times, as well as the clusters variability.
#' @return dependent.censoring returns an object of class "dcensoring" containing the fitted model.
#' An object of class "dcensoring" is a list containing at least the following components:
#' \itemize{
#'   \item \code{param_est} a vector containing estimated parameters (dependency parameter, regression coefficients associated with failure times, regression coefficients associated with dependent censoring times, and time distribution parameters (Weibull or piecewise exponential)).
#'   \item \code{stde} a vector containing the estimated standard errors of the estimated parameters vector.
#'   \item \code{crit} a vector containing the information criteria, Akaike's information criterion (AIC), Baysian information criterion (BIC), Hannan–Quinn information criterion (HQ).
#'   \item \code{pvalue} p-value of the estimated parameters vector.
#'   \item \code{n} number of observations in the dataset.
#'   \item \code{p} number of covariates associated with failure times (event of interest times).
#'   \item \code{q} number of covariates associated with dependent censoring times (informative censoring times or competitive risk times).
#'   \item \code{formula} formula used in the function call.
#'   \item \code{terms} the terms object used, containing the covariates associated with the failure times and with the dependent censoring times.
#'   \item \code{labels1} labels of the covariates associated with failure times.
#'   \item \code{labels2} labels of the covariates associated with dependent censoring times.
#'   \item \code{risco_a_T} a vector containing the cumulative baseline hazar of failure times.
#'   \item \code{risco_a_C} a vector containing the cumulative baseline hazar of dependent censoring times.
#'   \item \code{bi} a matrix containing the generated fragilities, one of the outputs of the function dependent.censoring, in which the individuals are in the rows and the Monte Carlo replicas in the columns.
#'   \item \code{X_T} a matrix of variables associated with failure times.
#'   \item \code{X_C} a matrix of variables associated with dependent censoring times.
#'   \item \code{time} a vector of the observable times.
#' }
#' @examples
#' \dontrun{
#' library(DepCens)
#' delta_t <- ifelse(KidneyMimic$cens==1,1,0)
#' delta_c <- ifelse(KidneyMimic$cens==2,1,0)
#' fit <- dependent.censoring(formula = time ~ x1 | x3, data=KidneyMimic, delta_t=delta_t,
#'                           delta_c=delta_c, ident=KidneyMimic$ident, approach = "mep")
#' summary_dc(fit)
#'}
dependent.censoring <- function(formula, data, delta_t, delta_c, ident, approach = c("weibull", "mep"), Num_intervals = 10){

  approach <- match.arg(approach)

  switch(approach,
         "weibull" = model_Weibull_dep(formula=formula, data=data, delta_t=delta_t, delta_c=delta_c, ident=ident),
         "mep" = model_MEP_dep(formula=formula, data=data, delta_t=delta_t, delta_c=delta_c, ident=ident, Num_intervals = 10))
}
