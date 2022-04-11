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
#'
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
