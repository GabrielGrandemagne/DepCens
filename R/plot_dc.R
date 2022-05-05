#' Plot the survival function
#'
#' @aliases plot_dc
#' @export
#' @description this graph helps to visualize the survival function.
#' @param object an object of the class "dcensoring".
#' @param scenario which defines the scenario in the graph (t: failure times, c: dependent censoring times, or both).
#' @details sometimes the smoothness treatment used (lowess function) results in a non-monotonous survival function.
#' @return a survival function graph of the fitted model.
#'
#' @examples
#' \dontrun{
#' fit <- dependent.censoring(formula = time ~ x1 | x3, data=KidneyMimic, delta_t=KidneyMimic$delta_t,
#'                           delta_c=KidneyMimic$delta_c, ident=KidneyMimic$ident, approach = "mep")
#' plot.dc(fit, scenario = "both")
#'}
#'
plot_dc <- function(object, scenario  = c("t", "c", "both")){

  scenario  <- match.arg(scenario )

  bmax <- object$bmax
  #Caso MEP
  if (is.null(bmax) == FALSE){
    switch(scenario ,
           "t" = invisible(plot.mep.t(object)),
           "c" = invisible(plot.mep.c(object)),
           "both"= invisible(plot.mep.t(object) + plot.mep.c(object)))
  }

  #Caso Weibull
  else{
    switch(scenario ,
           "t" = invisible(plot.weibull.t(object)),
           "c" = invisible(plot.weibull.c(object)),
           "both"= invisible(plot.weibull.t(object) + plot.weibull.c(object)))
  }
}
