#---------------------------------------------

#' Print the summary output
#'
#' @export
#' @param object an object of the class "dcensoring".
#' @param ... further arguments passed to or from other methods.
#' @return a summary of the fitted model.
#'
#' @examples
#' \dontrun{
#' fit <- dependent.censoring(formula = time ~ x1 | x3, data=KidneyMimic, delta_t=KidneyMimic$delta_t,
#'                           delta_c=KidneyMimic$delta_c, ident=KidneyMimic$ident, approach = "mep")
#' summary_dc(fit1)
#'}
summary_dc <- function(object, ...){

  bmax <- object$bmax

  if (is.null(bmax) == FALSE){
    p <- object$p
    q <- object$q
    Xlabels <- object$labels2
    Zlabels <- object$labels1

    cat("MEP APPROACH:")
    cat("\n")
    cat("\n")
    cat("Name\t","Estimate\t","Std. Error\t","CI INF \t","CI SUP\t","\n")
    cat("Alpha", format(object[[1]][(p+bmax+q+1)],nsmall=6), format(object[[2]][(p+bmax+q+1)], nsmall=6), format((object[[1]][(p+bmax+q+1)] - 1.96*object[[1]][(p+bmax+q+1)]),nsmall=6), format((object[[1]][(p+bmax+q+1)] + 1.96*object[[1]][(p+bmax+q+1)]), nsmall=6), sep = "\t", "\n" )
    cat("Sigma", format(object[[1]][(p+bmax+q+bmax+2)],nsmall=6), format(object[[1]][(p+bmax+q+bmax+2)], nsmall=6), format((object[[1]][(p+bmax+q+bmax+2)] - 1.96*object[[1]][(p+bmax+q+bmax+2)]),nsmall=6), format((object[[1]][(p+bmax+q+bmax+2)] + 1.96*object[[1]][(p+bmax+q+bmax+2)]),nsmall=6), sep = "\t", "\n")
    cat("\n")
    cat("Coefficients T:\n")
    cat("\n")
    cat("Name\t","Estimate\t","Std. Error\t","CI INF \t","CI SUP\t","\n")
    for (i in 1:p){
      cat(substr(sub(".*\\$","", Xlabels[i]), 1, 6), format(object[[1]][c(i)], nsmall=6),format(object[[1]][c(i)],nsmall=6),format((object[[1]][c(i)] - 1.96*object[[1]][c(i)]),nsmall=6),format((object[[1]][c(i)] + 1.96*object[[1]][c(i)]),nsmall=6),sep="\t","\n")
    }
    cat("\n")
    cat("Coefficients C:\n")
    cat("\n")
    cat("Name\t","Estimate\t","Std. Error\t","CI INF \t","CI SUP\t","\n")
    for (j in 1:q){
      cat(substr(sub(".*\\$","", Zlabels[j]), 1, 6), format(object[[1]][c(p+bmax+j)],nsmall=6),format(object[[1]][c(p+bmax+j)],nsmall=6),format((object[[1]][c(p+bmax+j)] - 1.96*object[[1]][c(p+bmax+j)]),nsmall=6),format((object[[1]][c(p+bmax+j)] + 1.96*object[[1]][c(p+bmax+j)]),nsmall=6),sep="\t","\n")
    }
  }
  else{
    p <- object$p
    q <- object$q
    Xlabels <- object$labels2
    Zlabels <- object$labels1

    cat("WEIBULL APPROACH:")
    cat("\n")
    cat("\n")
    cat("Name\t","Estimate\t","Std. Error\t","CI INF \t","CI SUP\t","\n")
    cat("Alpha",format(object[[1]][(5+p+q)],nsmall=6), format(object[[2]][(5+p+q)],nsmall=6), format((object[[1]][c(5+p+q)] - 1.96*object[[2]][c(5+p+q)]),nsmall=6),format((object[[1]][c(5+p+q)] + 1.96*object[[2]][c(5+p+q)]),nsmall=6), sep = "\t", "\n" )
    cat("Sigma",format(object[[1]][(6+p+q)],nsmall=6), format(object[[2]][(6+p+q)],nsmall=6), format((object[[1]][c(6+p+q)] - 1.96*object[[2]][c(6+p+q)]),nsmall=6),format((object[[1]][c(6+p+q)] + 1.96*object[[2]][c(6+p+q)]),nsmall=6), sep = "\t", "\n" )
    cat("\n")
    cat("Coefficients T:\n")
    cat("\n")
    cat("Name\t","Estimate\t","Std. Error\t","CI INF \t","CI SUP\t","\n")
    for (i in 1:p){
      cat(substr(sub(".*\\$","", Xlabels[i]), 1, 6),format(object[[1]][c(2+i)],nsmall=6),format(object[[2]][c(2+i)],nsmall=6),format((object[[1]][c(2+i)] - 1.96*object[[2]][c(2+i)]),nsmall=6),format((object[[1]][c(2+i)] + 1.96*object[[2]][c(2+i)]),nsmall=6), sep = "\t", "\n")
    }
    cat("\n")
    cat("Coefficients C:\n")
    cat("\n")
    cat("Name\t","Estimate\t","Std. Error\t","CI INF \t","CI SUP\t","\n")
    for (j in 1:q){
      cat(substr(sub(".*\\$","", Zlabels[j]), 1, 6),format(object[[1]][c(4+p+j)],nsmall=6), format(object[[2]][c(4+p+j)],nsmall=6),format((object[[1]][c(4+p+j)] - 1.96*object[[2]][c(4+p+j)]),nsmall=6),format((object[[1]][c(4+p+j)] + 1.96*object[[2]][c(4+p+j)]),nsmall=6), sep = "\t","\n")
    }
  }
}
