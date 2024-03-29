% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary_dc.R
\name{summary_dc}
\alias{summary_dc}
\title{Print the summary output}
\usage{
summary_dc(object, ...)
}
\arguments{
\item{object}{an object of the class "dcensoring".}

\item{...}{further arguments passed to or from other methods.}
}
\value{
a summary of the fitted model.
}
\description{
Print the summary output
}
\examples{
\donttest{
fit <- dependent.censoring(formula = time ~ x1 | x3, data=KidneyMimic, delta_t=KidneyMimic$delta_t,
                          delta_c=KidneyMimic$delta_c, ident=KidneyMimic$ident, dist = "mep")
summary_dc(fit)
}
}
