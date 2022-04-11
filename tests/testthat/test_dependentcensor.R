library(DepCens)

# MEP Approach
fit <- dependent.censoring(formula = time ~ x1 | x3, data=KidneyMimic, delta_t=KidneyMimic$delta_t,
                           delta_c=KidneyMimic$delta_c, ident=KidneyMimic$ident, approach = "mep")

# Weibull Approach
fit <- dependent.censoring(formula = time ~ x1 | x3 + x1, data=KidneyMimic, delta_t=KidneyMimic$delta_t,
                           delta_c=KidneyMimic$delta_c, ident=KidneyMimic$ident, approach = "weibull")