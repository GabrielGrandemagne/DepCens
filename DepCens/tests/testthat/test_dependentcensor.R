library(DepCens)

# MEP Approach
fit <- dependent.censoring(formula = time ~ x1 | x1, data=KidneyMimic, delta_t=KidneyMimic$delta_t,
                           delta_c=KidneyMimic$delta_c, ident=KidneyMimic$ident, dist = "mep")

# Weibull Approach
fit <- dependent.censoring(formula = time ~ x1 + x2 | x1 + x2, data=KidneyMimic, delta_t=KidneyMimic$delta_t,
                           delta_c=KidneyMimic$delta_c, ident=KidneyMimic$ident, dist = "weibull")
