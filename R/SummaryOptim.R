###################################################################################################
# Function: summary.optim
#
# Make summary table of model fit from optim(): Estimate, Std.Error, z.value, p.value.
#
#
# ARGUMENTS:
#
# optim.fit: Fitted model object from optim().
#
# VALUE:
#
# ds.out: Table with 4 columns: cbind(Estimate, Std.Error, z.value, p.value).
#
# Author: DVN
# Senturk D, Dalrymple DS, Mu Y, Nguyen DV (2014) Weighted hurdle regression method for joint
# modeling of cardiovascular events likelihood and rate in the U.S. dialysis population.
#
# 2.8.14
###################################################################################################
summary.optim <- function ( optim.fit, N=100 ){

  Estimate = optim.fit$par
  #cov.mat = solve(optim.fit$hessian)
  cov.mat= 2*optim.fit$value/(N - length(Estimate)) * solve(optim.fit$hessian)
  Std.Error = sqrt( diag(cov.mat)	)
  z.value = Estimate/Std.Error
  p.value = 2*( 1-pnorm(abs(z.value)) )

  ds.out = cbind(Estimate, Std.Error, z.value, p.value)
  colnames(ds.out) <- c('Estimate', 'Std Error', 'z value', 'p value')

  ds.out
}
###################################################################################################


#library(numDeriv)
#' @import numDeriv
summary.optim2 <- function ( optim.fit,f, N=100 ){

  Estimate = optim.fit$par
  #cov.mat = solve(optim.fit$hessian)
  cov.mat= 2*optim.fit$value/(N - length(Estimate)) * solve(hessian(f,optim.fit$par))
  Std.Error = sqrt( diag(cov.mat)	)
  z.value = Estimate/Std.Error
  p.value = 2*( 1-pnorm(abs(z.value)) )

  ds.out = cbind(Estimate, Std.Error, z.value, p.value)
  colnames(ds.out) <- c('Estimate', 'Std Error', 'z value', 'p value')

  ds.out
}
