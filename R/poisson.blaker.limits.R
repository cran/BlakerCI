poisson.blaker.limits <- function(x,level=.95,tol=1e-10,...) {
  if (x < 0) stop("Parameter x = ",x, " wrong!")
  if (level <= 0 || level >= 1) stop("Confidence level ",level," out of (0, 1)!")
  if (tol <= 0) stop("Numerical tolerance ",tol," nonpositive!")
  lower <- poisson.blaker.lower.limit(x,level,tol,...)
  upper <- poisson.blaker.upper.limit(x,level,tol,...)
  return(c(lower,upper))
}
