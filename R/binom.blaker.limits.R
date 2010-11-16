binom.blaker.limits <- function(x,n,level=.95,tol=1e-10) {
  lower <- binom.blaker.lower.limit(x,n,level,tol)
  upper <- 1 - binom.blaker.lower.limit(n-x,n,level,tol)
  return(c(lower,upper))
}
