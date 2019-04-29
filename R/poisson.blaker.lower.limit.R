poisson.blaker.lower.limit <- function(x,level,tol=1e-10,maxiter=100) {
  if (x <=0) return(0)
  if (x>0) {
    alpha <- 1-level
#   Clopper-Pearson limit (CPL)
    lower <- qgamma(alpha/2,x,1)
    p1 <- ppois(x-1,lower,lower.tail=FALSE)
    q1.cp <- qpois(p1,lower)-1
    upper <- x
    iter <- 0
    while (upper-lower >= tol) {
      iter <- iter+1
      if (iter > maxiter) {
        warning("Tolerance limit of ",tol, 
                             " not attained after ",maxiter, 
                             " iterations for x = ",x)
        break
      }
      mid <- (lower+upper)/2
      p1 <- ppois(x-1,mid,lower.tail=FALSE)
#   Blaker's limit is below the midpoint if either
#   (i)  acceptability at mid > alpha, or
#   (ii) acceptability function has a discontinuity between
#        the midpoint and CPL (first test).
      if (p1 >= ppois(q1.cp+1,mid) || p1 + ppois(q1.cp,mid) > alpha) {
        upper <- mid
      }
      else {
        lower <- mid
      }
    }
  return(lower)
  }
}  
