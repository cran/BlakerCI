poisson.blaker.upper.limit <- function(x,level,tol=1e-10,maxiter=100) {
    alpha <- 1-level
#   Clopper-Pearson limit (CPL)
    upper <- qgamma(1-alpha/2,x+1,1)
    p1 <- ppois(x,upper,lower.tail=TRUE)
    q1.cp <- qpois(1-p1,upper)+1
    lower <- x
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
      p1 <- ppois(x,mid,lower.tail=TRUE)
#   Blaker's limit is below the midpoint if either
#   (i)  acceptability at mid > alpha (NEW!! orig: >=), or
#   (ii) acceptability function has a discontinuity between
#        the midpoint and CPL (first test).
      if (p1 >= ppois(q1.cp-2,mid,lower.tail=FALSE) || p1 + ppois(q1.cp-1,mid,lower.tail=FALSE) > alpha) {
        lower <- mid
      }
      else {
        upper <- mid
      }
    }
  return(upper)

}  
