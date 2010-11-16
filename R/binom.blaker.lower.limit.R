binom.blaker.lower.limit <- function(x,n,level,tol=1e-10) {
  if (x <=0) return(0)
  if (x>0) {
    alpha <- 1-level
#   Clopper-Pearson limit (CPL)
    lower <- qbeta(alpha/2,x,n-x+1)
    p1 <- 1 - pbinom(x-1,n,lower)
    q1.cp <- qbinom(p1,n,lower)-1
    a1 <- p1 + pbinom(q1.cp,n,lower)
#   Blaker's limit = CPL?
    if (a1 >= alpha) {
      return(lower)
    }
    upper <- x/n
    while (upper-lower >= tol) {
      mid <- (lower+upper)/2
      p1 <- 1 - pbinom(x-1,n,mid)
      q1 <- qbinom(p1,n,mid)-1
#   Blaker's limit is below the midpoint if either
#   (i)  acceptance at mid >= alpha, or
#   (ii) acceptance function has a discontinuity between
#        the midpoint and CPL (test based on q1).
      if (q1 > q1.cp || p1 + pbinom(q1,n,mid) >= alpha) {
        upper <- mid
      }
      else {
        lower <- mid
      }
    }
  return(lower)
  }
}  
