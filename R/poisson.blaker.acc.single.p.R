poisson.blaker.acc.single.p <- function(x,p,type="orig",acc.tol=1e-10,output="acc",maxiter=100) {

# "Ordinary" acceptability at p.
if (p <= x) {
  p1.p <- ppois(x-1,p,lower.tail=FALSE)
  q1.p <- qpois(p1.p,p)-1
  a1.p <- min(p1.p + ppois(q1.p,p),1)
} else {
  p1.p <- ppois(x,p,lower.tail=TRUE)
  q1.p <- qpois(1-p1.p,p)+1
  a1.p <- min(p1.p + ppois(q1.p-1,p,lower.tail=FALSE),1)
}

# "Unimodalization"
  if (type == "unimod" && q1.p >= 0) { 
    if (p <= x) { 
      upper <- p
      a1.upp <- a1.p
      lower <- 0
      a1.low <- 1   
#
      iter <- 0
#   In 1.0-4, ... >= acc.tol replaced with ... > acc.tol
      while (a1.low > a1.p && a1.low - a1.upp > acc.tol) { 
        iter <- iter + 1
        if (iter > maxiter) {  
          warning("Convergence not attained after ",maxiter, 
                             " iterations for x = ",x,", p = ",p,
                             ", and acceptability tolerance limit of ",acc.tol)
          break
        } 
        mid <- (lower+upper)/2
        p1.mid <- ppois(x-1,mid,lower.tail=FALSE)
        p2.mid <- ppois(q1.p,mid)
        a1.mid <- p1.mid + p2.mid
        if (p1.mid < p2.mid) { 
          lower <- mid
          a1.low <- a1.mid
        } 
        else { 
          upper <- mid
          a1.upp <- a1.mid
        } 
      } 
    } else {   
## Unimodalisation for p > x
      lower <- p
      a1.low <- a1.p
      upper <- q1.p 
      a1.upp <- 1   

      iter <- 0
#   In 1.0-4, ... 
      while (a1.upp > a1.p && a1.upp - a1.low > acc.tol) { 
        iter <- iter + 1
        if (iter > maxiter) {  
          warning("Convergence not attained after ",maxiter, 
                             " iterations for x = ",x,", p = ",p,
                             ", and acceptability tolerance limit of ",acc.tol)
          break
        } 
        mid <- (lower+upper)/2
        p1.mid <- ppois(x,mid,lower.tail=TRUE)
        p2.mid <- ppois(q1.p-1,mid,lower.tail=FALSE)
        a1.mid <- p1.mid + p2.mid
        if (p1.mid < p2.mid) { 
          upper <- mid
          a1.upp <- a1.mid
        } 
        else { 
          lower <- mid
          a1.low <- a1.mid
        } 
      } 
 
    } 
    a1.p <- max(a1.p,a1.low)
  } 


  if (output == "acc") {
    return(a1.p)
  }
  else if (output == "both") {
    return(c(a1.p,q1.p))
  }
  else if (output == "q1") {
    return(q1.p)
  }


}
