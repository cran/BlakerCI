poisson.blaker.acc <- function(x,p,type=c("orig","unimod"),acc.tol=1e-10,...) {
  if (x < 0) stop("Parameter x = ",x, " wrong!")
  if (acc.tol <= 0) stop("Numerical tolerance ",acc.tol," nonpositive!")
  type <- match.arg(type)
  if (type != "orig" && type != "unimod") stop("Unknown type ",type,"!")
  m <- length(p)
  if (m < 1) stop("Empty vector p!")
  if (m < 2) {
    acc <- poisson.blaker.acc.single.p(x,p,type=type,acc.tol=acc.tol,...)
    return(acc)
  }
  else {
    if (max(p[2:m]-p[1:(m-1)]) <= 0) stop("Vector p not increasing!")
#   First, regardless of type ("orig"/"unimod"),
#   calculate "ordinary" acceptabilities.
    aq <- sapply(p,poisson.blaker.acc.single.p,x=x,acc.tol=acc.tol,output="both",...=...)
    acc <- aq[1,]
    if (type == "orig") {
      return(acc)
    }
    if (type == "unimod") {
      q1 <- aq[2,]
      p.hat <- x
      ind1 <- p <= p.hat
#     In each interval of continuity of the acceptability function,
#     "highlight" the leftmost (below x) or rightmost (above x)
#     p element.
      ind <- (ind1 & (q1 > c(-Inf,q1[1:(m-1)]))) |
             (!ind1 & (q1 < c(q1[2:m],Inf)))
#     Calculate the "unimodalized" version of the acceptability function
#     just at the "highlighted" points, and leave the rest
#     to cummax().
#     (The amount of slow iterative calculations is thus minimized.)
      acc[ind] <- sapply(p[ind],poisson.blaker.acc.single.p,x=x,type=type,acc.tol=acc.tol,...=...)
      m1 <- length(which(ind1))
      if (m1 > 1) {
        acc[1:m1] <- cummax(acc[1:m1])
      }
      if (m1 < m - 1) {
        acc[m:(m1+1)] <- cummax(acc[m:(m1+1)])
      }
      return(acc)
    }
  }
}
