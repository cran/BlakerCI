binom.blaker.acc.single.p <- function(x,n,p,type="orig",acc.tol=1e-10,output="acc") {
# Reduce calculations to p <= x/n and search, thus,
# below just for the LEFT interval limits.
  if (p > x/n) {
    x <- n - x
    p <- 1 - p
  }
# "Ordinary" acceptance at p.
  p1.p <- 1 - pbinom(x-1,n,p)
  q1.p <- qbinom(p1.p,n,p)-1
  a1.p <- min(p1.p + pbinom(q1.p,n,p),1)
# "Unimodalization"
  if (type == "unimod" && q1.p >= 0) {
    upper <- p
    a1.upp <- a1.p
    lower <- 0
    a1.low <- 1
    q1.low <- -1
#   Search for the first discontinuity point of the
#   acceptance function to the left of p.
#   Continue as long as there is a chance to find there
#   a higher acceptance value than at p.
    while (a1.low > a1.p && a1.low - a1.upp >= acc.tol) {
      mid <- (lower+upper)/2
      p1.mid <- 1 - pbinom(x-1,n,mid)
      q1.mid <- qbinom(p1.mid,n,mid)-1
      a1.mid <- p1.mid + pbinom(q1.p,n,mid)
      if (q1.mid < q1.p) {
        lower <- mid
        a1.low <- a1.mid
        q1.low <- q1.mid
      }
      else {
        upper <- mid
        a1.upp <- a1.mid
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
