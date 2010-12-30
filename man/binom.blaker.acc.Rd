\name{binom.blaker.acc}
\alias{binom.blaker.acc}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Blaker's binomial acceptability function, optionally
unimodalized.
}

\description{
Calculates values of the acceptability function for the binomial
distribution (see function \code{acceptbin} in Blaker (2000))
in a sequence of points (for, e.g., plotting purposes).
The acceptability function may optionally be 
\sQuote{unimodalized}, i.e. replaced with the smallest
greater or equal unimodal function. 
}
\usage{
binom.blaker.acc(x, n, p, type = c("orig", "unimod"),
      acc.tol = 1e-10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
     number of successes.
}
  \item{n}{
     number of trials.
}
  \item{p}{
     vector (length 1 allowed) of hypothesized binomial parameters
     (between 0 and 1).
     In case of more than one point, a growing sequence required. 
}
  \item{type}{
     for \code{type = "orig"}, original acceptability function
     calculated.
     For \code{type = "unimod"}, smallest unimodal function
     greater or equal to the acceptability function calculated
     instead. 
}
  \item{acc.tol}{
     numerical tolerance (relevant only for \code{type = "unimod"}).
}
}

\details{
For \code{type = "orig"}, essentially the same is calculated
as -- for single points -- by \code{acceptbin} function
from Blaker (2000).

Single values of the \sQuote{unimodalized} acceptability function 
(for \code{type = "unimod"}) are computed by an iterative
numerical algorithm implemented in internal function
\cr
\code{\link{binom.blaker.acc.single.p}}.
The function cited is called just once in each of the intervals
where the acceptability function is continuous
(namely in the leftmost one of those points of \code{p} 
that fall into the interval when dealing with points
below \code{x/n}, and the rightmost one when above 
\code{x/n}). The rest is done by function
\code{cummax}.  
This is considerably faster than calling 
\code{\link{binom.blaker.acc.single.p}}
for every point of \code{p}.
Note that applying \code{cummax} directly to
a vector of unmodified acceptability values
is even faster and provides a unimodal output;
it may, nevertheless, lack accuracy (see Examples). 
  
}

\value{
Vector of acceptability values (with or without unimodalization)
in points of \code{p}.
}
\references{
Blaker, H. (2000) Confidence curves and improved exact confidence 
   intervals for discrete distributions.
   Canadian Journal of Statistics 28: 783-798.\cr
   (Corrigenda: Canadian Journal of Statistics 29: 681.)

Fay, M. P. (2010). Confidence intervals that Match Fisher's exact 
   and Blaker's exact tests.
   Biostatistics. 11:373-374.

Fay, M.P. (2010a). Two-sided Exact Tests and Matching Confidence 
   Intervals for Discrete Data. 
   R Journal 2(1): 53-58.

Klaschka, J. (2010). On calculation of Blaker's binomial 
   confidence limits. COMPSTAT presentation.\cr
   \url{http://www-roc.inria.fr/axis/COMPSTAT2010/slides/slides_417_SP8.pdf}
   \cr
   Extended version: \url{www.cs.cas.cz/~klaschka/c10/417_ext.pdf}
}
\author{
Jan Klaschka \email{klaschka@cs.cas.cz}
}
\note{
%%  ~~further notes~~
Inspired by M.P. Fay (2010a), mentioning
\sQuote{unavoidable inconsistencies} between tests with non-unimodal
acceptability functions and confidence intervals derived from them.
When the acceptability functions are unimodalized
and the test modified accordingly (i.e. p-values slightly increased
in some cases), a perfectly matching test-CI pair is obtained. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
p <- seq(0,1,length=1001)
acc <- binom.blaker.acc(3,10,p)
acc1 <- binom.blaker.acc(3,10,p,type="unimod")
## The two functions look the same at first glance.
plot(p,acc,type="l")
lines(p,acc1,col="red")
legend(x=.7,y=.8,c("orig","unimod"),col=c("black","red"),lwd=1)

## There is, nevertheless, a difference.
plot(p,acc1-acc,type="l")

## Focussing the difference about p=0.4:
p <- seq(.395,.405,length=1001)
acc <- binom.blaker.acc(3,10,p)
acc1 <- binom.blaker.acc(3,10,p,type="unimod")
plot(p,acc,type="l",ylim=c(.749,.7495))
lines(p,acc1,col="red")
legend(x=.402,y=.7494,c("orig","unimod"),col=c("black","red"),lwd=1)

## Difference between type="unimod" and mere applying
## cummax to values obtained via type="orig":
p <- seq(0,1,length=1001)
x <- 59
n <- 355
## Upper confidence limit (at 0.95 level) is slightly above 0.209:
binom.blaker.limits(x,n) ## [1] 0.1300807 0.2090809
## Unmodified acceptability value fall below 0.05 at p = .209
## left to the limit (so that the null hypothesis p = .209
## would be rejected despite the fact that p lies within
## the confidence interval):
acc <- binom.blaker.acc(59,355,p)
rbind(p,acc)[,207:211]
##           [,1]       [,2]       [,3]       [,4]       [,5]
## p   0.20600000 0.20700000 0.20800000 0.20900000 0.21000000
## acc 0.06606867 0.05759836 0.05014189 0.04999082 0.04330283
##
## Modified acceptability is above 0.05 at p = 0.05 (so that
## hypothesis p = 0.05 is not rejected by the modified test):
acc1 <- binom.blaker.acc(59,355,p,type="unimod")
rbind(p,acc1)[,207:211]
##            [,1]       [,2]       [,3]       [,4]       [,5]
## p    0.20600000 0.20700000 0.20800000 0.20900000 0.21000000
## acc1 0.06608755 0.05759836 0.05014189 0.05000009 0.04331409
##
## Applying cummax to unmodified acceptabilities guarantees unimodality 
## but lacks accuracy, leaving the value at p = 0.209 below 0.05:
m <- max(which(p <= 59/355))
acc2 <- acc
acc2[1:m] <- cummax(acc2[1:m])
acc2[1001:(m+1)] <- cummax(acc2[1001:(m+1)])
rbind(p,acc2)[,207:211]
##            [,1]       [,2]       [,3]       [,4]       [,5]
## p    0.20600000 0.20700000 0.20800000 0.20900000 0.21000000
## acc2 0.06606867 0.05759836 0.05014189 0.04999082 0.04330283

}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
