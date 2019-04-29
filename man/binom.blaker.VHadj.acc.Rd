\name{binom.blaker.VHadj.acc}
\alias{binom.blaker.VHadj.acc}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Blaker's binomial acceptability function with Vos-Hudson adjustment.
%, optionally unimodalized.
}

\description{
Calculates  values of the Vos-Hudson adjusted acceptability function
in a sequence of points (for, e.g., plotting purposes).
The adjusted acceptability function may optionally be 
\dQuote{unimodalized}, i.e. replaced with the smallest
greater or equal unimodal function. 
}
\usage{
binom.blaker.VHadj.acc(x, n, p, type = c("orig", "unimod"),
      acc.tol = 1e-10, nmax=n+1000,int.eps=1e-12, ...)
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
     In case of more than one point, an increasing sequence required. 
}
  \item{type}{
     for \code{type = "orig"}, Vos-Hudson adjustment applied to 
     original acceptability function.
     For \code{type = "unimod"}, smallest unimodal function
     greater or equal to the adjusted acceptability function.
}
  \item{acc.tol}{
     numerical tolerance (relevant only for \code{type = "unimod"}).
}
  \item{nmax}{
      Pairs \code{(y, m)} of number of trials and number of successes 
      are allowed to contribute to the Vos-Hudson adjustment for only
      \code{m} up to \code{nmax}. Warning is returned when greater
      numbers of trials are suspected to have influence.
}
  \item{int.eps}{
      Maximum expected error of machine representation of integers
      calculated from reals via multiplication and division.
      (Used in order to round numbers correctly if they happen
      to be integer, e. g.
      \code{ceiling(xx - int.eps)} is calculated instead of 
      \code{ceiling(xx)}.)
} 
  \item{...}{
     additional arguments to be passed to \code{binom.blaker.acc.single.p};
     in fact, just \code{maxiter} (see \code{BlakerCI-internal}). 
}
}

\details{
The relationship between the adjusted acceptability function 
and the adjusted confidence intervals 
(see \code{\link{binom.blaker.VHadj.limits}})
is the same as between the unadjusted acceptability function
and confidence interval (see \code{\link{binom.blaker.acc}},
\code{\link{binom.blaker.limits}}): The confidence interval is the
convex hull of the set of those points where the function
exceeds 1 - confidence level. 
}

\value{
Vector of Vos-Hudson adjusted acceptability values (with or without unimodalization) in points of \code{p}.
}
%\references{
%Blaker, H. (2000) Confidence curves and improved exact confidence 
%   intervals for discrete distributions.
%   Canadian Journal of Statistics 28: 783-798.\cr
%   (Corrigenda: Canadian Journal of Statistics 29: 681.)

%Fay, M. P. (2010). Confidence intervals that Match Fisher's exact 
%   and Blaker's exact tests.
%   Biostatistics. 11:373-374.

%Fay, M.P. (2010a). Two-sided Exact Tests and Matching Confidence 
%   Intervals for Discrete Data. 
%   R Journal 2(1): 53-58.

%Klaschka, J. (2010). On calculation of Blaker's binomial 
%   confidence limits. COMPSTAT presentation.\cr
%   \url{http://www-roc.inria.fr/axis/COMPSTAT2010/slides/slides_417_SP8.pdf}
%   \cr
%   Extended version: \url{www.cs.cas.cz/~klaschka/c10/417_ext.pdf}
%}
\author{
Jan Klaschka \email{klaschka@cs.cas.cz}
}
%%\note{
%%  ~~further notes~~
%%}
\section{Warning}{
(1) Comparing output of the function with that of 
\code{\link{binom.blaker.acc}} cannot answer positively the question 
whether the unadjusted and adjusted functions are identical 
on an interval (but, up to the numerical accuracy, in the points 
of \code{p} only).

(2) The Warning section of the \code{\link{binom.blaker.VHadj.limits}}
documentation is relevant here, as well.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
p <- seq(0,1,length=10001)
acc.adj <- binom.blaker.VHadj.acc(6,13,p)
acc <- binom.blaker.acc(6,13,p)

plot(p,acc.adj,type="l",col="red",ylab="acceptability"
  ,main=paste("Vos-Hudson adjustment of acceptability function"
  ,"for 6 successes in 13 trials"
  , sep="\n")
)
lines(p,acc,type="l")
legend(x=.7,y=.8,c("unadjusted","adjustment"),col=c("black","red"),lwd=1)

## Plot of differences between the unadjusted and adjusted
## acceptability functions reveals some adjustment details
## hardly visible in the previous graph. 

plot(p,acc.adj-acc,type="l",ylab="acceptability difference")

## The narrow peak near 0.215 is close to the 
## Blaker's lower 0.95 confidence limit.
## 
## Focussing on the neighbourhood of 0.215:

p <- seq(0.21,0.22,length=1001)
acc.adj <- binom.blaker.VHadj.acc(6,13,p)
acc <- binom.blaker.acc(6,13,p)

plot(p,acc.adj,type="l",col="red",ylab="acceptability"
  ,main=paste("A detail of Vos-Hudson adjustment of acceptability function"
  ,"for 6 successes in 13 trials"
  ,sep="\n")
  ,ylim=c(0.02,0.09)
)
lines(p,acc,type="l")
legend(x=.210,y=.08,c("unadjusted","adjustment"),col=c("black","red"),lwd=1)

## The above adjustment results from the fact that, though
## 15 > 13 and 7/15 > 6/13, the acceptability function
## for 7 successes in 15 trials is greater that that for 6 successes
## in 13 trials on a short interval:

acc.7.15 <- binom.blaker.acc(7,15,p)
plot(p,acc,type="l",ylab="acceptability"
  ,main=paste("A detail of acceptability functions"
  ,sep="\n")
  ,ylim=c(0.02,0.09)
)
lines(p,acc.7.15,type="l",col="green")
legend(x=.210,y=.08,c("6 / 13","7 / 15"),col=c("black","green")
  ,title="succ / trials",lwd=1)

## The adjustment shifts the point where the 0.05 level is exceeded,
## i. e. the Blaker's lower 0.95 confidence limit, from 0.2158 to 0.2150. 
## (Compare with Examples in binom.blaker.VHadj.limits section.)
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
