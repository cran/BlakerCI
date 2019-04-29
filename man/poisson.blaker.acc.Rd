\name{poisson.blaker.acc}
\alias{poisson.blaker.acc}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Blaker's Poisson acceptability function, optionally
unimodalized.
}

\description{
Calculates values of the acceptability function for the Poisson
distribution (see Blaker (2000))
in a sequence of points (for, e.g., plotting purposes).
The acceptability function may optionally be 
\dQuote{unimodalized}, i.e. replaced with the smallest
greater or equal unimodal function. 
}
\usage{
poisson.blaker.acc(x, p, type = c("orig", "unimod"),
      acc.tol = 1e-10, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
     number of events.
}
  \item{p}{
     vector (length 1 allowed) of hypothesized Poisson parameters.
     In case of more than one point, an increasing sequence required. 
}
  \item{type}{
     for \code{type = "orig"}, original acceptability function
     calculated.\cr
     For \code{type = "unimod"}, smallest unimodal function
     greater or equal to the acceptability function calculated
     instead. 
}
  \item{acc.tol}{
     numerical tolerance (relevant only for \code{type = "unimod"}).
}
  \item{...}{
     additional arguments to be passed to 
     \code{\link{poisson.blaker.acc.single.p}};
     in fact, just \code{maxiter} (see \code{\link{BlakerCI-internal}}). 
}
}

\details{
Single values of the \dQuote{unimodalized} acceptability function 
(for \code{type = "unimod"}) are computed by an iterative
numerical algorithm implemented in internal function
\cr
\code{poisson.blaker.acc.single.p}.
The function cited is called just once in each of the intervals
where the acceptability function is continuous
(namely in the leftmost one of those points of \code{p} 
that fall into the interval when dealing with points
below \code{x}, and the rightmost one when above 
\code{x}). The rest is done by function
\code{cummax}.  
This is considerably faster than calling 
\code{poisson.blaker.acc.single.p}
for every point of \code{p}.
Note that applying \code{cummax} directly to
a vector of unmodified acceptability values
is even faster and provides a unimodal output;
it may, nevertheless, lack accuracy. 
  
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

%Fay, M. P. (2010). Confidence intervals that Match Fisher's exact 
%   and Blaker's exact tests.
%   Biostatistics. 11:373-374.

%Fay, M.P. (2010a). Two-sided Exact Tests and Matching Confidence 
Fay, M.P. (2010). Two-sided Exact Tests and Matching Confidence 
   Intervals for Discrete Data. 
   R Journal 2(1): 53-58.

%Klaschka, J. (2010). On calculation of Blaker's binomial 
%   confidence limits. COMPSTAT presentation.\cr
%   \url{http://www-roc.inria.fr/axis/COMPSTAT2010/slides/slides_417_SP8.pdf}
%   \cr
%   Extended version: \url{www.cs.cas.cz/~klaschka/c10/417_ext.pdf}
}
\author{
Jan Klaschka \email{klaschka@cs.cas.cz}
}
\note{
%%  ~~further notes~~
Inspired by M.P. Fay (2010), mentioning
\dQuote{unavoidable inconsistencies} between tests with non-unimodal
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
p <- seq(0,10,length=1001)
acc <- poisson.blaker.acc(3,p)
acc1 <- poisson.blaker.acc(3,p,type="unimod")
plot(p,acc,type="l")
lines(p,acc1,col="red")
legend(x=7,y=.8,c("orig","unimod"),col=c("black","red"),lwd=1)

## The two lines -- the unimodalized and original acceptabilities -- 
## look almost the same but some small differences are slightly 
## visible.

## They can be seen better this way:
plot(p,acc1-acc,type="l")

## Focussing on one of them:
p <- seq(5.05,5.6,length=1001)
acc <- poisson.blaker.acc(3,p)
acc1 <- poisson.blaker.acc(3,p,type="unimod")
plot(p,acc,type="l",ylim=c(.391,.396))
lines(p,acc1,col="red")
legend(x=5.4,y=.395,c("orig","unimod"),col=c("black","red"),lwd=1)

}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
