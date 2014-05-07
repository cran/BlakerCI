\name{binom.blaker.limits}
\alias{binom.blaker.limits}

\title{
  Blaker's binomial confidence limits
}

\description{
  Fast and accurate calculation of Blaker's binomial confidence limits.
}
\usage{
binom.blaker.limits(x, n, level = 0.95, tol = 1e-10, \dots)
}

\arguments{
  \item{x}{
   number of successes.
}
  \item{n}{
   number of trials.
}
  \item{level}{
   confidence level.
}
  \item{tol}{
   numerical tolerance.
}
  \item{\dots}{
   additional arguments to be passed to \code{binom.blaker.lower.limit};
   in fact, just \code{maxiter} (see \code{\link{BlakerCI-internal}}). 
}

}
\details{
   Note that the Blaker's \code{(1 - alpha)} confidence interval 
   is the convex hull of the set \code{C} of those points 
   where the acceptability function (Blaker (2000)) exceeds 
   level \code{alpha}. The original numerical algorithm from
   Blaker (2000) is prone, when \code{C} is a union
   of disjoint intervals, to skipping a short interval
   and finding inaccurate over-liberal confidence limits. 

   Function \code{binom.blaker.limits} is, by contrast,
   immune from such failures and yields always as 
   its result the whole confidence interval (Klaschka (2010)).

   Package \code{exactci} by M. P. Fay includes another algorithm
   that calculates Blaker's binomial confidence limits
   (see user-level function \code{binom.exact} and internal function 
   \code{exactbinomCI}). 
   It is more sophisticated than the original Blaker's one,
   but considerably slower and sometimes less accurate 
   than that of \code{binom.blaker.limits}.
}

\value{
   Length 2 vector -- the lower and upper confidence limits.
}
\references{
%% ~put references to the literature/web site here ~
Blaker, H. (2000) Confidence curves and improved exact confidence 
   intervals for discrete distributions.
   Canadian Journal of Statistics 28: 783-798.\cr
   (Corrigenda: Canadian Journal of Statistics 29: 681.)

%Fay, M. P. (2010). Confidence intervals that Match Fisher's exact 
%   and Blaker's exact tests.
%   Biostatistics. 11:373-374.

%Fay, M.P. (2010a). Two-sided Exact Tests and Matching Confidence 
%   Intervals for Discrete Data. 
%   R Journal 2(1): 53-58.

Klaschka, J. (2010). BlakerCI: An algorithm and R package
   for the Blaker's binomial confidence limits calculation.
   Technical report No. 1099, Institute of Computer Science,
   Academy of Sciences of the Czech Republic.\cr
\url{http://www.library.sk/i2/content.csg.cls?ictx=cav&repo=crepo1&key=99673653505}
}

\author{
Jan Klaschka \email{klaschka@cs.cas.cz}
}

\note{
%%  ~~further notes~~
Earlier 2010 algorithm versions were designed independently of 
(though already existing) M.P. Fay's packages \code{exact2x2}
and \code{exactci}.
Some later modifications, however, have been inspired
by Fay's programs. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\tabular{ll}{
\code{exactci:exactbinomci} \tab
      One of the options yields Blaker's
      limits. The algorithm is more \cr \tab
      sophisticated than
      the original Blaker's one. \cr
%\code{propCIs:blakerci, binGroup:binBlaker}
\code{propCIs:blakerci} \tab
      Implementation of the original algorithm from
      Blaker (2000). \cr
\code{binGroup:binBlaker} \tab
      Another implementation of the same algorithm.
}
}

\examples{
binom.blaker.limits(3,10) # [1] 0.08726443 0.61941066

## Example of a failure of the original algorithm:
## Requires PropCIs package.
## Tolerance 1e-4 - default in the Blaker's paper.
\dontrun{
blakerci(29,99,conf.level=0.95,tolerance=1e-4) ## [1] 0.2096386 0.3923087
## The correct upper limit should be 0.3929\dots,
## as demonstrated:
## (1) By the same function with a smaller tolerance:
blakerci(29,99,conf.level=0.95,tolerance=1e-7) ## [1] 0.2097022 0.3929079
## (2) By binom.blaker.limits 
##     (default confidence limit 0.95, default tolerance 1e-10):
binom.blaker.limits(29,99) ## [1] 0.2097022 0.3929079
## (3) By exactbinomCI function from package exactci
##     (default confidence level, default tolerance): 
exactbinomCI(29,99,tsmethod="blaker")[1:2] ## [1] 0.2097 0.3929
## The same function, smaller tolerance:
exactbinomCI(29,99,tsmethod="blaker",tol=1e-8)[1:2] 
                                                ## [1] 0.2097022 0.3929079

## Another example of a failure of the original algorithm 
## with even as small tolerance as 1e-6:
blakerci(59,355,conf.level=0.95,tolerance=1e-4) ## [1] 0.1299899 0.2085809
blakerci(59,355,conf.level=0.95,tolerance=1e-5) ## [1] 0.1300799 0.2085409
blakerci(59,355,conf.level=0.95,tolerance=1e-6) ## [1] 0.1300799 0.2085349
## Only for tolerance = 1e-7 the result is satisfactory
## and in agreement with binom.blaker.limits:
blakerci(59,355,conf.level=0.95,tolerance=1e-7) ## [1] 0.1300807 0.2090809
binom.blaker.limits(59,355)                     ## [1] 0.1300807 0.2090809
}



}



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
