\name{binom.blaker.limits}
\alias{binom.blaker.limits}

\title{
  Blaker's binomial confidence limits
}

\description{
  Fast and accurate calculation of Blaker's binomial confidence limits.
}
\usage{
binom.blaker.limits(x, n, level = 0.95, tol = 1e-10)
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

}
\details{
   Note that the Blaker's \code{(1 - alpha)} confidence interval 
   is the convex hull of the set \code{C} of those points 
   where the acceptance function (Blaker (2000)) exceeds 
   level \code{alpha}. The original numerical algorithm from
   Blaker (2000, 2001) is prone, when \code{C} is a union
   of disjoint intervals, to skipping a short interval
   and finding inaccurate over-optimistic confidence limits. 

   Function \code{\link{binom.blaker.limits}} is, by contrast,
   immune from such failures and yields always as 
   its result the whole confidence interval (Klaschka (2010)).
   The same is true for \code{exactci} package (user-level 
   function \code{binom.exact} and internal function 
   \code{exactbinomCI}). What makes \code{\link{binom.blaker.limits}} 
   still worth posting is the fact that it is simpler
   and calculates the Blaker's limits considerably faster, 
   especially for large \code{n}.
   (This may become important when, e.g., computing large tables
   of confidence limits.)
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

Fay, M. P. (2010). Confidence intervals that Match Fisher's exact 
   and Blaker's exact tests.
   Biostatistics. 11:373-374.

Fay, M.P. (2010a). Two-sided Exact Tests and Matching Confidence 
   Intervals for Discrete Data. 
   R Journal 2(1): 53-58.

Klaschka, J. (2010). On calculation of Blaker's binomial 
   confidence limits. COMPSTAT presentation (extended version).\cr
   \url{www.cs.cas.cz/~klaschka/c10/417_ext.pdf}
}

\author{
Jan Klaschka \email{klaschka@cs.cas.cz}
}

\note{
%%  ~~further notes~~
An early version written as a modification of the program
from Blaker (2000, 2001). At present, elements of original
functions \code{acceptbin} and \code{acceptinterval} are 
\sQuote{dissolved} within \code{\link{binom.blaker.lower.limit}} 
(internal function) and, in fact, throughout the package. 

Designed independently of M.P. Fay's works (otherwise
they would have been mentioned in \cite{Klaschka (2010)}). 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\tabular{ll}{
\code{exactci:exactbinomci} \tab
      one of the options yields similarly accurate Blaker's
      limits. \cr
\code{propCIs:blakerci, binGroup:binBlaker} \tab
      implementations of the original algorithm from
      Blaker (2000, 2001).
} 
}

\examples{
binom.blaker.limits(3,10) # [1] 0.08726443 0.61941066

## Example of a failure of the original algorithm:
## Requires PropCIs package.
## Tolerance 1e-4 -- default in the Blaker's paper.
\dontrun{
blakerci(29,99,conf.level=0.95,tolerance=1e-4) ## [1] 0.2096386 0.3923087
}
## The correct upper limit should be 0.3929\dots,
## as demonstrated:
## (1) By the same function with a smaller tolerance:
\dontrun{
blakerci(29,99,conf.level=0.95,tolerance=1e-7) ## [1] 0.2097022 0.3929079
} 
## (2) By binom.blaker.limits 
##     (default confidence limit 0.95, default tolerance 1e-10):
\dontrun{
binom.blaker.limits(29,99) ## [1] 0.2097022 0.3929079
}
## (3) By exactbinomCI function from package exactci
##     (default confidence level, default tolerance): 
\dontrun{
exactbinomCI(29,99,tsmethod="blaker")[1:2] ## [1] 0.2097 0.3929
}
## The same function, smaller tolerance:
\dontrun{
exactbinomCI(29,99,tsmethod="blaker",tol=1e-8)[1:2] ## [1] 0.2097022 0.3929079
}

## Another example of a failure of the original algorithm 
## with even as small tolerance as 1e-6:
\dontrun{
blakerci(59,355,conf.level=0.95,tolerance=1e-4) ## [1] 0.1299899 0.2085809
blakerci(59,355,conf.level=0.95,tolerance=1e-5) ## [1] 0.1300799 0.2085409
blakerci(59,355,conf.level=0.95,tolerance=1e-6) ## [1] 0.1300799 0.2085349
}
## Only for tolerance = 1e-7 the result is satisfactory
## and in agreement with binom.blaker.limits:
\dontrun{
blakerci(59,355,conf.level=0.95,tolerance=1e-7) ## [1] 0.1300807 0.2090809
binom.blaker.limits(59,355)                     ## [1] 0.1300807 0.2090809
}



}



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
