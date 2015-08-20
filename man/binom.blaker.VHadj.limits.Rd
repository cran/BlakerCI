\name{binom.blaker.VHadj.limits}
\alias{binom.blaker.VHadj.limits}

\title{
   Vos-Hudson adjustment of Blaker's binomial confidence limits
}

\description{
  Blaker's binomial confidence limits adjusted so that logical
  inconsistencies criticized by Vos and Hudson (2008) are avoided.
}
\usage{
binom.blaker.VHadj.limits(x, n, level = 0.95, tol = 1e-10, ...)
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
  \item{...}{
   additional arguments to be passed to \code{binom.blaker.VHadj.lower.limit}:
   \code{maxiter}, \code{nmax}, \code{int.eps}
  (see \code{\link{BlakerCI-internal}}). 
}

}
%\details{
%}

\value{
   Length 2 vector -- the lower and upper (adjusted) confidence    limits.
}
\references{
%% ~put references to the literature/web site here ~
Vos, P. W. & Hudson, S. (2008).
    Problems with binomial two-sided tests and the associated
    confidence intervals.
    Australian & New Zealand Journal of Statistics
    50(1): 81-89.

Lecoutre, B. & Poitevineau, J. (2014). New results for computing
   Blaker's exact confidence interval limits for usual
   one-parameter discrete distributions.
   Communications in Statistics - Simulation and Computation, 
   \url{http://dx.doi.org/10.1080/03610918.2014.911900}.
}

\author{
Jan Klaschka \email{klaschka@cs.cas.cz}
}

\note{
%%  ~~further notes~~
Vos & Hudson (2008) gave examples of mutually contradictory
inferences yielded by some binomial tests and confidence
intervals, including the Blaker's confidence interval.
Their objections may be interpreted as follows:
When the number of trials is increased so that the success
proportion increases (decreases) or remains the same,
the lower (upper) confidence limit at the same confidence level
should not decrease (increase).\cr
The adjustment implemented in \code{binom.blaker.VHadj.limits}
replaces the lower (upper) Blaker's confidence limit
for \code{x} successes in \code{n} trials
with the infimum (supremum) of the Blaker's lower (upper)
confidence limits over such pairs \code{y}, \code{m}
that \code{m} is not less that \code{n}, and 
\code{y/m} is not less (greater) than \code{x/n}.

Note that Lecoutre & Poitevineau (2014), refering to the 
criticism by
Vos & Hudson, proposed a modification of the Blaker's
confidence limits. 
Their adjustment, however, eliminates only a subset of
\dQuote{discrepancies} treated by \code{binom.blaker.VHadj.limits},
namely nonmonotonicities of upper (lower) Blaker's confidence
bounds in the number of trials when the number of successes
(failures) remains the same.
}

\section{Warning}{
The stopping rule used is not fully justified: 

The Clopper-Pearson \code{1 - alpha} confidence bounds
for \code{x} successes in \code{n} trials 
may be expressed as \code{qbeta(alpha/2,x,n-x+1)} and
\code{qbeta(1-alpha/2,x+1,n-x)},
and can be generalized this way to real 
(i. e. not only integer) values of \code{x}.

The stopping rule used in \code{binom.blaker.VHadj.limits}
relies on the hypothesis that
the generalized lower (upper) Clopper-Pearson confidence bounds 
grow (decrease) whenever the number of trials grows, 
and the proportion of successes grows (decreases) or 
remains unchanged (with obvious exceptions in extremes).

Though I firmly trust the hypothesis, I can prove it, so far, 
just for integer numbers of successes (i. e. for \dQuote{ordinary} 
Clopper-Pearson confidence bounds, not the generalized ones), 
and lack a general proof.
Should the hypothesis be invalid, the stopping rule
implemented in \code{binom.blaker.VHadj.limits}
would be incorrect, and the process of modifying
the Blaker's confidence bounds could be incomplete
in some cases. 




}

%% ~Make other sections like Warning with %\section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%\tabular{ll}{
%\code{exactci:exactbinomci} \tab
%      One of the options yields Blaker's
%      limits. The algorithm is more 
%      sophisti- \cr \tab cated than
%      the original Blaker's one. \cr
%\code{propCIs:blakerci, binGroup:binBlaker}
%\code{propCIs:blakerci} \tab
%      Implementation of the original algorithm from
%      Blaker (2000). \cr
%\code{binGroup:binBlaker} \tab
%      Another implementation of the same algorithm.
%}
%}

\examples{
binom.blaker.VHadj.limits(6,13) # [1] 0.2150187 0.7395922

## Note that the lower limit differs from the 
## unadjusted version:

binom.blaker.limits(6,13)       # [1] 0.2158050 0.7395922

## The (unadjusted) lower limit was replaced with the
## Blaker's lower limit (both unadjusted and adjusted) 
## assigned to 7 successes in 15 trials:

binom.blaker.limits(7,15)       # [1] 0.2150187 0.7096627
binom.blaker.VHadj.limits(7,15) # [1] 0.2150187 0.7096627

## The adjustment avoids a contradiction between
## inferences corresponding to
## 6 successes in 13 trials, and 7 successess in 15 trials: 
## Though the latter situation means a higher succes proportion
## in a higher number of trials, it is assigned a smaller
## (unadjusted) Blaker's 95\% lower confidence limit. 
}







% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
