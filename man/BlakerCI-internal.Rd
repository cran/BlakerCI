\name{BlakerCI-internal}
\alias{BlakerCI-internal}
\alias{binom.blaker.lower.limit}
\alias{binom.blaker.acc.single.p}
\alias{binom.blaker.VHadj.lower.limit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal functions, not expected to be called by the user.
}
\description{
Calculation of the lower Blaker's confidence limit
as defined by Blaker (\code{binom.blaker.lower.limit}),
or with so called Vos-Hudson adjustment (\code{binom.blaker.VHadj.lower.limit});
a single acceptability value, optionally \dQuote{unimodalized}
(\code{binom.blaker.acc.single.p}).
}
\usage{
binom.blaker.lower.limit(x, n, level, tol = 1e-10, maxiter=100)
binom.blaker.VHadj.lower.limit(x,n,level,tol=1e-10,maxiter=100,
      nmax=n+1000,int.eps=1e-10)
binom.blaker.acc.single.p(x, n, p, type = "orig", acc.tol = 1e-10,
      output = "acc", maxiter=100)
}
%- maybe also 'usage' for other objects documented here.
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
      numerical tolerance (for the confidence limit).
}
  \item{p}{
      point (single value between 0 and 1) where to calculate the acceptability.
}
  \item{type}{
      \code{"orig"}, or \code{"unimod"} -- either unmodified, or unimodalized 
      acceptability \cr
      (see \code{\link{binom.blaker.acc}}).      
}
  \item{acc.tol}{
      numerical tolerance (for the acceptability values when \code{type = "unimod"}). 
}
  \item{output}{
      the acceptability value output (\code{output = "acc"}, the default)
      or, instead, an auxiliary integer-valued parameter \code{q1}, used
      for testing whether points belong to the same continuous segment 
      of the acceptability function
      (\code{output = "q1"}), or both (\code{output = "both"}).
}
  \item{maxiter}{
      Maximum number of interval halving iterations during the search
      for a confidence limit (\code{binom.blaker.lower.limit}) or
      a discontinuity point of the acceptability function 
      (\code{binom.blaker.acc.single.p} with \code{type = "unimod"}). 
      When the required accuracy is not reached in \code{maxiter} steps -- 
      typically when too small \code{tol} or \code{acc.tol}
      exceeds capabilities of machine arithmetic -- last step's result
      is returned with warning.
}

  \item{nmax}{
      Pairs \code{(y, m)} of number of trials and number of successes 
      are allowed to contribute to the Vos-Hudson adjustment for only
      \code{m} up to \code{nmax}. Warning is returned when greater
      numbers of trials are suspect to have influence.
}
  \item{int.eps}{
      Maximum expected error of machine representation of integers
      calculated via multiplication and division from reals.
      (Used in order to round numbers correctly if they happen
      to be integer, e. g.
      \code{ceiling(xx - int.eps)} is calculated instead of 
      \code{ceiling(xx)}.)
}

}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
For \code{binom.blaker.lower.limit} and
\code{binom.blaker.VHadj.lower.limit}, a single number -- 
the lower confidence limit. \cr
For \code{binom.blaker.acc.single.p} -- depending on the
\code{output} parameter -- a single acceptability value, or a single
auxiliary integer, or both.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Jan Klaschka \email{klaschka@cs.cas.cz}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{binom.blaker.limits}}, \code{\link{binom.blaker.VHadj.limits}},\cr
\code{\link{binom.blaker.acc}}, \code{\link{binom.blaker.VHadj.acc}}.
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
